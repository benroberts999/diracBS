#pragma once
#include "Angular/Angular_tables.hpp"
#include "qip/Check.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <string>

namespace UnitTest {

namespace helper {

// This loops through every possibly C^k and three-j symbol (up to given maximum
// j), and calculates {C^k, 3j} symbol multiple ways (using lookup table, +
// basic formulas etc.). Checks all these against each other, and returns the
// eps of the worst comparison
double ck_loop(const Angular::Ck_ab &Ck, Angular::Ck_ab &Ck_m);

// This does the same, but for the 6j symbols
double sj_loop(const Angular::SixJ &sj, Angular::SixJ &sj_m);

} // namespace helper

//******************************************************************************
//! Unit tests for angular functions/classes (threeJ symbols, lookup tables etc)
bool Angular(std::ostream &obuff) {
  bool pass = true;

  {
    // Maximum value of 2*j (for initial run)
    const int max2j_1 = 5;
    // Form C^k lookup tables. One used statically, one dynamically re-sized
    // ("dynamic" one uses the _mutable lookup functions, that calclate the
    // angular factor if it doesn't exist already)
    Angular::Ck_ab Ck(max2j_1);
    Angular::Ck_ab Ck_dynamic(0);

    Angular::SixJ sj(max2j_1);
    Angular::SixJ sj_dynamic(0);

    const auto eps = helper::ck_loop(Ck, Ck_dynamic);
    pass &= qip::check_value(&obuff, "AngularTables Ck", eps, 0.0, 1.0e-14);

    const auto eps6 = helper::sj_loop(sj, sj_dynamic);
    pass &= qip::check_value(&obuff, "AngularTables 6j", eps6, 0.0, 1.0e-14);

    // Test the "extend" capability (extend tables to larger j values:)
    const int max2j_2 = 15;

    // "Extend" the tables:
    Ck.fill(max2j_2);
    sj.fill(max2j_2);

    const auto eps2 = helper::ck_loop(Ck, Ck_dynamic);
    pass &= qip::check_value(&obuff, "AngularTables Ck - extend", eps2, 0.0,
                             1.0e-13);

    const auto eps62 = helper::sj_loop(sj, sj_dynamic);
    pass &= qip::check_value(&obuff, "AngularTables 6j - extend", eps62, 0.0,
                             1.0e-14);
  }

  return pass;
}

} // namespace UnitTest

//******************************************************************************

double UnitTest::helper::ck_loop(const Angular::Ck_ab &Ck,
                                 Angular::Ck_ab &Ck_m) {
  const auto max2j = Ck.max_tj();
  double max = 0.0;

  for (int kia = 0;; ++kia) {
    const auto tja = Angular::twojFromIndex(kia);
    if (tja > max2j)
      break;
    const auto ka = Angular::kappaFromIndex(kia);
    for (int kib = 0;; ++kib) {
      const auto tjb = Angular::twojFromIndex(kib);
      if (tjb > max2j)
        break;
      const auto kb = Angular::kappaFromIndex(kib);

      // loop through all k multipolarities:
      for (int k = 0; k <= max2j; ++k) {

        const auto c1 = Ck.get_Ckab(k, ka, kb);
        const auto c1_m = Ck_m.get_Ckab_mutable(k, ka, kb);
        const auto c2 = Angular::Ck_kk(k, ka, kb);
        const auto sign = ((tja + 1) / 2) % 2 == 0 ? 1 : -1;
        const auto c3 = Angular::parity(Angular::l_k(ka), Angular::l_k(kb), k) *
                        Angular::threej_2(tja, tjb, 2 * k, -1, 1, 0) * sign *
                        std::sqrt((tja + 1) * (tjb + 1));

        const auto tc1 = Ck.get_tildeCkab(k, ka, kb);
        const auto tc2 =
            Angular::parity(Angular::l_k(ka), Angular::l_k(kb), k) *
            Angular::threej_2(tja, tjb, 2 * k, -1, 1, 0) *
            std::sqrt((tja + 1) * (tjb + 1));

        const auto tj1 = Ck.get_3jkab(k, ka, kb);
        const auto tj2 = Angular::threej_2(tja, tjb, 2 * k, -1, 1, 0);
        const auto tj3 = gsl_sf_coupling_3j(tja, tjb, 2 * k, -1, 1, 0);

        const auto l1 = Ck.get_Lambdakab(k, ka, kb);
        const auto l2 =
            tj1 * tj1 * Angular::parity(Angular::l_k(ka), Angular::l_k(kb), k);

        max = qip::max_abs(max, c1 - c2, c1_m - c1, c2 - c3, tc1 - tc2,
                           tj1 - tj2, tj2 - tj3, l1 - l2);
      }
    }
  }
  return max;
}

double UnitTest::helper::sj_loop(const Angular::SixJ &sj, Angular::SixJ &sj_m) {
  const auto max2j = sj.max_tj();

  double max = 0.0;

  // Loop through all 2j's:
  for (int tja = 1; tja <= max2j; tja += 2) {
    for (int tjb = 1; tjb <= max2j; tjb += 2) {
      for (int tjc = 1; tjc <= max2j; tjc += 2) {
        for (int tjd = 1; tjd <= max2j; tjd += 2) {

          // loop through all k multipolarities:
          for (int k = 0; k <= max2j; ++k) {
            for (int l = 0; l <= max2j; ++l) {
              // Calculate 6j symbol, multiple ways:
              const auto s1 = sj.get_6j(tja, tjb, tjc, tjd, k, l);
              const auto s1_m = sj_m.get_6j_mutable(tja, tjb, tjc, tjd, k, l);
              const auto s2 = Angular::sixj_2(tja, tjb, 2 * k, tjc, tjd, 2 * l);
              const auto s3 =
                  gsl_sf_coupling_6j(tja, tjb, 2 * k, tjc, tjd, 2 * l);
              max = qip::max_abs(max, s1 - s1_m, s1 - s2, s2 - s3);
            }
          }
        }
      }
    }
  }
  return max;
}
