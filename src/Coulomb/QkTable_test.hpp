#pragma once
#include "Coulomb/QkTable.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"

namespace UnitTest {

//******************************************************************************
//******************************************************************************
//! Unit tests for Coulomb integrals (y^k_ab, R^k_abcd, lookup tables etc).
//! Also: tests quadrature integation method
bool QkTable(std::ostream &) {
  bool pass = true;

  using namespace Coulomb;

  Coulomb::QkTable<float, Symmetry::Rk> qk{};

  Wavefunction wf({1000, 1.0e-5, 50.0, 10.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.hartreeFockCore("HartreeFock", 0.0, "[Ne]");
  wf.formBasis({"15spdf", 30, 7, 1.0e-5, 1.0e-6, 30.0, false});

  // Form the Coulomb lookup tables:
  const Coulomb::YkTable yk(wf.rgrid, &wf.basis);

  std::cout << "Filling ... " << std::flush;
  {
    IO::ChronoTimer t("\nFill");
    qk.fill(yk);
  }
  std::cout << "done\n\n" << std::flush;

  qk.count();
  std::cout << "\n\n";

  // Why is updating so slow??? Should be just as fast as adding?????
  // No. Because map starts out very small, so it is quick to check if element
  // exists already. As map gets larger, this is slower. On 'update' the map is
  // the maximum size right from the start

  // std::cout << "Updating ... " <<
  // std::flush;
  // {
  //   IO::ChronoTimer t("\nUpdate");
  //   qk.update(yk);
  // }
  // std::cout << "done\n" << std::flush;

  {
    IO::ChronoTimer t("Direct calc");
    double sum1 = 0.0;
    for (const auto &a : wf.basis) {
      for (const auto &b : wf.basis) {
        for (const auto &c : wf.basis) {
          for (const auto &d : wf.basis) {
            const auto [kmin, kmax] = Coulomb::k_minmax_Q(a, b, c, d);
            for (int k = kmin; k <= kmax; k += 2) {
              const auto yk_bd = yk.ptr_yk_ab(k, b, d);
              if (yk_bd == nullptr)
                continue;
              sum1 += Coulomb::Qk_abcd(a, b, c, d, k, *yk_bd, yk.Ck());
            }
          }
        }
      }
    }
    std::cout << "sum1=" << sum1 << "\n" << std::flush;
  }
  std::cout << "\n";

  {
    IO::ChronoTimer t("Use table");
    double sum2 = 0.0;
    for (const auto &a : wf.basis) {
      for (const auto &b : wf.basis) {
        for (const auto &c : wf.basis) {
          for (const auto &d : wf.basis) {
            const auto Qk = qk.Qk(a, b, c, d);
            if (Qk == nullptr)
              continue;
            for (const auto q : *Qk) {
              sum2 += double(q);
            }
          }
        }
      }
    }
    std::cout << "sum2=" << sum2 << "\n" << std::flush;
  }
  std::cout << "\n";

  {
    IO::ChronoTimer t("Use table, in-order access");
    double sum3 = 0.0;
    for (const auto &[key, value] : qk) {
      for (const auto q : value) {
        sum3 += double(q);
      }
    }
    std::cout << "sum3=" << sum3 << "\n" << std::flush;
  }
  std::cout << "\n";

  return pass;
}

} // namespace UnitTest
