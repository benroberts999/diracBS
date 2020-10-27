#pragma once
#include "DiracOperator/Operators.hpp"
#include "MixedStates.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/Format.hpp"
#include <iomanip>
#include <string>

namespace UnitTest {

//******************************************************************************
namespace helper {

struct Case {
  double eps;
  std::string name = "";
  double w = 0.0;
};

// Compares: <B||dA> to <B||h||A> / (e_A - e_B + w)
// Where dA is solution to: (H - e_A - w)|dA> = -(h - de_A)|A>
// Calculates above comparison for each pair of valence states (A,B), for which
// the matrix element (and denominator) is non-zero. Also loops over different w
// (frequency) values.
// Returns the best and worst cases
inline std::pair<Case, Case> MS_loops(const Wavefunction &wf,
                                      const DiracOperator::TensorOperator *h);

} // namespace helper

//******************************************************************************
//******************************************************************************
//! Unit tests Mixed States (TDHF method, solving non-local DE)
bool MixedStates(std::ostream &obuff) {
  bool passQ = true;

  // Create wavefunction object
  Wavefunction wf({10000, 1.0e-7, 200.0, 0.33 * 200.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.hartreeFockCore("HartreeFock", 0.0, "[Xe]");
  wf.hartreeFockValence("7sp5d4f");

  // Define E1, E2, pnc, M1, and hyperfine operators:
  const auto hE1 = DiracOperator::E1(*wf.rgrid);
  const auto hE2 = DiracOperator::Ek(*wf.rgrid, 2);
  const auto c = Nuclear::c_hdr_formula_rrms_t(wf.get_rrms());
  const auto hpnc = DiracOperator::PNCnsi(c, Nuclear::default_t, *wf.rgrid);
  const auto hM1 = DiracOperator::M1(*wf.rgrid, wf.alpha, 0.0);
  // Use "spherical ball" model for hyperfine (doesn't work with zero-size)
  const auto hhfs = DiracOperator::Hyperfine(
      1.0, 1.0, std::sqrt(5.0 / 3) * wf.get_rrms(), *wf.rgrid,
      DiracOperator::Hyperfine::sphericalBall_F());

  const std::vector<const DiracOperator::TensorOperator *> hs{&hE1, &hE2, &hpnc,
                                                              &hhfs, &hM1};

  // For each operator
  for (const auto h : hs) {

    const auto [best, worst] = helper::MS_loops(wf, h);

    const auto wbest = qip::fstring("%.2g", best.w);
    const auto wworst = qip::fstring("%.2g", worst.w);

    passQ &= qip::check_value(&obuff,
                              "MS:" + h->name() + " worst (" + worst.name +
                                  " w=" + wworst + ")",
                              worst.eps, 0.0, 4e-5);
    passQ &= qip::check_value(
        &obuff, "MS:" + h->name() + " best (" + best.name + " w=" + wbest + ")",
        best.eps, 0.0, 6e-7);
  }

  return passQ;
}

} // namespace UnitTest

//******************************************************************************
inline std::pair<UnitTest::helper::Case, UnitTest::helper::Case>
UnitTest::helper::MS_loops(const Wavefunction &wf,
                           const DiracOperator::TensorOperator *h) {
  // compare the best and worst! (of each operator)
  Case best{999.0};
  Case worst{0.0};

  const auto omega_mults = std::vector{0.0, 0.25};

#pragma omp parallel for
  for (auto i = 0ul; i < wf.valence.size(); ++i) {
    const auto &Fv = wf.valence[i];
    const auto vl = wf.get_Vlocal(Fv.l());
    for (const auto &Fm : wf.valence) {

      // Only do for cases that make sense:
      if (Fm == Fv)
        continue; // gives 0/0
      if (h->isZero(Fv.k, Fm.k))
        continue;
      if (Fm.k != Fv.k && std::abs(Fm.n - Fv.n) != 0)
        continue;

      const auto h_mv = h->reducedME(Fm, Fv);

      // Don't do in cases where ME is extremely small. NOTE: This requires
      // a good choice of units, since ME may be small due to small
      // coeficient, which should not count as 'small'
      if (std::abs(h_mv) < 1.0e-8) {
        continue;
      }

      // Form 'rhs': (h - de_A)|A>
      auto hFv = h->reduced_rhs(Fm.k, Fv);
      if (Fm.k == Fv.k) {
        auto de = h->reducedME(Fv, Fv);
        hFv -= de * Fv;
      }

      // loop over few frequencies:
      for (const auto &w_mult : omega_mults) {
        const auto w = std::abs(Fv.en * w_mult);

        const auto dFv =
            HF::solveMixedState(Fm.k, Fv, w, vl, wf.alpha, wf.core, hFv);

        const auto lhs = Fm * dFv;
        const auto rhs = h_mv / (Fv.en - Fm.en + w);
        const auto eps = (lhs - rhs) / (lhs + rhs);

// find the best and worst case:
#pragma omp critical(find_max)
        {
          if (std::abs(eps) > std::abs(worst.eps)) {
            worst.eps = eps;
            worst.name = Fm.shortSymbol() + "|" + Fv.shortSymbol();
            worst.w = w;
          }
          if (std::abs(eps) < std::abs(best.eps)) {
            best.eps = eps;
            best.name = Fm.shortSymbol() + "|" + Fv.shortSymbol();
            best.w = w;
          }
        }
      }
    }
  }
  return {best, worst};
}
