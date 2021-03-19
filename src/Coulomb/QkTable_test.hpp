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
bool QkTable(std::ostream &obuff) {
  bool pass = true;

  using namespace Coulomb;

  Coulomb::QkTable<float, Symmetry::Qk> qk{};

  Wavefunction wf({1000, 1.0e-5, 50.0, 10.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.hartreeFockCore("HartreeFock", 0.0, "[Ne]");
  wf.formBasis({"20spdfg", 30, 7, 1.0e-5, 1.0e-6, 30.0, false});

  // Form the Coulomb lookup tables:

  // YkTable stores Hartee Y-functions Y_ab(r)
  // These save much time when calculating Q^k coeficients
  const Coulomb::YkTable yk(wf.rgrid, &wf.basis);

  std::cout << "Filling ... " << std::flush;
  {
    IO::ChronoTimer t("Fill");
    qk.fill(yk);
  }
  std::cout << "done\n\n" << std::flush;

  qk.count();
  std::cout << "\n\n";

  return true;

  // Compare the speed of using Qk lookup table vs. direct calculation
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
    // XXX This is *faster* when there are no symmetries
    // (Filling is slower, as expected)
    // Only possible reason is NormalOrder?? But, on profiling, doesn't account
    // for the difference??
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
    std::cout << "sum2=" << sum2 << " (should = sum1)\n" << std::flush;
  }
  std::cout
      << "Why is 'use' _slower_ when accounting for symmetry? NormalOrder?\n";
  std::cout << "\n";

  // If we just need to access each Q, we can take advantage of map, and get
  // O(1) lookup by just going through each element. In practise, I doubt this
  // is ever actually useful..
  {
    IO::ChronoTimer t("Use table, in-order access");
    double sum3 = 0.0;
    for (const auto &[key, value] : qk) {
      for (const auto q : value) {
        sum3 += double(q);
      }
    }
    std::cout << "sum3=" << sum3 << " (does not = sum1)\n" << std::flush;
    // nb: sum3 not equal to sum1, since sum1 sums over many equivilant Q^k.
    // If symmetry = none, then should be equal
  }
  std::cout << "\n";

  // Here begins actual unit tests. Ensure the Q (and P,W,R) values
  // stores/derived from table match those calculated from scratch (noting that
  // stored in Qk as float, but calc'd as double, so may be rounding errors)

  {
    // Test the table, including working out which Qk each one is
    // i.e., use stored index's to reconstruct which orbitals used, use them
    float max_dev = 0.0;
    for (const auto &[key, value] : qk) {
      const auto [a, b, c, d] = key;
      const auto match_index = [](auto i) {
        return [=](const auto &f) { return f.nk_index() == i; };
      };
      const auto &Fa =
          *std::find_if(wf.basis.cbegin(), wf.basis.cend(), match_index(a));
      const auto &Fb =
          *std::find_if(wf.basis.cbegin(), wf.basis.cend(), match_index(b));
      const auto &Fc =
          *std::find_if(wf.basis.cbegin(), wf.basis.cend(), match_index(c));
      const auto &Fd =
          *std::find_if(wf.basis.cbegin(), wf.basis.cend(), match_index(d));
      const auto [kmin, kmax] = Coulomb::k_minmax_Q(Fa, Fb, Fc, Fd);
      auto qk_it = value.cbegin();
      for (int k = kmin; k <= kmax; k += 2) {
        const auto qk0 = Coulomb::Qk_abcd(Fa, Fb, Fc, Fd, k);
        const auto dev = std::abs(float(*qk_it) - float(qk0));
        if (dev > max_dev)
          max_dev = dev;
        ++qk_it;
      }
    }
    pass &= qip::check_value(&obuff, "QkTable", max_dev, 0.0f, 1.0e-5f);
  }

  {
    // Test table versions of Q,P,W,R vs. 'from scratch' ones.
    // (only for subset, since otherwise v. slow)
    float max_devQR = 0.0;
    float max_devPW = 0.0;
    for (const auto &a : wf.basis) {
      for (auto ib = 0ul; ib < wf.basis.size(); ib += 2) {
        const auto &b = wf.basis[ib];
        for (auto ic = 1ul; ic < wf.basis.size(); ic += 3) {
          const auto &c = wf.basis[ic];
          for (auto id = 0ul; id < wf.basis.size(); id += 4) {
            const auto &d = wf.basis[id];
            const auto [kmin, kmax] = Coulomb::k_minmax_Q(a, b, c, d);
            for (int k = 0; k < kmax + 3; ++k) {
              const auto q1 = qk.Q(k, a, b, c, d);
              const auto p1 = qk.P(k, a, b, c, d);
              const auto w1 = qk.W(k, a, b, c, d);
              const auto r1 = qk.R(k, a, b, c, d);
              const auto q0 = Coulomb::Qk_abcd(a, b, c, d, k);
              const auto p0 = Coulomb::Pk_abcd(a, b, c, d, k);
              const auto w0 = Coulomb::Wk_abcd(a, b, c, d, k);
              const auto r0 = Coulomb::Rk_abcd(a, b, c, d, k);
              const auto devQ = std::abs(float(q1) - float(q0));
              const auto devP = std::abs(float(p1) - float(p0));
              const auto devW = std::abs(float(w1) - float(w0));
              const auto devR = std::abs(float(r1) - float(r0));
              const auto devQR = std::max({devQ, devR});
              const auto devPW = std::max({devP, devW});
              if (devQR > max_devQR)
                max_devQR = devQR;
              if (devPW > max_devPW)
                max_devPW = devPW;
            }
          }
        }
      }
    }
    pass &= qip::check_value(&obuff, "QkTable: Q,R", max_devQR, 0.0f, 5.0e-6f);
    pass &= qip::check_value(&obuff, "QkTable: P,W", max_devPW, 0.0f, 5.0e-6f);
  }

  return pass;
}

} // namespace UnitTest
