#include "ExternalField.hpp"
#include "Angular/Angular.hpp"
#include "Angular/Wigner_369j.hpp"
#include "Dirac/DiracOperator.hpp"
#include "Dirac/DiracSpinor.hpp"
#include "Dirac/Wavefunction.hpp"
#include "HF/HartreeFockClass.hpp"
#include "IO/ChronoTimer.hpp"
#include <algorithm>
#include <vector>
//
#include "Dirac/BSplineBasis.hpp" // XXX add to make

//******************************************************************************
ExternalField::ExternalField(const DiracOperator *const h,
                             const std::vector<DiracSpinor> &core,
                             const std::vector<double> &vl, const double alpha,
                             const double omega)
    : m_h(h), p_core(&core), m_vl(vl), m_alpha(alpha), m_omega(omega), //
      static_fieldQ(std::abs(omega) > 0 ? false : false),              //
      m_rank(h->rank()), m_pi(h->parity()), m_imag(h->imaginaryQ())    //
// w>0 typically. Allowed to be -ve for tests?
{
  ChronoTimer timer("ExternalField (construct)");
  m_X.resize(core.size());
  for (auto ic = 0u; ic < core.size(); ic++) {
    const auto &phic = core[ic];
    const auto pi_ch = phic.parity() * m_pi;
    const auto tj_c = phic.twoj();
    const auto tjmin_tmp = tj_c - 2 * m_rank;
    const auto tjmin = tjmin_tmp < 1 ? 1 : tjmin_tmp;
    const auto tjmax = tj_c + 2 * m_rank;
    for (int tj = tjmin; tj <= tjmax; tj += 2) {
      const auto l_minus = (tj - 1) / 2;
      const auto pi_chla = Wigner::parity_l(l_minus) * pi_ch;
      const auto l = (pi_chla == 1) ? l_minus : l_minus + 1;
      // XXX Add option to restrict l = lc ?
      const auto kappa = Wigner::kappa_twojl(tj, l);
      m_X[ic].emplace_back(0, kappa, *(phic.p_rgrid));
    }
  }
  if (!static_fieldQ)
    m_Y = m_X;

  // for (const auto &phic : core)

  auto max_twoj = [](const DiracSpinor &a, const DiracSpinor &b) {
    return a.twoj() < b.twoj();
  };

  auto max_tj_core = std::max_element(core.cbegin(), core.cend(), max_twoj);
  auto max_tj_dPsi = max_tj_core->twoj() + 2 * m_rank;
  m_6j.fill(m_rank, max_tj_dPsi);

  // for (const auto &phic : core) {
  //   auto &dPsi = get_dPsis(phic, dPsiType::X);
  //   std::cout << phic.symbol() << " : ";
  //   for (const auto &dX : dPsi) {
  //     std::cout << dX.symbol() << " ";
  //   }
  //   std::cout << "\n";
  // }
}
//******************************************************************************
std::size_t ExternalField::core_index(const DiracSpinor &phic) {
  // XXX Note: no bounds checking here!
  return static_cast<std::size_t>(
      std::find(p_core->cbegin(), p_core->cend(), phic) - p_core->cbegin());
}

//******************************************************************************
std::vector<DiracSpinor> &ExternalField::get_dPsis(const DiracSpinor &phic,
                                                   dPsiType XorY) {
  auto index = core_index(phic);
  return XorY == dPsiType::X ? m_X[index] : m_Y[index];
}

//******************************************************************************
const DiracSpinor &ExternalField::get_dPsi_x(const DiracSpinor &phic,
                                             dPsiType XorY, const int kappa_x) {
  const auto &dPsis = get_dPsis(phic, XorY);
  auto match_kappa_x = [=](const DiracSpinor &phi) { return phi.k == kappa_x; };
  return *std::find_if(dPsis.cbegin(), dPsis.cend(), match_kappa_x);
}

//******************************************************************************
void ExternalField::solve_TDHFcore() {
  // XXX THIS should have omega inside! ?
  //
  ChronoTimer timer("solve_TDHFcore");

#pragma omp parallel for
  for (auto ic = 0ul; ic < p_core->size(); ic++) {
    const auto &phic = (*p_core)[ic];
    auto &dPsis_X = m_X[ic];
    for (auto &Xx : dPsis_X) {
      const auto hPsic = m_h->reduced_rhs(Xx.k, phic); // always same!?
      HartreeFock::solveMixedState(Xx, phic, m_omega, m_vl, m_alpha, *p_core,
                                   hPsic);
    }
    if (!static_fieldQ) {
      auto &dPsis_Y = m_Y[ic];
      for (auto &Yx : dPsis_Y) {
        // XXX rhs? or lhs? lhs. rhs, but with conjugate!?
        auto hPsic = m_h->reduced_rhs(Yx.k, phic); // always same!?
        if (m_h->imaginaryQ())
          hPsic *= -1.0;
        HartreeFock::solveMixedState(Yx, phic, -m_omega, m_vl, m_alpha, *p_core,
                                     hPsic);
      }
    }
  }

  // for (const auto &phi_b : *p_core) {
  //   auto &X_betas = get_dPsis(phi_b, dPsiType::X);
  //   std::cout << "r " << phi_b.symbol() << " ";
  //   for (const auto &phi_beta : X_betas) {
  //     std::cout << phi_beta.symbol() << " ";
  //   }
  //   std::cout << "\n";
  //   for (int ir = 0; ir < phi_b.p_rgrid->num_points; ir++) {
  //     std::cout << phi_b.p_rgrid->r[ir] << " " << phi_b.f[ir] << " ";
  //     for (const auto &phi_beta : X_betas) {
  //       std::cout << phi_beta.f[ir] << " ";
  //     }
  //     std::cout << "\n";
  //   }
  //   std::cin.get();
  // }

  // std::cout << "x\n" << std::flush;
  // for (const auto &phi_b : *p_core) {
  //   auto &X_betas = get_dPsis(phi_b, dPsiType::X);
  //   std::cout << phi_b.symbol() << ":\n";
  //   for (const auto &phi_beta : X_betas) {
  //     auto dv = dV_ab(phi_beta, phi_b);
  //     std::cout << "  " << phi_beta.symbol() << " " << dv << "\n";
  //   }
  // }
  // std::cout << "x\n" << std::flush;

  // for (const auto &phic : *p_core) {
  //   auto &dPsi = get_dPsis(phic, dPsiType::X);
  //   for (const auto &Xx : dPsi) {
  //     std::cout << phic.symbol() << " " << Xx.symbol() << " " << Xx * Xx << "
  //     "
  //               << Xx * phic << "\n";
  //   }
  // }
  //
}

//******************************************************************************
void ExternalField::solve_TDHFcore_matrix(const Wavefunction &wf) {
  static bool first = true;
  //

  ChronoTimer timer("solve_TDHFcore_matrix");

  const std::size_t nspl = 30;
  const std::size_t kspl = 7;
  const double rmin = 1.0e-4; //?
  const double rmax = 50.0;   //?

  // XXX Do not need to construct basis each time! Just do once (each kappa)

  // A := H - (e-w)S
  // b := -h -dV +deS_c

  auto a_damp = first ? 0.0 : 0.5;

  auto max_ki = 0;
  for (const auto &Px : m_X) {
    for (const auto &x : Px) {
      auto ki = Angular::indexFromKappa(x.k);
      if (ki > max_ki)
        max_ki = ki;
    }
  }
  std::vector<std::vector<DiracSpinor>> basis_kappa;
  basis_kappa.reserve(max_ki + 1);
  for (int ki = 0; ki <= max_ki; ki++) {
    auto k = Angular::kappaFromIndex(ki);
    basis_kappa.push_back(SplineBasis::form_spline_basis(
        k, nspl, kspl, rmin, rmax, *((*p_core)[0].p_rgrid), m_alpha));
  }

  // for (const auto &phic : *p_core) {
  auto tmp_X = m_X;
  // #pragma omp parallel for
  for (std::size_t ic = 0; ic < p_core->size(); ic++) {
    const auto &phic = (*p_core)[ic];
    // auto &dPsi = get_dPsis(phic, dPsiType::X);
    auto &dPsi = tmp_X[ic];
    auto de = m_h->reducedME(phic, phic); // no dV ?
    for (auto &Xx : dPsi) {
      // auto basis = SplineBasis::form_spline_basis(Xx.k, nspl, kspl, rmin,
      // rmax,
      //                                             *(phic.p_rgrid), m_alpha);
      auto ki = Angular::indexFromKappa(Xx.k);
      const auto &basis = basis_kappa[ki];

      LinAlg::Vector bi((int)basis.size());
      for (int i = 0; i < (int)basis.size(); ++i) {
        const auto &xi = basis[i];
        // fill LHS vector, b
        auto hi = m_h->reducedME(xi, phic);
        auto dV = first ? 0.0 : dV_ab(xi, phic);
        auto deS = (xi.k == phic.k) ? de * (xi * phic) : 0.0; // reduced? OK?
        bi[i] = -hi - dV + deS;
      }
      auto [Hij, Sij] = SplineBasis::fill_Hamiltonian_matrix(basis, wf);
      auto Aij = Hij - (phic.en - m_omega) * Sij;

      // auto tm = std::min(phic.twoj(), phic.twoj());
      // auto rme_factor = m_h->rme3js(Xx.twoj(), phic.twoj(), tm, 0); //??
      // // Aij *= (1.0 / rme_factor);
      // bi *= rme_factor;
      // std::cout << rme_factor << "\n";
      // auto f = ((Xx.twoj() + 1.0) * (phic.twoj() + 1.0));
      // std::cout << f * rme_factor << "\n";
      // bi *= rme_factor;
      // if (rme_factor == 0)
      //   continue;

      auto c = LinAlg::solve_Axeqb(Aij, bi);
      // c *= rme_factor;

      Xx.scale(a_damp); // or 0.5, for DAMP!
      for (int i = 0; i < (int)basis.size(); ++i) {
        Xx += (1.0 - a_damp) * c[i] * basis[i];
      }
    }
  }

  if (!static_fieldQ) {
    // for (const auto &phic : *p_core) {
    //   auto &dPsi = get_dPsis(phic, dPsiType::Y);
    auto tmp_Y = m_Y;
    // #pragma omp parallel for
    for (std::size_t ic = 0; ic < p_core->size(); ic++) {
      const auto &phic = (*p_core)[ic];
      // auto &dPsi = get_dPsis(phic, dPsiType::X);
      auto &dPsi = tmp_Y[ic];
      auto de = m_h->reducedME(phic, phic); // no dV ?
      for (auto &Yy : dPsi) {
        // auto basis = SplineBasis::form_spline_basis(Xx.k, nspl, kspl, rmin,
        // rmax, *(phic.p_rgrid), m_alpha);
        auto ki = Angular::indexFromKappa(Yy.k);
        const auto &basis = basis_kappa[ki];
        LinAlg::Vector bi((int)basis.size());
        for (int i = 0; i < (int)basis.size(); ++i) {
          const auto &yi = basis[i];
          // fill LHS vector, b
          // auto hi = m_h->reducedME(phic, yi); // ?? XXX
          auto hi = m_h->reducedME(yi, phic); // ?? XXX
          auto dV = first ? 0.0 : dV_ab_Y(yi, phic);
          // auto dV = dV_ab_Y(phic, yi);
          auto deS = (yi.k == phic.k) ? de * (yi * phic) : 0.0; // reduced? OK?
          bi[i] = -hi - dV + deS;
        }
        auto [Hij, Sij] = SplineBasis::fill_Hamiltonian_matrix(basis, wf);
        auto Aij = Hij - (phic.en + m_omega) * Sij;

        // auto tm = std::min(phic.twoj(), phic.twoj());
        // auto rme_factor = m_h->rme3js(phic.twoj(), Yy.twoj(), tm, 0); //??
        // // // Aij *= (1.0 / rme_factor);
        // bi *= rme_factor;
        // std::cout << rme_factor << "\n";
        // auto f = ((Yy.twoj() + 1.0) * (phic.twoj() + 1.0));
        // std::cout << f * rme_factor << "\n";
        // bi *= rme_factor;
        // if (rme_factor == 0)
        //   continue;

        auto c = LinAlg::solve_Axeqb(Aij, bi);
        // c *= rme_factor;

        Yy.scale(a_damp); // or 0.5, for DAMP!
        for (int i = 0; i < (int)basis.size(); ++i) {
          Yy += (1.0 - a_damp) * c[i] * basis[i]; // XXX -ve
        }
      }
    }
    m_Y = tmp_Y;
  }
  // std::cin.get();
  m_X = tmp_X;
  first = false;
}

//******************************************************************************
// does it matter if a or b is in the core?
double ExternalField::dV_ab(const DiracSpinor &phi_alpha,
                            const DiracSpinor &phi_a) {

  // XXX Take advantage of internal y_bd - symmetry save time? No bother, this
  // is very fast! [solve_TDHFcore is ~2000x slower]
  // XXX Am I accounting for 3j symbol twice?? [using projections of dPsi??]

  // ChronoTimer timer("dV_ab");

  const auto k = m_h->rank();
  const auto m1tkp1 = Wigner::evenQ(k + 1) ? 1 : -1;
  const auto tkp1 = double(2 * k + 1);

  const auto tjalpha = phi_alpha.twoj();
  const auto tja = phi_a.twoj();

  double rme_sum_dirX = 0.0;
  double rme_sum_excX = 0.0;
#pragma omp parallel for
  for (auto ic = 0u; ic < p_core->size(); ic++) {
    const auto &phi_b = (*p_core)[ic];
    const auto tjb = phi_b.twoj();
    const auto &X_betas = get_dPsis(phi_b, dPsiType::X);
    double rme_sum_dirX_c = 0.0;
    double rme_sum_excX_c = 0.0;
    for (const auto &phi_beta : X_betas) {
      const auto tjbeta = phi_beta.twoj();
      const auto m1jaljbe = Wigner::evenQ_2(tjalpha + tjbeta) ? 1 : -1;
      const auto Qkabcd =
          Coulomb::Qk_abcd_any(phi_alpha, phi_b, phi_a, phi_beta, k);
      rme_sum_dirX_c += m1jaljbe * Qkabcd;
      // exchange part:
      const auto bmc = std::abs(tjalpha - tjbeta);
      const auto bpc = tjalpha + tjbeta;
      const auto amd = std::abs(tja - tjb);
      const auto apd = tja + tjb;
      const auto l_min = std::max(amd, bmc) / 2;
      const auto l_max = std::min(apd, bpc) / 2;
      for (int l = l_min; l <= l_max; ++l) {
        const auto m1jaljb = Wigner::evenQ_2(tjalpha + tjb) ? 1 : -1;
        const auto sixj = m_6j.get_6j(tja, tjalpha, tjbeta, tjb, k, l);
        if (sixj == 0)
          continue;
        const auto Qlabcd =
            Coulomb::Qk_abcd_any(phi_alpha, phi_a, phi_beta, phi_b, k);
        rme_sum_excX_c += m1jaljb * sixj * Qlabcd;
      }
    }
#pragma omp critical(sum_X_core)
    {
      rme_sum_dirX += rme_sum_dirX_c;
      rme_sum_excX += rme_sum_excX_c;
    }
  }

  double rme_sum_dirY = 0.0;
  double rme_sum_excY = 0.0;
  if (!static_fieldQ) {
#pragma omp parallel for
    for (auto ic = 0u; ic < p_core->size(); ic++) {
      const auto &phi_b = (*p_core)[ic];
      const auto tjb = phi_b.twoj();
      const auto &Y_betas = get_dPsis(phi_b, dPsiType::Y);
      double rme_sum_dirY_c = 0.0;
      double rme_sum_excY_c = 0.0;
      // std::cout << phi_b.symbol() << "\n";
      for (const auto &phi_beta : Y_betas) {
        const auto tjbeta = phi_beta.twoj();
        const auto m1jaljbe = Wigner::evenQ_2(tjalpha + tjbeta) ? 1 : -1;
        const auto Qkabcd =
            Coulomb::Qk_abcd_any(phi_alpha, phi_b, phi_a, phi_beta, k);
        // std::cout << " " << phi_beta.symbol() << " " << Qkabcd << "\n";
        rme_sum_dirY_c += m1jaljbe * Qkabcd;
        // exchange part:
        const auto amd = std::abs(tja - tjbeta);
        const auto apd = tja + tjbeta;
        const auto bmc = std::abs(tjalpha - tjb);
        const auto bpc = tjalpha + tjb;
        const auto l_min = std::max(amd, bmc) / 2;
        const auto l_max = std::min(apd, bpc) / 2;
        for (int l = l_min; l <= l_max; ++l) {
          const auto m1jaljb = Wigner::evenQ_2(tjbeta + tjb) ? 1 : -1;
          const auto sixj = m_6j.get_6j(tja, tjalpha, tjb, tjbeta, k, l);
          if (sixj == 0)
            continue;
          const auto Qlabcd =
              Coulomb::Qk_abcd_any(phi_alpha, phi_a, phi_b, phi_beta, k);
          rme_sum_excY_c += m1jaljb * sixj * Qlabcd;
        }
      }
#pragma omp critical(sum_Y_core)
      {
        rme_sum_dirY += rme_sum_dirY_c;
        rme_sum_excY += rme_sum_excY_c;
      }
    }
  } else {
    rme_sum_dirY = rme_sum_dirX;
    rme_sum_excY = rme_sum_excX;
  }
  auto rme_sum_dir = rme_sum_dirX + rme_sum_dirY; // Works w/ - ??
  auto rme_sum_exc = rme_sum_excX + rme_sum_excY;
  rme_sum_dir *= m1tkp1 / tkp1;
  rme_sum_exc *= m1tkp1;

  // std::cout << rme_sum_dirX << " " << rme_sum_dirY << " / " << rme_sum_excX
  //           << " " << rme_sum_excY << "\n";

  return rme_sum_dir + rme_sum_exc;
}

//******************************************************************************
// does it matter if a or b is in the core?
double ExternalField::dV_ab_Y(const DiracSpinor &phi_alpha,
                              const DiracSpinor &phi_a) {

  // XXX Take advantage of internal y_bd - symmetry save time? No bother, this
  // is very fast! [solve_TDHFcore is ~2000x slower]
  // XXX Am I accounting for 3j symbol twice?? [using projections of dPsi??]

  // ChronoTimer timer("dV_ab");

  const auto k = m_h->rank();
  const auto m1tkp1 = Wigner::evenQ(k + 1) ? 1 : -1;
  const auto tkp1 = double(2 * k + 1);

  const auto tjalpha = phi_alpha.twoj();
  const auto tja = phi_a.twoj();

  double rme_sum_dirX = 0.0;
  double rme_sum_excX = 0.0;
#pragma omp parallel for
  for (auto ic = 0u; ic < p_core->size(); ic++) {
    const auto &phi_b = (*p_core)[ic];
    const auto tjb = phi_b.twoj();
    const auto &X_betas = get_dPsis(phi_b, dPsiType::Y);
    double rme_sum_dirX_c = 0.0;
    double rme_sum_excX_c = 0.0;
    for (const auto &phi_beta : X_betas) {
      const auto tjbeta = phi_beta.twoj();
      const auto m1jaljbe = Wigner::evenQ_2(tjalpha + tjbeta) ? 1 : -1;
      const auto Qkabcd =
          Coulomb::Qk_abcd_any(phi_alpha, phi_b, phi_a, phi_beta, k);
      rme_sum_dirX_c += m1jaljbe * Qkabcd;
      // exchange part:
      const auto bmc = std::abs(tjalpha - tjbeta);
      const auto bpc = tjalpha + tjbeta;
      const auto amd = std::abs(tja - tjb);
      const auto apd = tja + tjb;
      const auto l_min = std::max(amd, bmc) / 2;
      const auto l_max = std::min(apd, bpc) / 2;
      for (int l = l_min; l <= l_max; ++l) {
        const auto m1jaljb = Wigner::evenQ_2(tjalpha + tjb) ? 1 : -1;
        const auto sixj = m_6j.get_6j(tja, tjalpha, tjbeta, tjb, k, l);
        if (sixj == 0)
          continue;
        const auto Qlabcd =
            Coulomb::Qk_abcd_any(phi_alpha, phi_a, phi_beta, phi_b, k);
        rme_sum_excX_c += m1jaljb * sixj * Qlabcd;
      }
    }
#pragma omp critical(sum_X_core)
    {
      rme_sum_dirX += rme_sum_dirX_c;
      rme_sum_excX += rme_sum_excX_c;
    }
  }

  double rme_sum_dirY = 0.0;
  double rme_sum_excY = 0.0;
  if (!static_fieldQ) {
#pragma omp parallel for
    for (auto ic = 0u; ic < p_core->size(); ic++) {
      const auto &phi_b = (*p_core)[ic];
      const auto tjb = phi_b.twoj();
      const auto &Y_betas = get_dPsis(phi_b, dPsiType::X);
      double rme_sum_dirY_c = 0.0;
      double rme_sum_excY_c = 0.0;
      // std::cout << phi_b.symbol() << "\n";
      for (const auto &phi_beta : Y_betas) {
        const auto tjbeta = phi_beta.twoj();
        const auto m1jaljbe = Wigner::evenQ_2(tjalpha + tjbeta) ? 1 : -1;
        const auto Qkabcd =
            Coulomb::Qk_abcd_any(phi_alpha, phi_b, phi_a, phi_beta, k);
        // std::cout << " " << phi_beta.symbol() << " " << Qkabcd << "\n";
        rme_sum_dirY_c += m1jaljbe * Qkabcd;
        // exchange part:
        const auto amd = std::abs(tja - tjbeta);
        const auto apd = tja + tjbeta;
        const auto bmc = std::abs(tjalpha - tjb);
        const auto bpc = tjalpha + tjb;
        const auto l_min = std::max(amd, bmc) / 2;
        const auto l_max = std::min(apd, bpc) / 2;
        for (int l = l_min; l <= l_max; ++l) {
          const auto m1jaljb = Wigner::evenQ_2(tjbeta + tjb) ? 1 : -1;
          const auto sixj = m_6j.get_6j(tja, tjalpha, tjb, tjbeta, k, l);
          if (sixj == 0)
            continue;
          const auto Qlabcd =
              Coulomb::Qk_abcd_any(phi_alpha, phi_a, phi_b, phi_beta, k);
          rme_sum_excY_c += m1jaljb * sixj * Qlabcd;
        }
      }
#pragma omp critical(sum_Y_core)
      {
        rme_sum_dirY += rme_sum_dirY_c;
        rme_sum_excY += rme_sum_excY_c;
      }
    }
  } else {
    rme_sum_dirY = rme_sum_dirX;
    rme_sum_excY = rme_sum_excX;
  }
  auto rme_sum_dir = rme_sum_dirX + rme_sum_dirY; // Works w/ - ??
  auto rme_sum_exc = rme_sum_excX + rme_sum_excY;
  rme_sum_dir *= m1tkp1 / tkp1;
  rme_sum_exc *= m1tkp1;

  // std::cout << rme_sum_dirX << " " << rme_sum_dirY << " / " << rme_sum_excX
  //           << " " << rme_sum_excY << "\n";

  return rme_sum_dir + rme_sum_exc;
}

//******************************************************************************
DiracSpinor ExternalField::dV_ab_rhs(const DiracSpinor &phia,
                                     const DiracSpinor &phib) {
  return (phia * phib) * phib; // XX just stand-in
}
