#include "DiagramRPA.hpp"
#include "Angular/Angular_369j.hpp"
#include "Angular/Angular_tables.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/YkTable.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/safeProfiler.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

namespace MBPT {

//******************************************************************************
DiagramRPA::DiagramRPA(const DiracOperator::TensorOperator *const h,
                       const std::vector<DiracSpinor> &basis,
                       const std::vector<DiracSpinor> &core)
    : m_k(h->rank()), m_pi(h->parity()), m_imag(h->imaginaryQ()) {

  // Set up basis:
  for (const auto &Fi : basis) {
    const bool inCore = std::find(cbegin(core), cend(core), Fi) != cend(core);
    if (inCore) {
      holes.push_back(Fi);
    } else {
      excited.push_back(Fi);
    }
  }

  // Calc t0 (and setup t) [RPA MEs for hole-excited]
  setup_ts(h);
  fill_W_matrix(h);
}

//******************************************************************************
DiagramRPA::DiagramRPA(const DiracOperator::TensorOperator *const h,
                       const DiagramRPA *const drpa)
    : m_k(h->rank()), m_pi(h->parity()), m_imag(h->imaginaryQ()) {
  //
  if (m_k != drpa->m_k || m_pi != drpa->m_pi) {
    std::cerr << "\nFAIL21 in DiagramRPA: Cannot use 'eat' constructor for "
                 "different rank/parity operators!\n";
    std::abort();
  }

  // Set up basis:
  holes = drpa->holes;
  excited = drpa->excited;

  setup_ts(h);

  // "eat" W matrices from other rpa
  Wanmb = drpa->Wanmb;
  Wabmn = drpa->Wabmn;
  Wmnab = drpa->Wmnab;
  Wmban = drpa->Wmban;
}

//******************************************************************************
void DiagramRPA::fill_W_matrix(const DiracOperator::TensorOperator *const h) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  if (holes.empty() || excited.empty()) {
    std::cout << "\nWARNING 64 in DiagramRPA: no basis! RPA will be zero\n";
    return;
  }
  // IO::ChronoTimer sw("fill_W_matrix");

  const auto maxtj_c =
      std::max_element(holes.cbegin(), holes.cend(), DiracSpinor::comp_j)
          ->twoj();
  const auto maxtj_e =
      std::max_element(excited.cbegin(), excited.cend(), DiracSpinor::comp_j)
          ->twoj();
  const auto maxtj = std::max(maxtj_c, maxtj_e);

  const Coulomb::YkTable Yhe(holes.front().p_rgrid, &holes, &excited);
  const Coulomb::YkTable Yee(holes.front().p_rgrid, &excited);
  const Coulomb::YkTable Yhh(holes.front().p_rgrid, &holes);
  const auto &Ck = maxtj_e > maxtj_c ? Yee.Ck() : Yhh.Ck();
  const Angular::SixJ sj(maxtj, maxtj);

  // RPA: store W Coulomb integrals (used only for Core RPA its)
  std::cout << "Filling RPA Diagram matrix .. " << std::flush;
  Wanmb.resize(holes.size());
  Wabmn.resize(holes.size());
#pragma omp parallel for
  for (std::size_t i = 0; i < holes.size(); i++) {
    const auto &Fa = holes[i];
    auto &Wa_nmb = Wanmb[i];
    auto &Wa_bmn = Wabmn[i];
    Wa_nmb.reserve(excited.size());
    Wa_bmn.reserve(excited.size());
    for (const auto &Fn : excited) {
      auto &Wan_mb = Wa_nmb.emplace_back();
      auto &Wab_mn = Wa_bmn.emplace_back();
      Wan_mb.reserve(excited.size());
      Wab_mn.reserve(excited.size());
      for (const auto &Fm : excited) {
        auto &Wanm_b = Wan_mb.emplace_back();
        auto &Wabm_n = Wab_mn.emplace_back();
        Wanm_b.reserve(holes.size());
        Wabm_n.reserve(holes.size());
        for (const auto &Fb : holes) {
          if (h->isZero(Fb.k, Fn.k)) {
            Wanm_b.emplace_back(0.0);
            Wabm_n.emplace_back(0.0);
            continue;
          }
          const auto yknb = Yhe.ptr_yk_ab(m_k, Fb, Fn);
          const auto &ybm = Yhe.get_y_ab(Fb, Fm);
          const auto &ynm = Yee.get_y_ab(Fn, Fm);
          const auto xQ =
              yknb ? Coulomb::Qk_abcd(Fa, Fn, Fm, Fb, m_k, *yknb, Ck) : 0.0;
          const auto xP = Coulomb::Pk_abcd(Fa, Fn, Fm, Fb, m_k, ynm, Ck, sj);
          const auto yQ =
              yknb ? Coulomb::Qk_abcd(Fa, Fb, Fm, Fn, m_k, *yknb, Ck) : 0.0;
          const auto yP = Coulomb::Pk_abcd(Fa, Fb, Fm, Fn, m_k, ybm, Ck, sj);
          Wanm_b.push_back(xQ + xP);
          Wabm_n.push_back(yQ + yP);
        }
      }
    }
  }
  Wmnab.resize(excited.size());
  Wmban.resize(excited.size());
  std::cout << "." << std::flush;
#pragma omp parallel for
  for (std::size_t i = 0; i < excited.size(); i++) {
    const auto &Fm = excited[i];
    auto &Wa_nmb = Wmnab[i];
    auto &Wa_bmn = Wmban[i];
    Wa_nmb.reserve(excited.size());
    Wa_bmn.reserve(excited.size());
    for (const auto &Fn : excited) {
      auto &Wan_mb = Wa_nmb.emplace_back();
      auto &Wab_mn = Wa_bmn.emplace_back();
      Wan_mb.reserve(holes.size());
      Wab_mn.reserve(holes.size());
      for (const auto &Fa : holes) {
        auto &Wanm_b = Wan_mb.emplace_back();
        auto &Wabm_n = Wab_mn.emplace_back();
        Wanm_b.reserve(holes.size());
        Wabm_n.reserve(holes.size());
        for (const auto &Fb : holes) {
          if (h->isZero(Fb.k, Fn.k)) {
            Wanm_b.emplace_back(0.0);
            Wabm_n.emplace_back(0.0);
            continue;
          }
          const auto yknb = Yhe.ptr_yk_ab(m_k, Fb, Fn);
          const auto &yna = Yhe.get_y_ab(Fa, Fn);
          const auto &yba = Yhh.get_y_ab(Fb, Fa);
          const auto xQ =
              yknb ? Coulomb::Qk_abcd(Fm, Fn, Fa, Fb, m_k, *yknb, Ck) : 0.0;
          const auto xP = Coulomb::Pk_abcd(Fm, Fn, Fa, Fb, m_k, yna, Ck, sj);
          const auto yQ =
              yknb ? Coulomb::Qk_abcd(Fm, Fb, Fa, Fn, m_k, *yknb, Ck) : 0.0;
          const auto yP = Coulomb::Pk_abcd(Fm, Fb, Fa, Fn, m_k, yba, Ck, sj);
          Wanm_b.push_back(xQ + xP);
          Wabm_n.push_back(yQ + yP);
        }
      }
    }
  }
  std::cout << " done.\n" << std::flush;
}

//******************************************************************************
void DiagramRPA::setup_ts(const DiracOperator::TensorOperator *const h) {
  if (holes.empty() || excited.empty())
    return;

  // Calc t0 (and setup t)
  for (const auto &Fa : holes) {
    std::vector<double> t0a_m;
    for (const auto &Fm : excited) {
      t0a_m.push_back(h->reducedME(Fa, Fm));
    }
    t0am.push_back(t0a_m);
  }
  for (const auto &Fm : excited) {
    std::vector<double> t0m_a;
    for (const auto &Fa : holes) {
      t0m_a.push_back(h->reducedME(Fm, Fa));
    }
    t0ma.push_back(t0m_a);
  }
  clear_tam();
}
//******************************************************************************
void DiagramRPA::clear_tam() {
  tam = t0am;
  tma = t0ma;
}

//******************************************************************************
double DiagramRPA::dV(const DiracSpinor &Fw, const DiracSpinor &Fv,
                      const bool first_order) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  if (holes.empty() || excited.empty())
    return 0.0;

  const auto orderOK = Fv.en <= Fw.en;
  const auto &Fi = orderOK ? Fv : Fw;
  const auto &Ff = orderOK ? Fw : Fv;
  const auto f = (1.0 / (2 * m_k + 1));

  std::vector<double> sum_a(holes.size());
#pragma omp parallel for
  for (std::size_t ia = 0; ia < holes.size(); ia++) {
    const auto &Fa = holes[ia];
    const auto s1 = ((Fa.twoj() - Ff.twoj() + 2 * m_k) % 4 == 0) ? 1 : -1;
    for (std::size_t im = 0; im < excited.size(); im++) {
      const auto &Fm = excited[im];
      if (t0am[ia][im] == 0.0)
        continue;
      const auto s2 = ((Fa.twoj() - Fm.twoj()) % 4 == 0) ? 1 : -1;
      const auto Wwmva = Coulomb::Wk_abcd(Ff, Fm, Fi, Fa, m_k);
      const auto Wwavm = Coulomb::Wk_abcd(Ff, Fa, Fi, Fm, m_k);
      const auto ttam = first_order ? t0am[ia][im] : tam[ia][im];
      const auto ttma = first_order ? t0ma[im][ia] : tma[im][ia];
      const auto A = ttam * Wwmva / (Fa.en - Fm.en - m_omega);
      const auto B = Wwavm * ttma / (Fa.en - Fm.en + m_omega);
      sum_a[ia] += s1 * (A + s2 * B);
    }
  }
  return f * std::accumulate(begin(sum_a), end(sum_a), 0.0);
}

//******************************************************************************
void DiagramRPA::rpa_core(const double omega, const bool print) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  m_omega = std::abs(omega);
  // IO::ChronoTimer sw("rpa_core");

  if (holes.empty() || excited.empty())
    return;

  if (print) {
    printf("RPA(D) (w=%.3f): .. \r", m_omega);
    std::cout << std::flush;
  }
  int it = 1;
  auto eps = 0.0;
  const auto f = (1.0 / (2 * m_k + 1));
  for (; it <= max_its; it++) {
    std::vector<double> eps_m(excited.size()); //"thread-safe" eps..?
    // XXX "Small" race condition in here???
#pragma omp parallel for
    for (std::size_t im = 0; im < excited.size(); im++) {
      const auto &Fm = excited[im];
      double eps_worst_a = 0.0;
      for (std::size_t ia = 0; ia < holes.size(); ia++) {
        const auto &Fa = holes[ia];

        double sum_am = 0;
        double sum_ma = 0;

        // Can replace this with dV?? NO. 1) it calcs W (since valence)
        // 2) not thread safe
        for (std::size_t ib = 0; ib < holes.size(); ib++) {
          const auto &Fb = holes[ib];
          const auto s1 = ((Fb.twoj() - Fa.twoj() + 2 * m_k) % 4 == 0) ? 1 : -1;
          const auto s3 = ((Fb.twoj() - Fm.twoj() + 2 * m_k) % 4 == 0) ? 1 : -1;
          for (std::size_t in = 0; in < excited.size(); in++) {
            const auto &Fn = excited[in];

            const auto tbn = tam[ib][in];
            if (tbn == 0.0)
              continue;
            const auto tnb = tma[in][ib];
            const auto tdem = tbn / (Fb.en - Fn.en - m_omega);
            const auto s2 = ((Fb.twoj() - Fn.twoj()) % 4 == 0) ? 1 : -1;
            const auto stdep = s2 * tnb / (Fb.en - Fn.en + m_omega);
            const auto A = tdem * Wanmb[ia][in][im][ib];
            const auto B = stdep * Wabmn[ia][in][im][ib];
            const auto C = tdem * Wmnab[im][in][ia][ib];
            const auto D = stdep * Wmban[im][in][ia][ib];
            sum_am += s1 * (A + B);
            sum_ma += s3 * (C + D);
          }
        }

        const auto prev = tam[ia][im];
        // 0.5 factor is for damping. f*sum is dV
        tam[ia][im] = 0.5 * (tam[ia][im] + t0am[ia][im] + f * sum_am);
        tma[im][ia] = 0.5 * (tma[im][ia] + t0ma[im][ia] + f * sum_ma);
        const auto delta = std::abs((tam[ia][im] - prev) / tam[ia][im]);
        if (delta > eps_worst_a)
          eps_worst_a = delta;
      } // a (holes)
      eps_m[im] = eps_worst_a;
    } // m (excited)
    // XXX "small" race condition somewhere regarding eps??
    // The itteraion it converges on always seems to be the same..
    // but the value for eps printed changes slightly each run???
    eps = *std::max_element(cbegin(eps_m), cend(eps_m));
    if (eps < eps_targ)
      break;
    if (print && it % 15 == 0) {
      printf("RPA(D) (w=%.3f): %2i %.1e \r", m_omega, it, eps);
      std::cout << std::flush;
    }
  } // its
  if (print) {
    printf("RPA(D) (w=%.3f): %2i %.1e\n", m_omega, it, eps);
  }
  m_core_eps = eps;
}

} // namespace MBPT
