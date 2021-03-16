#pragma once
#include "Coulomb.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "YkTable.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <map>

//! Coulomb
namespace Coulomb {

/*!
@brief
@details



*/

enum class Symmetry { Rk, Wk, none };
// XXX Better to have a base class (with 'none')
// And then specialised derived classes Wk and Rk ?

// Rk symmetry:
// {abcd} = cbad = adcb = cdab = badc = bcda = dabc = dcba
// Normal ordering: i = min(a,b,c,d) is first index
// Two options (second and fourth may be swapped): choose 2nd to be smallest
// Wk symmetry:
// {abcd} = badc = cdab = dcba

template <typename Real = float, Symmetry symmetry = Symmetry::none> //
class QkTable {
  using Index = uint16_t;
  using IndexSet = std::array<Index, 4>;

private:
  std::map<IndexSet, std::vector<Real>> m_data{};

public:
  QkTable() {}

  //! Gives arrow access to all map functions
  auto operator-> () { return &m_data; }

  //! Fill all non-zero Qk integrals (nb: often don't need all)
  void fill(const YkTable &yk) {
    assert(yk.a_is_b());
    const auto &orbs = yk.get_a();
    // nb: can do this more efficiently [without needing .find()] if we know
    // which symmetry a <{a,c,b}, and b<d
    for (const auto &a : orbs) {
      for (const auto &b : orbs) {
        for (const auto &c : orbs) {
          for (const auto &d : orbs) {
            addQ(a, b, c, d, yk);
          }
        }
      }
    }
  }

  //! Update all non-zero Qk integrals (nb: often don't need all)
  void update(const YkTable &yk) {
    assert(yk.a_is_b());
    const auto &orbs = yk.get_a();
    // nb: can do this more efficiently [without needing .find()] if we know
    // which symmetry a <{a,c,b}, and b<d
    for (const auto &a : orbs) {
      for (const auto &b : orbs) {
        for (const auto &c : orbs) {
          for (const auto &d : orbs) {
            updateQ(a, b, c, d, yk);
          }
        }
      }
    }
  }

  //! Add a specific new Q (for each k) to map. If exists, does nothing
  void addQ(const DiracSpinor &a, const DiracSpinor &b, const DiracSpinor &c,
            const DiracSpinor &d, const YkTable &yk) {
    const auto [kmin, kmax] = Coulomb::k_minmax_Q(a, b, c, d);
    if (kmax > kmin)
      return;
    const auto [it, new_element] =
        m_data.insert({NormalOrder(a, b, c, d), std::vector<Real>{}});
    // If already exists, do nothing
    if (!new_element)
      return;
    assert(it->empty());
    const auto num_ks = std::size_t((kmax - kmin) / 2 + 1);
    it->reserve(num_ks);
    for (int k = kmin; k <= kmax; k += 2) {
      const auto yk_bd = yk.ptr_yk_ab(k, b, d);
      assert(yk_bd != nullptr);
      it->push_back(Qk_abcd(a, b, c, d, k, *yk_bd, yk.Ck()));
    }
  }

  //! Updates a specific new Q (for each k) in map. If absent, does nothing!
  void updateQ(const DiracSpinor &a, const DiracSpinor &b, const DiracSpinor &c,
               const DiracSpinor &d, const YkTable &yk) {
    const auto map_it = m_data.find(NormalOrder(a, b, c, d));
    if (map_it == m_data.end())
      return;
    const auto [kmin, kmax] = Coulomb::k_minmax_Q(a, b, c, d);
    auto vec_it = map_it->begin();
    for (int k = kmin; k <= kmax; k += 2) {
      assert(vec_it != map_it->end());
      const auto yk_bd = yk.ptr_yk_ab(k, b, d);
      assert(yk_bd != nullptr);
      *vec_it = Qk_abcd(a, b, c, d, k, *yk_bd, yk.Ck());
      vec_it++;
    }
  }

  //! Return pointer stored Qk (each k). If absent, nullptr
  const std::vector<Real> *Qk(const DiracSpinor &a, const DiracSpinor &b,
                              const DiracSpinor &c,
                              const DiracSpinor &d) const {
    const auto it = m_data.find(NormalOrder(a, b, c, d));
    return (it == m_data.end()) ? nullptr : &*it;
  }

  //! Returns specific Q^k_abcd. If absent, zero
  Real Q(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const {
    const auto qk = Qk(a, b, c, d);
    if (qk == nullptr)
      return Real{0.0};
    const auto [kmin, kmax] = Coulomb::k_minmax_Q(a, b, c, d);
    if ((k < kmin) || (k > kmax) || (k % 2 != kmin % 2))
      Real{0.0};
    const auto index = std::size_t((k - kmin) / 2);
    assert(qk->size() > index);
    return (*qk)[index];
  }

  //! Uses stored Q's to return R. Zero if absent
  Real R(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const {
    const auto tQk = Q(k, a, b, c, d);
    const auto s = Angular::neg1pow(k);
    const auto tCkac = Angular::tildeCk_kk(k, a.k, c.k);
    const auto tCkbd = Angular::tildeCk_kk(k, b.k, d.k);
    return tQk / (s * tCkac * tCkbd);
  }

  //! Uses stored Q's to return P. Zero if absent
  Real P(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const {
    Real Pk_abcd{0.0};
    const auto Ql_abdc = Qk(a, b, d, c); // exchange!
    if (Ql_abdc == nullptr)
      return Pk_abcd;
    const auto [lmin, lmax] = Coulomb::k_minmax_Q(a, b, d, c);
    int l = lmin;
    for (auto &q_abdc : *Ql_abdc) {
      assert(l <= lmax);
      const auto sixj = Coulomb::sixj(a, c, k, b, d, l);
      Pk_abcd += sixj * q_abdc;
      l += 2;
    }
    Pk_abcd *= (2 * k + 1);
    return Pk_abcd;
  }

  //! Uses stored Q's to return W. Zero if absent
  Real W(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const {
    return Q(k, a, b, c, d) + P(k, a, b, c, d);
  }

public:
  IndexSet NormalOrder(const DiracSpinor &a, const DiracSpinor &b,
                       const DiracSpinor &c, const DiracSpinor &d) {
    return NormalOrder(a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index());
  }

  IndexSet NormalOrder(int a, int b, int c, int d) {
    return NormalOrder(Index(a), Index(b), Index(c), Index(d));
  }

  IndexSet NormalOrder(const IndexSet &abcd) {
    const auto [a, b, c, d] = abcd;
    return NormalOrder(a, b, c, d);
  }

  IndexSet NormalOrder(Index a, Index b, Index c, Index d) {
    if constexpr (symmetry == Symmetry::Rk) {
      return NormalOrder_Rk(a, b, c, d);
    } else if constexpr (symmetry == Symmetry::Wk) {
      return NormalOrder_Wk(a, b, c, d);
    } else if constexpr (symmetry == Symmetry::none) {
      return NormalOrder_none(a, b, c, d);
    } else {
      assert(false);
    }
  }

private:
  IndexSet NormalOrder_Rk(Index a, Index b, Index c, Index d) {
    // put smallest first
    const auto min = std::min({a, b, c, d});
    if (min == a) {
      // options are abcd, and adcb
      return (b < d) ? std::array{a, b, c, d} : std::array{a, d, c, b};
    } else if (min == b) {
      // options are badc, and bcda
      return (a < c) ? std::array{b, a, d, c} : std::array{b, c, d, a};
    } else if (min == c) {
      // options are cbad, and cdab
      return (b < d) ? std::array{c, b, a, d} : std::array{c, d, a, b};
    } else if (min == d) {
      // options are dabc, and dcba
      return (a < c) ? std::array{d, a, b, c} : std::array{d, c, b, a};
    }
    assert(false);
    // Check - what happens when 'min' is not unique?
    // Must still be fine
  }

  IndexSet NormalOrder_Wk(Index a, Index b, Index c, Index d) {
    // put smallest first
    const auto min = std::min({a, b, c, d});
    if (min == a) {
      return {a, b, c, d};
    } else if (min == b) {
      return {b, a, d, c};
    } else if (min == c) {
      return {c, d, a, b};
    } else if (min == d) {
      return {d, c, b, a};
    }
    assert(false);
  }

  IndexSet NormalOrder_none(Index a, Index b, Index c, Index d) {
    return {a, b, c, d};
  }
};

} // namespace Coulomb
