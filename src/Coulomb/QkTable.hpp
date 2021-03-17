#pragma once
#include "Coulomb.hpp"
#include "IO/SafeProfiler.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "YkTable.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <map>

namespace Coulomb {

/*!
@brief
@details


*/

//! Symmetry (state index order) for tables
enum class Symmetry { Qk, Wk, none };
// XXX Better to have a base class (with 'none') ?
// And then specialised derived classes Wk and Rk ?

// Qk symmetry:
// {abcd} = cbad = adcb = cdab = badc = bcda = dabc = dcba
// Normal ordering: i = min(a,b,c,d) is first index
// Two options (second and fourth may be swapped): choose 2nd to be smallest
// Wk symmetry:
// {abcd} = badc = cdab = dcba

//! Class to store Q^k integrals, and similar, accounting for symmetry.
//! @details Mostly, a wrapper for std::map
template <typename Real = float, Symmetry symmetry = Symmetry::none>
class QkTable {

  using Index = uint16_t;
  using IndexSet = std::array<Index, 4>;
#if defined(_OPENMP)
  static constexpr bool use_omp = true;
#else
  static constexpr bool use_omp = false;
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

private:
  std::map<IndexSet, std::vector<Real>> m_data{};

public:
  // 'rule of zero' for now seems fine
  // QkTable() {}

  //! Gives arrow access to all underlying map functions
  auto operator-> () { return &m_data; }

  //! Provide these for ranged for loops etc.
  auto begin() { return m_data.begin(); }
  auto cbegin() const { return m_data.cbegin(); }
  auto end() { return m_data.end(); }
  auto cend() const { return m_data.cend(); }

  //! Fill all non-zero Qk integrals (nb: often don't need all)
  //! @details Note: You can also just 'add' indevedidual terms using ->insert()
  void fill(const YkTable &yk) {
    if constexpr (use_omp && symmetry == Symmetry::Qk) {
      fill_Rk(yk);
    } else if constexpr (use_omp && symmetry == Symmetry::Wk) {
      fill_Wk(yk);
    } else if constexpr (use_omp) {
      fill_noSymm(yk);
    } else {
      fill_series(yk);
    }
  }

  //! tmp
  void count() const {
    std::size_t num = 0;
    for (const auto &[key, value] : m_data) {
      num += value.size();
      // None of the vectors should be empty - just a consistancy check
      if (value.size() == 0)
        std::cout << "**\n";
    }
    std::cout << "Qk table contains: " << num << " non-zero elements\n";
  }

  //! Add a specific new Q (for each k) to map. If exists, does nothing
  void addQ(const DiracSpinor &a, const DiracSpinor &b, const DiracSpinor &c,
            const DiracSpinor &d, const YkTable &yk) {
    [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
    // min/max multipolarity (accounts for parity):
    const auto [kmin, kmax] = Coulomb::k_minmax_Q(a, b, c, d);
    if (kmin > kmax)
      return;
    const auto [it, new_element] =
        m_data.insert({NormalOrder(a, b, c, d), std::vector<Real>{}});
    // If already exists, do nothing
    if (!new_element)
      return;
    auto &[key, value] = *it;
    assert(value.empty());
    const auto num_ks = std::size_t((kmax - kmin) / 2 + 1);
    value.reserve(num_ks);
    for (int k = kmin; k <= kmax; k += 2) {
      const auto yk_bd = yk.ptr_yk_ab(k, b, d);
      assert(yk_bd != nullptr);
      value.push_back(Real(Qk_abcd(a, b, c, d, k, *yk_bd, yk.Ck())));
    }
  }

  //! Updates a specific new Q (for each k) in map. If absent, does nothing!
  void updateQ(const DiracSpinor &a, const DiracSpinor &b, const DiracSpinor &c,
               const DiracSpinor &d, const YkTable &yk) {
    const auto map_it = m_data.find(NormalOrder(a, b, c, d));
    if (map_it == m_data.end())
      return;
    auto &[key, value] = *map_it;
    const auto [kmin, kmax] = Coulomb::k_minmax_Q(a, b, c, d);
    auto vec_it = value.begin();
    for (int k = kmin; k <= kmax; k += 2) {
      assert(vec_it != value.end());
      const auto yk_bd = yk.ptr_yk_ab(k, b, d);
      assert(yk_bd != nullptr);
      *vec_it = Real(Qk_abcd(a, b, c, d, k, *yk_bd, yk.Ck()));
      vec_it++;
    }
  }

  //! Return (const) pointer stored Qk (each k). If absent, nullptr
  const std::vector<Real> *Qk(const DiracSpinor &a, const DiracSpinor &b,
                              const DiracSpinor &c,
                              const DiracSpinor &d) const {
    [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
    const auto Norm_abcd = NormalOrder(a, b, c, d);
    {
      [[maybe_unused]] auto sp2 = IO::Profile::safeProfiler(__func__, "find");
      const auto it = m_data.find(Norm_abcd);
      if ((it == m_data.cend()))
        return nullptr;
      const auto &[key, value] = *it;
      return &value;
    }
  }

  //! Returns specific Q^k_abcd. If absent, zero
  Real Q(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const {
    [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
    const auto qk = Qk(a, b, c, d);
    if (qk == nullptr)
      return Real{0.0};
    const auto [kmin, kmax] = Coulomb::k_minmax_Q(a, b, c, d);
    if ((k < kmin) || (k > kmax) || (k % 2 != kmin % 2))
      return Real{0.0};
    const auto index = std::size_t((k - kmin) / 2);
    assert(qk->size() > index);
    return (*qk)[index];
  }

  //! Uses stored Q's to return R (see Coulomb for definition). Zero if absent
  Real R(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const {
    [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
    const auto tQk = Q(k, a, b, c, d);
    const auto s = Angular::neg1pow(k);
    const auto tCkac = Angular::tildeCk_kk(k, a.k, c.k);
    const auto tCkbd = Angular::tildeCk_kk(k, b.k, d.k);
    return tQk / Real(s * tCkac * tCkbd);
  }

  //! Uses stored Q's to return P. Zero if absent
  Real P(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const {
    [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
    Real Pk_abcd{0.0};
    const auto Ql_abdc = Qk(a, b, d, c); // exchange!
    if (Ql_abdc == nullptr)
      return Pk_abcd;
    const auto [lmin, lmax] = Coulomb::k_minmax_Q(a, b, d, c);
    int l = lmin;
    for (auto &q_abdc : *Ql_abdc) {
      assert(l <= lmax);
      const auto sixj = Real(Coulomb::sixj(a, c, k, b, d, l));
      Pk_abcd += sixj * q_abdc;
      l += 2;
    }
    Pk_abcd *= Real(2 * k + 1);
    return Pk_abcd;
  }

  //! Uses stored Q's to return W (W = Q + P). Zero if absent
  Real W(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const {
    return Q(k, a, b, c, d) + P(k, a, b, c, d);
  }

public:
  //! Returns index set {a,b,c,d} in "Normal order" (depends on Symmetry)
  IndexSet NormalOrder(const DiracSpinor &a, const DiracSpinor &b,
                       const DiracSpinor &c, const DiracSpinor &d) const {
    return NormalOrder(a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index());
  }

  IndexSet NormalOrder(int a, int b, int c, int d) const {
    return NormalOrder(Index(a), Index(b), Index(c), Index(d));
  }

  IndexSet NormalOrder(const IndexSet &abcd) const {
    const auto [a, b, c, d] = abcd;
    return NormalOrder(a, b, c, d);
  }

  IndexSet NormalOrder(Index a, Index b, Index c, Index d) const {
    if constexpr (symmetry == Symmetry::Qk) {
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
  IndexSet NormalOrder_Rk(Index a, Index b, Index c, Index d) const {
    [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
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

  IndexSet NormalOrder_Wk(Index a, Index b, Index c, Index d) const {
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

  IndexSet NormalOrder_none(Index a, Index b, Index c, Index d) const {
    [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
    return {a, b, c, d};
  }

  //****************************************************************************
  void fill_Rk(const YkTable &yk) {
    m_data.clear();
    assert(yk.a_is_b());
    const auto &basis = yk.get_a();

    const auto a_size = basis.size();
    std::vector<QkTable<Real, symmetry>> maps(a_size);

// To ensure no collisions in below merges, enforce symmetry here
// This is only for efficiency, no change to logic without this
#pragma omp parallel for
    for (auto i = 0ul; i < a_size; ++i) {
      auto &tmap = maps[i];
      const auto &a = basis[i];
      for (const auto &b : basis) {
        if (a > b)
          continue;
        for (const auto &c : basis) {
          if (a > c)
            continue;
          for (const auto &d : basis) {
            if (a > d || b > d)
              continue;
            tmap.addQ(a, b, c, d, yk);
          }
        }
      }
    }

    for (auto &map : maps) {
      m_data.merge(map.m_data);
      map->clear();
    }
  }

  //----------------------------------------------------------------------------
  void fill_Wk(const YkTable &yk) {
    m_data.clear();
    assert(yk.a_is_b());
    const auto &basis = yk.get_a();

    const auto a_size = basis.size();
    std::vector<QkTable<Real, symmetry>> maps(a_size);

// ensure no collisions in below merges, by enforcing symmetry here
#pragma omp parallel for
    for (auto i = 0ul; i < a_size; ++i) {
      auto &tmap = maps[i];
      const auto &a = basis[i];
      for (const auto &b : basis) {
        if (a > b)
          continue;
        for (const auto &c : basis) {
          if (a > c)
            continue;
          for (const auto &d : basis) {
            if (a > d)
              continue;
            tmap.addQ(a, b, c, d, yk);
          }
        }
      }
    }

    for (auto &map : maps) {
      m_data.merge(map.m_data);
      map->clear();
    }
  }

  //----------------------------------------------------------------------------
  void fill_noSymm(const YkTable &yk) {
    m_data.clear();
    assert(yk.a_is_b());
    const auto &basis = yk.get_a();

    const auto a_size = basis.size();
    std::vector<QkTable<Real, symmetry>> maps(a_size);

#pragma omp parallel for
    for (auto i = 0ul; i < a_size; ++i) {
      auto &tmap = maps[i];
      const auto &a = basis[i];
      for (const auto &b : basis) {
        for (const auto &c : basis) {
          for (const auto &d : basis) {
            tmap.addQ(a, b, c, d, yk);
          }
        }
      }
    }

    for (auto &map : maps) {
      m_data.merge(map.m_data);
      map->clear();
    }
  }

  //----------------------------------------------------------------------------
  void fill_series(const YkTable &yk) {
    m_data.clear();
    assert(yk.a_is_b());
    const auto &basis = yk.get_a();

    if constexpr (symmetry == Symmetry::Qk) {
      for (const auto &a : basis) {
        for (const auto &b : basis) {
          if (a > b)
            continue;
          for (const auto &c : basis) {
            if (a > c)
              continue;
            for (const auto &d : basis) {
              if (a > d || b > d)
                continue;
              continue;
              addQ(a, b, c, d, yk);
            }
          }
        }
      }
    } else if constexpr (symmetry == Symmetry::Wk) {
      for (const auto &a : basis) {
        for (const auto &b : basis) {
          if (a > b)
            continue;
          for (const auto &c : basis) {
            if (a > c)
              continue;
            for (const auto &d : basis) {
              if (a > d)
                continue;
              continue;
              addQ(a, b, c, d, yk);
            }
          }
        }
      }
    } else {
      for (const auto &a : basis) {
        for (const auto &b : basis) {
          for (const auto &c : basis) {
            for (const auto &d : basis) {
              addQ(a, b, c, d, yk);
            }
          }
        }
      }
    }
  }
}; // namespace Coulomb

} // namespace Coulomb
