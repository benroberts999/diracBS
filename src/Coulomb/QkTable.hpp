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

  // Gives arrow access to all map functions
  auto operator-> () { return &m_data; }

  //! Fill all non-zero Qk integrals (nb: often don't need all)
  void fill(const YkTable &yk);
  //! Update all non-zero Qk integrals (nb: often don't need all)
  void update(const YkTable &yk);

  //! Add a specific new Q (for each k) to map. If exists, does nothing
  void addQ(const DiracSpinor &a, const DiracSpinor &b, const DiracSpinor &c,
            const DiracSpinor &d, const YkTable &yk);
  //! Updates a specific new Q (for each k) in map. If abscent, does nothing!
  void updateQ(const DiracSpinor &a, const DiracSpinor &b, const DiracSpinor &c,
               const DiracSpinor &d, const YkTable &yk);

  // These assume zero of not found

  //! Return stored Qk (each k). If abscent, empty
  const std::vector<Real> &Qk(const DiracSpinor &a, const DiracSpinor &b,
                              const DiracSpinor &c, const DiracSpinor &d);
  //! Returns specific Q^k_abcd. If abscent, zero
  Real Q(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d);
  //! Uses stored Q's to return R. Zero if abscent
  Real R(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d);
  //! Uses stored Q's to return P. Zero if abscent
  Real P(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d);
  //! Uses stored Q's to return W. Zero if abscent
  Real W(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d);

  // Functions to construct each Q? Use Coulomb..
  // Qk(a,b,c,d);
  // Q(k,a,b,c,d);
  // R(k,a,b,c,d);
  // W(k,a,b,c,d);
  // P(k,a,b,c,d);

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
