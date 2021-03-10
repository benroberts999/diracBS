#pragma once
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
  std::map<IndexSet, std::vector<Real>> m_data;

public:
  QkTable();

  std::pair<IndexSet, std::vector<Real>> *get(IndexSet abcd) {
    auto it = m_data.find(NormalOrder_Rk(abcd));
    if (it != m_data.end())
      return &*it;
  }

private:
  IndexSet NormalOrder(const IndexSet &abcd) {
    const auto [a, b, c, d] = abcd;
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
