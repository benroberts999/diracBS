#pragma once
#include "Coulomb/QkTable.hpp"
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

  return pass;
}

} // namespace UnitTest
