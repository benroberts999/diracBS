#pragma once
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <tuple>
#include <utility>
#include <vector>

//! Defines SqMatrix, Vector, and linear-algebra solvers (incl Eigensystems)
namespace LinAlg {

//******************************************************************************
//! Basic Square matrix class of constant construct-time size
class SqMatrix {

public:
  const std::size_t n = 0;
  gsl_matrix *m = nullptr;

public:
  SqMatrix() {}
  SqMatrix(std::size_t in_n);
  // SqMatrix(const std::initializer_list<double> &l);
  ~SqMatrix();
  SqMatrix(const SqMatrix &matrix);
  SqMatrix &operator=(const SqMatrix &other);

public:
  //! Constructs a diagonal unit matrix (identity)
  void make_identity();
  //! Sets all elements to zero
  void zero();
  //! All values with |Mij|<value are set to zero
  void clip_low(double value);
  //! All values with |Mij|>value are set to +/-value
  void clip_high(double value);

  //! Forces Matrix to be symmetric by: Mij -> (Mij + Mji)/2
  void enforce_symmetric(); // change to "symmetrise"?
  //! Returns largest elemt of matrix [Mij - Mji]; zero if symmetric
  [[nodiscard]] double check_symmetric();
  //! Prints Matrix to screen (for tests)
  void print();

  //! Returns the transpose of matrix: not destructive
  [[nodiscard]] SqMatrix transpose() const;
  //! Determinate via LU decomp. Note: expensive.
  [[nodiscard]] double determinant() const;
  //! Inverts the matrix: nb: destructive
  void invert();
  //! Returns the inverce of matrix: not destructive
  [[nodiscard]] SqMatrix inverse() const;

  double *operator[](int i) const;
  double *operator[](std::size_t i) const;
  friend SqMatrix operator*(const SqMatrix &lhs, const SqMatrix &rhs);
  SqMatrix &operator+=(const SqMatrix rhs);
  friend SqMatrix operator+(SqMatrix lhs, const SqMatrix &rhs);
  SqMatrix &operator-=(const SqMatrix rhs);
  friend SqMatrix operator-(SqMatrix lhs, const SqMatrix &rhs);
  SqMatrix &operator*=(const double x);
  friend SqMatrix operator*(const double x, SqMatrix rhs);

  void mult_elements_by(const SqMatrix &rhs);
  friend SqMatrix mult_elements(SqMatrix lhs, const SqMatrix &rhs);
};

//******************************************************************************
template <typename T> class Complex {
public:
  T re = 0;
  T im = 0;
  Complex<T> conj() const { return {re, -im}; }
  //! norm2 = re^2 + im^2, no sqrt (ruins T)
  T norm2() const { return re * re + im * im; }
  //! Mult two complex:
  friend Complex<T> operator*(const Complex<T> &a, const Complex<T> &b) {
    return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re};
  }
  Complex<T> &operator*=(const Complex<T> x) {
    *this = (*this) * x;
    return *this;
  }
  //! Mult by const:
  Complex<T> &operator*=(const T x) {
    this->re *= x;
    this->im *= x;
    return *this;
  }
  friend Complex<T> operator*(Complex<T> a, const T &x) { return a *= x; }
  friend Complex<T> operator*(const T &x, Complex<T> a) { return a *= x; }
  //! Add/subtrac:
  Complex<T> &operator+=(const Complex<T> x) {
    this->re += x.re;
    this->im += x.im;
    return *this;
  }
  friend Complex<T> operator+(Complex<T> a, const Complex<T> &b) {
    return a += b;
  }
  Complex<T> &operator-=(const Complex<T> x) {
    this->re -= x.re;
    this->im -= x.im;
    return *this;
  }
  friend Complex<T> operator-(Complex<T> a, const Complex<T> &b) {
    return a -= b;
  }
};

//******************************************************************************
//! Basic Complex Square matrix class of constant construct-time size
class ComplexSqMatrix {

public:
  const std::size_t n = 0;
  gsl_matrix_complex *m = nullptr;

public:
  ComplexSqMatrix() {}
  ComplexSqMatrix(std::size_t in_n);
  ~ComplexSqMatrix();
  ComplexSqMatrix(const ComplexSqMatrix &other);
  ComplexSqMatrix &operator=(const ComplexSqMatrix &other);

public:
  //! Constructs a diagonal unit matrix
  void make_identity();
  //! Sets all elements to zero
  void zero();

  //! Returns the transpose of matrix: not destructive
  [[nodiscard]] ComplexSqMatrix transpose() const;
  //! Inverts the matrix: nb: destructive
  void invert();
  //! Returns the inverce of matrix: not destructive
  [[nodiscard]] ComplexSqMatrix inverse() const;

  //! make a ComplexSqMatrix from a SqMatrix: C = x*mR, x is complex
  static ComplexSqMatrix make_complex(const Complex<double> &x,
                                      const SqMatrix &mR);

  //! For testing??
  Complex<double> get_copy(std::size_t i, std::size_t j) const;

  //! Get the real part (copy) of the complex matrix
  [[nodiscard]] SqMatrix real() const;
  //! Get the imaginary part (copy) of the complex matrix
  [[nodiscard]] SqMatrix imaginary() const;
  //! multiply elements in place Aij -> Aij*Bij
  void mult_elements_by(const ComplexSqMatrix &rhs);
  //! multiply elements Aij = Bij*Cij
  friend ComplexSqMatrix mult_elements(ComplexSqMatrix lhs,
                                       const ComplexSqMatrix &rhs);

  //! Multiply elements by constant
  ComplexSqMatrix &operator*=(const Complex<double> &x);
  friend ComplexSqMatrix operator*(const Complex<double> &x,
                                   ComplexSqMatrix rhs);
  friend ComplexSqMatrix operator*(ComplexSqMatrix rhs,
                                   const Complex<double> &x);

  //! Matrix multiplication of 2 C-mats
  friend ComplexSqMatrix operator*(const ComplexSqMatrix &x,
                                   const ComplexSqMatrix &y);

  //! Add + subtract Complex matrices
  ComplexSqMatrix &operator+=(const ComplexSqMatrix &rhs);
  ComplexSqMatrix &operator-=(const ComplexSqMatrix &rhs);
  friend ComplexSqMatrix operator+(ComplexSqMatrix lhs,
                                   const ComplexSqMatrix &rhs);
  friend ComplexSqMatrix operator-(ComplexSqMatrix lhs,
                                   const ComplexSqMatrix &rhs);
};

//******************************************************************************
//! Basic vector class of constant construct-time size
class Vector {
public:
  const std::size_t n;
  gsl_vector *vec;

public:
  Vector(const std::size_t in_n);
  template <typename T> Vector(const std::initializer_list<T> &l);
  Vector &operator=(const Vector &other);
  Vector(const Vector &other);
  ~Vector();

  void clip_low(const double value);
  void clip_high(const double value);
  void print();

  double &operator[](int i) const;
  double &operator[](std::size_t i) const;
  Vector &operator+=(const Vector rhs);
  friend Vector operator+(Vector lhs, const Vector &rhs);
  Vector &operator-=(const Vector rhs);
  friend Vector operator-(Vector lhs, const Vector &rhs);
  Vector &operator*=(const double x);
  friend Vector operator*(const double x, Vector rhs);
  friend Vector operator*(const SqMatrix &Aij, const Vector &bj);

  friend double inner_product(const Vector &a, const Vector &b);
  friend double operator*(const Vector &a, const Vector &b);
  friend SqMatrix outer_product(const Vector &a, const Vector &b);
};

//******************************************************************************
//******************************************************************************

//! Solves Matrix equationL A*x = b for x
Vector solve_Axeqb(const SqMatrix &Am, const Vector &b);

//------------------------------------------------------------------------------
//! @brief Solves Av = ev for eigenvalues e and eigenvectors v
//! for Real Symmetric Matrices
//! @details Eigensystems: Note: A and B are destroyed! Don't use afterwards
//! (Can't avoid this without needless copy). Optionally sorts e and v by e
[[nodiscard]] std::pair<Vector, SqMatrix>
realSymmetricEigensystem(SqMatrix &A, bool sort = true);

//! @brief Solves Av = eBv for eigenvalues e and eigenvectors v for Real
//! Generalized Symmetric-Definite Eigensystems
[[nodiscard]] std::pair<Vector, SqMatrix>
realSymmetricEigensystem(SqMatrix &A, SqMatrix &B, bool sort = true);

//! @briefSolves for Av = ev for Real Nonsymmetric Matrices.
//! @details e and v are complex; returned as {real, imag} (seperate
//! vectors/Matrix)
[[nodiscard]] std::tuple<Vector, Vector, SqMatrix, SqMatrix>
realNonSymmetricEigensystem(SqMatrix &A, bool sort = true);

//! @brief Solves Av = eBv for eigenvalues e and eigenvectors v for Real
//! Generalized Non-Symmetric-Definite Eigensystems.
[[nodiscard]] std::tuple<Vector, Vector, SqMatrix, SqMatrix>
realNonSymmetricEigensystem(SqMatrix &A, SqMatrix &B, bool sort = true);

} // namespace LinAlg
