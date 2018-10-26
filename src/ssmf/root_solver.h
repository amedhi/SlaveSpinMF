/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-19 11:24:03
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-04-19 11:26:10
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef ROOT_SOLVER_H
#define ROOT_SOLVER_H
#include <iostream>
#include <memory>
#include <vector>
#include <Eigen/Dense>
#include "../lattice/matrix.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

namespace root {

class gsl_solver 
{
public:
  gsl_solver();
  ~gsl_solver() {}
  void allocate(const std::size_t& n);
  //void set_problem(int (*f) (const gsl_vector * x, void * params, gsl_vector * f),
  //  const std::size_t& n, void * params);
  //void find_root(const std::vector<double>& x_init);
  int find_root(void* params, int (*func) (const gsl_vector* x, void* params, gsl_vector* fvec),
    std::vector<double>& xvec, const double& eps_ftol=1.0E-6);
private:
  std::size_t ndim_{0};
  gsl_multiroot_function func_;
  std::unique_ptr<gsl_vector> xvec_{nullptr};
  std::unique_ptr<const gsl_multiroot_fsolver_type> solver_type_{nullptr};
  std::unique_ptr<gsl_multiroot_fsolver> solver_{nullptr};
};


class RootSolver
{
public:
  RootSolver() { fx_.resize(ndim_); J_.resize(ndim_,ndim_); }
  ~RootSolver() {}
  void init(const std::size_t& n) { ndim_=n; fx_.resize(ndim_); J_.resize(ndim_,ndim_); }
  template<class F>
  void solve(F func, const RealVector& x0);
private:
  std::size_t ndim_{1};
  RealVector fx_;
  RealMatrix J_;
  void line_search(const RealVector& x_old, const RealVector& fx_old, 
    const double& fnorm, const RealVector& grad_f, const RealVector &pvec, 
    const RealVector& x, const RealVector& fx);
};

template<class F>
void RootSolver::solve(F func, const RealVector& x0)
{
  double fnorm;
  RealVector grad_f(ndim_);
  RealVector x(ndim_);
  RealVector x_old(ndim_);
  RealVector fx_old(ndim_);
  RealVector pvec(ndim_);

  bool need_derivative = true;
  int max_iter = 1;
  x = x0;
  for (int iter=0; iter<max_iter; ++iter) {
    func(x, fx_, J_, need_derivative);
    fnorm = 0.5*fx_.squaredNorm();
    // calculate Grad(fx) for the line search
    grad_f = J_.transpose() * fx_;
    // store x and fx
    x_old = x;
    fx_old = fx_;
    // solve the linear equations
    pvec = -J_.colPivHouseholderQr().solve(fx_);
    std::cout << "pvec = "<<pvec.transpose()<<"\n";
    line_search(x_old, fx_old, fnorm, grad_f, pvec, x, fx_);
  }

}



} // end namespace root
#endif
