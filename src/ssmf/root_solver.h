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
  RootSolver() { init(1); }
  ~RootSolver() {}
  void init(const std::size_t& n) 
    { 
      ndim_=n; 
      x_.resize(ndim_);
      fx_.resize(ndim_); 
      Jmat_.resize(ndim_,ndim_); 
      x_old_.resize(ndim_);
      fx_old_.resize(ndim_);
      fgrad_.resize(ndim_);
      pvec_.resize(ndim_);
    }
  template<class F>
  int solve(F func, RealVector& x);
private:
  std::size_t ndim_{1};
  RealVector x_;
  RealVector fx_;
  RealMatrix Jmat_;
  RealVector x_old_;
  RealVector fx_old_;
  RealVector fgrad_;
  RealVector pvec_;
  double fnorm_;
  double fnorm_old_;

  int MAXITER_{100};
  double STPMX_{100.0};
  double LFUNC_{1.0E-6};
  double ftol_{1.0E-9};
  double min_tol_{1.0E-9};
  double xtol_{1.0E-9};
  double stpmax_{1.0};

  template<class F>
  int line_search(F func, const RealVector& x_old, const double& fnorm_old, 
    const RealVector& fgrad, const RealVector &pvec, const double& stpmax,
    RealVector& x, RealVector& fx, double& fnorm);
};

template<class F>
int RootSolver::solve(F func, RealVector& x0)
{
  // Test for initial guess being root
  x_ = x0;
  func(x_,fx_,Jmat_,false);
  if (fx_.array().abs().maxCoeff() < 0.01*ftol_) return 0;
  // stpmax for line searches
  stpmax_ = STPMX_*std::max(x_.norm(),double(ndim_));
  // start iteration
  for (int iter=0; iter<MAXITER_; ++iter) {
    func(x_,fx_,Jmat_,true);
    // calculate Grad(fx) for the line search
    fgrad_ = Jmat_.transpose()*fx_;
    // the norm function
    fnorm_ = 0.5*fx_.squaredNorm();
    // store x, fx, fnorm
    x_old_ = x_;
    fx_old_ = fx_;
    fnorm_old_ = fnorm_;
    // Search direction by solving the linear equations 
    pvec_ = -Jmat_.colPivHouseholderQr().solve(fx_);
    // do line search
    int info = line_search(func,x_old_,fnorm_old_,fgrad_,pvec_,stpmax_,x_,fx_,fnorm_);
    if (info == 2) return info;
    // line search returns 'x_' , 'fx_' & 'fnorm_'
    // std::cout << "fnorm = " << fnorm_ << "\n";
    // check for convergence on function values
    if (fx_.array().abs().maxCoeff() < ftol_) {
      x0 = x_;
      return 0;
    } 
    // check for grad_f zero case
    if (info == 1) {
      double test = 0.0;
      double den = std::max(fnorm_,0.5*ndim_);
      for (int i=0; i<ndim_; ++i) {
        double tmp = std::abs(fgrad_(i))*std::max(std::abs(x_(i)),1.0)/den;
        if (test < tmp) test = tmp;
      }
      x0 = x_;
      return test < min_tol_ ? 1 : 0;
    }
    // check for convergence on del_x
    double dxmax = 0.0;
    for (int i=0; i<ndim_; ++i) {
      double dx = std::abs(x_(i)-x_old_(i))/std::max(std::abs(x_(i)),1.0);
      if (dxmax < dx) dxmax = dx;
    }
    if (dxmax < xtol_) {
      x0 = x_;
      return 0;
    }
  }
  std::cout << "** RootSolver::solve: iteration exceeded\n";
  //getchar();
  return 1;
}

template<class F>
int RootSolver::line_search(F func, const RealVector& x_old, const double& fnorm_old, 
  const RealVector& fgrad, const RealVector& pvec, const double& stpmax,
  RealVector& x, RealVector& fx, double& fnorm)
{
  double lambda, lambda2, lambda_tmp, lambda_min;
  double fnorm2, slope;

  // scale if attempted step is too big
  RealVector spvec = pvec;
  double stplen = pvec.norm();
  if (stplen > stpmax) {
    double scale = stpmax/stplen;
    spvec *= scale;
  }

  // slope
  slope = fgrad.dot(spvec);
  if (slope >= 0.0) {
    std::cout << "** RootSolver::line_search: Roundoff problem\n";
    return 2;
  }

  // lambda_min
  double d = 0.0;
  for (int i=0; i<ndim_; ++i) {
    double p = std::abs(spvec(i))/std::max(std::abs(x_old(i)),1.0);
    if (d < p) d = p;
  }
  lambda_min = xtol_/d;

  // Always try full Newton step first
  lambda = 1.0;
  while (true) {
    x = x_old + lambda * spvec;
    //std::cout << "x = " << x.transpose() << "\n"; getchar();
    func(x,fx,Jmat_,false);
    fnorm = 0.5*fx.squaredNorm();
    if (lambda < lambda_min) {
      // Convergence on \Delta x. For zero finding the calling 
      // program should verify its convergence.
      x = x_old;
      return 1;
    }
    else if (fnorm <= fnorm_old + LFUNC_*lambda*slope) {
      return 0;
    }
    else {
      // Backtrack
      if (lambda == 1.0) {
        // first time 
        lambda_tmp = -slope/(2.0*(fnorm - fnorm_old - slope));
      }
      else {
        // subsequent backtracks
        double rhs1 = fnorm-fnorm_old-lambda*slope;
        double rhs2 = fnorm2-fnorm_old-lambda2*slope;
        double a = (rhs1/(lambda*lambda)-rhs2/(lambda2*lambda2))/(lambda-lambda2);
        double b = (-lambda2*rhs1/(lambda*lambda)+lambda*rhs2/(lambda2*lambda2))
                    /(lambda-lambda2);
        if (std::abs(a)<1.0E-15) {
          lambda_tmp = -slope/(2.0*b);
        }
        else {
          double disc = b*b - 3.0*a*slope;
          if (disc < 0.0) lambda_tmp = 0.5*lambda;
          else if (b <= 0.0) lambda_tmp = (-b+std::sqrt(disc))/(3.0*a);
          else lambda_tmp = -slope/(b+std::sqrt(disc));
        }
        // lambda <= 0.5*lambda1
        if (lambda_tmp > 0.5*lambda) lambda_tmp = 0.5*lambda;
      }
    }
    lambda2 = lambda;
    fnorm2 = fnorm;
    // lambda => 0.1*lambda1
    lambda = std::max(lambda_tmp, 0.1*lambda);
  }

}




} // end namespace root
#endif
