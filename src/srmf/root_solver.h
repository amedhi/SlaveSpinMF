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
  void find_root(void* params, int (*func) (const gsl_vector* x, void* params, gsl_vector* fvec),
  const std::vector<double>& x_init);
private:
  std::size_t ndim_{0};
  gsl_multiroot_function func_;
  std::unique_ptr<gsl_vector> xvec_{nullptr};
  std::unique_ptr<const gsl_multiroot_fsolver_type> solver_type_{nullptr};
  std::unique_ptr<gsl_multiroot_fsolver> solver_{nullptr};
};




} // end namespace srmf
#endif
