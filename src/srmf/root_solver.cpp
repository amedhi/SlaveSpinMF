/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@macbook
* @Date:   2018-09-18 17:33:36
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-09-20 23:34:27
*----------------------------------------------------------------------------*/
#include <cassert>
#include "root_solver.h"

namespace root {

gsl_solver::gsl_solver(void) 
{
  solver_type_.reset(gsl_multiroot_fsolver_hybrids);
}

void gsl_solver::allocate(const std::size_t& n)
{
  ndim_ = n;
  xvec_.reset(gsl_vector_alloc(ndim_));
  solver_.reset(gsl_multiroot_fsolver_alloc(solver_type_.get(), ndim_));
}

/*void gsl_solver::set_problem(int (*func) (const gsl_vector * x, void * params, gsl_vector * f),
    const std::size_t& n, void * params)
{
  ndim_ = n;
  func_ = {func, ndim_, params};
  xvec_.reset(gsl_vector_alloc(ndim_));
  solver_.reset(gsl_multiroot_fsolver_alloc(solver_type_.get(), ndim_));
  gsl_multiroot_fsolver_set(solver_.get(),&func_,xvec_.get());
}*/

void gsl_solver::find_root(void* params, int (*func) (const gsl_vector* x, void* params, gsl_vector* fvec),
  const std::vector<double>& x_init)
{
  assert(x_init.size() == ndim_);
  func_ = {func, ndim_, params};
  for (int i=0; i<ndim_; ++i) gsl_vector_set(xvec_.get(), i, x_init[i]);
  //for (int i=0; i<ndim_; ++i) std::cout << gsl_vector_get(xvec_.get(), i) << "\n";
  gsl_multiroot_fsolver_set(solver_.get(), &func_, xvec_.get());
  int status = GSL_CONTINUE;
  int iter = 0;
  while (status==GSL_CONTINUE && iter < 50) {
    int status = gsl_multiroot_fsolver_iterate(solver_.get());
    ++iter;
  }
  std::cout << "gsl_solver::find_root: status = " << gsl_strerror(status) << "\n";
}


//void gsl_solver::find_root(const std::vector<double>& x_init)
//{
//  assert(x_init.size() == ndim_);
//  for (int i=0; i<ndim_; ++i) gsl_vector_set(xvec_.get(), i, x_init[i]);
//}


} // end namespace
