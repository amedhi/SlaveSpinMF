/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@macbook
* @Date:   2018-09-18 17:33:36
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-10-26 10:32:00
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

int gsl_solver::find_root(void* params, int (*func) (const gsl_vector* x, void* params, gsl_vector* fvec),
  std::vector<double>& xvec, const double& eps_ftol)
{
  assert(xvec.size() == ndim_);
  func_ = {func, ndim_, params};
  for (int i=0; i<ndim_; ++i) gsl_vector_set(xvec_.get(), i, xvec[i]);
  //for (int i=0; i<ndim_; ++i) std::cout << gsl_vector_get(xvec_.get(), i) << "\n";
  gsl_multiroot_fsolver_set(solver_.get(), &func_, xvec_.get());
  int status = GSL_CONTINUE;
  int iter = 0;
  while (status==GSL_CONTINUE && iter < 50) {
    //std::cout << "gsl iter = " << iter+1 << "\n";
    status = gsl_multiroot_fsolver_iterate(solver_.get());
    status = gsl_multiroot_test_residual(solver_.get()->f, eps_ftol);
    ++iter;
  }
  if (status != GSL_SUCCESS) {
    std::cout << ">> gsl_solver::find_root: iteration exited with status '" 
    << gsl_strerror(status)<< "'\n";
  }
  // solution
  for (int i=0; i<ndim_; ++i) xvec[i] = gsl_vector_get(solver_.get()->x,i);
  //gsl_multiroot_fsolver_free(solver_.get());
  //gsl_vector_free(xvec_.get());
  return status;
}



//void gsl_solver::find_root(const std::vector<double>& x_init)
//{
//  assert(x_init.size() == ndim_);
//  for (int i=0; i<ndim_; ++i) gsl_vector_set(xvec_.get(), i, x_init[i]);
//}


} // end namespace
