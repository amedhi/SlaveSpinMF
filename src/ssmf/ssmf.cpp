/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-22 13:00:28
*----------------------------------------------------------------------------*/
#include <iomanip>
#include <fstream>
#include "ssmf.h"

namespace srmf {


SRMF::SRMF(const input::Parameters& inputs) 
  : graph_(inputs) 
  //, blochbasis_(graph_)
  , model_(inputs, graph_.lattice())
  //, sr_parms_(inputs, graph_,model_)
  , mf_params_(inputs, graph_,model_)
  , spinon_model_(inputs,model_,graph_,mf_params_)
  , boson_model_(inputs,model_,graph_,mf_params_)
{
  // For solving for LM parameter equation
  fx_dim_ = graph_.lattice().num_basis_orbitals();
  x_vec_.resize(fx_dim_);
  fx_vec_.resize(fx_dim_);
  solver_.allocate(fx_dim_);
}

int gsl_problem_eqn1(const gsl_vector* x, void* parms, gsl_vector* f)
{
  SRMF * pThis = ((class SRMF *) parms);
  for (int i=0; i<pThis->fx_dim_; ++i) {
    pThis->x_vec_[i] = gsl_vector_get(x,i);
  }
  int status = pThis->spinon_energy_eqn(pThis->x_vec_, pThis->fx_vec_);
  for (int i=0; i<pThis->fx_dim_; ++i) {
    gsl_vector_set(f, i, pThis->fx_vec_[i]);
  }
  if (status ==0 ) return GSL_SUCCESS;
  else return GSL_FAILURE;
}

int SRMF::spinon_energy_eqn(const std::vector<double>& x, std::vector<double>& fx)
{
  spinon_model_.set_shifted_en(x);
  selconsistent_solve();
  auto lambda = mf_params_.site(0).lm_params();
  std::cout<<"x(0) = "<< x[0] << " x(1) ="<<x[1]<<"\n";
  std::cout<<"lambda = "<< lambda.transpose()<<"\n"; getchar();
  for (int i=0; i<fx_dim_; ++i) {
    fx[i] = x[i] - lambda[i] - spinon_model_.orbital_en()[i];
    std::cout << "fx["<<i<<"] = "<< fx[i] << "\n";
  }
  return 0;
}

int SRMF::run(const input::Parameters& inputs) 
{
  spinon_model_.update(inputs);
  boson_model_.update(spinon_model_);

  //mf_params_.
  boson_ke_.resize(mf_params_.num_bonds());
  spinon_ke_.resize(mf_params_.num_bonds());
  diff_boson_ke_.resize(mf_params_.num_bonds());
  diff_spinon_ke_.resize(mf_params_.num_bonds());
  boson_ke_norm_.resize(mf_params_.num_bonds());
  spinon_ke_norm_.resize(mf_params_.num_bonds());

  for (int i=0; i<mf_params_.num_bonds(); ++i) {
    int rows = mf_params_.bond(i).spinon_ke(0).rows();
    int cols = mf_params_.bond(i).spinon_ke(0).cols();
    boson_ke_[i].resize(rows, cols);
    spinon_ke_[i].resize(rows, cols);
    diff_boson_ke_[i].resize(rows, cols);
    diff_spinon_ke_[i].resize(rows, cols);
  }

  mf_params_.init_params();
  for (int i=0; i<mf_params_.num_bonds(); ++i) {
    spinon_ke_[i] = mf_params_.bond(i).spinon_ke(0);
    boson_ke_[i] = mf_params_.bond(i).boson_ke(0);
  }

  // find spinon shifted local energies for U=0;
  /*
  double U = spinon_model_.get_parameter_value("U");
  if (std::abs(U)<1.0E-12) {
    for (int i=0; i<fx_dim_; ++i) {
      x_vec_[i] = spinon_model_.orbital_en()[i];
    }
    solver_.find_root(this, &gsl_problem_eqn1, x_vec_, lm_ftol_);
    spinon_model_.set_shifted_en(x_vec_);
  }
  std::cout << "shifted_e0 =" << x_vec_[0] << "\n";
  */
  // solve 
  selconsistent_solve();
  return 0;
}

int SRMF::selconsistent_solve(void) 
{
  int max_sb_iter = 10;
  bool converged = false;
  for (int iter=0; iter<max_sb_iter; ++iter) {
    spinon_model_.solve(graph_,mf_params_);
    boson_model_.solve(mf_params_);
    // check convergence
    for (int i=0; i<mf_params_.num_bonds(); ++i) {
      diff_spinon_ke_[i] = spinon_ke_[i] - mf_params_.bond(i).spinon_ke(0);
      diff_boson_ke_[i] = boson_ke_[i] - mf_params_.bond(i).boson_ke(0);
      spinon_ke_norm_[i] = diff_spinon_ke_[i].abs2().maxCoeff();
      boson_ke_norm_[i] = diff_boson_ke_[i].abs2().maxCoeff();
    }
    double sp_norm = spinon_ke_norm_.maxCoeff();
    double sb_norm = boson_ke_norm_.maxCoeff();

    //std::cout<<"ssmf iter="<<iter+1<<", norm=("<<sp_norm<<","<<sb_norm<<")\n";
    if (sp_norm<1.0E-6 && sb_norm<1.0E-6) {
      converged = true;
      break;
    }

    // continue
    for (int i=0; i<mf_params_.num_bonds(); ++i) {
      spinon_ke_[i] = mf_params_.bond(i).spinon_ke(0);
      boson_ke_[i] = mf_params_.bond(i).boson_ke(0);
    }
  }
  //if (converged) std::cout<<"ssmf converged!\n";

  std::cout<<spinon_model_.get_parameter_value("U")<<"  ";
  for (int i=0; i<mf_params_.num_sites(); ++i) {
    //std::cout<<"Z["<<i<<"] = "<< sr_parms_.site(i).qp_weights().transpose()<<"\n";
    std::cout<<mf_params_.site(i).qp_weights().transpose()<<"\n";
    std::cout<<mf_params_.site(i).lm_params().transpose()<<"\n";
  }
  return 0;
}


void SRMF::print_copyright(std::ostream& os)
{
  std::cout << "SRMF\n";
}


} // end namespace diag
