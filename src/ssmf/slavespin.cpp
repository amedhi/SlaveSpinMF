/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-19 11:24:03
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-06-19 22:59:00
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "slavespin.h"
#include <stdexcept>
#include <string>
#include <cassert>
#include <boost/algorithm/string.hpp>
#include <boost/math/tools/roots.hpp>


namespace srmf {

SlaveSpin::SlaveSpin(const input::Parameters& inputs, const model::Hamiltonian& model, const lattice::LatticeGraph& graph, const MF_Params& mf_params)
  //: rotor_graph_(graph)
{
  // SlaveSpin lattice has only one original lattice unit cell 
  num_sites_ = mf_params.num_sites();
  num_bonds_ = mf_params.num_bonds();

  // Assuming all sites have same 'site_dim'.
  spin_orbitals_ = mf_params.site(0).spin_orbitals(); // including 'UP' & 'DN' spins
  site_dim_ = spin_orbitals_.size(); 
  /*
    'spin_orbitals' Note:  If NO SOC, the assumption is that
    the first half of the spin-orbitals correspond to UP-spins 
    and the second half correspond to DOWN-spins.
    If SOC, then the spin-orbitals will follow the same ordering
    as in the site basis in the given Hamiltonian.
  */

  // storages
  // Lagrange multipliers
  lm_params_.resize(num_sites_);
  for (auto& elem : lm_params_) elem = realArray1D::Ones(site_dim_);
  //lm_params_noint_.resize(num_sites_);
  //for (auto& elem : lm_params_noint_) elem = realArray1D::Zero(site_dim_);
  // gauge factors, c for the operator O+ = (S- + cS+)
  gauge_factors_.resize(num_sites_);
  for (auto& elem : gauge_factors_) elem = realArray1D::Ones(site_dim_);
  // order parameters: <O+>
  site_order_params_.resize(num_sites_);
  for (auto& elem : site_order_params_) elem = cmplArray1D::Ones(site_dim_);
  // <O+O->
  site_avg_OplusMinus_.resize(num_sites_);
  for (auto& elem : site_avg_OplusMinus_) elem = cmplArray1D::Ones(site_dim_);
  // Quasi-particle weights
  qp_weights_.resize(num_sites_);
  for (auto& elem : qp_weights_) elem = realArray1D::Ones(site_dim_);
  // Renormalized site couplings
  renorm_site_couplings_.resize(num_sites_);
  for (auto& elem : renorm_site_couplings_) elem = cmplArray1D::Ones(site_dim_);
  // Renormalized bond couplings
  renorm_bond_couplings_.resize(num_bonds_);
  for (auto& elem : renorm_bond_couplings_) 
    elem = cmplArray2D::Ones(site_dim_,site_dim_);
  // <Sz> values
  Sz_avg_.resize(num_sites_);
  for (auto& elem : Sz_avg_) elem.resize(site_dim_);

  spinon_density_.resize(num_sites_);
  for (auto& elem : spinon_density_) elem.resize(site_dim_);

  // boson model parameters  
  switch (model.id2()) {
    case model::model_id2::HUBBARD:
      modelparams_.set_id(HUBBARD);
      modelparams_.update_U(model.get_parameter_value("U"));
      break;
    case model::model_id2::BHZ:
      modelparams_.set_id(BHZ);
      //modelparams_.update_e0(model.get_parameter_value("e0"));
      modelparams_.update_U(model.get_parameter_value("U"));
      break;
    case model::model_id2::PYROCHLORE:
      modelparams_.set_id(PYROCHLORE);
      modelparams_.update_U(model.get_parameter_value("U"));
      modelparams_.update_J(model.get_parameter_value("J"));
      modelparams_.update_lambda(model.get_parameter_value("lambda"));
      break;
    default: 
      throw std::range_error("*error: SlaveSpin: 'model' not implemented"); 
  }

  // SlaveSpin clusters
  std::string theory_name = inputs.set_value("theory", "Z2");
  boost::to_upper(theory_name);
  if (theory_name=="Z2") theory_ = theory_t::Z2;
  else if (theory_name=="U1") theory_ = theory_t::U1;
  else throw std::range_error("*error: SlaveSpin: invalid 'theory' name"); 
  // whether solve only for single site
  solve_single_site_ = inputs.set_value("solve_single_site", false);

  make_clusters(mf_params);
  for (auto& cluster : clusters_) 
    cluster.init_hamiltonian(modelparams_,gauge_factors_,lm_params_,renorm_site_couplings_);

  // For solving for LM parameter equation
  fx_dim_ = num_sites_*site_dim_;
  x_vec_.resize(fx_dim_);
  fx_vec_.resize(fx_dim_);
  gsl_solver_.allocate(fx_dim_);
}

void SlaveSpin::update(const model::Hamiltonian& model)
{
  switch (modelparams_.id()) {
    case HUBBARD:
      modelparams_.update_U(model.get_parameter_value("U"));
      break;
    case BHZ:
      //modelparams_.update_e0(model.get_parameter_value("e0"));
      modelparams_.update_U(model.get_parameter_value("U"));
      break;
    case PYROCHLORE:
      modelparams_.update_U(model.get_parameter_value("U"));
      modelparams_.update_J(model.get_parameter_value("J"));
      modelparams_.update_lambda(model.get_parameter_value("lambda"));
      break;
    default: break;
  }
  for (auto& cluster : clusters_) {
    cluster.update_interaction_matrix(modelparams_);
  }
}

void SlaveSpin::solve(MF_Params& mf_params) 
{
  // set constrained density
  for (int i=0; i<num_sites_; ++i) {
    spinon_density_[i] = mf_params.site(i).spinon_density();
  }
  for (auto& cluster : clusters_) {
    cluster.set_spinon_density(spinon_density_);
  }

  // gauge factors
  if (!gauge_factors_set_) {
    for (int site=0; site<num_sites_; ++site) {
      for (auto alpha : spin_orbitals_) {
        double nf = spinon_density_[site][alpha];
        if (theory_==theory_t::Z2) {
          //delta_ = 0.0;
          gauge_factors_[site][alpha]=1.0/std::sqrt((delta_+nf)*(1.0-nf+delta_))-1.0;
        }
        else {// U1 Theory
          gauge_factors_[site][alpha]=1.0/std::sqrt((nf+delta_)*(1.0-nf+delta_)); 
        }
        //std::cout << "nf["<<site<<"]["<<alpha<<"] = "<<nf<<"\n";
        //std::cout << "c["<<site<<"]["<<alpha<<"] = "<<gauge_factors_[site][alpha]<<"\n";
        //getchar();
      }
    }
    //gauge_factors_set_ = true;
  }

  // set SOC couplings
  if (modelparams_.id()==PYROCHLORE) {
    //if (!gauge_factors_solved_) {
    //  solve_gauge_factors(mf_params);
    //  gauge_factors_solved_ = true;
    //}
    for (auto& cluster : clusters_) {
      cluster.update_soc_matrix(mf_params, gauge_factors_);
    }
  }

  // renormalized bond couplings
  set_bond_couplings(mf_params);

  // --------Self Consistent solution-------
  self_consistent_solve(mf_params);

  // QP weights
  for (int site=0; site<num_sites_; ++site) {
    qp_weights_[site] = site_order_params_[site].abs2();
    //std::cout<<"Z["<<site<<"] = "<< qp_weights_[site].transpose()<<"\n";
    //std::cout<<"<Z+>["<<site<<"] = "<< site_order_params_[site].transpose()<<"\n";
  }
  // update boson parameters
  for (int i=0; i<num_sites_; ++i) {
    mf_params.site(i).lm_params() = lm_params_[i];
    mf_params.site(i).qp_weights() = qp_weights_[i];
    //std::cout<<"lambda["<<i<<"] = "<< lm_params_[i].transpose()<<"\n";
    //std::cout<<"Z["<<i<<"] = "<< qp_weights_[i].transpose()<<"\n";
  }
  if (modelparams_.id()==PYROCHLORE) {
    set_renormalized_soc(mf_params);
  }
  update_bond_order_params(mf_params);
  update_renorm_site_potential(mf_params);

  /*
  // Solve for 'lambda' for the corresponding 'non-interacting' problem
  double U = 0.0;
  for (auto& cluster : clusters_) {
    cluster.update_parameters(U);
  }
  self_consistent_solve(srparams);
  for (int i=0; i<num_sites_; ++i) {
    srparams.site(i).lm_params_noint() = lm_params_[i];
    std::cout<<"lambda0["<<i<<"] = "<< lm_params_[i].transpose()<<"\n";
  }
  // restore interactions 
  for (auto& cluster : clusters_) {
    cluster.update_parameters(U_);
  }
  //getchar();
  */

}

void SlaveSpin::solve_gauge_factors(const MF_Params& mf_params)
{
  // set interaction term zero
  double U = modelparams_.get_U();
  double J = modelparams_.get_J();
  modelparams_.update_U(0.0);
  modelparams_.update_J(0.0);
  for (auto& cluster : clusters_) {
    cluster.update_interaction_matrix(modelparams_);
  }
  set_bond_couplings(mf_params);

  realArray1D a(site_dim_);
  realArray1D b(site_dim_);
  realArray1D c(site_dim_);
  realArray1D func_a(site_dim_);
  realArray1D func_b(site_dim_);
  realArray1D func_c(site_dim_);
  int max_iter = 50;
  for (int site=0; site<clusters_.size(); ++site) {
    a = realArray1D::Constant(site_dim_,1.0); 
    b = realArray1D::Constant(site_dim_,10.0);;
    gauge_factors_[site] = a;
    func_a = gauge_factors_func(mf_params, site);
    gauge_factors_[site] = b;
    func_b = gauge_factors_func(mf_params, site);

    for (const auto& alpha : spin_orbitals_) {
      // bisection
      if (func_a(alpha)*func_b(alpha)>=0.0) {
        throw std::range_error("SlaveSpin::solve_gauge_factors: roots not bracketed\n");
      }
      int iter;
      for (iter=0; iter<max_iter; ++iter) {
        double c_alpha = 0.5*(a[alpha] + b[alpha]);
        gauge_factors_[site][alpha] = c_alpha;
        func_c = gauge_factors_func(mf_params, site);
        if (std::abs(func_c[alpha])<1.0E-6 || 0.5*std::abs(a[alpha]-b[alpha])<1.0E-6) break;
        if (func_a[alpha] * func_c[alpha]<0.0) {
          b[alpha] = c_alpha;
          func_b[alpha] = func_c[alpha];
        }
        else {
          a[alpha] = c_alpha;
          func_a[alpha] = func_c[alpha];
        }
      }
      if (iter >= max_iter) {
        std::cout << "** SlaveSpin::solve_gauge_factors: iteration exceeded\n";
      }
      //std::cout << "func_a = " << func_a.transpose() << "\n";
      //std::cout << "func_b = " << func_b.transpose() << "\n";
      //std::cout << "gauge_factors_ = " << gauge_factors_[site].transpose() << "\n";
      //getchar();
    }
  }
  // restore interaction term zero
  modelparams_.update_U(U);
  modelparams_.update_J(J);
  for (auto& cluster : clusters_) {
    cluster.update_interaction_matrix(modelparams_);
  }
}

realArray1D SlaveSpin::gauge_factors_func(const MF_Params& mf_params, const int& i)
{
  cmpl_siteparms_t trial_order_params(num_sites_);
  for (auto& elem : trial_order_params) elem = cmplArray1D::Constant(site_dim_,1.0);
  cmplArray1D order_params_diff(site_dim_);
  int max_iter = 100;
  bool converged = false;

  clusters_[i].update_soc_matrix(mf_params, gauge_factors_);
  for (int iter=0; iter<max_iter; ++iter) {
    set_site_couplings(mf_params, trial_order_params);
    clusters_[i].update_hamiltonian(gauge_factors_,renorm_site_couplings_);
    clusters_[i].solve_lm_params(lm_params_);
    clusters_[i].update_hamiltonian(lm_params_);
    clusters_[i].solve_hamiltonian();
    clusters_[i].get_avg_Ominus(gauge_factors_,site_order_params_);
    order_params_diff = site_order_params_[i]-trial_order_params[i];
    double norm = order_params_diff.abs2().maxCoeff();
    if (norm<1.0E-6) {
      converged = true; 
      break;
    } 
    trial_order_params[i] = site_order_params_[i];
  }
  // qp_weights should be 1
  qp_weights_[i] = site_order_params_[i].abs2();
  return qp_weights_[i]-realArray1D::Constant(site_dim_,1.0);
}


void SlaveSpin::self_consistent_solve(const MF_Params& mf_params)
{
  // trial 'qp_weights'
  cmpl_siteparms_t trial_order_params(num_sites_);
  for (auto& elem : trial_order_params) elem = cmplArray1D::Constant(site_dim_,1.0);
  cmplArray1D order_params_diff(num_sites_*site_dim_);

  int max_iter = 100;
  bool converged = false;
  bool print_progress = false;
  for (int iter=0; iter<max_iter; ++iter) {
    // update 'renormalized site couplings' for the new 'order parameters'
    set_site_couplings(mf_params, trial_order_params);

    for (auto& cluster : clusters_) {
      cluster.update_hamiltonian(gauge_factors_,renorm_site_couplings_);
    }
    // solve for new LM-parameters to satisfy the 'slave-spin constraint'
    update_lm_params();
    /*for (int site=0; site<num_sites_; ++site) {
      std::cout<<"lambda["<<site<<"] = "<< lm_params_[site].transpose()<<"\n";
    }
    getchar();
    */
    // update cluster hamiltonians for new LM-parameters
    for (auto& cluster : clusters_) {
      cluster.update_hamiltonian(lm_params_);
      /*for (int i=0; i<cluster.hamiltonian_matrix().rows(); ++i) {
        std::cout << "h("<<i<<","<<i<<") = "<<cluster.hamiltonian_matrix()(i,i) << "\n";
      }
      getchar();*/
    }
    // calculate new site op-s
    update_site_order_params();
    /*
    for (int site=0; site<num_sites_; ++site) {
      std::cout << "<Zplus>["<<site<<"] = " << site_order_params_[site].transpose() << "\n";
    }
    getchar();
    */
    

    // check convergence
    int i = 0;
    for (int site=0; site<num_sites_; ++site) {
      for (auto& alpha: spin_orbitals_) {
        order_params_diff(i) = site_order_params_[site][alpha]-trial_order_params[site][alpha];
        ++i; 
      }
    }
    double norm = order_params_diff.abs2().maxCoeff();
    if (print_progress) {
      std::cout << "boson iter = " << iter+1 << ", norm = " << norm << "\n";
    }
    if (norm<1.0E-6) {
      converged = true;
      break;
    } 
    // continue
    for (int site=0; site<num_sites_; ++site) {
      trial_order_params[site] = site_order_params_[site];
    }
  }
  if (converged) {
    if(print_progress) std::cout<<"Bosons converged!\n";
  } 
}

void SlaveSpin::set_bond_couplings(const MF_Params& mf_params) 
{
  // Assuming spinorbitals includes both UP & DOWN spins
  int num_orb = static_cast<int>(site_dim_)/2;
  assert(2*num_orb == static_cast<int>(site_dim_));
  for (int i=0; i<num_bonds_; ++i) {
    renorm_bond_couplings_[i] = mf_params.bond(i).spinon_renormed_cc(0);
    // renorm_bond_couplings_[i].setOnes();
    //std::cout << "bond_field["<<i<<"]=\n"<< renorm_bond_couplings_[i] << "\n"; getchar();
  }
}

void SlaveSpin::set_site_couplings(const MF_Params& mf_params, 
  const cmpl_siteparms_t& site_order_params) 
{
  // Renormalizing field on each site due to hopping terms from the 'bath'
  if (cluster_type_==cluster_t::SITE) {
    for (auto& elem : renorm_site_couplings_) elem.setZero();
    for (int i=0; i<mf_params.num_bonds(); ++i) {
      auto s = mf_params.bond(i).src();
      auto t = mf_params.bond(i).tgt();

      // for source site
      auto tchi = renorm_bond_couplings_[i].matrix();
      // tchi_phi = Matrix(tchi) * Vector(phi)
      auto tchi_phi = tchi * site_order_params[t].matrix(); 
      renorm_site_couplings_[s] += tchi_phi.array(); 

      // for target site
      auto tchi_t = renorm_bond_couplings_[i].matrix().adjoint();
      auto tchi_phi_t = tchi_t * site_order_params[s].matrix(); 
      renorm_site_couplings_[t] += tchi_phi_t.array(); 
    }
    // take only the real part
    /*
    for (int i=0; i<num_sites_; ++i) {
      for (int m=0; m<renorm_site_couplings_[i].rows(); ++m) {
        for (int n=0; n<renorm_site_couplings_[i].cols(); ++n) {
          renorm_site_couplings_[i](m,n) = std::real(renorm_site_couplings_[i](m,n));
        }
      }
    }*/


    /*
    int id = 0;
    for (const auto& bond : mf_params.bonds()) {
      auto s = bond.src();
      auto t = bond.tgt();

      // partial sum over 'orbital' index of the neighbour site
      auto tchi_phi = renorm_bond_couplings_[id++].rowwise()*site_order_params[t].transpose();
      auto tchi_phi_sum = tchi_phi.rowwise().sum();
      //std::cout << tchi_phi_sum.transpose() << "\n"; getchar();

      // sum over all neighbouring sites
      renorm_site_couplings_[s] += tchi_phi_sum; // * site_order_params[t];
      renorm_site_couplings_[t] += tchi_phi_sum.conjugate(); 
    }
    */
    /* 
    for (int i=0; i<num_sites_; ++i) {
      std::cout << "site field ["<<i<<"] = " << renorm_site_couplings_[i].transpose() << "\n"; 
    }
    getchar();
    */  
  }
  else if (cluster_type_ == cluster_t::BOND) {
  }
  else if (cluster_type_ == cluster_t::CELL) {
  }
  else {
  }
}

void SlaveSpin::set_renormalized_soc(MF_Params& mf_params)
{
  if (cluster_type_==cluster_t::SITE) {
    cmplArray2D flip_ampl(site_dim_,site_dim_);
    for (int i=0; i<mf_params.num_sites(); ++i) {
      for (int m=0; m<site_dim_; ++m) {
        for (int n=0; n<site_dim_; ++n) {
          flip_ampl(m,n) = site_order_params_[i](m) * std::conj(site_order_params_[i](n));
        }
      }
      mf_params.site(i).boson_flip_ampl() = flip_ampl;
      //std::cout << "flip_ampl["<<i<<"] =\n" << flip_ampl << "\n";
      //getchar();
      mf_params.site(i).set_boson_renormalization();
    }
  }
  else {
  }
}

void SlaveSpin::update_bond_order_params(MF_Params& mf_params) 
{
  // Renormalizing field on each site due to hopping terms from the 'bath'
  cmplArray2D bond_op(site_dim_,site_dim_);
  if (cluster_type_==cluster_t::SITE) {
    for (int i=0; i<mf_params.num_bonds(); ++i) {
      auto s = mf_params.bond(i).src();
      auto t = mf_params.bond(i).tgt();
      for (int m=0; m<site_dim_; ++m) {
        for (int n=0; n<site_dim_; ++n) {
          bond_op(m,n) = site_order_params_[s](m) * std::conj(site_order_params_[t](n));
        }
      }
      mf_params.bond(i).boson_ke(0) = bond_op;
      //std::cout << "bond_op["<<i<<"] =\n" << bond_op << "\n";
      mf_params.bond(i).set_boson_renormalization();
    }
  }
  else if (cluster_type_ == cluster_t::BOND) {
  }
  else if (cluster_type_ == cluster_t::CELL) {
  }
  else {
  }
}

void SlaveSpin::update_renorm_site_potential(MF_Params& mf_params) 
{
  if (theory_==Z2) {
    for (int site=0; site<num_sites_; ++site) {
      mf_params.site(site).lm_params() = lm_params_[site];
    }
  }
  else { // U1 Theory
    cmplArray2D ee = cmplArray2D::Zero(site_dim_,site_dim_);
    for (int i=0; i<mf_params.num_bonds(); ++i) {
      ee += mf_params.bond(i).spinon_renormed_cc(0) * 
            mf_params.bond(i).boson_ke(0);
    }
    //std::cout << ee << "\n\n"; getchar();
    auto ebar = 4.0*ee.rowwise().sum().real();
    // gauge factors
    realArray1D mu(site_dim_);
    for (int site=0; site<num_sites_; ++site) {
      for (auto alpha : spin_orbitals_) {
        double nf = spinon_density_[site][alpha];
        double eta = (2.0*nf-1)/(4.0*nf*(1.0-nf));
        mu[alpha] = 2*ebar[alpha]*eta;
        //std::cout << "mu["<<site<<"]["<<alpha<<"] = "<<mu[alpha]<<"\n";
        //std::cout << "lambda["<<site<<"]["<<alpha<<"] = "<<lm_params_[site][alpha]<<"\n";
      }
      //std::cout << "nf["<<site<<"] = "<<spinon_density_[site].transpose() <<"\n";
      //std::cout << "lambda["<<site<<"] = "<<srparams.site(site).lm_params().transpose() <<"\n";
      //std::cout << "mu_gen["<<site<<"] = "<<mu.transpose() <<"\n";
      mf_params.site(site).lm_params() = (lm_params_[site]-mu);
      //std::cout << "lambda["<<site<<"] = "<<srparams.site(site).lm_params().transpose() <<"\n";
    }
  }
  //getchar();
  /*
  real_siteparms_t Splus_avg;
  Splus_avg.resize(num_sites_);
  for (auto& elem : Splus_avg) elem = realArray1D::Zero(site_dim_);
  for (auto& cluster : clusters_) {
    cluster.get_avg_Splus(Splus_avg);
  }
  */
}

int gsl_problem_equation(const gsl_vector* x, void* parms, gsl_vector* f)
{
  SlaveSpin * pThis = ((class SlaveSpin *) parms);
  for (int i=0; i<pThis->fx_dim_; ++i) {
    pThis->x_vec_[i] = gsl_vector_get(x,i);
  }
  int status = pThis->constraint_equation(pThis->x_vec_, pThis->fx_vec_);
  for (int i=0; i<pThis->fx_dim_; ++i) {
    gsl_vector_set(f, i, pThis->fx_vec_[i]);
  }
  if (status ==0 ) return GSL_SUCCESS;
  else return GSL_FAILURE;
}

int SlaveSpin::constraint_equation(const std::vector<double>& x, std::vector<double>& fx)
{
  // read the new LM-parameters
  int i=0;
  for (int site=0; site<num_sites_; ++site) {
    for (auto& alpha: spin_orbitals_) {
      lm_params_[site][alpha] = x[i];
      ++i; 
      //std::cout << "lambda = " << lm_params_[site][alpha] << "\n";
    }
  }
  //getchar();
  // calculate <Sz> for the new parameters
  for (auto& cluster : clusters_) {
      cluster.update_hamiltonian(lm_params_);
      cluster.solve_hamiltonian();
      cluster.get_avg_Sz(Sz_avg_);
      //cluster.groundstate_dLambda(); 
      // groundstate derivative
      //std::cout << "dX/dLambda =\n"<< cluster.groundstate_dLambda() << "\n";
      //getchar();
  }
  // LHS of the constraint equation: fx = (<Sz> + 1/2) - n_f
  i=0;
  for (int site=0; site<num_sites_; ++site) {
    for (auto& alpha: spin_orbitals_) {
      fx[i] = (Sz_avg_[site][alpha]+0.5) - spinon_density_[site][alpha];
      ++i; 
      //std::cout << "fx = " << fx[i] << "\n";
    }
  }
  //getchar();
  return 0;
}

void SlaveSpin::update_lm_params(void)
{
  if (solve_single_site_) {
    clusters_[0].solve_lm_params(lm_params_);
    for (int i=1; i<clusters_.size(); ++i) {
      lm_params_[i] = lm_params_[0];
    }
  }
  else {
    for (auto& cluster : clusters_) {
      cluster.solve_lm_params(lm_params_);
    }
  }

  //for (unsigned site=0; site<num_sites_; ++site) {
  //  std::cout << "lambda["<<site<<"] = " << lm_params_[site].transpose() << "\n";
  //}
  /*
  //solver_.set_problem(&lagrange_eqn, 2, this);
  //std::vector<double> x0(num_sites_*site_dim_, 1.0); 
  for (auto& x : x_vec_) x = 1.0;
  solver_.find_root(this, &gsl_problem_equation, x_vec_, lm_ftol_);
  // solved 'LM' parameters
  int i = 0;
  for (unsigned site=0; site<num_sites_; ++site) {
    for (auto& alpha: spin_orbitals_) {
      lm_params_[site][alpha] = x_vec_[i];
      std::cout << "lambda_alpha["<<i<<"] = " << x_vec_[i] << "\n";
      ++i; 
    }
  }
  getchar();
  */
}

void SlaveSpin::update_site_order_params(void)
{
  // calculate <O+> 
  if (theory_==Z2) {
    if (solve_single_site_) {
      clusters_[0].solve_hamiltonian();
      clusters_[0].get_avg_Ominus(gauge_factors_,site_order_params_);
      for (int i=1; i<clusters_.size(); ++i) {
        site_order_params_[i] = site_order_params_[0];
      }
    }
    else {
      for (auto& cluster : clusters_) {
        cluster.solve_hamiltonian();
        cluster.get_avg_Ominus(gauge_factors_,site_order_params_);
      }
    }
  }
  else { // U1 Theory
    if (solve_single_site_) {
      clusters_[0].solve_hamiltonian();
      clusters_[0].get_avg_Zminus(gauge_factors_,site_order_params_);
      for (int i=1; i<clusters_.size(); ++i) {
        site_order_params_[i] = site_order_params_[0];
      }
    }
    else {
      for (auto& cluster : clusters_) {
        cluster.solve_hamiltonian();
        cluster.get_avg_Zminus(gauge_factors_,site_order_params_);
      }
    }
  }
  /*
  for (unsigned site=0; site<num_sites_; ++site) {
    std::cout << "<Oplus>["<<site<<"] = " << site_order_params_[site].transpose() << "\n";
  }
  getchar();
  */
}

void SlaveSpin::make_clusters(const MF_Params& mf_params)
{
  cluster_type_ = cluster_t::SITE;
  //bonds_ = mf_params.bonds();
  clusters_.clear();
  switch (cluster_type_) {
    case cluster_t::SITE:
      for (unsigned i=0; i<num_sites_; ++i) {
        clusters_.push_back({cluster_t::SITE,theory_,i,mf_params.site(i).spin_orbitals()});
      }
      break;
    case cluster_t::BOND:
      break;
    case cluster_t::CELL:
      break;
  }
}

void Cluster::init_hamiltonian(const ModelParams& p, const real_siteparms_t& gauge_factors, 
  const real_siteparms_t& lm_params, const cmpl_siteparms_t& site_fields)
{
  if (type_!=cluster_t::SITE) {
    throw std::range_error("Cluster::init_hamiltonian: defined only for 'SITE' cluster\n");
  }
  interaction_elems_.setZero();
  lagrange_elems_.setZero();
  soc_mat_.setZero();
  hmatrix_.setZero();
}

// update 'interaction terms'
void Cluster::update_interaction_matrix(const ModelParams& p)
{
  int site_id = 0;
  SlaveSpinBasis::idx_t i, j;
  // interaction matrix elements
  if (p.id()==HUBBARD || p.id()==BHZ) {
    double U_half = 0.5*p.get_U();
    double Sz;
    for (i=0; i<basis_dim_; ++i) {
      double total_Sz = 0.0;
      for (auto& alpha : spin_orbitals_) {
        std::tie(Sz,j) = basis_.apply_Sz(site_id,alpha,i);
        total_Sz += Sz;
      }
      interaction_elems_(i) = U_half*total_Sz*total_Sz;
      //std::cout << "H_U=" << interaction_elems_(i)  << "\n";
    }
  }
  else if (p.id()==PYROCHLORE) {
    // |basis> = |0U, 0D, 1U, 1D, 2U, 2D>
    double U = p.get_U();
    double J = p.get_J();
    //double J = 0.1*U;
    double Uprime = U - 2.0*J; 
    double U_one_body = -0.5*(14.0*U - 37.0*J);
    double Sz, Sz_up, Sz_dn;
    for (i=0; i<basis_dim_; ++i) {
      // U term
      double sum = 0.0;
      for (int m=0; m<spin_orbitals_.size()-1; m+=2) {
        std::tie(Sz_up,j) = basis_.apply_Sz(site_id,m,i);
        std::tie(Sz_dn,j) = basis_.apply_Sz(site_id,m+1,i);
        sum += Sz_up*Sz_dn;
      }
      interaction_elems_(i) = U*sum;
      //std::cout << "H_U=" << interaction_elems_(i)  << "\n";

      // Uprime term
      sum = 0.0;
      for (int m=0; m<spin_orbitals_.size()-1; m+=2) {
        std::tie(Sz_up,j) = basis_.apply_Sz(site_id,m,i);
        for (int n=0; n<spin_orbitals_.size()-1; n+=2) {
          if (m == n) continue;
          std::tie(Sz_dn,j) = basis_.apply_Sz(site_id,n+1,i);
          sum += Sz_up*Sz_dn;
        }
      }
      interaction_elems_(i) += Uprime*sum;

      // U-J term
      sum = 0.0;
      double Sz2_up, Sz2_dn;
      for (int m=0; m<spin_orbitals_.size()-1; m+=2) {
        std::tie(Sz_up,j) = basis_.apply_Sz(site_id,m,i);
        std::tie(Sz_dn,j) = basis_.apply_Sz(site_id,m+1,i);
        for (int n=0; n<m; n+=2) {
          std::tie(Sz2_up,j) = basis_.apply_Sz(site_id,n,i);
          std::tie(Sz2_dn,j) = basis_.apply_Sz(site_id,n+1,i);
          //std::cout << "m, n = " << m/2+1 << "  " << n/2+1 << "\n"; 
          sum += Sz_up*Sz2_up + Sz_dn*Sz2_dn;
        }
      }
      interaction_elems_(i) += (Uprime-J)*sum;
      //getchar();

      // One body terms
      ///*
      double total_Sz = 0.0;
      for (auto& alpha : spin_orbitals_) {
        std::tie(Sz,j) = basis_.apply_Sz(site_id,alpha,i);
        total_Sz += Sz;
      }
      interaction_elems_(i) += U_one_body*total_Sz;
      //*/
      //std::cout << "H_U=" << interaction_elems_(i)  << "\n";
    }
  }
  else {
    throw std::range_error("Cluster::init_hamiltonian: undefined model\n");
  }
}

void Cluster::update_soc_matrix(const MF_Params& mfp, 
  const real_siteparms_t& gauge_factors)
{
  SlaveSpinBasis::idx_t i,j,k;
  double mat_elem, mat_elem2;
  // Spin-Orbit coupling term
  soc_mat_.setZero();
  int site_id = 0;
  if (theory_==theory_t::Z2) {
    for (j=0; j<basis_dim_; ++j) {
      for (auto& beta : spin_orbitals_) {
        double c = gauge_factors[site_][beta];
        std::tie(mat_elem,k) = basis_.apply_Ominus(c,site_id,beta,j);
        if (k == basis_.null_idx()) continue;
        for (auto& alpha : spin_orbitals_) {
          double c2 = gauge_factors[site_][alpha];
          std::tie(mat_elem2,i) = basis_.apply_Oplus(c2,site_id,alpha,k);
          if (i != basis_.null_idx()) {
            soc_mat_(i,j) += mat_elem*mat_elem2*
              mfp.site(site_).spinon_renormed_soc()(alpha,beta);
            // hc term taken care by the loops over 'alpha' & 'beta'
          }
        }
      }
    }
  } // end Z2 theory
  else {// U1 Theory
    for (j=0; j<basis_dim_; ++j) {
      for (auto& beta : spin_orbitals_) {
        double c = gauge_factors[site_][beta];
        std::tie(mat_elem,k) = basis_.apply_Zminus(c,site_id,beta,j);
        if (k == basis_.null_idx()) continue;
        for (auto& alpha : spin_orbitals_) {
          double c2 = gauge_factors[site_][alpha];
          std::tie(mat_elem2,i) = basis_.apply_Zplus(c2,site_id,alpha,k);
          if (i != basis_.null_idx()) {
            soc_mat_(i,j) += mat_elem*mat_elem2*
              mfp.site(site_).spinon_renormed_soc()(alpha,beta);
            // hc term taken care by the loops over 'alpha' & 'beta'
          }
        }
      }
    }
  } // end U1 Theory
  //std::cout << "soc_mat = " << soc_mat_ << "\n";
  //getchar();
}

// update for new 'site_couplings'
void Cluster::update_hamiltonian(const real_siteparms_t& gauge_factors,
  const cmpl_siteparms_t& new_site_couplings)
{
  int site_id = 0;
  //std::cout << "rotor::Cluster::update_hamiltonian: ---Problem-------\n";
  SlaveSpinBasis::idx_t i, j;
  double mat_elem;
  hmatrix_.setZero();
  if (theory_==theory_t::Z2) {
    // Hopping term (in main Hamiltonian)
    for (j=0; j<basis_dim_; ++j) {
      for (auto& alpha : spin_orbitals_) {
        double c = gauge_factors[site_][alpha]; // 'site_' here is not mistake
        std::tie(mat_elem,i) = basis_.apply_Oplus(c,site_id,alpha,j);
        if (i != basis_.null_idx()) {
          auto term = mat_elem * new_site_couplings[site_][alpha];
          hmatrix_(i,j) += term;
          hmatrix_(j,i) += std::conj(term);
          //std::cout << "H["<<i<<","<<j<<"]=" << term << "\n"; getchar();
        }
      }
    }
  }
  else { // U1 Theory
    for (j=0; j<basis_dim_; ++j) {
      for (auto& alpha : spin_orbitals_) {
        double c = gauge_factors[site_][alpha];
        std::tie(mat_elem,i) = basis_.apply_Zplus(c,site_id,alpha,j);
        if (i != basis_.null_idx()) {
          auto term = mat_elem * new_site_couplings[site_][alpha];
          hmatrix_(i,j) += term;
          hmatrix_(j,i) += std::conj(term);
        //std::cout<<"elem("<<i<<","<<j<<") = "<<term<<"\n";
        //hmatrix_(i,j) = std::conj(term);
        //hmatrix_(j,i) += std::conj(term);
        //std::cout << "H["<<i<<","<<j<<"]=" << term << "\n"; getchar();
        }
      }
    }
  }
  // add SOC part
  hmatrix_ += soc_mat_;
  // diagonal elements
  for (i=0; i<basis_dim_; ++i) {
    hmatrix_(i,i) = interaction_elems_(i)+lagrange_elems_(i);
    //std::cout << "hmatrix("<<i<<","<<i<<") = "<<hmatrix_(i,i)<<"\n"; getchar();
  }
}

// update for new 'lm parameters'
void Cluster::update_hamiltonian(const real_siteparms_t& new_lm_params)
{
  int site_id = 0;
  // lagrange fields
  SlaveSpinBasis::idx_t i, j;
  double Sz;
  for (i=0; i<basis_dim_; ++i) {
    double sum = 0.0;
    for (auto& alpha : spin_orbitals_) {
      std::tie(Sz,j) = basis_.apply_Sz(site_id,alpha,i);
      sum += (Sz+0.5) * new_lm_params[site_][alpha];
    }
    lagrange_elems_(i) = sum;
    //std::cout << "lagrange_elems["<<i<<"]="<< sum << "\n"; getchar();
  }

  for (i=0; i<basis_dim_; ++i) {
    hmatrix_(i,i) = interaction_elems_(i)+lagrange_elems_(i);
  }
  //std::cout<<"H (LM-update) =\n" << hmatrix_ << "\n"; getchar();
  // final matrix 
}

void Cluster::solve_hamiltonian(void) const
{
  //std::cout << "hmatrix =\n" << hmatrix_ << "\n"; getchar();
  eigen_solver_.compute(hmatrix_);
  groundstate_ = eigen_solver_.eigenvectors().col(0);
  //std::cout << "Eigenvalues = " << eigen_solver_.eigenvalues().transpose() << "\n";
  //std::cout << "Eigenvector = " << groundstate_.transpose() << "\n";
  //getchar();
}

void Cluster::set_spinon_density(const real_siteparms_t& spinon_density) 
{
  spinon_density_[site_] = spinon_density[site_];
}

// update for new 'lm parameters'
void Cluster::solve_lm_params(real_siteparms_t& lm_params)
{
  RealVector lambda(total_spinorbitals_);
  for (int i=0; i<total_spinorbitals_; ++i) lambda[i] = 1.0;
  root_solver_.solve([this](const RealVector& x, RealVector& fx, 
    RealMatrix& J, const bool& need_derivative) 
    {return lambda_equation(x,fx,J,need_derivative);}, lambda);
  // read the new LM-parameters
  int i=0;
  for (auto& alpha: spin_orbitals_) {
    lm_params[site_][alpha] = lambda[i];
    ++i; 
  }
  //std::cout<<"lambda["<<site_<<"] = "<< lm_params[site_].transpose()<<"\n";
  //getchar();
}

int Cluster::lambda_equation(const RealVector& lambda, RealVector& func, 
  RealMatrix& Jmat, const bool& need_derivative) 
{
  int site_id = 0;
  // read the new LM-parameters
  int n=0;
  for (auto& alpha: spin_orbitals_) {
    lm_params_[site_][alpha] = lambda[n];
    ++n; 
  }
  update_hamiltonian(lm_params_);
  solve_hamiltonian();
  // average 'Sz'
  RealVector Sz_avg = RealVector::Zero(spin_orbitals_.size());
  SlaveSpinBasis::idx_t i, j;
  double Sz;
  for (auto& alpha : spin_orbitals_) {
    for (i=0; i<basis_dim_; ++i) {
      std::tie(Sz,j) = basis_.apply_Sz(site_id,alpha,i);
      Sz_avg(alpha) += Sz * std::norm(groundstate_(i));
    }
  }
  //std::cout << "<Sz> = "<<Sz_avg.transpose() << "\n"; //getchar();
  // function values
  for (auto& alpha : spin_orbitals_) {
    func(alpha) = (Sz_avg(alpha)+0.5) - spinon_density_[site_][alpha];
    //std::cout << "func = " << func(alpha) << "\n";
  }
  //std::cout << "sp density = " << spinon_density_[site_].transpose() << "\n";
  //getchar();
  if (!need_derivative) return 0;

  // First compute derivative of the groundstate 
  // (See D.V. Murthy and R.T. Haftka, on "Eigenvector Derivative")

  // Derivative of H-matrix wrt 'lambda'-s: are diagonal matrices 
  for (i=0; i<basis_dim_; ++i) {
    for (auto& alpha : spin_orbitals_) {
      std::tie(Sz,j) = basis_.apply_Sz(site_id,alpha,i);
      // column 'alpha' stores the diagonal elements of dH/d\lambda_alpha
      H_dLambda_(i,alpha) = Sz + 0.5;
    }
  }
  //std::cout << "H_dLambda_=\n" << H_dLambda_ << "\n";
  //getchar();

  cmplVector Cvec(basis_dim_);
  cmplVector xvec(basis_dim_);
  for (auto& alpha : spin_orbitals_) {
    // dH * groundstate
    for (int i=0; i<basis_dim_; ++i) {
      xvec(i) = H_dLambda_(i,alpha) * groundstate_(i);
    }
    // coefficients c_{0j}
    for (int j=1; j<basis_dim_; ++j) {
      auto xc = eigen_solver_.eigenvectors().col(j).transpose().conjugate() 
              * xvec;
      Cvec(j) = xc(0,0)/(eigen_solver_.eigenvalues()(0)-eigen_solver_.eigenvalues()(j));
    }
    Cvec(0) = 0.0;
    cmplVector dX = cmplVector::Zero(basis_dim_);
    for (int i=1; i<basis_dim_; ++i) {
      dX += Cvec(i) * eigen_solver_.eigenvectors().col(i);
    }
    // column 'alpha' stores the groundstate derivative wrt \lambda_\alpha
    groundstate_dLambda_.col(alpha) = dX;
  }
  //std::cout << "dLambda=\n" << groundstate_dLambda_ << "\n";
  //getchar();

  // Jacobian
  //std::cout << "Cvec = "<< Cvec.transpose() << "\n";
  //std::cout << "DX =\n"<< groundstate_dLambda_ << "\n";
  // d<S^z_i>/d\lambda_j
  for (auto& alpha : spin_orbitals_) {
    for (auto& beta : spin_orbitals_) {
      double sum = 0.0;
      for (i=0; i<basis_dim_; ++i) {
        std::tie(Sz,j) = basis_.apply_Sz(site_id,alpha,i);
        sum += Sz*std::real(std::conj(groundstate_(i))*groundstate_dLambda_(i,beta));
      }
      Jmat(alpha, beta) = 2.0*sum;
      //std::cout << "J["<<alpha<<","<<beta<<"] = "<< Jmat(alpha, beta) << "\n";
    }
  }
  return 0;
}

int Cluster::rosenbrock_f(const RealVector& x, RealVector& fx, 
  RealMatrix& dfx, const bool& need_derivative) 
{
  double a = 1.0;
  double b = 10.0;
  fx(0) = a*(1.0-x(0));
  fx(1) = b*(x(1)-x(0)*x(0));
  if (!need_derivative) return 0;
  dfx(0,0) = -a;
  dfx(0,1) = 0.0;
  dfx(1,0) = -2.0*b*x(0);
  dfx(1,1) = b;
  return 0;
}

const ComplexMatrix& Cluster::groundstate_dLambda(void) 
{
  int site_id = 0;
  /*RealVector x0(2);
  root_solver_.solve([this](const RealVector& x, RealVector& fx, 
    RealMatrix& J, const bool& need_derivative) 
    {return lambda_equation(x, fx, J, need_derivative);}, x0);
  */
  root_solver_.init(2);
  RealVector x(2);
  x(0) = 10.0; x(1) = -50.0;
  root_solver_.solve([this](const RealVector& x, RealVector& fx, 
    RealMatrix& J, const bool& need_derivative) 
    {return rosenbrock_f(x,fx,J,need_derivative);}, x);
  std::cout << "sol = " << x.transpose() << "\n";
  std::cout << "---Testing rosenbrock_f----------\n";
  exit(0);

  // See D.V. Murthy and R.T. Haftka, on "Eigenvector Derivative".
  // Derivative of H-matrix wrt 'lambda'-s: are diagonal matrix 
  SlaveSpinBasis::idx_t i, j;
  double Sz;
  for (i=0; i<basis_dim_; ++i) {
    for (auto& alpha : spin_orbitals_) {
      std::tie(Sz,j) = basis_.apply_Sz(site_id,alpha,i);
      // column 'alpha' stores the diagonal elements of dH/d\lambda_alpha
      H_dLambda_(i,alpha) = Sz + 0.5;
    }
  }

  cmplVector Cvec(basis_dim_);
  cmplVector xvec(basis_dim_);
  for (auto& alpha : spin_orbitals_) {
    for (int j=1; j<basis_dim_; ++j) {
      // dH * groundstate
      for (int i=0; i<basis_dim_; ++i) {
        xvec(i) = H_dLambda_(i,alpha) * groundstate_(i);
      }
      auto xc = eigen_solver_.eigenvectors().col(j).transpose().conjugate() 
              * xvec;
      Cvec(j) = xc(0,0)/(eigen_solver_.eigenvalues()(0)-eigen_solver_.eigenvalues()(j));
    }
    Cvec(0) = 0.0;
    cmplVector dX = cmplVector::Zero(basis_dim_);
    for (int i=1; i<basis_dim_; ++i) {
      dX += Cvec(i) * eigen_solver_.eigenvectors().col(i);
    }
    // column 'alpha' stores the groundstate derivative wrt \lambda_\alpha
    groundstate_dLambda_.col(alpha) = dX;
  }

  //std::cout << "Cvec = "<< Cvec.transpose() << "\n";
  //std::cout << "DX =\n"<< groundstate_dLambda_ << "\n";

  // d<S^z_i>/d\lambda_j
  for (auto& alpha : spin_orbitals_) {
    for (auto& beta : spin_orbitals_) {
      double sum = 0.0;
      for (i=0; i<basis_dim_; ++i) {
        std::tie(Sz,j) = basis_.apply_Sz(site_id,alpha,i);
        sum += 2.0*Sz*std::real(std::conj(groundstate_(i))*groundstate_dLambda_(i,beta));
      }
      std::cout << "J["<<alpha<<","<<beta<<"] = "<< sum << "\n";
    }
  }
  std::cout << "\n\n";
  getchar();

  return groundstate_dLambda_;
}

void Cluster::get_avg_Sz(real_siteparms_t& Sz_avg) const
{
  int site_id = 0;
  for (auto& alpha : spin_orbitals_) Sz_avg[site_][alpha] = 0.0;
  SlaveSpinBasis::idx_t i, j;
  double Sz;
  for (i=0; i<basis_dim_; ++i) {
    for (auto& alpha : spin_orbitals_) {
      std::tie(Sz,j) = basis_.apply_Sz(site_id,alpha,i);
      Sz_avg[site_][alpha] += Sz * std::norm(groundstate_(i));
    }
  }
  //std::cout << "<Sz>["<<site_<<"] = "<< Sz_avg[site_].transpose() << "\n";
}


void Cluster::get_avg_Splus(real_siteparms_t& Splus_avg) const
{
  int site_id = 0;
  for (auto& alpha : spin_orbitals_) Splus_avg[site_][alpha] = 0.0;
  SlaveSpinBasis::idx_t i, j;
  double mat_elem;
  for (i=0; i<basis_dim_; ++i) {
    auto c_i = groundstate_(i);
    for (auto& alpha : spin_orbitals_) {
      std::tie(mat_elem,j) = basis_.apply_Splus(site_id,alpha,i);
      if (j != basis_.null_idx()) {
        auto c_j = groundstate_(j);
        Splus_avg[site_][alpha] += std::real(mat_elem*std::conj(c_j)*c_i);
      }
    }
  }
  std::cout << "<Splus>["<<site_<<"] = "<< Splus_avg[site_].transpose() << "\n";
  //std::cout << "groundstate = "<< groundstate_.transpose() << "\n";
}

void Cluster::get_avg_Zminus(const real_siteparms_t& gauge_factors, 
  cmpl_siteparms_t& Zplus_avg) const
{
  int site_id = 0;
  Zplus_avg[site_].setZero(); 
  SlaveSpinBasis::idx_t i, j;
  double mat_elem;
  for (auto& alpha : spin_orbitals_) {
    double c = gauge_factors[site_][alpha];
    for (i=0; i<basis_dim_; ++i) {
      auto c_i = groundstate_(i);
      std::tie(mat_elem,j) = basis_.apply_Zminus(c,site_id,alpha,i);
      if (j != basis_.null_idx()) {
        auto c_j = groundstate_(j);
        Zplus_avg[site_][alpha] += mat_elem*std::conj(c_j)*c_i;
        /*if (alpha==0) {
          std::cout << mat_elem << "\n";
          std::cout << std::conj(c_j)*c_i << "\n";
        }*/
      }
    }
  }
  //std::cout << "groundstate = "<< groundstate_ << "\n";
  //std::cout << "<Zplus>["<<site_<<"] = "<< Oplus_avg[site_].transpose() << "\n";
  //getchar();
}

void Cluster::get_avg_Ominus(const real_siteparms_t& gauge_factors, 
  cmpl_siteparms_t& Ominus_avg) const
{
  int site_id = 0;
  Ominus_avg[site_].setZero(); 
  SlaveSpinBasis::idx_t i, j;
  double mat_elem;
  for (auto& alpha : spin_orbitals_) {
    double c = gauge_factors[site_][alpha];
    for (i=0; i<basis_dim_; ++i) {
      auto c_i = groundstate_(i);
      std::tie(mat_elem,j) = basis_.apply_Ominus(c,site_id,alpha,i);
      if (j != basis_.null_idx()) {
        auto c_j = groundstate_(j);
        Ominus_avg[site_][alpha] += mat_elem*std::conj(c_j)*c_i;
      }
    }
  }
  //std::cout << "groundstate = "<< groundstate_ << "\n";
  //std::cout << "<Zplus>["<<site_<<"] = "<< Ominus_avg[site_].transpose() << "\n";
  //getchar();
}


void Cluster::get_avg_Oplus_Ominus(const real_siteparms_t& gauge_factors, 
  cmpl_siteparms_t& Opm_avg) const
{
  int site_id = 0;
  Opm_avg[site_].setZero(); 
  SlaveSpinBasis::idx_t i, j, k;
  double e1, e2;
  for (auto& alpha : spin_orbitals_) {
    double c = gauge_factors[site_][alpha];
    std::complex<double> sum = 0.0;
    for (i=0; i<basis_dim_; ++i) {
      auto c_i = groundstate_(i);
      std::tie(e1,j) = basis_.apply_Ominus(c,site_id,alpha,i);
      if (j != basis_.null_idx()) {
        std::tie(e2,k) = basis_.apply_Oplus(c,site_id,alpha,j);
        if (k != basis_.null_idx()) {
          auto c_k = groundstate_(k);
          sum += e1*e2*std::conj(c_k)*c_i;
          std::cout << "e1="<<e1<<" e2="<<e2<<"\n"; getchar();
        }
      }
    }
    Opm_avg[site_][alpha] = sum;
  }
  //std::cout << "<Oplus>["<<site_<<"] = "<< Oplus_avg[site_].transpose() << "\n";
  //std::cout << "groundstate = "<< groundstate_.transpose() << "\n";
  //getchar();
}


} // end namespace srmf
