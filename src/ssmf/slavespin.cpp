/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-19 11:24:03
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-10-25 22:37:19
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "slavespin.h"
#include <stdexcept>
#include <string>
#include <cassert>
#include <boost/algorithm/string.hpp>
#include <boost/math/tools/roots.hpp>


namespace srmf {

SlaveSpin::SlaveSpin(const input::Parameters& inputs, const model::Hamiltonian& model, 
  const lattice::LatticeGraph& graph, const SB_Params& srparams)
  //: rotor_graph_(graph)
{
  // SlaveSpin lattice has only one original lattice unit cell 
  num_sites_ = srparams.num_sites();
  num_bonds_ = srparams.num_bonds();

  // Assuming all sites have same 'site_dim'.
  spin_orbitals_ = srparams.site(0).spin_orbitals(); // including 'UP' & 'DN' spins
  site_dim_ = spin_orbitals_.size(); 

  // boson model parameters  
  double U = model.get_parameter_value("U");

  // storages
  // Lagrange multipliers
  lm_params_.resize(num_sites_);
  for (auto& elem : lm_params_) elem = realArray1D::Ones(site_dim_);
  // gauge factors, c for the operator O+ = (S- + cS+)
  gauge_factors_.resize(num_sites_);
  for (auto& elem : gauge_factors_) elem = realArray1D::Ones(site_dim_);
  // order parameters: <O+>
  site_order_params_.resize(num_sites_);
  for (auto& elem : site_order_params_) elem = cmplArray1D::Ones(site_dim_);
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


  // SlaveSpin clusters
  make_clusters(srparams);
  for (auto& cluster : clusters_) 
    cluster.init_hamiltonian(U,gauge_factors_,lm_params_,renorm_site_couplings_);

  // For solving for LM parameter equation
  fx_dim_ = num_sites_*site_dim_;
  x_vec_.resize(fx_dim_);
  fx_vec_.resize(fx_dim_);
  solver_.allocate(fx_dim_);

  // 'slave boson' basis
  //----------Assuming 'SITE' cluster-------------
  /*
  sites_per_cluster_ = 1;
  ssbasis_.construct(sites_per_cluster_,site_dim_);
  basis_dim_ = ssbasis_.dim();
  make_clusters(srparams);
  cluster_hams_.resize(clusters_.size());
  for (auto& mat : cluster_hams_) {
    mat = ComplexMatrix::Zero(basis_dim_, basis_dim_);
  }
  */
  /*
  using namespace model;
  std::string name;
  double defval;
  CouplingConstant cc;
  rotor_model_.init(graph.lattice());
  rotor_model_.add_parameter(name="U", defval=U);
  rotor_model_.add_parameter(name="mu", defval=0.0);
  //--Operators are directly implemeted in cluster hamiltonian
  //rotor_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
  //rotor_model_.add_siteterm(name="e^i_theta", cc="1", op::ciup_dag());
  //rotor_model_.add_siteterm(name="hubbard", cc="U", op::hubbard_int());
  // finalize model
  rotor_model_.finalize(graph.lattice());
  */
}

void SlaveSpin::update(const model::Hamiltonian& model)
{
  double U = model.get_parameter_value("U");
  for (auto& cluster : clusters_) {
    cluster.update_parameters(U);
  }
}

void SlaveSpin::solve(SB_Params& srparams) 
{
  // set constrained density
  //realArray1D tmp(4);
  //tmp(0)=0.121894; tmp(1)=0.87810; tmp(2)=0.121894; tmp(3)=0.87810;
  for (int i=0; i<num_sites_; ++i) {
    spinon_density_[i] = srparams.site(i).spinon_density();
    //spinon_density_[i] = tmp;
  }
  // gauge factors
  for (int site=0; site<num_sites_; ++site) {
    for (auto alpha : spin_orbitals_) {
      double nf = spinon_density_[site][alpha];
      gauge_factors_[site][alpha] = 1.0/std::sqrt(nf*(1.0-nf))-1.0;
      //std::cout << "nf["<<site<<"]["<<alpha<<"] = "<<nf<<"\n";
      //std::cout << "c["<<site<<"]["<<alpha<<"] = "<<gauge_factors_[site][alpha]<<"\n";
    }
  }

  // renormalized bond couplings
  set_bond_couplings(srparams);

  // trial 'qp_weights'
  cmpl_siteparms_t trial_order_params(num_sites_);
  for (auto& elem : trial_order_params) elem = cmplArray1D::Constant(site_dim_,1.0);

  cmplArray1D order_params_diff(num_sites_*site_dim_);

  int max_iter = 100;
  bool converged = false;
  bool print_progress = false;
  for (int iter=0; iter<max_iter; ++iter) {
    // update 'renormalized site couplings' for the new 'order parameters'
    set_site_couplings(srparams, trial_order_params);
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
    }
    // calculate new site op-s
    update_site_order_params();
    // check convergence
    int i = 0;
    for (unsigned site=0; site<num_sites_; ++site) {
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
    for (unsigned site=0; site<num_sites_; ++site) {
      trial_order_params[site] = site_order_params_[site];
    }
  }
  if (converged) {
    if(print_progress) std::cout<<"Bosons converged!\n";
  } 
  // QP weights
  for (int site=0; site<num_sites_; ++site) {
    qp_weights_[site] = site_order_params_[site].abs2();
    //std::cout<<"Z["<<site<<"] = "<< qp_weights_[site].transpose()<<"\n";
  }

  // update boson parameters
  for (int i=0; i<num_sites_; ++i) {
    srparams.site(i).lm_params() = lm_params_[i];
    srparams.site(i).qp_weights() = qp_weights_[i];
    //std::cout<<"lambda["<<i<<"] = "<< lm_params_[i].transpose()<<"\n";
    //std::cout<<"Z["<<i<<"] = "<< qp_weights_[i].transpose()<<"\n";
  }
  update_bond_order_params(srparams);
}

void SlaveSpin::update_bond_order_params(SB_Params& srparams) 
{
  // Renormalizing field on each site due to hopping terms from the 'bath'
  cmplArray2D bond_op(site_dim_,site_dim_);
  if (cluster_type_==cluster_t::SITE) {
    for (int i=0; i<srparams.num_bonds(); ++i) {
      auto s = srparams.bond(i).src();
      auto t = srparams.bond(i).tgt();
      for (int m=0; m<site_dim_; ++m) {
        for (int n=0; n<site_dim_; ++n) {
          bond_op(m,n) = site_order_params_[s](m) * std::conj(site_order_params_[t](n));
        }
      }
      srparams.bond(i).boson_ke() = bond_op;
      //std::cout << "bond_op["<<i<<"] =\n" << bond_op << "\n";
    }
  }
  else if (cluster_type_ == cluster_t::BOND) {
  }
  else if (cluster_type_ == cluster_t::CELL) {
  }
  else {
  }
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
  for (unsigned site=0; site<num_sites_; ++site) {
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
      cluster.groundstate_dLambda(); 
      // groundstate derivative
      //std::cout << "dX/dLambda =\n"<< cluster.groundstate_dLambda() << "\n";
      //getchar();
  }
  // LHS of the constraint equation: fx = (<Sz> + 1/2) - n_f
  i=0;
  for (unsigned site=0; site<num_sites_; ++site) {
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
  //solver_.set_problem(&lagrange_eqn, 2, this);
  //std::vector<double> x0(num_sites_*site_dim_, 1.0); 
  for (auto& x : x_vec_) x = 1.0;
  solver_.find_root(this, &gsl_problem_equation, x_vec_, lm_ftol_);
  // solved 'LM' parameters
  int i = 0;
  for (unsigned site=0; site<num_sites_; ++site) {
    for (auto& alpha: spin_orbitals_) {
      lm_params_[site][alpha] = x_vec_[i];
      //std::cout << "lambda["<<i<<"] = " << x_vec_[i] << "\n";
      ++i; 
    }
  }
  //getchar();
}

void SlaveSpin::update_site_order_params(void)
{
  // calculate <O+> 
  for (auto& cluster : clusters_) {
      cluster.solve_hamiltonian();
      cluster.get_avg_Oplus(gauge_factors_,site_order_params_);
  }
  /*for (unsigned site=0; site<num_sites_; ++site) {
    std::cout << "<Oplus>["<<site<<"] = " << site_order_params_[site].transpose() << "\n";
  }*/
}

void SlaveSpin::set_bond_couplings(const SB_Params& srparams) 
{
  // Assuming spinorbitals includes both UP & DOWN spins
  int num_orb = static_cast<int>(site_dim_)/2;
  assert(2*num_orb == static_cast<int>(site_dim_));

  for (int i=0; i<num_bonds_; ++i) {
    /*renorm_bond_couplings_[i].setZero();
    auto bfield_spinUP = srparams.bond(i).spinon_ke() * srparams.bond(i).term_cc(0);
    // first diagonal block for spin-UP
    renorm_bond_couplings_[i].block(0,0,num_orb,num_orb) = bfield_spinUP;
    // second diagonal block for spin-DN
    // assuming spin-UP and DN symmetry
    renorm_bond_couplings_[i].block(num_orb,num_orb,num_orb,num_orb) = bfield_spinUP;
    */
    renorm_bond_couplings_[i] = srparams.bond(i).spinon_renormed_cc(0);
    //std::cout << "bond_field["<<i<<"]=\n"<< renorm_bond_couplings_[i] << "\n"; getchar();
  }
}

void SlaveSpin::set_site_couplings(const SB_Params& srparams, 
  const cmpl_siteparms_t& site_order_params) 
{
  // Renormalizing field on each site due to hopping terms from the 'bath'
  if (cluster_type_==cluster_t::SITE) {
    for (auto& elem : renorm_site_couplings_) elem.setZero();
    int id = 0;
    for (const auto& bond : srparams.bonds()) {
      auto s = bond.src();
      auto t = bond.tgt();

      // partial sum over 'orbital' index of the neighbour site
      //auto tchi_sum = renorm_bond_couplings_[id++].rowwise().sum();
      auto tchi_phi = renorm_bond_couplings_[id++].rowwise()*site_order_params[s].transpose();
      auto tchi_phi_sum = tchi_phi.rowwise().sum();
      /*std::cout << "chi=\n"<<renorm_bond_couplings_[id-1]<<"\n";
      std::cout << "phi="<<site_order_params[s].transpose()<<"\n";
      std::cout << "tchi=\n"<<tchi_phi<<"\n";
      std::cout << "tchi_sum=\n"<<tchi_phi_sum.transpose()<<"\n\n\n";
      getchar();
      */

      // sum over all neighbouring sites
      renorm_site_couplings_[s] += tchi_phi_sum; // * site_order_params[t];
      renorm_site_couplings_[t] += tchi_phi_sum.conjugate(); // * site_order_params[s];
    }
    /*for (int i=0; i<num_sites_; ++i) {
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

void SlaveSpin::make_clusters(const SB_Params& srparams)
{
  cluster_type_ = cluster_t::SITE;
  bonds_ = srparams.bonds();
  //site_links_ = srparams.site_links();
  clusters_.clear();
  switch (cluster_type_) {
    case cluster_t::SITE:
      for (unsigned i=0; i<num_sites_; ++i) {
        clusters_.push_back({cluster_t::SITE,i,srparams.site(i).spin_orbitals()});
      }
      break;
    case cluster_t::BOND:
      break;
    case cluster_t::CELL:
      break;
  }
}

void Cluster::init_hamiltonian(const double& U, const real_siteparms_t& gauge_factors,
  const real_siteparms_t& lm_params, const cmpl_siteparms_t& site_fields)
{
  //std::cout << "SlaveSpin:: applying Zplus instead of Splus\n";
  hmatrix_.setZero();
  // interaction matrix elements
  double U_half = 0.5*U;
  double Sz;
  SlaveSpinBasis::idx_t i, j;
  for (i=0; i<basis_dim_; ++i) {
    double total_Sz = 0.0;
    for (auto& alpha : spin_orbitals_) {
      std::tie(Sz,j) = basis_.apply_Sz(site_,alpha,i);
      total_Sz += Sz;
    }
    interaction_elems_(i) = U_half*total_Sz*total_Sz;
    //std::cout << "H_U=" << interaction_elems_(i)  << "\n";
  }

  // LM parameter term
  for (i=0; i<basis_dim_; ++i) {
    double sum = 0.0;
    for (auto& alpha : spin_orbitals_) {
      std::tie(Sz,j) = basis_.apply_Sz(site_,alpha,i);
      sum += (Sz+0.5) * lm_params[site_][alpha];
    }
    //std::cout << "lambda=" << sum << "\n";
    lagrange_elems_(i) = sum;
  }

  // site operator
  // operator O+
  double mat_elem;
  for (j=0; j<basis_dim_; ++j) {
    for (auto& alpha : spin_orbitals_) {
      double c = gauge_factors[site_][alpha];
      std::tie(mat_elem,i) = basis_.apply_Oplus(c,site_,alpha,j);
      if (i != basis_.null_idx()) {
        auto term = mat_elem * site_fields[site_][alpha];
        hmatrix_(i,j) += term;
        hmatrix_(j,i) += std::conj(term);
        //std::cout << "H["<<i<<","<<"j]=" << term << "\n"; getchar();
      }
    }
  }

  // final matrix 
  for (i=0; i<basis_dim_; ++i) {
    hmatrix_(i,i) = interaction_elems_(i) + lagrange_elems_(i);
  }
  //std::cout << "cluster_mat =\n" << hmatrix_ << "\n";
}

void Cluster::update_parameters(const double& U)
{
  // interaction matrix elements
  double U_half = 0.5*U;
  double mat_elem;
  SlaveSpinBasis::idx_t i, j;
  for (i=0; i<basis_dim_; ++i) {
    double total_Sz = 0.0;
    for (auto& alpha : spin_orbitals_) {
      std::tie(mat_elem,j) = basis_.apply_Sz(site_,alpha,i);
      total_Sz += mat_elem;
    }
    interaction_elems_(i) = U_half*total_Sz*total_Sz;
    //std::cout << "total_Sz=" << total_Sz  << "\n";
    //std::cout << "H_U=" << interaction_elems_(i)  << "\n";
  }
  //getchar();
}

// update for new 'site_couplings'
void Cluster::update_hamiltonian(const real_siteparms_t& gauge_factors,
  const cmpl_siteparms_t& new_site_couplings)
{
  //std::cout << "rotor::Cluster::update_hamiltonian: ---Problem-------\n";
  SlaveSpinBasis::idx_t i, j;
  double mat_elem;
  hmatrix_.setZero();
  for (j=0; j<basis_dim_; ++j) {
    for (auto& alpha : spin_orbitals_) {
      double c = gauge_factors[site_][alpha];
      std::tie(mat_elem,i) = basis_.apply_Oplus(c,site_,alpha,j);
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
  // diagonal elements
  for (i=0; i<basis_dim_; ++i) {
    hmatrix_(i,i) = interaction_elems_(i) + lagrange_elems_(i);
  }
  //std::cout<<"H (site-coupling-update) =\n" << hmatrix_ << "\n"; getchar();
  //getchar();
}

// update for new 'lm parameters'
void Cluster::update_hamiltonian(const real_siteparms_t& new_lm_params)
{
  // lagrange fields
  SlaveSpinBasis::idx_t i, j;
  double Sz;
  for (i=0; i<basis_dim_; ++i) {
    double sum = 0.0;
    for (auto& alpha : spin_orbitals_) {
      std::tie(Sz,j) = basis_.apply_Sz(site_,alpha,i);
      sum += (Sz+0.5) * new_lm_params[site_][alpha];
    }
    //std::cout << "lambda=" << sum << "\n";
    lagrange_elems_(i) = sum;
  }
  for (i=0; i<basis_dim_; ++i) {
    hmatrix_(i,i) = interaction_elems_(i) + lagrange_elems_(i);
  }
  //std::cout<<"H (LM-update) =\n" << hmatrix_ << "\n"; getchar();
}

void Cluster::solve_hamiltonian(void) const
{
  eigen_solver_.compute(hmatrix_);
  groundstate_ = eigen_solver_.eigenvectors().col(0);
  //std::cout << "Eigenvalues = " << eigen_solver_.eigenvalues().transpose() << "\n";
  //std::cout << "Eigenvector = " << groundstate_.transpose() << "\n";
  //getchar();
}

int Cluster::lambda_equation(const RealVector& lambda, RealVector& f_lambda, 
  RealMatrix& df_lambda, const bool& need_derivative) 
{
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
  /*RealVector x0(2);
  root_solver_.solve([this](const RealVector& x, RealVector& fx, 
    RealMatrix& J, const bool& need_derivative) 
    {return lambda_equation(x, fx, J, need_derivative);}, x0);
  */
  root_solver_.init(2);
  RealVector x0(2);
  x0(0) = -10.0; x0(1) = -5.0;
  root_solver_.solve([this](const RealVector& x, RealVector& fx, 
    RealMatrix& J, const bool& need_derivative) 
    {return rosenbrock_f(x,fx,J,need_derivative);}, x0);
  std::cout << "---Testing rosenbrock_f----------\n";
  exit(0);



  // See D.V. Murthy and R.T. Haftka, on "Eigenvector Derivative".

  // Derivative of H-matrix wrt 'lambda'-s: are diagonal matrix 
  SlaveSpinBasis::idx_t i, j;
  double Sz;
  for (i=0; i<basis_dim_; ++i) {
    for (auto& alpha : spin_orbitals_) {
      std::tie(Sz,j) = basis_.apply_Sz(site_,alpha,i);
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
        std::tie(Sz,j) = basis_.apply_Sz(site_,alpha,i);
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
  for (auto& alpha : spin_orbitals_) Sz_avg[site_][alpha] = 0.0;
  SlaveSpinBasis::idx_t i, j;
  double Sz;
  for (i=0; i<basis_dim_; ++i) {
    for (auto& alpha : spin_orbitals_) {
      std::tie(Sz,j) = basis_.apply_Sz(site_,alpha,i);
      Sz_avg[site_][alpha] += Sz * std::norm(groundstate_(i));
    }
  }
  //std::cout << "<Sz>["<<site_<<"] = "<< Sz_avg[site_].transpose() << "\n";
}


void Cluster::get_avg_Splus(real_siteparms_t& Splus_avg) const
{
  for (auto& alpha : spin_orbitals_) Splus_avg[site_][alpha] = 0.0;
  SlaveSpinBasis::idx_t i, j;
  double mat_elem;
  for (i=0; i<basis_dim_; ++i) {
    auto c_i = groundstate_(i);
    for (auto& alpha : spin_orbitals_) {
      std::tie(mat_elem,j) = basis_.apply_Oplus(1.0,site_,alpha,i);
      if (j != basis_.null_idx()) {
        auto c_j = groundstate_(j);
        Splus_avg[site_][alpha] += std::real(mat_elem*std::conj(c_j)*c_i);
      }
    }
  }
  std::cout << "<Oplus>["<<site_<<"] = "<< Splus_avg[site_].transpose() << "\n";
  //std::cout << "groundstate = "<< groundstate_.transpose() << "\n";
}

void Cluster::get_avg_Oplus(const real_siteparms_t& gauge_factors, 
  cmpl_siteparms_t& Oplus_avg) const
{
  Oplus_avg[site_].setZero(); 
  SlaveSpinBasis::idx_t i, j;
  double mat_elem;
  for (i=0; i<basis_dim_; ++i) {
    auto c_i = groundstate_(i);
    for (auto& alpha : spin_orbitals_) {
      double c = gauge_factors[site_][alpha];
      std::tie(mat_elem,j) = basis_.apply_Ominus(c,site_,alpha,i);
      if (j != basis_.null_idx()) {
        auto c_j = groundstate_(j);
        Oplus_avg[site_][alpha] += mat_elem*std::conj(c_j)*c_i;
        //std::cout << "mat_elem = "<< mat_elem << "\n";
        //std::cout << "cj*ci["<<alpha<<"] = "<< std::conj(c_j)*c_i << "\n";
      }
    }
  }
  //std::cout << "<Oplus>["<<site_<<"] = "<< Oplus_avg[site_].transpose() << "\n";
  //std::cout << "groundstate = "<< groundstate_.transpose() << "\n";
  //getchar();
}


} // end namespace srmf
