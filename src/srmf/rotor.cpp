/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-19 11:24:03
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-09-24 16:00:59
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "rotor.h"
#include <stdexcept>
#include <string>
#include <cassert>
#include <boost/algorithm/string.hpp>
#include <boost/math/tools/roots.hpp>


namespace srmf {


Rotor::Rotor(const input::Parameters& inputs, const model::Hamiltonian& model, 
  const lattice::LatticeGraph& graph, const SR_Params& srparams)
  //: rotor_graph_(graph)
{
  std::cout << "----Rotor::Rotor-------\n";
  // Rotor lattice has only one original lattice unit cell 
  num_sites_ = srparams.num_sites();
  num_bonds_ = srparams.num_bonds();

  // Assuming all sites have same 'site_dim'.
  spin_orbitals_ = srparams.site(0).spin_orbitals(); // including 'UP' & 'DN' spins
  site_dim_ = spin_orbitals_.size(); 

  // boson model parameters  
  double U = model.get_parameter_value("U");
  U_half_ = 0.5 * U;

  // storages
  // Lagrange multipliers
  lm_params_.resize(num_sites_);
  for (auto& elem : lm_params_) elem = realArray1D::Ones(site_dim_);
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


  // Rotor clusters
  make_clusters(srparams);
  for (auto& cluster : clusters_) 
    cluster.init_hamiltonian(U,lm_params_,renorm_site_couplings_);

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

void Rotor::solve(SR_Params& srparams) 
{
  // set constrained density
  for (int i=0; i<num_sites_; ++i) {
    spinon_density_[i] = srparams.site(i).spinon_density();
  }

  // renormalized bond couplings
  set_bond_couplings(srparams);

  // trial 'qp_weights'
  real_siteparms_t trial_qp_weights(num_sites_);
  for (auto& elem : trial_qp_weights) elem = realArray1D::Constant(site_dim_,0.5);

  realArray1D qp_weights_diff(num_sites_*site_dim_);

  int max_iter = 100;
  for (int iter=0; iter<max_iter; ++iter) {
    // update 'renormalized site couplings' for the new 'qp_weights'
    set_site_couplings(srparams, trial_qp_weights);
    for (auto& cluster : clusters_) {
      cluster.update_hamiltonian(renorm_site_couplings_);
    }
    // solve for new LM-parameters to satisfy the 'slave-spin constraint'
    update_lm_params();
    // update cluster hamiltonians for new LM-parameters
    for (auto& cluster : clusters_) {
      cluster.update_hamiltonian(lm_params_);
    }
    // calculate new QP weights
    update_qp_weights();
    // check convergence
    int i = 0;
    for (unsigned site=0; site<num_sites_; ++site) {
      for (auto& alpha: spin_orbitals_) {
        qp_weights_diff(i) = qp_weights_[site][alpha]-trial_qp_weights[site][alpha];
        ++i; 
      }
    }
    double norm = qp_weights_diff.abs2().maxCoeff();
    std::cout << "iter = " << iter+1 << ", norm = " << norm << "\n";
    if (norm<1.0E-6) break;
    // continue
    for (unsigned site=0; site<num_sites_; ++site) {
      trial_qp_weights[site] = qp_weights_[site];
    }


    //update_with_phi(trial_phi_);
    /*
    double mu = solve_for_mu();
    update_with_mu(mu);
    solve_clusters();
    eval_site_phi();
    diff_phi_ = site_phi_ - trial_phi_;
    double norm = diff_phi_.abs2().maxCoeff();
    std::cout << "iter = " << iter+1 << ", norm = " << norm << "\n";
    if (norm<1.0E-8) break;
    trial_phi_ = site_phi_;
    */
  }
  // bond ke parameters
  //eval_bond_ke();
}

int gsl_problem_equation(const gsl_vector* x, void* parms, gsl_vector* f)
{
  Rotor * pThis = ((class Rotor *) parms);
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

int Rotor::constraint_equation(const std::vector<double>& x, std::vector<double>& fx)
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
  // calculate <Sz> for the new parameters
  for (auto& cluster : clusters_) {
      cluster.update_hamiltonian(lm_params_);
      cluster.solve_hamiltonian();
      cluster.get_avg_Sz(Sz_avg_);
  }
  // LHS of the constraint equation: fx = (<Sz> + 1/2) - n_f
  i=0;
  for (unsigned site=0; site<num_sites_; ++site) {
    for (auto& alpha: spin_orbitals_) {
      fx[i] = (0.5 + Sz_avg_[site][alpha]) - spinon_density_[site][alpha];
      ++i; 
      //std::cout << "fx = " << fx[i] << "\n";
    }
  }
  //getchar();
  return 0;
}

void Rotor::update_lm_params(void)
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
      std::cout << "lambda["<<i<<"] = " << x_vec_[i] << "\n";
      ++i; 
    }
  }
}

void Rotor::update_qp_weights(void)
{
  // calculate Z=<Sz> 
  for (auto& cluster : clusters_) {
      cluster.solve_hamiltonian();
      cluster.get_avg_Splus(qp_weights_);
  }
  for (unsigned site=0; site<num_sites_; ++site) {
    std::cout << "Z["<<site<<"] = " << qp_weights_[site].transpose() << "\n";
  }
}

void Rotor::set_bond_couplings(const SR_Params& srparams) 
{
  // Assuming spinorbitals includes both UP & DOWN spins
  int num_orb = static_cast<int>(site_dim_)/2;
  assert(2*num_orb == static_cast<int>(site_dim_));

  for (int i=0; i<num_bonds_; ++i) {
    renorm_bond_couplings_[i].setZero();
    auto bfield_spinUP = srparams.bond(i).spinon_ke() * srparams.bond(i).term_cc(0);
    // first diagonal block for spin-UP
    renorm_bond_couplings_[i].block(0,0,num_orb,num_orb) = bfield_spinUP;
    // second diagonal block for spin-DN
    // assuming spin-UP and DN symmetry
    renorm_bond_couplings_[i].block(num_orb,num_orb,num_orb,num_orb) = bfield_spinUP;
    //std::cout << "bond_field["<<i<<"]="<< renorm_bond_couplings_[i] << "\n";
  }

  /*
  // Renormalizing field on each site due to hopping from the 'bath'
  ArrayXcd tchi_field(site_dim_);
  if (cluster_type_==cluster_t::SITE) {
    for (int i=0; i<num_sites_; ++i) {
      ArrayXcd tchi_sum = ArrayXcd::Zero(num_orb);
      for (const auto& b: srparams.site(i).connected_bonds()) {
        // t * chi product for UP-spins
        auto tchi_UP = srparams.bond(b).spinon_ke() * srparams.bond(b).term_cc(0);
        // partial sum over 'orbital' index of the neighbour site
        auto tchi_orbital_sum = tchi_UP.rowwise().sum();
        // sum over all neighbouring sites
        if (srparams.site(i).is_outgoing_bond(b)) {
          tchi_sum += tchi_orbital_sum;
        }
        else {
          tchi_sum += tchi_orbital_sum.conjugate();
        }
      }
      // Assuming spin-UP and DN symmetry
      for (int j=0; j<num_orb; ++j) tchi_field(j) = tchi_sum(j); 
      for (int j=num_orb; j<site_dim_; ++j) tchi_field(j) = tchi_sum(j-num_orb); 
      renorm_site_couplings_[i] = tchi_field;
      //std::cout << "site field ["<<i<<"] = " << tchi_field.transpose() << "\n"; 
    }
  }
  else if (cluster_type_ == cluster_t::BOND) {
  }
  else if (cluster_type_ == cluster_t::CELL) {
  }
  else {
  }
  */
}

void Rotor::set_site_couplings(const SR_Params& srparams, 
  const real_siteparms_t& site_qp_weights) 
{
  // Renormalizing field on each site due to hopping terms from the 'bath'
  if (cluster_type_==cluster_t::SITE) {
    for (auto& elem : renorm_site_couplings_) elem.setZero();
    int id = 0;
    for (const auto& bond : srparams.bonds()) {
      auto s = bond.src();
      auto t = bond.tgt();

      // partial sum over 'orbital' index of the neighbour site
      auto tchi_sum = renorm_bond_couplings_[id++].rowwise().sum();

      // sum over all neighbouring sites
      renorm_site_couplings_[s] += tchi_sum * site_qp_weights[t];
      renorm_site_couplings_[t] += tchi_sum.conjugate() * site_qp_weights[s];
    }
    //for (int i=0; i<num_sites_; ++i) {
    //  std::cout << "site field ["<<i<<"] = " << renorm_site_couplings_[i].transpose() << "\n"; 
    //}
  }
  else if (cluster_type_ == cluster_t::BOND) {
  }
  else if (cluster_type_ == cluster_t::CELL) {
  }
  else {
  }
}

void Rotor::make_clusters(const SR_Params& srparams)
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

void Cluster::init_hamiltonian(const double& U, const real_siteparms_t& lagrange_fields,
    const cmpl_siteparms_t& site_fields)
{
  //std::cout << "Rotor:: applying Zplus instead of Splus\n";
  hmatrix_.setZero();
  // interaction matrix elements
  double U_half = 0.5 * U;
  double mat_elem;
  SlaveSpinBasis::idx_t i, j;
  for (i=0; i<basis_dim_; ++i) {
    double total_Sz = 0.0;
    for (auto& alpha : spin_orbitals_) {
      std::tie(mat_elem,j) = basis_.apply_Sz(site_,alpha,i);
      total_Sz += mat_elem;
    }
    interaction_elems_(i) = U_half * total_Sz  * total_Sz;
    //std::cout << "H_U=" << interaction_elems_(i)  << "\n";
  }

  // lagrange fields
  double Sz;
  for (i=0; i<basis_dim_; ++i) {
    double sum = 0.0;
    for (auto& alpha : spin_orbitals_) {
      std::tie(Sz,j) = basis_.apply_Sz(site_,alpha,i);
      sum += (Sz+0.5) * lagrange_fields[site_][alpha];
    }
    //std::cout << "lambda=" << sum << "\n";
    lagrange_elems_(i) = sum;
  }

  // site operator
  // operator S+
  for (i=0; i<basis_dim_; ++i) {
    for (auto& alpha : spin_orbitals_) {
      std::tie(mat_elem,j) = basis_.apply_Splus(site_,alpha,i);
      //std::tie(mat_elem,j) = basis_.apply_Zplus(site_,alpha,i);
      if (j != basis_.null_idx()) {
        auto term = mat_elem * site_fields[site_][alpha];
        hmatrix_(i,j) = term;
        hmatrix_(j,i) = std::conj(term);
      }
    }
  }

  // final matrix 
  for (i=0; i<basis_dim_; ++i) {
    hmatrix_(i,i) = interaction_elems_(i) + lagrange_elems_(i);
  }
  //std::cout << "cluster_mat =\n" << hmatrix_ << "\n";
}

// update for new 'site_couplings'
void Cluster::update_hamiltonian(const cmpl_siteparms_t& new_site_couplings)
{
  // site operator
  // operator S+
  double mat_elem;
  SlaveSpinBasis::idx_t i, j;
  for (i=0; i<basis_dim_; ++i) {
    for (auto& alpha : spin_orbitals_) {
      std::tie(mat_elem,j) = basis_.apply_Splus(site_,alpha,i);
      //std::tie(mat_elem,j) = basis_.apply_Zplus(site_,alpha,i);
      if (j != basis_.null_idx()) {
        auto term = mat_elem * new_site_couplings[site_][alpha];
        hmatrix_(i,j) = term;
        hmatrix_(j,i) = std::conj(term);
      }
    }
  }
  //std::cout<<"H (site-coupling-update) =\n" << hmatrix_ << "\n"; getchar();
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
      sum += (Sz + 0.5) * new_lm_params[site_][alpha];
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
      std::tie(mat_elem,j) = basis_.apply_Splus(site_,alpha,i);
      if (j != basis_.null_idx()) {
        auto c_j = groundstate_(j);
        Splus_avg[site_][alpha] += std::real(mat_elem*std::conj(c_j)*c_i);
      }
    }
  }
  std::cout << "<Zplus>["<<site_<<"] = "<< Splus_avg[site_].transpose() << "\n";
  //std::cout << "groundstate = "<< groundstate_.transpose() << "\n";
}


} // end namespace srmf
