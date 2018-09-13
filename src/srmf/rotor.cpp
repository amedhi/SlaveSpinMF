/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-19 11:24:03
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-09-13 15:55:00
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "rotor.h"
#include <stdexcept>
#include <string>
#include <cassert>
#include <boost/algorithm/string.hpp>
#include <boost/math/tools/roots.hpp>

namespace srmf {


/*
Cluster::Cluster(const cluster_t& type, const int& id, const rotor_basis& basis)
  : type_{type}, id_{id}, dim_{basis.dim()}
{
  switch (type_) {
    case cluster_t::SITE: num_sites_ = 1; break;
    case cluster_t::BOND: num_sites_ = 2; break;
    case cluster_t::CELL: num_sites_ = 2; break;
  }
  ham_.resize(dim_,dim_);
}
void Cluster::solve()
{
}
*/


Rotor::Rotor(const input::Parameters& inputs, const model::Hamiltonian& model, 
  const lattice::LatticeGraph& graph, const SR_Params& srparams)
  //: rotor_graph_(graph)
{
  std::cout << "----Rotor::Rotor-------\n";
  // Rotor lattice has only one original lattice unit cell 
  num_sites_ = srparams.num_sites();
  num_bonds_ = srparams.num_bonds();
  //sites_ = srparams.sites();
  //bonds_ = srparams.bonds();

  // slave spin basis
  //----------Assuming all sites have same 'site_dim'-------------
  spin_orbitals_ = srparams.site(0).spin_orbitals(); // including 'UP' & 'DN' spins
  site_dim_ = spin_orbitals_.size(); 

  // 'slave boson' basis
  //----------Assuming 'SITE' cluster-------------
  sites_per_cluster_ = 1;
  ssbasis_.construct(sites_per_cluster_,site_dim_);
  basis_dim_ = ssbasis_.dim();

  make_clusters(srparams);
  cluster_hams_.resize(clusters_.size());
  for (auto& mat : cluster_hams_) {
    mat = ComplexMatrix::Zero(basis_dim_, basis_dim_);
  }

  // Rotor model 
  using namespace model;
  std::string name;
  double defval;
  CouplingConstant cc;
  double U = model.get_parameter_value("U");
  U_half_ = 0.5 * U;
  rotor_model_.init(graph.lattice());
  rotor_model_.add_parameter(name="U", defval=U);
  rotor_model_.add_parameter(name="mu", defval=0.0);
  //--Operators are directly implemeted in "init_matrix_matrix"
  /*rotor_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
  rotor_model_.add_siteterm(name="e^i_theta", cc="1", op::ciup_dag());
  rotor_model_.add_siteterm(name="hubbard", cc="U", op::hubbard_int());
  */
  // finalize model
  rotor_model_.finalize(graph.lattice());

  // storages
  renorm_site_fields_.resize(num_sites_);
  for (auto& elem : renorm_site_fields_) elem = ComplexArray1D::Ones(site_dim_);
  site_lagrange_fields_.resize(num_sites_);
  for (auto& elem : site_lagrange_fields_) elem = ArrayXd::Ones(site_dim_);

  site_qp_weights_.resize(num_sites_);
  for (auto& elem : site_qp_weights_) elem.resize(site_dim_);
  renorm_bond_fields_.resize(num_bonds_);
  for (auto& elem : renorm_bond_fields_) 
    elem.resize(site_dim_,site_dim_);


  bond_tchi_.resize(num_bonds_);
  bond_ke_.resize(num_bonds_);
  site_density_.resize(num_sites_);
  site_mu_.resize(num_sites_);
  site_phi_.resize(num_sites_);
  trial_phi_.resize(num_sites_);
  diff_phi_.resize(num_sites_);
  site_mfp_.resize(num_sites_);
  cluster_groundstate_.resize(dim_, clusters_.size());
  for (int i=0; i<num_sites_; ++i) {
    site_mu_[i] = 0.0;
    site_phi_[i] = 1.0;
  }
  init_matrix_elems(srparams);
  //construct_cluster_hams();
  for (auto& cluster : cluster_list_) 
    cluster.init_hamiltonian(U,site_lagrange_fields_,renorm_site_fields_);
}

void Rotor::solve(SR_Params& srparams) 
{
  // set constrained density
  for (int i=0; i<num_sites_; ++i) site_density_(i) = 1.0-srparams.spinon_density(i);
  constrained_density_ = site_density_.sum()/num_sites_;
  std::cout << "cons_density = " << constrained_density_ << "\n";

  // renormalizing parameters from spinon sector
  set_renormalized_bond_fields(srparams);

  std::vector<ComplexArray1D> trial_qp_weights(num_sites_);
  for (auto& elem : trial_qp_weights) elem = ComplexArray1D::Ones(site_dim_);

  //set_renormalized_site_fields(srparams, trial_qp_weights);

  int max_iter = 1;
  for (int iter=0; iter<max_iter; ++iter) {
    set_renormalized_site_fields(srparams, trial_qp_weights);

    //for (auto& cluster : cluster_list_) {
    //   cluster.update_siteop_elems(renorm_site_fields_);
    //}

    //solve_lagrange_fields()
    //update_chemical_potential();

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

void Rotor::solve_lagrange_fields(void)
{
}

void Rotor::set_renormalized_bond_fields(const SR_Params& srparams) 
{
  // Assuming spinorbitals includes both UP & DOWN spins
  int num_orb = static_cast<int>(site_dim_)/2;
  assert(2*num_orb == static_cast<int>(site_dim_));

  for (int i=0; i<num_bonds_; ++i) {
    renorm_bond_fields_[i].setZero();
    auto bfield_spinUP = srparams.bond(i).spinon_ke() * srparams.bond(i).term_cc(0);
    // first diagonal block for spin-UP
    renorm_bond_fields_[i].block(0,0,num_orb,num_orb) = bfield_spinUP;
    // second diagonal block for spin-DN
    // assuming spin-UP and DN symmetry
    renorm_bond_fields_[i].block(num_orb,num_orb,num_orb,num_orb) = bfield_spinUP;
    //std::cout << "bond_field["<<i<<"]="<< renorm_bond_fields_[i] << "\n";
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
      renorm_site_fields_[i] = tchi_field;
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

void Rotor::set_renormalized_site_fields(const SR_Params& srparams, 
  const std::vector<ComplexArray1D>& site_qp_weights) 
{
  // Renormalizing field on each site due to hopping terms from the 'bath'
  if (cluster_type_==cluster_t::SITE) {
    for (auto& elem : renorm_site_fields_) elem.setZero();
    int id = 0;
    for (const auto& bond : srparams.bonds()) {
      auto s = bond.src();
      auto t = bond.tgt();

      // partial sum over 'orbital' index of the neighbour site
      auto tchi_sum = renorm_bond_fields_[id++].rowwise().sum();

      // sum over all neighbouring sites
      renorm_site_fields_[s] += tchi_sum * site_qp_weights[t];
      renorm_site_fields_[t] += tchi_sum.conjugate() * site_qp_weights[s];
    }
    //for (int i=0; i<num_sites_; ++i) {
    //  std::cout << "site field ["<<i<<"] = " << renorm_site_fields_[i].transpose() << "\n"; 
    //}
  }
  else if (cluster_type_ == cluster_t::BOND) {
  }
  else if (cluster_type_ == cluster_t::CELL) {
  }
  else {
  }
}

double Rotor::solve_for_mu(void) 
{
  double guess = 0.0;
  double factor = 2.0;
  const boost::uintmax_t maxit = 20; 
  boost::uintmax_t it = maxit;      
  bool is_rising = true;
  boost::math::tools::eps_tolerance<double> tol(3);
  std::pair<double,double> r = boost::math::tools::bracket_and_solve_root(
    [this](double mu) { return avg_particle_density_eqn(mu); },
    guess, factor, is_rising, tol, it);
  return r.first + 0.5 * (r.second-r.first);
} 

double Rotor::avg_particle_density_eqn(const double& mu) 
{
  // set chemical potential
  update_with_mu(mu);
  solve_clusters();
  eval_particle_density(); 
  return avg_density()-constrained_density_;
}

void Rotor::make_clusters(const SR_Params& srparams)
{
  bonds_ = srparams.bonds();
  //site_links_ = srparams.site_links();
  clusters_.clear();
  cluster_list_.clear();
  switch (cluster_type_) {
    case cluster_t::SITE:
      sites_per_cluster_ = 1;
      clusters_.resize(num_sites_);
      for (unsigned i=0; i<num_sites_; ++i) clusters_[i].push_back(i);

      for (unsigned i=0; i<num_sites_; ++i) {
        cluster_list_.push_back({cluster_t::SITE,i,srparams.site(i).spin_orbitals()});
      }



      break;
    case cluster_t::BOND:
      sites_per_cluster_ = 2;
      clusters_.resize(bonds_.size());
      for (int i=0; i<bonds_.size(); ++i) clusters_[i].push_back(i);
      break;
    case cluster_t::CELL:
      sites_per_cluster_ = num_sites_;
      clusters_.resize(1);
      for (int i=0; i<bonds_.size(); ++i) clusters_[0].push_back(i);
      break;
  }
}

void Rotor::solve_clusters(void)
{
  int i = 0;
  for (auto& mat : cluster_hams_) {
    // solve the hamiltonian
    eigen_solver_.compute(mat);
    //std::cout << "gndstate = " << eigen_solver_.eigenvectors().col(0) << "\n";
    //getchar();
    // store ground state
    cluster_groundstate_.col(i) = eigen_solver_.eigenvectors().col(0);
    i++;
  }
}

void Rotor::update_with_mu(const double& new_mu) 
{
  for (auto& mu : site_mu_) mu = new_mu;
  if (cluster_type_==cluster_t::SITE) {
    int i = 0;
    for (auto& mat : cluster_hams_) {
      double mu = site_mu_[i];
      for (const auto& elem: diagonal_elems_) {
        int n = elem.row();
        double ntheta = elem.value();
        mat(n,n) = (U_half_*ntheta-mu)*ntheta;
      }
      i++;
    }
  }
  else if (cluster_type_==cluster_t::BOND) {
    throw std::runtime_error("Rotor::update_with_mu: undefined 'cluster_type'");
  }
  else if (cluster_type_==cluster_t::CELL) {
    throw std::runtime_error("Rotor::update_with_mu: undefined 'cluster_type'");
  }
  else {
    throw std::runtime_error("Rotor::update_with_mu: undefined 'cluster_type'");
  }
}

void Rotor::update_with_phi(const ArrayXcd& new_phi)
{
  site_phi_ = new_phi;
  if (cluster_type_==cluster_t::SITE) {
    int i = 0;
    for (auto& mat : cluster_hams_) {
      auto matrix_elem = site_phi_[i] * site_mfp_[i];
      for (const auto& elem: siteop_elems_) {
        int m = elem.row();
        int n = elem.col();
        mat(m,n) = -matrix_elem;
        mat(n,m) = -std::conj(matrix_elem);
      }
      i++;
    }
  }
  else if (cluster_type_==cluster_t::BOND) {
    throw std::runtime_error("Rotor::update_with_phi: undefined 'cluster_type'");
  }
  else if (cluster_type_==cluster_t::CELL) {
    throw std::runtime_error("Rotor::update_with_phi: undefined 'cluster_type'");
  }
  else {
    throw std::runtime_error("Rotor::update_with_phi: undefined 'cluster_type'");
  }
}

void Rotor::eval_particle_density(void)
{
  double ntheta, sum;
  switch (cluster_type_) {
    case cluster_t::SITE:
      for (int i=0; i<clusters_.size(); ++i) {
        sum = 0.0;
        for (auto n=0; n<dim_; ++n) {
          ntheta = basis_.apply_ni(n);
          sum += ntheta * std::norm(cluster_groundstate_(n,i));
        }
        site_density_[i] = sum;
        //std::cout << "site density = " << i << ": " << sum << "\n";
      }
      break;
    case cluster_t::BOND:
      break;
    case cluster_t::CELL:
      break;
  }
}

void Rotor::eval_site_phi(void)
{
  rotor_basis::idx_t j; // i
  std::complex<double> sum;
  switch (cluster_type_) {
    case cluster_t::SITE:
      for (int n=0; n<clusters_.size(); ++n) {
        sum = 0.0;
        for (j=0; j<dim_-1; ++j) {
          auto c_j = cluster_groundstate_(j,n);
          /*i = basis_.apply_cidag(j);
          if (i != basis_.null_idx()) {
            auto c_i = cluster_groundstate_(i,n);
            sum += std::conj(c_i)*c_j;
          }*/
          auto c_i = cluster_groundstate_(j+1,n);
          sum += std::conj(c_i)*c_j;
        }
        site_phi_[n] = sum;
        std::cout << "site phi = " << n << ": " << sum << "\n";
      }
      break;
    case cluster_t::BOND:
      break;
    case cluster_t::CELL:
      break;
  }
}

void Rotor::eval_bond_ke(void)
{
  if (cluster_type_==cluster_t::SITE) {
    int n = 0;
    for (auto& bond : bonds_) {
      auto phi_i = site_phi_[bond.src()];
      auto phi_j = site_phi_[bond.tgt()];
      bond_ke_[n++] = std::conj(phi_i) * phi_j;
    }
  }
  else if (cluster_type_ == cluster_t::BOND) {
  }
  else if (cluster_type_ == cluster_t::CELL) {
  }
  else {
  }
}

void Rotor::construct_cluster_hams(void)
{
  if (cluster_type_==cluster_t::SITE) {
    int s = 0;
    for (auto& mat : cluster_hams_) {
      // interaction operator
      for (const auto& elem: interaction_elems_) {
        int n = elem.row();
        double total_Sz = elem.value();
        mat(n,n) = U_half_ * total_Sz  * total_Sz;
        //std::cout << "mat(n,n) = " << mat(n,n) << "\n";
      }

      // lagrange fields
      //for (i=0; i<basis_dim_; ++i) {
      //  double total_Sz = 0.0;
      //  for (auto& alpha : spin_orbitals_) {
      //    std::tie(mat_elem,j) = ssbasis_.apply_Sz(s,alpha,i);
      //    total_Sz += mat_elem;
      //  }
      //  //std::cout << "total_Sz=" << total_Sz << "\n";
      //  interaction_elems_.push_back({i,i,total_Sz});
      //}


      // site operators
      auto coupling = site_phi_[s] * site_mfp_[s];
      for (const auto& elem: siteop_elems_) {
        int m = elem.row();
        int n = elem.col();
        mat(m,n) = -coupling * elem.value();
        mat(n,m) = std::conj(mat(m,n));
      }
      s++;
      // check
      // std::cout << "ham = " << mat << "\n";
    }
  }
  else if (cluster_type_==cluster_t::BOND) {
    throw std::runtime_error("Rotor::construct_cluster_hams: undefined 'cluster_type'");
  }
  else if (cluster_type_==cluster_t::CELL) {
    throw std::runtime_error("Rotor::construct_cluster_hams: undefined 'cluster_type'");
  }
  else {
    throw std::runtime_error("Rotor::construct_cluster_hams: undefined 'cluster_type'");
  }
}

void Rotor::init_matrix_elems(const SR_Params& srparams) 
{
  // non-zero matrix elements of various operators in the Hamiltonian
  SlaveSpinBasis::idx_t i, j;
  double mat_elem;
  unsigned site = 0; 
  // interaction matrix elements
  interaction_elems_.clear();
  for (i=0; i<basis_dim_; ++i) {
    double total_Sz = 0.0;
    for (auto& alpha : spin_orbitals_) {
      std::tie(mat_elem,j) = ssbasis_.apply_Sz(site,alpha,i);
      total_Sz += mat_elem;
    }
    //std::cout << "total_Sz=" << total_Sz << "\n";
    interaction_elems_.push_back({i,i,total_Sz});
  }

  siteop_elems_.clear();
  // operator S+
  for (i=0; i<basis_dim_; ++i) {
    for (auto& alpha : spin_orbitals_) {
      std::tie(mat_elem,j) = ssbasis_.apply_Splus(site,alpha,i);
      if (j != ssbasis_.null_idx()) {
        siteop_elems_.push_back({i,j,mat_elem});
        //std::cout << "siteop elem ("<<i<<","<<j<<")="<<mat_elem<<"\n";
      }
    }
  }

/*
  // non-zero matrix elements of various operators in the Hamiltonian
  rotor_basis::idx_t i, j;
  int mat_elem;
  unsigned site = 0;
  diagonal_elems_.clear();
  for (i=0; i<dim_; ++i) {
    std::tie(mat_elem, j) = basis_.apply_ni(i, site);
    diagonal_elems_.push_back({i,i,mat_elem});
  }
  //for (auto& elem : siteop_elems_) std::cout << elem.row() << " " << elem.col() << "\n";
  bondop_elems_.clear();
  */
}

void Cluster::init_hamiltonian(const double& U, const std::vector<ArrayXd>& lagrange_fields,
    const std::vector<ArrayXcd>& site_fields)
{
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
    //std::cout << "total_Sz=" << total_Sz << "\n";
    interaction_elems_(i) = U_half * total_Sz  * total_Sz;
  }

  // lagrange fields
  double Sz;
  for (i=0; i<basis_dim_; ++i) {
    double sum = 0.0;
    for (auto& alpha : spin_orbitals_) {
      std::tie(Sz,j) = basis_.apply_Sz(site_,alpha,i);
      sum += Sz * lagrange_fields[site_][alpha];
    }
    //std::cout << "lambda=" << sum << "\n";
    lagrange_elems_(i) = sum;
  }

  // site operator
  // operator S+
  for (i=0; i<basis_dim_; ++i) {
    for (auto& alpha : spin_orbitals_) {
      std::tie(mat_elem,j) = basis_.apply_Splus(site_,alpha,i);
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
  std::cout << "cluster_mat = " << hmatrix_ << "\n";
}


} // end namespace srmf
