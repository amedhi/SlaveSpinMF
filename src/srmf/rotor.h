/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-19 11:24:03
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-04-19 11:26:10
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef ROTOR_H
#define ROTOR_H

#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <tuple>
#include <stdexcept>
//#include <Eigen/Sparse>
#include "../scheduler/task.h"
#include "../lattice/graph.h"
#include "../lattice/matrix.h"
//#include "../model/quantum_op.h"
#include "../model/model.h"
#include "srparams.h"
#include "rbasis_states.h"
#include "root_solver.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
//#include "./blochbasis.h"

namespace srmf {

enum cluster_t {SITE, BOND, CELL};

class MatrixElem
{
public:
  MatrixElem() {}
  MatrixElem(const int& row, const int& col, const double& val) 
    : row_{row}, col_{col}, value_{val} {}
  ~MatrixElem() {}
  void change_value(const double& val) { value_=val; }
  const int& row(void) const { return row_; }
  const int& col(void) const { return col_; }
  const double& value(void) const { return value_; }
private:
  int row_{0};
  int col_{0};
  double value_{0.0};
};

class Cluster 
{
public:
  Cluster() {}
  Cluster(const cluster_t& type, const unsigned& site, 
    const std::vector<unsigned>& spin_orbitals) 
  : type_{type}, site_{site}, spin_orbitals_{spin_orbitals}
  {
    basis_.construct(1, spin_orbitals_.size());
    basis_dim_ = basis_.dim();
    interaction_elems_.resize(basis_dim_);
    lagrange_elems_.resize(basis_dim_);
    hmatrix_.resize(basis_dim_, basis_dim_);
    groundstate_.resize(basis_dim_);
  }
  ~Cluster() {}
  void init_hamiltonian(const double& U, const real_siteparms_t& gauge_factors, 
    const real_siteparms_t& lagrange_fields, const cmpl_siteparms_t& site_fields);
  void update_hamiltonian(const real_siteparms_t& new_lm_params);
  void update_hamiltonian(const real_siteparms_t& gauge_factors, const cmpl_siteparms_t& new_site_couplings);
  //void get_groundstate(ComplexVector& eigvec) const;
  void solve_hamiltonian(void) const;
  const ComplexVector& groundstate(void) const { return groundstate_; }
  void get_avg_Sz(real_siteparms_t& Sz_avg) const;
  void get_avg_Splus(real_siteparms_t& Splus_avg) const;
  void get_avg_Oplus(const real_siteparms_t& gauge_factors, cmpl_siteparms_t& order_params) const;
private:
  cluster_t type_{cluster_t::SITE};
  unsigned site_{0};
  unsigned basis_dim_{0};
  std::vector<unsigned> spin_orbitals_; 
  SlaveSpinBasis basis_;

  RealVector interaction_elems_; 
  RealVector lagrange_elems_; 
  ComplexMatrix hmatrix_;
  mutable ComplexVector groundstate_;
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> eigen_solver_;
};


class Rotor 
{
public:
  //Rotor() {}
  Rotor(const input::Parameters& inputs, const model::Hamiltonian& model, 
    const lattice::LatticeGraph& graph, const SR_Params& srparams);
  ~Rotor() {}
  void solve(SR_Params& srparams);
  //int init(const lattice::Lattice& lattice) override;
  //int finalize(const lattice::LatticeGraph& graph);
  //void update(const input::Parameters& inputs);
  friend int gsl_problem_equation(const gsl_vector* x, void* parms, gsl_vector* f);
private:
  //using LatticeGraph = lattice::LatticeGraph;
  using Model = model::Hamiltonian;
  //Model rotor_model_;
  unsigned num_sites_;
  unsigned num_bonds_;
  std::vector<sr_site> sites_; 
  std::vector<sr_bond> bonds_; 

  // clusters
  unsigned site_dim_;
  std::vector<unsigned> spin_orbitals_; // for a single site
  cluster_t cluster_type_;
  std::vector<Cluster> clusters_;

  // gsl solver
  double lm_ftol_{1.0E-8};
  unsigned fx_dim_;
  std::vector<double> x_vec_;
  std::vector<double> fx_vec_;
  root::gsl_solver solver_;

  // site & bond parameters
  double U_half_{0.0};
  real_siteparms_t lm_params_;
  real_siteparms_t qp_weights_;
  real_siteparms_t gauge_factors_;
  cmpl_siteparms_t site_order_params_;
  cmpl_siteparms_t renorm_site_couplings_;
  cmpl_bondparms_t renorm_bond_couplings_;
  real_siteparms_t spinon_density_;

  real_siteparms_t Sz_avg_;

  ComplexMatrix cluster_groundstate_;
  std::vector<MatrixElem> interaction_elems_; 
  std::vector<MatrixElem> diagonal_elems_; 
  std::vector<MatrixElem> siteop_elems_; 
  std::vector<MatrixElem> bondop_elems_; 
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> eigen_solver_;
  /*
  Eigen::SparseMatrix<double> diag_terms_;
  Eigen::SparseMatrix<std::complex<double>> offdiag_siteterms_;
  Eigen::SparseMatrix<std::complex<double>> offdiag_bondterms_;
  //Eigen::SparseMatrix<double> work_;
  //ComplexMatrix psi_work2_;
  */
  void make_clusters(const SR_Params& srparams);
  void init_matrix_elems(const SR_Params& srparams);
  void set_bond_couplings(const SR_Params& srparams);
  void set_site_couplings(const SR_Params& srparams, 
    const cmpl_siteparms_t& site_order_params);
  void set_site_fields(void);
  void update_lm_params(void);
  void update_order_params(void);
  int constraint_equation(const std::vector<double>& x, std::vector<double>& fx);
  double solve_for_mu(void);
  void solve_clusters(void);
  void eval_particle_density(void);
  void eval_site_phi(void);
  void eval_bond_ke(void);
  void construct_cluster_hams(void);
  void update_with_mu(const double& new_mu);
  void update_with_phi(const cmplArray1D& new_phi);
  void solve_number_density(double mu);
  double avg_particle_density_eqn(const double& mu); 
};




} // end namespace srmf

#endif
