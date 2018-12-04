/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-19 11:24:03
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-04-19 11:26:10
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef SLAVESPIN_H
#define SLAVESPIN_H

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
#include "sb_params.h"
#include "boson_basis.h"
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

class ModelParams
{
public:
  ModelParams() {}
  ~ModelParams() {}
  void update_e0(const realArray1D& e0) { e0_=e0; }
  void update_U(const double& U) { U_=U; }
  void update_U1(const double& U) { U1_=U; hunde_coupling_=true; }
  void update_hunde_J(const double& J) { hunde_J_=J; hunde_coupling_=true; }
  const realArray1D& get_e0(void) const { return e0_; }
  const double& get_U(void) const { return U_; }
  const double& get_U1(void) const { return U1_; }
  const double& get_hunde_J(void) const { return hunde_J_; }
  bool have_hunde_coupling(void) const { return hunde_coupling_; }
private:
  bool hunde_coupling_{false};
  double U_{0.0};
  double U1_{0.0}; // inter-orbital
  double hunde_J_{0.0};
  realArray1D e0_; // orbital energy
};

class Cluster 
{
public:
  Cluster() {}
  Cluster(const cluster_t& type, const unsigned& site, 
    const std::vector<unsigned>& spin_orbitals) 
  : type_{type}, site_{site}, spin_orbitals_{spin_orbitals}
  {
    total_spinorbitals_ = spin_orbitals_.size();
    basis_.construct(1, spin_orbitals_.size());
    basis_dim_ = basis_.dim();
    lm_params_.resize(1);
    lm_params_[0].resize(spin_orbitals_.size());
    spinon_density_.resize(1);
    spinon_density_[0].resize(spin_orbitals_.size());
    interaction_elems_.resize(basis_dim_);
    lagrange_elems_.resize(basis_dim_);
    //orbital_en_elems_.resize(basis_dim_);
    hmatrix_.resize(basis_dim_, basis_dim_);
    groundstate_.resize(basis_dim_);
    H_dLambda_.resize(basis_dim_,total_spinorbitals_);
    groundstate_dLambda_.resize(basis_dim_,total_spinorbitals_);
    root_solver_.init(total_spinorbitals_);
  }
  ~Cluster() {}
  void init_hamiltonian(const ModelParams& p, const real_siteparms_t& gauge_factors, 
    const real_siteparms_t& lagrange_fields, const cmpl_siteparms_t& renorm_site_fields);
  void set_spinon_density(const real_siteparms_t& spinon_density);
  void solve_lm_params(real_siteparms_t& lm_params);
  void update_hamiltonian(const ModelParams& p);
  void update_hamiltonian(const real_siteparms_t& new_lm_params);
  void update_hamiltonian(const real_siteparms_t& gauge_factors, const cmpl_siteparms_t& new_site_couplings);
  //void get_groundstate(ComplexVector& eigvec) const;
  void solve_hamiltonian(void) const;
  const ComplexVector& groundstate(void) const { return groundstate_; }
  const ComplexMatrix& groundstate_dLambda(void); 
  void get_avg_Sz(real_siteparms_t& Sz_avg) const;
  void get_avg_Splus(real_siteparms_t& Splus_avg) const;
  void get_avg_Zplus(const real_siteparms_t& gauge_factors, cmpl_siteparms_t& order_params) const;
  void get_avg_Oplus_Ominus(const real_siteparms_t& gauge_factors, cmpl_siteparms_t& Opm_avg) const;
private:
  cluster_t type_{cluster_t::SITE};
  unsigned site_{0};
  unsigned basis_dim_{0};
  unsigned total_spinorbitals_{0};
  std::vector<unsigned> spin_orbitals_; 
  SlaveSpinBasis basis_;
  root::RootSolver root_solver_;

  real_siteparms_t lm_params_;
  real_siteparms_t spinon_density_;
  RealVector interaction_elems_; 
  RealVector lagrange_elems_; 
  //RealVector orbital_en_elems_; 
  //std::vector<MatrixElem> orbital_en_elems_;
  ComplexMatrix hmatrix_;
  ComplexMatrix H_dLambda_;
  ComplexMatrix groundstate_dLambda_;
  mutable ComplexVector groundstate_;
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> eigen_solver_;

  int lambda_equation(const RealVector& lambda, RealVector& f_lambda, 
    RealMatrix& df_lambda, const bool& need_derivative);
  int rosenbrock_f(const RealVector& x, RealVector& fx, 
    RealMatrix& dfx, const bool& need_derivative);
};


class SlaveSpin 
{
public:
  //SlaveSpin() {}
  SlaveSpin(const input::Parameters& inputs, const model::Hamiltonian& model, 
    const lattice::LatticeGraph& graph, const SB_Params& srparams);
  ~SlaveSpin() {}
  void update(const model::Hamiltonian& model);
  void solve(SB_Params& srparams);
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
  std::vector<sb_site> sites_; 
  std::vector<sb_bond> bonds_; 
  ModelParams modelparams_;
  double delta_{1.0E-4};

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
  real_siteparms_t lm_params_;
  //real_siteparms_t lm_params_noint_;
  real_siteparms_t qp_weights_;
  real_siteparms_t gauge_factors_;
  real_siteparms_t spinon_density_;
  cmpl_siteparms_t site_order_params_;
  cmpl_siteparms_t site_avg_OplusMinus_;
  cmpl_siteparms_t renorm_site_couplings_;
  cmpl_bondparms_t renorm_bond_couplings_;

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
  void self_consistent_solve(const SB_Params& srparams);
  void make_clusters(const SB_Params& srparams);
  void init_matrix_elems(const SB_Params& srparams);
  void set_bond_couplings(const SB_Params& srparams);
  void set_site_couplings(const SB_Params& srparams, 
    const cmpl_siteparms_t& site_order_params);
  void set_site_fields(void);
  void update_lm_params(void);
  void update_site_order_params(void);
  void update_bond_order_params(SB_Params& srparams);
  void update_renorm_site_potential(SB_Params& srparams);
  int constraint_equation(const std::vector<double>& x, std::vector<double>& fx);
  double solve_for_mu(void);
  void solve_clusters(void);
  void eval_particle_density(void);
  void eval_site_phi(void);
  void construct_cluster_hams(void);
  void update_with_mu(const double& new_mu);
  void update_with_phi(const cmplArray1D& new_phi);
  void solve_number_density(double mu);
  double avg_particle_density_eqn(const double& mu); 
};




} // end namespace srmf

#endif
