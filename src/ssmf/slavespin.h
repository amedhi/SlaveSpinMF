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
//#include "sb_params.h"
#include "mf_params.h"
#include "boson_basis.h"
#include "root_solver.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
//#include "./blochbasis.h"

namespace ssmf {

enum cluster_t {SITE, BOND, CELL};
enum theory_t {Z2, U1};
enum model_id {HUBBARD, HUBBARD_NBAND, BHZ, PYROCHLORE, HUBBARD_SQIRIDATE};

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
  void set_id(const model_id& id) { id_=id; }
  void update_e0(const realArray1D& e0) { e0_=e0; }
  void update_U(const double& U) { U_=U; }
  void update_U1(const double& U) { U1_=U; hunde_coupling_=true; }
  void update_J(const double& J) { hunde_J_=J; hunde_coupling_=true; }
  void update_lambda(const double& lambda) { SO_lambda_=lambda; }
  void set_J_relative(const bool& yesno) { J_is_relative_=yesno; }
  void set_spinflip_terms(const bool& yesno ) { have_spinflip_terms_=yesno; }
  const model_id& id(void) const { return id_; }
  const realArray1D& get_e0(void) const { return e0_; }
  const double& get_U(void) const { return U_; }
  const double& get_U1(void) const { return U1_; }
  const double& get_J(void) const { return hunde_J_; }
  const double& get_lambda(void) const { return SO_lambda_; }
  bool have_hunde_coupling(void) const { return hunde_coupling_; }
  bool J_is_relative(void) const { return J_is_relative_; }
  bool have_spinflip_terms(void) const { return have_spinflip_terms_; }
private:
  model_id id_;
  bool hunde_coupling_{false};
  bool J_is_relative_{false};
  bool have_spinflip_terms_{true};
  double U_{0.0};
  double U1_{0.0}; // inter-orbital
  double hunde_J_{0.0};
  double SO_lambda_{0.0};
  realArray1D e0_; // orbital energy
};

class Cluster 
{
public:
  Cluster() {}
  Cluster(const cluster_t& type, const theory_t& theory, 
    const int& site, const std::vector<int>& spin_orbitals) 
  : type_{type}, theory_{theory}, site_{site}, spin_orbitals_{spin_orbitals}
  {
    end_site_id_ = site_+1;
    total_spinorbitals_ = spin_orbitals_.size();
    basis_.construct(1, spin_orbitals_.size());
    basis_dim_ = basis_.dim();
    //std::cout << total_spinorbitals_ << " " << basis_dim_ << "\n";
    lm_params_.resize(end_site_id_);
    for (auto& elem : lm_params_) elem.resize(spin_orbitals_.size());
    gauge_factors_.resize(end_site_id_);
    for (auto& elem : gauge_factors_) elem.resize(spin_orbitals_.size());
    site_couplings_.resize(end_site_id_);
    for (auto& elem : site_couplings_) {
      elem.resize(spin_orbitals_.size());
      elem = cmplArray1D::Ones(spin_orbitals_.size());
    }

    spinon_density_.resize(end_site_id_);
    for (auto& elem : spinon_density_) elem.resize(spin_orbitals_.size());
    interaction_elems_.resize(basis_dim_);
    lagrange_elems_.resize(basis_dim_);
    //orbital_en_elems_.resize(basis_dim_);
    soc_couplings_.resize(end_site_id_);
    for (auto& elem : soc_couplings_) {
      elem = cmplArray2D::Zero(spin_orbitals_.size(),spin_orbitals_.size());
    }
    soc_mat_.resize(basis_dim_,basis_dim_);
    hmatrix_.resize(basis_dim_, basis_dim_);
    interaction_mat_.resize(basis_dim_,basis_dim_);
    groundstate_.resize(basis_dim_);
    H_dLambda_.resize(basis_dim_,total_spinorbitals_);
    groundstate_dLambda_.resize(basis_dim_,total_spinorbitals_);
    root_solver_.init(total_spinorbitals_);
  }
  ~Cluster() {}
  void reset_mfp(void) {
    for (auto& elem : lm_params_) elem.setZero();
    for (auto& elem : gauge_factors_) elem.setConstant(1.0);
    for (auto& elem : site_couplings_) {
      elem = cmplArray1D::Ones(spin_orbitals_.size());
    }
  }
  void init_hamiltonian(const ModelParams& p, const real_siteparms_t& gauge_factors, 
    const real_siteparms_t& lagrange_fields, const cmpl_siteparms_t& renorm_site_fields);
  void set_spinon_density(const real_siteparms_t& spinon_density);
  void set_site_couplings(const cmpl_siteparms_t& site_couplings);
  void set_soc_couplings(const MF_Params& mfp);
  void solve_lm_params(real_siteparms_t& lm_params);
  void update_soc_matrix(const real_siteparms_t& gauge_factors);
  void update_interaction_matrix(const ModelParams& p);
  void update_hamiltonian(const real_siteparms_t& new_lm_params);
  void update_hamiltonian(const real_siteparms_t& gauge_factors, const cmpl_siteparms_t& new_site_couplings);
  //void get_groundstate(ComplexVector& eigvec) const;
  void solve_hamiltonian(void) const;
  const ComplexMatrix& hamiltonian_matrix(void) const { return hmatrix_; }
  const ComplexVector& groundstate(void) const { return groundstate_; }
  const ComplexMatrix& groundstate_dLambda(void); 
  double get_interaction_energy(const ModelParams& p);
  void get_avg_Sz(real_siteparms_t& Sz_avg) const;
  void get_avg_Splus(real_siteparms_t& Splus_avg) const;
  void get_avg_Zminus(const real_siteparms_t& gauge_factors, cmpl_siteparms_t& order_params) const;
  void get_avg_Ominus(const real_siteparms_t& gauge_factors, cmpl_siteparms_t& order_params) const;
  void get_avg_OplusOminus(const real_siteparms_t& gauge_factors, cmplArray2D& Opm_avg) const;
  void get_avg_ZplusZminus(const real_siteparms_t& gauge_factors, cmplArray2D& Zpm_avg) const;
  friend int gsl_lambda_equation(const gsl_vector* x, void* parms, gsl_vector* f);
private:
  cluster_t type_{cluster_t::SITE};
  theory_t theory_{theory_t::Z2};
  int site_{0};
  int end_site_id_{0};
  bool real_amplitudes_{true};
  int basis_dim_{0};
  int total_spinorbitals_{0};
  std::vector<int> spin_orbitals_; 
  SlaveSpinBasis basis_;
  root::RootSolver root_solver_;

  real_siteparms_t lm_params_;
  real_siteparms_t gauge_factors_;
  cmpl_siteparms_t site_couplings_;
  real_siteparms_t spinon_density_;
  RealVector interaction_elems_; 
  RealVector lagrange_elems_; 
  //RealVector orbital_en_elems_; 

  cmpl_bondparms_t soc_couplings_; 
  ComplexMatrix soc_mat_;
  ComplexMatrix hmatrix_;
  ComplexMatrix interaction_mat_;
  ComplexMatrix H_dLambda_;
  ComplexMatrix groundstate_dLambda_;
  mutable ComplexVector groundstate_;
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> eigen_solver_;

  int lambda_equation(const RealVector& lambda, RealVector& f_lambda, 
    RealMatrix& df_lambda, const bool& need_derivative);
  int rosenbrock_f(const RealVector& x, RealVector& fx, 
    RealMatrix& dfx, const bool& need_derivative);
  //realArray1D gauge_and_lambda_func(const MF_Params& mf_params, const int& i)
  int gauge_and_lambda_func(const std::vector<double>& x, std::vector<double>& fx); 
};


class SlaveSpin 
{
public:
  //SlaveSpin() {}
  SlaveSpin(const input::Parameters& inputs, const model::Hamiltonian& model, 
    const lattice::LatticeGraph& graph, const MF_Params& mf_params);
  ~SlaveSpin() {}
  const std::string& info_str(void) const { return info_str_; }
  void update(const model::Hamiltonian& model);
  void solve(MF_Params& mf_params);
  double interaction_energy(void);
  //int init(const lattice::Lattice& lattice) override;
  //int finalize(const lattice::LatticeGraph& graph);
  //void update(const input::Parameters& inputs);
  friend int gsl_problem_equation(const gsl_vector* x, void* parms, gsl_vector* f);
private:
  theory_t theory_{theory_t::Z2};
  bool solve_single_site_{false};
  bool gauge_factors_set_{false};
  bool gauge_factors_solved_{false};
  bool set_fixed_gauge_{true};
  bool SO_coupling_{false};
  bool use_previous_iter_{false};
  double fixed_gauge_{1.0};
  std::string info_str_;
  //using LatticeGraph = lattice::LatticeGraph;
  using Model = model::Hamiltonian;
  //Model rotor_model_;
  int num_sites_;
  int num_bonds_;
  //std::vector<sb_site> sites_; 
  //std::vector<sb_bond> bonds_; 
  ModelParams modelparams_;
  double delta_{1.0E-4};

  // clusters
  int site_dim_;
  std::vector<int> spin_orbitals_; // for a single site
  cluster_t cluster_type_;
  std::vector<Cluster> clusters_;

  int max_iter_{500};
  double conv_tol_{1.0E-8};

  // gsl solver
  double gsl_ftol_{1.0E-12};
  unsigned fx_dim_;
  std::vector<double> x_vec_;
  std::vector<double> fx_vec_;
  root::gsl_solver gsl_solver_;
  int solve_cluster_{0};

  // site & bond parameters
  real_siteparms_t lm_params_;
  real_siteparms_t lm_params_noint_;
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
  void set_info_string(void); 
  void solve_gauge_factors(const MF_Params& mf_params);
  realArray1D gauge_factors_func(const MF_Params& mf_params, const int& i);
  void self_consistent_solve(const MF_Params& mf_params);
  void make_clusters(const MF_Params& mf_params);
  void init_matrix_elems(const MF_Params& mf_params);
  void set_bond_couplings(const MF_Params& mf_params);
  void set_lm_params_noint(const MF_Params& mf_params);
  void set_noninteracting_params(const MF_Params& mf_params);
  void set_site_couplings(const MF_Params& mf_params, 
    const cmpl_siteparms_t& site_order_params);
  void set_site_fields(void);
  void update_lm_params(void);
  void get_order_params(cmpl_siteparms_t& order_params);
  void set_renormalized_soc(MF_Params& mf_params);
  void update_bond_order_params(MF_Params& mf_params);
  void update_renorm_site_potential(MF_Params& mf_params);
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




} // end namespace ssmf

#endif
