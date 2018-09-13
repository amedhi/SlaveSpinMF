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
#include <Eigen/Sparse>
#include "../scheduler/task.h"
#include "../lattice/graph.h"
#include "../lattice/matrix.h"
#include "../model/quantum_op.h"
#include "../model/model.h"
#include "srparams.h"
#include "rbasis_states.h"
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
  }
  ~Cluster() {}
  void init_hamiltonian(const double& U, const std::vector<ArrayXd>& lagrange_fields,
    const std::vector<ArrayXcd>& site_fields);
  void update_siteop_elems(const std::vector<ArrayXcd>& site_fields);
  void update_lagrange_elems(const std::vector<ArrayXd>& lagrange_fields);
private:
  cluster_t type_{cluster_t::SITE};
  unsigned site_{0};
  unsigned basis_dim_{0};
  std::vector<unsigned> spin_orbitals_; 
  SlaveSpinBasis basis_;

  RealVector interaction_elems_; 
  RealVector lagrange_elems_; 
  ComplexMatrix hmatrix_;
};

/*
template <class T=double> 
struct f_of_mu 
{
  f_of_mu(const T& mu): mu_(mu) {};
private:
  T mu_;
};*/

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
  double avg_density(void) const { return site_density_.sum()/num_sites_; }
  //void update(const input::Parameters& inputs);
private:
  //using LatticeGraph = lattice::LatticeGraph;
  using Model = model::Hamiltonian;
  //using bond = sr_bond;
  //using links = SR_Params::links;
  //LatticeGraph rotor_graph_;
  Model rotor_model_;
  unsigned num_sites_;
  unsigned num_bonds_;
  std::vector<sr_site> sites_; 
  std::vector<sr_bond> bonds_; 

  // slave spin hamiltonian
  SlaveSpinBasis ssbasis_;
  unsigned site_dim_;
  unsigned basis_dim_;
  unsigned dim_;
  std::vector<unsigned> spin_orbitals_; // for a single site

  // clusters
  cluster_t cluster_type_;
  //std::vector<links> site_links_; 
  std::vector<Cluster> cluster_list_;
  std::vector<std::vector<unsigned> > clusters_;
  unsigned sites_per_cluster_{1};
  // site & bond parameters
  double U_half_{0.0};
  double constrained_density_{0.0};
  double avg_density_{0.0};
  ArrayXd site_density_;
  std::vector<double> site_mu_;
  ArrayXcd site_phi_;
  ArrayXcd site_mfp_;
  std::vector<ArrayXd> site_lagrange_fields_;
  std::vector<ComplexArray1D> site_qp_weights_;
  std::vector<ComplexArray> renorm_bond_fields_;
  std::vector<ArrayXcd> renorm_site_fields_;
  ArrayXcd bond_tchi_;
  ArrayXcd bond_ke_;


  ArrayXcd trial_phi_;
  ArrayXcd diff_phi_;

  // hamiltonian matrix
  rotor_basis basis_;
  std::vector<ComplexMatrix> cluster_hams_;
  std::vector<ComplexMatrix> cluster_hams_siteop_;
  std::vector<ArrayXd> cluster_hams_lag_;

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
  void set_renormalized_bond_fields(const SR_Params& srparams);
  void set_renormalized_site_fields(const SR_Params& srparams, 
    const std::vector<ComplexArray1D>& site_qp_weights);
  void set_site_fields(void);
  void solve_lagrange_fields(void);
  double solve_for_mu(void);
  void solve_clusters(void);
  void eval_particle_density(void);
  void eval_site_phi(void);
  void eval_bond_ke(void);
  void construct_cluster_hams(void);
  void update_with_mu(const double& new_mu);
  void update_with_phi(const ArrayXcd& new_phi);
  void solve_number_density(double mu);
  double avg_particle_density_eqn(const double& mu); 
};




} // end namespace srmf

#endif
