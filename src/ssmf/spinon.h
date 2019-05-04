/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 14:51:12
* Last Modified by:   amedhi
* Last Modified time: 2017-06-11 16:49:03
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef SPINON_H
#define SPINON_H

#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <stdexcept>
#include "../scheduler/task.h"
#include "../lattice/graph.h"
#include "../lattice/matrix.h"
#include "../lattice/blochbasis.h"
#include "../model/quantum_op.h"
#include "../model/model.h"
//#include "sb_params.h"
#include "mf_params.h"
//#include "./blochbasis.h"

constexpr std::complex<double> ii(void) { return std::complex<double>{0.0,static_cast<double>(1.0)}; }

namespace srmf {

class UnitcellTerm
{
public:
  UnitcellTerm() {}
  ~UnitcellTerm() {}
  void build_bondterm(const model::HamiltonianTerm& sterm, const lattice::LatticeGraph& graph);
  void build_siteterm(const model::HamiltonianTerm& sterm, const lattice::LatticeGraph& graph);
  void eval_coupling_constant(const model::ModelParams& cvals, const model::ModelParams& pvals);
  void update_bondterm_cc(const int& term_id, const MF_Params& mf_params);
  void update_siteterm_cc(const MF_Params& mf_params);
  const unsigned& num_out_bonds(void) const { return num_out_bonds_; } 
  const Vector3d& bond_vector(const int& i) const { return bond_vectors_[i]; }
  const ComplexMatrix& coeff_matrix(const int& i=0) const { return coeff_matrices_[i]; }
  //const double& coupling(const unsigned& site_type) const; 
  const model::op::quantum_op& qn_operator(void) const { return op_; }
private:
  model::op::quantum_op op_;
  unsigned num_out_bonds_;
  unsigned num_basis_sites_;
  unsigned dim_;
  std::vector<ComplexMatrix> coeff_matrices_;
  std::vector<strMatrix> expr_matrices_;
  std::vector<Vector3d> bond_vectors_;
};

class Spinon : public model::Hamiltonian
{
public:
  Spinon() {}
  Spinon(const input::Parameters& inputs, const model::Hamiltonian& model, 
    const lattice::LatticeGraph& graph, const MF_Params& mf_params);
  ~Spinon() {}
  int init(const lattice::Lattice& lattice) override;
  int finalize(const lattice::LatticeGraph& graph);
  void solve(const lattice::LatticeGraph& graph, MF_Params& mf_params);
  void update(const input::Parameters& inputs);
  void update_terms(void) override;
  //const realArray1D& orbital_en(void) const { return orbital_en_; }
  void set_shifted_en(const std::vector<double>& shifted_e0) 
    { for (int i=0; i<kblock_dim_; ++i) orbital_en_shifted_[i]=shifted_e0[i]; }
  void update_site_parameter(const std::string& pname, const double& pvalue);
  void construct_kspace_block(const MF_Params& mf_params, const Vector3d& kvec);
  const ComplexMatrix& quadratic_spinup_block(void) const { return quadratic_block_up_; }
  const ComplexMatrix& pairing_part(void) const { return pairing_block_; }
  const realArray1D& orbital_en(void) const override { return orbital_en_; }
private:
  using Model = model::Hamiltonian;
  std::vector<UnitcellTerm> usite_terms_;
  std::vector<UnitcellTerm> ubond_terms_;
  basis::BlochBasis blochbasis_;
  bool have_TP_symmetry_{true};
  bool SO_coupling_{false};
  int spin_multiply_{2};
  int num_sites_{0};
  int num_unitcells_{0};
  int num_bonds_{0};
  int num_kpoints_{0};
  int num_symm_kpoints_{0};
  int kblock_dim_{0};
  int num_basis_sites_{0};
  // matrices in kspace representation
  ComplexMatrix quadratic_block_up_;
  ComplexMatrix quadratic_block_dn_;
  ComplexMatrix pairing_block_;
  realArray1D orbital_en_;
  realArray1D orbital_en_shifted_;
  ComplexMatrix work; //, work2;
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_k_up_;
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_k_dn_;

  // ground state (Fermi-Sea) representations
  int num_total_states_{0};
  int num_spins_{0};
  int num_upspins_{0};
  int num_dnspins_{0};
  bool metallic_{false};
  double hole_doping_{0.0};
  double last_hole_doping_{10.39}; // unlikely input
  double band_filling_{1.0};
  // fermi level
  int num_fill_particles_{0};
  int smear_func_order_{4};
  double fermi_energy_;
  double smear_width_;
  bool degeneracy_warning_{false};
  std::vector<std::pair<int,int>> qn_list_; // list of (k,n)
  std::vector<double> ek_list_;
  std::vector<int> idx_;
  struct kshell_t {int k; int nmin; int nmax; realArray1D smear_wt; };
  std::vector<kshell_t> kshells_up_;
  std::vector<kshell_t> kshells_dn_;

  void build_unitcell_terms(const lattice::LatticeGraph& graph);
  void update_unitcell_terms(void);
  void set_particle_num(const input::Parameters& inputs);
  void construct_groundstate(const MF_Params& mf_params);
  void construct_groundstate_v2(const MF_Params& mf_params);
  void compute_averages(const lattice::LatticeGraph& graph, MF_Params& mf_params);
  double MethfesselPaxton_func(const int& N, const double& x);
  double MarzariVenderbilt_smear(const double& x);
  double particle_constraint_eqn(const double& mu);
};


} // end namespace srmf

#endif