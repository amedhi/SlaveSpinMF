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
#include "srparams.h"
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
  const unsigned& num_out_bonds(void) const { return num_out_bonds_; } 
  const Vector3d& bond_vector(const unsigned& i) const { return bond_vectors_[i]; }
  const ComplexMatrix& coeff_matrix(const unsigned& i=0) const { return coeff_matrices_[i]; }
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
    const lattice::LatticeGraph& graph, const SR_Params& srparams);
  ~Spinon() {}
  int init(const lattice::Lattice& lattice) override;
  int finalize(const lattice::LatticeGraph& graph);
  void solve(const lattice::LatticeGraph& graph, SR_Params& srparams);
  void update(const input::Parameters& inputs);
  void update_terms(void) override;
  void update_site_parameter(const std::string& pname, const double& pvalue);
  void construct_kspace_block(const Vector3d& kvec);
  const ComplexMatrix& quadratic_spinup_block(void) const { return quadratic_block_up_; }
  const ComplexMatrix& pairing_part(void) const { return pairing_block_; }
private:
  using Model = model::Hamiltonian;
  std::vector<UnitcellTerm> usite_terms_;
  std::vector<UnitcellTerm> ubond_terms_;
  basis::BlochBasis blochbasis_;
  // matrices in kspace representation
  unsigned num_sites_{0};
  unsigned num_bonds_{0};
  unsigned num_kpoints_{0};
  unsigned kblock_dim_{0};
  unsigned num_basis_sites_{0};
  //unsigned dim_{0};
  bool have_TP_symmetry_{true};
  ComplexMatrix quadratic_block_up_;
  ComplexMatrix quadratic_block_dn_;
  ComplexMatrix pairing_block_;
  ComplexMatrix work; //, work2;
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_k_up_;
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_k_dn_;

  // ground state (Fermi-Sea) representations
  unsigned num_total_states_{0};
  unsigned num_spins_{0};
  unsigned num_upspins_{0};
  unsigned num_dnspins_{0};
  double hole_doping_{0.0};
  double last_hole_doping_{10.39}; // unlikely input
  double band_filling_{1.0};
  double fermi_energy_;
  struct kshell_t {unsigned k; unsigned nmin; unsigned nmax;};
  std::vector<kshell_t> kshells_up_;
  std::vector<kshell_t> kshells_dn_;

  void build_unitcell_terms(const lattice::LatticeGraph& graph);
  void update_unitcell_terms(void);
  void set_particle_num(const input::Parameters& inputs);
  void construct_groundstate(void);
  void compute_averages(const lattice::LatticeGraph& graph, SR_Params& srparams);
};


} // end namespace srmf

#endif