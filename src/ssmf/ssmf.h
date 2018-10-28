/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef DIAG_H
#define DIAG_H

#include <iostream>
#include <Eigen/Eigenvalues>
#include "../scheduler/worker.h"
#include "../lattice/constants.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
//#include "../lattice/blochbasis.h"
#include "../model/model.h"
#include "sb_params.h"
#include "slavespin.h"
#include "spinon.h"

namespace srmf {

class SRMF : public scheduler::Worker
{
public:
  SRMF(const input::Parameters& parms); 
  ~SRMF() {}
  int start(const input::Parameters& parms) override { return 0; }
  int run(const input::Parameters& parms) override;
  void finish(void) override {} 
  void dostep(void) override {} 
  void halt(void) override {} 
  static void print_copyright(std::ostream& os);
  friend int gsl_problem_eqn1(const gsl_vector* x, void* parms, gsl_vector* f);
private:
  lattice::LatticeGraph graph_;
  basis::BlochBasis blochbasis_;
  model::Hamiltonian model_;
  SB_Params sr_parms_;
  Spinon spinon_model_;

  // gsl solver
  //double lm_ftol_{1.0E-8};
  unsigned fx_dim_;
  std::vector<double> x_vec_;
  std::vector<double> fx_vec_;
  root::gsl_solver solver_;
  SlaveSpin boson_model_;
  //unsigned num_kpoints_{1};
  //unsigned kblock_dim_{1};
  //mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_k_up_;
  //std::vector<Vector3d> symm_line_;
  // outputs
  /*bool need_chern_number_{false};
  bool need_ebands_full_{false};
  bool need_ebands_symm_{false};
  bool need_band_gap_{false};*/
  // work arrays
  cmpl_bondparms_t boson_ke_;
  cmpl_bondparms_t spinon_ke_;
  cmpl_bondparms_t diff_boson_ke_;
  cmpl_bondparms_t diff_spinon_ke_;
  realArray1D boson_ke_norm_;
  realArray1D spinon_ke_norm_;

  //int compute_chern_number(void);
  //int compute_band_gap(void);
  int selconsistent_solve(void);
  int spinon_energy_eqn(const std::vector<double>& x, std::vector<double>& fx);
};


} // end namespace srmf

#endif
