/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef DIAG_H
#define DIAG_H

#include <iostream>
#include <memory>
#include <Eigen/Eigenvalues>
#include "../scheduler/worker.h"
#include "../lattice/constants.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "mf_params.h"
#include "slavespin.h"
#include "spinon.h"
#include "datafile.h"
#include <boost/filesystem.hpp>

namespace ssmf {

class SSMF : public scheduler::Worker
{
public:
  SSMF(const input::Parameters& parms); 
  ~SSMF() {}
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
  MF_Params mf_params_;
  Spinon spinon_model_;
  SlaveSpin boson_model_;
  bool diag_only_{false};
  double conv_tol_{1.0E-8};
  int max_ssmf_iter_{200};

  // output data
  std::vector<std::shared_ptr<file::DataFile>> mfp_files_;
  std::vector<std::shared_ptr<file::DataFile>> siteavg_files_;
  file::DataFile file_conv_data_;
  file::DataFile file_mfp_;
  file::DataFile file_sp_site_;
  file::DataFile file_sp_bond_;
  file::DataFile file_energy_;
  bool heading_printed_{false};

  // gsl solver
  //double lm_ftol_{1.0E-8};
  //unsigned fx_dim_;
  //std::vector<double> x_vec_;
  //std::vector<double> fx_vec_;
  //root::gsl_solver solver_;

  // work arrays
  cmpl_bondparms_t boson_ke_;
  cmpl_bondparms_t spinon_ke_;
  cmpl_bondparms_t diff_boson_ke_;
  cmpl_bondparms_t diff_spinon_ke_;
  realArray1D boson_ke_norm_;
  realArray1D spinon_ke_norm_;
  realArray1D qp_weights_norm_;

  mutable std::ostringstream info_str_;

  //int compute_chern_number(void);
  //int compute_band_gap(void);
  void print_output(void);
  void make_info_str(const input::Parameters& inputs);
  int selconsistent_solve(void);
  int spinon_energy_eqn(const std::vector<double>& x, std::vector<double>& fx);
};


} // end namespace ssmf

#endif
