/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-22 13:00:28
*----------------------------------------------------------------------------*/
#include <iomanip>
#include <fstream>
#include "srmf.h"

namespace srmf {

SRMF::SRMF(const input::Parameters& inputs) 
  : graph_(inputs) 
  //, blochbasis_(graph_)
  , model_(inputs, graph_.lattice())
  , sr_parms_(inputs, graph_,model_)
  , spinon_model_(inputs,model_,graph_,sr_parms_)
  , boson_model_(inputs,model_,graph_,sr_parms_)
{
}

int SRMF::run(const input::Parameters& inputs) 
{
  spinon_model_.update(inputs);
  boson_model_.update(spinon_model_);

  //sr_parms_.
  cmpl_bondparms_t boson_ke(sr_parms_.num_bonds());
  cmpl_bondparms_t spinon_ke(sr_parms_.num_bonds());
  cmpl_bondparms_t diff_boson_ke(sr_parms_.num_bonds());
  cmpl_bondparms_t diff_spinon_ke(sr_parms_.num_bonds());
  realArray1D boson_ke_norm(sr_parms_.num_bonds());
  realArray1D spinon_ke_norm(sr_parms_.num_bonds());

  sr_parms_.init_mf_params();
  for (int i=0; i<sr_parms_.num_bonds(); ++i) {
    //sr_parms_.bond(i).set_spinon_ke();
    //sr_parms_.bond(i).set_boson_ke();
    spinon_ke[i] = sr_parms_.bond(i).spinon_ke();
    boson_ke[i] = sr_parms_.bond(i).boson_ke();
  }

  int max_sb_iter = 10;
  bool converged = false;
  for (int iter=0; iter<max_sb_iter; ++iter) {
    spinon_model_.solve(graph_, sr_parms_);
    boson_model_.solve(sr_parms_);
    // check convergence
    for (int i=0; i<sr_parms_.num_bonds(); ++i) {
      diff_spinon_ke[i] = spinon_ke[i] - sr_parms_.bond(i).spinon_ke();
      diff_boson_ke[i] = boson_ke[i] - sr_parms_.bond(i).boson_ke();
      spinon_ke_norm[i] = diff_spinon_ke[i].abs2().maxCoeff();
      boson_ke_norm[i] = diff_boson_ke[i].abs2().maxCoeff();
    }
    double sp_norm = spinon_ke_norm.maxCoeff();
    double sb_norm = boson_ke_norm.maxCoeff();

    //std::cout<<"ssmf iter="<<iter+1<<", norm=("<<sp_norm<<","<<sb_norm<<")\n";
    if (sp_norm<1.0E-6 && sb_norm<1.0E-6) {
      converged = true;
      break;
    }

    // continue
    for (int i=0; i<sr_parms_.num_bonds(); ++i) {
      spinon_ke[i] = sr_parms_.bond(i).spinon_ke();
      boson_ke[i] = sr_parms_.bond(i).boson_ke();
    }
  }
  //if (converged) std::cout<<"ssmf converged!\n";

  std::cout<<spinon_model_.get_parameter_value("U")<<"  ";
  for (int i=0; i<sr_parms_.num_sites(); ++i) {
    //std::cout<<"Z["<<i<<"] = "<< sr_parms_.site(i).qp_weights().transpose()<<"\n";
    std::cout<<sr_parms_.site(i).qp_weights().transpose()<<"\n";
  }


  return 0;
}


void SRMF::print_copyright(std::ostream& os)
{
  std::cout << "SRMF\n";
}


} // end namespace diag
