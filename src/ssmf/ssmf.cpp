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
#include <boost/algorithm/string.hpp>
#include "ssmf.h"

namespace ssmf {


SSMF::SSMF(const input::Parameters& inputs) 
  : graph_(inputs) 
  , model_(inputs, graph_.lattice())
  , mf_params_(inputs, graph_,model_)
  , spinon_model_(inputs,model_,graph_,mf_params_)
  , boson_model_(inputs,model_,graph_,mf_params_)
{
  // For solving for LM parameter equation
  /*
  fx_dim_ = graph_.lattice().num_basis_orbitals();
  x_vec_.resize(fx_dim_);
  fx_vec_.resize(fx_dim_);
  solver_.allocate(fx_dim_);
  */
  std::string job = inputs.set_value("job", "SSMF");
  boost::trim(job);
  boost::to_upper(job);
  if (job == "DIAGONALIZATION") diag_only_ = true;
  diag_only_ = false;
  conv_tol_ = inputs.set_value("ssmf_conv_tol", 1.0E-6);

  std::string prefix = inputs.set_value("prefix", "results");
  prefix = "./"+prefix+"/";
  boost::filesystem::path prefix_dir(prefix);
  boost::filesystem::create_directory(prefix_dir);
  make_info_str(inputs);
  //file_conv_data_.init(prefix, "conv_data", info_str_.str());
  for (int i=0; i<mf_params_.num_sites(); ++i) {
    mfp_files_.push_back(std::make_shared<file::DataFile>());
    mfp_files_.back()->init(prefix, "mfp_site"+std::to_string(i), info_str_.str());
    siteavg_files_.push_back(std::make_shared<file::DataFile>());
    siteavg_files_.back()->init(prefix, "siteavg_site"+std::to_string(i), info_str_.str());
  }
  //file_mfp_.init(prefix, "mfp", info_str_.str());
  //file_sp_site_.init(prefix, "sp_site", info_str_.str());
  //file_sp_bond_.init(prefix, "sp_bond", info_str_.str());
  file_energy_.init(prefix, "energy", info_str_.str());
  //site_avg_.init("site_avg", heading);
  //bond_avg_.init("bond_avg", heading);
  spinon_model_.init_files(prefix, info_str_.str());
}


int SSMF::run(const input::Parameters& inputs) 
{
  spinon_model_.update(inputs);
  boson_model_.update(spinon_model_);

  //mf_params_.
  boson_ke_.resize(mf_params_.num_bonds());
  spinon_ke_.resize(mf_params_.num_bonds());
  diff_boson_ke_.resize(mf_params_.num_bonds());
  diff_spinon_ke_.resize(mf_params_.num_bonds());
  boson_ke_norm_.resize(mf_params_.num_bonds());
  spinon_ke_norm_.resize(mf_params_.num_bonds());
  qp_weights_norm_.resize(mf_params_.num_sites());

  for (int i=0; i<mf_params_.num_bonds(); ++i) {
    int rows = mf_params_.bond(i).spinon_ke(0).rows();
    int cols = mf_params_.bond(i).spinon_ke(0).cols();
    boson_ke_[i].resize(rows, cols);
    spinon_ke_[i].resize(rows, cols);
    diff_boson_ke_[i].resize(rows, cols);
    diff_spinon_ke_[i].resize(rows, cols);
  }

  mf_params_.init_params();
  for (int i=0; i<mf_params_.num_bonds(); ++i) {
    spinon_ke_[i] = mf_params_.bond(i).spinon_ke(0);
    boson_ke_[i] = mf_params_.bond(i).boson_ke(0);
  }

  // only 'diagonalization'
  if (diag_only_) {
    spinon_model_.solve(graph_,mf_params_);
    spinon_model_.print_output(mf_params_);
    return 0;
  }

  // find spinon shifted local energies for U=0;
  /*
  double U = spinon_model_.get_parameter_value("U");
  if (std::abs(U)<1.0E-12) {
    for (int i=0; i<fx_dim_; ++i) {
      x_vec_[i] = spinon_model_.orbital_en()[i];
    }
    solver_.find_root(this, &gsl_problem_eqn1, x_vec_, lm_ftol_);
    spinon_model_.set_shifted_en(x_vec_);
  }
  std::cout << "shifted_e0 =" << x_vec_[0] << "\n";
  */
  // solve 
  selconsistent_solve();
  //spinon_model_.print_results();
  //boson_model_.print_results();

  return 0;
}

int SSMF::selconsistent_solve(void) 
{
  bool converged = false;
  for (int iter=0; iter<max_ssmf_iter_; ++iter) {
    spinon_model_.solve(graph_,mf_params_);
    boson_model_.solve(mf_params_);
    // check convergence
    for (int i=0; i<mf_params_.num_bonds(); ++i) {
      diff_spinon_ke_[i] = spinon_ke_[i] - mf_params_.bond(i).spinon_ke(0);
      diff_boson_ke_[i] = boson_ke_[i] - mf_params_.bond(i).boson_ke(0);
      spinon_ke_norm_[i] = diff_spinon_ke_[i].abs2().maxCoeff();
      boson_ke_norm_[i] = diff_boson_ke_[i].abs2().maxCoeff();
      //std::cout << "spinon_ke_norm["<<i<<"] = "<<spinon_ke_norm_[i] << "\n";
    }
    //getchar();
    double sp_norm = spinon_ke_norm_.maxCoeff();
    double sb_norm = boson_ke_norm_.maxCoeff();

    std::cout<<"ssmf iter="<<iter+1<<", norm=("<<sp_norm<<","<<sb_norm<<")\n";
    if (sp_norm<conv_tol_ && sb_norm<conv_tol_) {
      converged = true;
      break;
    }
    // stop if all QP weights becomes zero
    for (int i=0; i<mf_params_.num_sites(); ++i) {
      //qp_weights_norm_[i] = mf_params_.site(i).qp_weights().minCoeff();
      qp_weights_norm_[i] = mf_params_.site(i).qp_weights().maxCoeff();
    }
    double z_norm = qp_weights_norm_.maxCoeff();
    if (z_norm<1.0E-6) {
      converged = true;
      break;
    }

    // continue
    for (int i=0; i<mf_params_.num_bonds(); ++i) {
      spinon_ke_[i] = mf_params_.bond(i).spinon_ke(0);
      boson_ke_[i] = mf_params_.bond(i).boson_ke(0);
    }
  }
  //if (converged) std::cout<<"ssmf converged!\n";

  for (int i=0; i<mf_params_.num_sites(); ++i) {
    std::cout<<spinon_model_.get_parameter_value("U")<<"  ";
    //std::cout<<"Z["<<i<<"] = "<< sr_parms_.site(i).qp_weights().transpose()<<"\n";
    std::cout<<mf_params_.site(i).qp_weights().transpose()<<"\n";
  }
  for (int i=0; i<mf_params_.num_sites(); ++i) {
    std::cout<<spinon_model_.get_parameter_value("U")<<"  ";
    std::cout<<mf_params_.site(i).lm_params().transpose()<<"\n";
  }
  double KE = spinon_model_.energy(mf_params_);
  double PE = boson_model_.interaction_energy();
  double E = KE+PE;
  std::cout << "KE = " << KE << "\n";
  std::cout << "PE = " << PE << "\n";
  std::cout << "TE = " << E << "\n";

  spinon_model_.print_output(mf_params_);
  print_output();

  return 0;
}

void SSMF::print_output(void)
{
  //---------------------Z and lambda-------------------------
  int s = -1;
  for (auto& file : mfp_files_) {
    s++;
    file->open();
    if (!heading_printed_) { 
      file->fs()<<std::left<< "   ";
      for (const auto& pname : spinon_model_.pnames()) {
        file->fs()<<std::left<<std::setw(15)<< pname;
      }
      for (int m=0; m<mf_params_.site(0).spin_orbitals().size(); ++m) {
        file->fs()<<std::left<<std::setw(15)<<"Z"+std::to_string(m);
      }
      for (int m=0; m<mf_params_.site(0).spin_orbitals().size(); ++m) {
        file->fs()<<std::left<<std::setw(15)<<"h"+std::to_string(m);
      }
      file->fs() << "\n";
      file->fs()<<"#"<< std::string(72, '-') << "\n";
    }
    file->fs()<<std::scientific<<std::uppercase<<std::setprecision(6);
    for (const auto& pval : spinon_model_.pvals()) {
      file->fs()<<std::right<<std::setw(15)<<pval;
    }
    for (int m=0; m<mf_params_.site(s).spin_orbitals().size(); ++m) {
      file->fs()<<std::setw(15)<<mf_params_.site(s).qp_weights()[m];
    }
    for (int m=0; m<mf_params_.site(s).spin_orbitals().size(); ++m) {
      file->fs()<<std::setw(15)<<mf_params_.site(s).lm_params()[m];
    }
    file->fs() << "\n";
    file->close();
  }

  s = -1;
  for (auto& file : siteavg_files_) {
    s++;
    file->open();
    if (!heading_printed_) { 
      file->fs()<<std::left<< "   ";
      for (const auto& pname : spinon_model_.pnames()) {
        file->fs()<<std::left<<std::setw(15)<< pname;
      }
      for (int m=0; m<mf_params_.site(0).spin_orbitals().size(); ++m) {
        file->fs()<<std::left<<std::setw(15)<<"n"+std::to_string(m);
      }
      file->fs() << "\n";
      file->fs()<<"#"<< std::string(72, '-') << "\n";
    }
    file->fs()<<std::scientific<<std::uppercase<<std::setprecision(6);
    for (const auto& pval : spinon_model_.pvals()) {
      file->fs()<<std::right<<std::setw(15)<<pval;
    }
    for (int m=0; m<mf_params_.site(s).spin_orbitals().size(); ++m) {
      file->fs()<<std::setw(15)<<mf_params_.site(s).spinon_density()[m];
    }
    file->fs() << "\n";
    file->close();
  }


  //---------------------energy-------------------------
  file_energy_.open();
  if (!heading_printed_) { 
    file_energy_.fs()<<std::left<< "   ";
    for (const auto& pname : spinon_model_.pnames()) {
      file_energy_.fs()<<std::left<<std::setw(15)<< pname;
    }
    file_energy_.fs()<<std::left<<std::setw(15)<< "TE";
    file_energy_.fs()<<std::left<<std::setw(15)<< "KE";
    file_energy_.fs()<<std::left<<std::setw(15)<< "PE";
    file_energy_.fs()<<"\n";
    file_energy_.fs()<<"#"<< std::string(72, '-') << "\n";
  }
  double KE = spinon_model_.energy(mf_params_);
  double PE = boson_model_.interaction_energy();
  file_energy_.fs()<<std::scientific<<std::uppercase<<std::setprecision(6);
  for (const auto& pval : spinon_model_.pvals()) {
    file_energy_.fs()<<std::right<<std::setw(15)<<pval;
  }
  file_energy_.fs()<<std::right<<std::setw(15)<<KE+PE;
  file_energy_.fs()<<std::right<<std::setw(15)<<KE;
  file_energy_.fs()<<std::right<<std::setw(15)<<PE;
  file_energy_.fs() << "\n";
  file_energy_.close();


  heading_printed_ = true;
}

void SSMF::make_info_str(const input::Parameters& inputs)
{
  info_str_.clear();
  print_copyright(info_str_);
  info_str_ << "# "<< inputs.job_id() <<"\n"; 
  info_str_ << spinon_model_.info_str(); 
  info_str_ << boson_model_.info_str(); 
  info_str_ << "#" << std::string(72, '-') << "\n";
  info_str_ << "# SSMF solution:\n"; 
  info_str_ << "# Max iterations = "<<max_ssmf_iter_<<"\n"; 
  info_str_ << std::scientific<<std::uppercase<<std::setprecision(2);
  info_str_ << "# Conv tolerance = "<<conv_tol_<<"\n"; 

  //info_str_ << config.info_str(); 
  //info_str_ << "# Samples = " << num_measure_steps_;
  //info_str_ << ", warmup = " << num_warmup_steps_;
  //info_str_ << ", min_interval = " << min_interval_;
  //info_str_ << ", max_interval = " << max_interval_ << "\n";
}

void SSMF::print_copyright(std::ostream& os)
{
  os << "#" << std::string(72,'-') << "\n";
  os << "#" << " Program: Slave Spinon Mean Field (SSMF) calculation\n";
  os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n";
  os << "#" << std::string(72,'-') << "\n";
}

/*
int gsl_problem_eqn1(const gsl_vector* x, void* parms, gsl_vector* f)
{
  SSMF * pThis = ((class SSMF *) parms);
  for (int i=0; i<pThis->fx_dim_; ++i) {
    pThis->x_vec_[i] = gsl_vector_get(x,i);
  }
  int status = pThis->spinon_energy_eqn(pThis->x_vec_, pThis->fx_vec_);
  for (int i=0; i<pThis->fx_dim_; ++i) {
    gsl_vector_set(f, i, pThis->fx_vec_[i]);
  }
  if (status ==0 ) return GSL_SUCCESS;
  else return GSL_FAILURE;
}

int SSMF::spinon_energy_eqn(const std::vector<double>& x, std::vector<double>& fx)
{
  spinon_model_.set_shifted_en(x);
  selconsistent_solve();
  auto lambda = mf_params_.site(0).lm_params();
  std::cout<<"x(0) = "<< x[0] << " x(1) ="<<x[1]<<"\n";
  std::cout<<"lambda = "<< lambda.transpose()<<"\n"; getchar();
  for (int i=0; i<fx_dim_; ++i) {
    fx[i] = x[i] - lambda[i] - spinon_model_.orbital_en()[i];
    std::cout << "fx["<<i<<"] = "<< fx[i] << "\n";
  }
  return 0;
}
*/


} // end namespace diag
