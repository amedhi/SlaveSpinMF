/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 18:54:09
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-17 12:12:07
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <numeric>
#include "./spinon.h"
#include <boost/algorithm/string.hpp>
#include <boost/math/tools/roots.hpp>
#include "../expression/complex_expression.h"
//#include "../expression/expression.h"

namespace ssmf {

Spinon::Spinon(const input::Parameters& inputs, const model::Hamiltonian& model, const lattice::LatticeGraph& graph, const MF_Params& mf_params)
: model::Hamiltonian(model)
, blochbasis_(graph)
{
  // model parameters
  if (have_parameter("U")) {
    modelp_names_.push_back("U");
    modelp_vals_.push_back(get_parameter_value("U"));
  }
  if (have_parameter("J")) {
    modelp_names_.push_back("J");
    modelp_vals_.push_back(get_parameter_value("J"));
  }
  if (have_parameter("lambda")) {
    modelp_names_.push_back("lambda");
    modelp_vals_.push_back(get_parameter_value("lambda"));
  }
  if (have_parameter("ext_field")) {
    modelp_names_.push_back("ext_field");
    modelp_vals_.push_back(get_parameter_value("ext_field"));
  }

  // add mean-field terms to the spinon model
  if (graph.lattice().id()==lattice::lattice_id::SQUARE_2BAND) {
  }
  if (graph.lattice().id()==lattice::lattice_id::PYROCHLORE_3D) {
  }
  int nowarn;
  no_spinon_lambda_ = inputs.set_value("no_spinon_lambda",true,nowarn);
  assume_fixed_groundstate_ = inputs.set_value("assume_fixed_groundstate", false, nowarn);
  have_TP_symmetry_ = model.have_TP_symmetry();
  SO_coupling_ = model.is_spinorbit_coupled();
  num_sites_ = graph.num_sites();
  num_bonds_ = graph.num_bonds();
  num_unitcells_ = graph.lattice().num_unitcells();
  num_kpoints_ = blochbasis_.num_kpoints(); 
  num_symm_kpoints_ = blochbasis_.num_symm_kpoints(); 
  kblock_dim_ = blochbasis_.subspace_dimension();
  num_basis_sites_ = graph.lattice().num_basis_sites();
  num_total_states_ = num_unitcells_ * kblock_dim_;
  if (num_sites_ != num_unitcells_*num_basis_sites_) {
    throw std::logic_error("Spinon::Spinon: Incorrect 'lattice'\n");
  }
  //dim_ = graph.lattice().num_basis_orbitals();
  quadratic_block_up_.resize(kblock_dim_,kblock_dim_);
  pairing_block_.resize(kblock_dim_,kblock_dim_);
  work.resize(kblock_dim_,kblock_dim_);
  build_unitcell_terms(graph);
  set_particle_num(inputs);
  set_info_string();
  iteration_zero_ = true;
}

int Spinon::init(const lattice::Lattice& lattice)
{
  return Model::init(lattice);
}

int Spinon::finalize(const lattice::LatticeGraph& graph)
{
  Model::finalize(graph.lattice());
  num_basis_sites_ = graph.lattice().num_basis_sites();
  //dim_ = graph.lattice().num_basis_orbitals();
  quadratic_block_up_.resize(kblock_dim_,kblock_dim_);
  pairing_block_.resize(kblock_dim_,kblock_dim_);
  work.resize(kblock_dim_,kblock_dim_);
  build_unitcell_terms(graph);
  set_info_string();
  iteration_zero_ = true;
  return 0;
}

void Spinon::set_info_string(void)
{
  std::ostringstream info_strm;
  info_strm.clear();
  info_strm << "#" << std::string(72, '-') << "\n";
  info_strm << "# Spinon Sector:\n";
  if (SO_coupling_) info_strm << "# SO coupling = YES\n";
  else info_strm << "# SO coupling = NO\n";
  if (no_spinon_lambda_) info_strm << "# Lagrange coupling = OFF\n";
  else info_strm << "# Lagrange coupling = ON\n";
  info_strm << "# Number of states = "<<num_total_states_<<"\n";
  info_strm << "# Number of particles = "<<num_spins_;
  info_strm << " (per site = "<<double(num_spins_)/num_sites_<<")\n";
  info_strm << "# Hole doping = "<<hole_doping_<<"\n";
  info_str_ = model::Hamiltonian::info_str()+info_strm.str();
}

void Spinon::update(const input::Parameters& inputs)
{
  //std::cout << "START>>>>>>>>>>\n";
  //std::cout << "<<<<<<<<<<<<END\n";
  Model::update_parameters(inputs);
  update_terms();
  for (int i=0; i<modelp_names_.size(); ++i) {
    modelp_vals_[i] = get_parameter_value(modelp_names_[i]);
  }
}

void Spinon::init_files(const std::string& prefix, const std::string& heading)
{
  // observables
  file_bands_.init(prefix, "bands", heading);
}

void Spinon::solve(const lattice::LatticeGraph& graph, MF_Params& mf_params)
{
  // renormalized bond couplings
  int term_id=0;
  for (auto& term : ubond_terms_) {
    if (term.qn_operator().is_quadratic() && term.qn_operator().spin_up()) {
      term.update_bondterm_cc(term_id++, mf_params);
    }
  }
  // renormalized soc couplings
  for (auto& term : usite_terms_) {
    if (term.qn_operator().id()==model::op_id::spin_flip) {
      term.update_siteterm_cc(mf_params);
    }
  }
  if (!assume_fixed_groundstate_) {
    construct_groundstate_v2(mf_params);
    //construct_groundstate(mf_params);
  }
  else {
    if (!groundstate_determined_) {
      construct_groundstate_v2(mf_params);
      //construct_groundstate(mf_params);
      groundstate_determined_ = true;
    }
  }
  compute_averages(graph,mf_params);
  //std::cout << "START>>>>>>>>>>\n";
  //std::cout << "END>>>>>>>>>>\n";
  //std::cout << "KE = " << energy(mf_params) << "\n";
  //std::exit(0);
}

void Spinon::compute_averages(const lattice::LatticeGraph& graph, MF_Params& mf_params)
{
  //Eigen::VectorXcd amplitude_vec1, amplitude_vec2;
  Eigen::Matrix<std::complex<double>,1,Eigen::Dynamic> amplitude_vec1, amplitude_vec2;
  for (int i=0; i<mf_params.num_sites(); ++i) {
    mf_params.site(i).spinon_density().setZero();
    mf_params.site(i).spinon_flip_ampl().setZero();
  }
  for (int i=0; i<mf_params.bonds().size(); ++i) {
    mf_params.bond(i).spinon_ke(0).setZero();
  }

  for (int i=0; i<kshells_up_.size(); ++i) {
    int k = kshells_up_[i].k;
    Vector3d kvec = blochbasis_.kvector(k);
    double symm_wt = blochbasis_.kweight(k);
    //std::cout << "Hack: Spinon:compute_averages\n";
    //kvec = Vector3d(-0.785, 0.785, 0.785);
    //std::cout << "symm_wt = " << symm_wt << "\n"; getchar();
    construct_kspace_block(mf_params, kvec);
    es_k_up_.compute(quadratic_spinup_block());

    // site density
    int nmin = kshells_up_[i].nmin;
    int nmax = kshells_up_[i].nmax;
    int nbands = nmax-nmin+1;
    amplitude_vec1.resize(nbands);
    amplitude_vec2.resize(nbands);
    //Eigen::VectorXcd eigvec_wt(nbands);
    for (int j=0; j<mf_params.num_sites(); ++j) {
      int site_dim = mf_params.site(j).dim();
      realArray1D n_avg(site_dim); 
      for (int m=0; m<site_dim; ++m) {
        auto ii = mf_params.site(j).state_indices()[m];
        double norm = 0.0;
        for (int band=nmin; band<=nmax; ++band) {
          norm += std::norm(es_k_up_.eigenvectors().row(ii)[band])
                  * kshells_up_[i].smear_wt(band);
        }
        n_avg(m) = norm;
      }
      //std::cout << n_avg.transpose(); getchar();
      mf_params.site(j).spinon_density() += symm_wt*n_avg;
    }

    // spin-flip term 
    if (SO_coupling_) {
      for (int j=0; j<mf_params.num_sites(); ++j) {
        int site_dim = mf_params.site(j).dim();
        cmplArray2D flip_ampl(site_dim,site_dim);
        for (int m=0; m<site_dim; ++m) {
          auto ii = mf_params.site(j).state_indices()[m];
          amplitude_vec1 = es_k_up_.eigenvectors().block(ii,nmin,1,nbands);
          for (int n=0; n<site_dim; ++n) {
            auto jj = mf_params.site(j).state_indices()[n];
            amplitude_vec2 = es_k_up_.eigenvectors().block(jj,nmin,1,nbands);
            //flip_ampl(m,n) = amplitude_vec1.dot(amplitude_vec2);
            std::complex<double> chi_sum = 0.0;
            for (int band=nmin; band<=nmax; ++band) {
              chi_sum += std::conj(amplitude_vec1(band)) * amplitude_vec2(band)
                       * kshells_up_[i].smear_wt(band);
            }
            flip_ampl(m,n) = chi_sum;
          }
        }
        mf_params.site(j).spinon_flip_ampl() += symm_wt*flip_ampl;
      }
    }

    // bond averages
    for (int j=0; j<mf_params.num_bonds(); ++j) {
      Vector3d delta = mf_params.bond(j).vector();
      std::complex<double> exp_kdotr = std::exp(ii()*kvec.dot(delta));
      int src = mf_params.bond(j).src();
      int tgt = mf_params.bond(j).tgt();
      int rows = mf_params.site(src).dim();
      int cols = mf_params.site(tgt).dim();
      //sr_bond bond(mf_params.bonds()[j]);
      //std::cout << "bond: src - tgt = "<<j<<": "<<src<<"-"<<tgt<< "\n";
      cmplArray2D ke_matrix=cmplArray2D::Zero(rows, cols);
      for (int m=0; m<rows; ++m) {
        auto ii = mf_params.site(src).state_indices()[m];
        amplitude_vec1 = es_k_up_.eigenvectors().block(ii,nmin,1,nbands);
        //std::cout << "m, ii = " << m << ", " << ii << "\n";
        //std::cout << "amplitude_vec1 = " << amplitude_vec1.transpose() << "\n"; getchar();
        for (int n=0; n<cols; ++n) {
          auto jj = mf_params.site(tgt).state_indices()[n];
          amplitude_vec2 = es_k_up_.eigenvectors().block(jj,nmin,1,nbands);
          //std::cout << "n, jj = " << n << ", " << jj << "\n";
          //getchar();
          //std::complex<double> chi_sum = amplitude_vec1.dot(amplitude_vec2);
          std::complex<double> chi_sum = 0.0;
          for (int band=nmin; band<=nmax; ++band) {
            chi_sum += std::conj(amplitude_vec1(band)) * amplitude_vec2(band)
                      * kshells_up_[i].smear_wt(band);
          }

          ke_matrix(m,n) = exp_kdotr * chi_sum;
          //std::complex<double> chi_sum(0.0);
          //for (unsigned band=nmin; band<=nmax; ++band) {
          //  chi_sum += exp_kdotr * std::conj(es_k_up_.eigenvectors().row(ii)[band]) * 
          //         es_k_up_.eigenvectors().row(jj)[band];
          //}
          //ke_matrix(m,n) = chi_sum;
        }
      }
      mf_params.bond(j).spinon_ke(0) += symm_wt*ke_matrix;
    }
  }

  // final site density 
  /*
  realArray1D n_sum(mf_params.site(0).dim());
  n_sum.setZero();
  for (int i=0; i<mf_params.num_sites(); ++i) {
    n_sum += mf_params.site(i).spinon_density();
  }
  for (int i=0; i<mf_params.num_sites(); ++i) {
    mf_params.site(i).spinon_density() = n_sum/mf_params.num_sites();
  }
  */

  for (int i=0; i<mf_params.num_sites(); ++i) {
    //mf_params.site(i).spinon_density() /= num_kpoints_;
    realArray1D n_avg = mf_params.site(i).spinon_density()/num_kpoints_;
    mf_params.site(i).spinon_density() = n_avg;
    // print
    //*
    std::cout<<"site-"<<i<<":"<<"\n";
    for (int m=0; m<mf_params.site(i).dim(); ++m) {
      std::cout << "<n_up>["<<m<<"] = " << mf_params.site(i).spinon_density()[m] << "\n";
    }
    std::cout << "\n";
    //*/
  }
  //getchar();
  // onsite energy
  double onsite_energy = 0.0;
  for (const auto& term : usite_terms_) {
    if (term.name() != "ExtField") {
      auto e0 = term.coeff_matrix().diagonal().real();
      for (int i=0; i<mf_params.num_sites(); ++i) {
        int site_dim = mf_params.site(i).dim();
        for (int m=0; m<site_dim; ++m) {
          auto ii = mf_params.site(i).state_indices()[m];
          onsite_energy += e0(ii)*mf_params.site(i).spinon_density()[m];
        }
      }
      onsite_energy /= mf_params.num_sites();
    }
  }

  // final spin-flip amplitudes
  double soc_energy = 0.0;
  if (SO_coupling_) {
    for (int i=0; i<mf_params.num_sites(); ++i) {
      cmplArray2D flip_ampl = mf_params.site(i).spinon_flip_ampl()/num_kpoints_;
      mf_params.site(i).spinon_flip_ampl() = flip_ampl;
      mf_params.site(i).set_spinon_renormalization();
      // print
      /*
      std::cout << std::setprecision(2);
      std::cout << std::scientific;
      std::cout<<"site-"<<i<<":"<<"\n";
      std::cout << "sf ampl = \n" << mf_params.site(i).spinon_flip_ampl() << "\n";
      std::cout << "\n";
      */
      flip_ampl = flip_ampl*mf_params.site(i).boson_renormed_soc();
      soc_energy += flip_ampl.real().sum();
    }
    soc_energy /= mf_params.num_sites();
  }

  // final KE average 
  double bond_en = 0.0;
  for (int i=0; i<mf_params.num_bonds(); ++i) {
    cmplArray2D ke_matrix = mf_params.bond(i).spinon_ke(0)/num_kpoints_;
    mf_params.bond(i).spinon_ke(0) = ke_matrix;
    mf_params.bond(i).set_spinon_renormalization();

    //bond_en += 2*mf_params.bond(i).spinon_renormed_cc(0).real().sum();
    ke_matrix = ke_matrix*mf_params.bond(i).boson_renormed_cc(0);
    bond_en += 2*ke_matrix.real().sum();
    //std::cout << "ke_matrix =\n" << ke_matrix << "\n"; getchar();
    // print
    /*
    std::cout<<"bond-"<<i<<":"<<"\n";
    for (int m=0; m<mf_params.bond(i).spinon_ke(0).rows(); ++m) {
      for (int n=0; n<mf_params.bond(i).spinon_ke(0).cols(); ++n) {
        std::cout << "t   ["<<m<<","<<n<<"] = " << mf_params.bond(i).term_coupling(0)(m,n) << "\n";
        std::cout << "chi ["<<m<<","<<n<<"] = " << mf_params.bond(i).spinon_ke(0)(m,n) << "\n";
        std::cout << "tchi["<<m<<","<<n<<"] = " << mf_params.bond(i).spinon_renormed_cc(0)(m,n) << "\n\n";
      }
    }
    std::cout << "\n";
    getchar();
    */
  }
  bond_en /= mf_params.num_bonds();
  total_energy_ = onsite_energy + bond_en + soc_energy;
  mf_params.ke_per_site() = total_energy_;
  /*
  std::cout << "bond en = " << bond_en << "\n";
  std::cout << "site en = " << onsite_energy << "\n";
  std::cout << "soc en = " << soc_energy << "\n";
  getchar();
  */
}

void Spinon::print_output(const MF_Params& mf_params)
{
  file_bands_.open();
  int k = 0;
  for (const auto& kvec : blochbasis_.symm_path_k()) {
    construct_kspace_block(mf_params, kvec);
    es_k_up_.compute(quadratic_spinup_block(), Eigen::EigenvaluesOnly);
    file_bands_.fs()<<std::setw(6)<<k++; 
    file_bands_.fs()<<std::setw(14)<<kvec(0)<<std::setw(14)<< kvec(1); 
    file_bands_.fs()<<std::setw(14)<<es_k_up_.eigenvalues().transpose()<<"\n";
  }
  file_bands_.fs()<<"\n"; 
  file_bands_.close();
}

void Spinon::construct_groundstate(const MF_Params& mf_params)
{
  //if (have_TP_symmetry_) {
    /* Has T.P (Time Reversal * Inversion) symmetry. 
       So we have e_k(up) = e_k(dn).
    */
  if (true) {
    std::vector<std::pair<unsigned,unsigned>> qn_list; // list of (k,n)
    std::vector<double> ek;
    for (unsigned k=0; k<num_kpoints_; ++k) {
      Vector3d kvec = blochbasis_.kvector(k);
      construct_kspace_block(mf_params, kvec);
      es_k_up_.compute(quadratic_spinup_block(), Eigen::EigenvaluesOnly);
      ek.insert(ek.end(),es_k_up_.eigenvalues().data(),
        es_k_up_.eigenvalues().data()+kblock_dim_);
      //std::cout << kvec.transpose() << " " << es_k_up_.eigenvalues() << "\n"; getchar();
      for (unsigned n=0; n<kblock_dim_; ++n) {
        qn_list.push_back({k, n});
      }
    }
    //std::sort(ek.begin(),ek.end());
    //for (const auto& e : ek) std::cout << e << "\n";

    // Indices in the original ek-array is sorted according to increasing ek
    std::vector<int> idx(ek.size());
    std::iota(idx.begin(),idx.end(),0);
    std::sort(idx.begin(),idx.end(),[&ek](const int& i1, const int& i2) 
      {return ek[i1]<ek[i2];});
    /*for (int i=0; i<ek.size(); ++i) {
      std::cout << i << "  " << idx[i] << "  " << ek[idx[i]] << "\n";
    }*/
    // mean energy
    /*double e0 = 0.0;
    for (int i=0; i<num_upspins_; ++i) {
      e0 += ek[idx[i]];
    }
    e0 = 2 * e0 / num_sites_;
    std::cout << "e0 = " << e0 << "\n"; getchar();
    */
    // check for degeneracy 
    int num_fill_particles = num_spins_;
    double degeneracy_tol = 1.0E-12;
    int top_filled_level = num_fill_particles-1;
    fermi_energy_ = ek[idx[top_filled_level]];
    int num_degen_states = 1;
    int num_valence_particle = 1;
    // look upward in energy
    for (int i=top_filled_level+1; i<ek.size(); ++i) {
      if (std::abs(fermi_energy_-ek[idx[i]])>degeneracy_tol) break;
      num_degen_states++;
    }
    // look downward in energy
    if (num_degen_states>1) {
      for (int i=top_filled_level-1; i>=0; --i) {
        if (std::abs(fermi_energy_ - ek[idx[i]])>degeneracy_tol) break;
        num_degen_states++;
        num_valence_particle++;
      }
      // warn
      degeneracy_warning_ = true;
      if (degeneracy_warning_) {
        std::cout << ">> warning: Spinon groundstate degeneracy: " << 2*num_valence_particle <<
          " particles in " << 2*num_degen_states << " states." << "\n";
      }
    }
    /* 
      Filled k-shells. A k-shell is a group of energy levels having same 
      value of quantum number k.
    */
    // find 'nmax' values of filled k-shells
    std::vector<int> shell_nmax(num_kpoints_);
    for (auto& elem : shell_nmax) elem = -1; // invalid default value
    int k, n;
    for (int i=0; i<num_fill_particles; ++i) {
      int state = idx[i]; 
      std::tie(k,n) = qn_list[state];
      if (shell_nmax[k] < n) shell_nmax[k] = n;
    }
    // store the filled k-shells
    kshells_up_.clear();
    for (int k=0; k<num_kpoints_; ++k) {
      int nmax = shell_nmax[k];
      if (nmax != -1) {
        realArray1D smear_wts = realArray1D::Ones(nmax+1);
        kshells_up_.push_back({k,0,nmax,smear_wts});
      }
    }
    /*
    for (unsigned k=0; k<kshells_up_.size(); ++k) {
      std::cout << kshells_up_[k].k << " " << kshells_up_[k].nmin << "  "
          << kshells_up_[k].nmax << "\n";
    }
    */
  }
  else {
    /* 
      No T.P (Time Reversal * Inversion) symmetry. 
      Need to consider both UP and DN states.
    */
    throw std::range_error("Spinon::construct_groundstate: case not implemented\n");
  }
}

void Spinon::construct_groundstate_v2(const MF_Params& mf_params)
{
  // Check 'MethfesselPaxton_func'
  /*double x=-3.0;
  for (int n=0; n<600; ++n) {
    double fx = MethfesselPaxton_func(4,x);
    std::cout << x << "   " << fx << "\n";
    x += 0.01;
  }
  exit(0);*/
  //std::ofstream fs("bands.txt");
  //if (have_TP_symmetry_) {
    /* Has T.P (Time Reversal * Inversion) symmetry. 
       So we have e_k(up) = e_k(dn).
    */
  if (true) {
    qn_list_.clear();
    ek_list_.clear();
    for (int k=0; k<num_symm_kpoints_; ++k) {
      Vector3d kvec = blochbasis_.kvector(k);
      construct_kspace_block(mf_params, kvec);
      es_k_up_.compute(quadratic_spinup_block(), Eigen::EigenvaluesOnly);
      //std::cout << quadratic_spinup_block() << "\n"; getchar();
      //std::cout << es_k_up_.eigenvalues().transpose() << "\n";
      //fs<<k<<" "<<kvec.transpose()<<" "<<es_k_up_.eigenvalues().transpose()<<"\n";
      ek_list_.insert(ek_list_.end(),es_k_up_.eigenvalues().data(),
        es_k_up_.eigenvalues().data()+kblock_dim_);
      //std::cout << kvec.transpose() << "\n" << es_k_up_.eigenvalues() << "\n"; getchar();
      for (int n=0; n<kblock_dim_; ++n) {
        qn_list_.push_back({k, n});
      }
    }
    //fs.close();
    //std::sort(ek.begin(),ek.end());
    //for (const auto& e : ek) std::cout << e << "\n";

    // Indices in the original ek-array is sorted according to increasing ek
    idx_.resize(ek_list_.size());
    std::iota(idx_.begin(),idx_.end(),0);
    std::sort(idx_.begin(),idx_.end(),[this](const int& i1, const int& i2) 
      {return ek_list_[i1]<ek_list_[i2];});
    /* 
    for (int i=0; i<ek_list_.size(); ++i) {
      int j = idx_[i];
      std::cout<<"(i, k, n, E) = "<<i<<" "<<qn_list_[j].first
      <<" "<<qn_list_[j].second<<" "<<ek_list_[j]<< "\n";
    }
    */
    
    // Check whether metallic
    metallic_ = false;
    num_fill_particles_ = num_spins_;
    //std::cout << "num_fill_particles = " << num_fill_particles_ << "\n";
    //getchar();

    // Smearing Parameters
    bandwidth_ = ek_list_[idx_.back()]-ek_list_[idx_.front()];
    smear_width_ = 0.4*bandwidth_/num_unitcells_;
    smear_func_order_ = 4;
    //std::cout << "BW = " << bandwidth << "\n";
    //std::cout << "W = " << smear_width_ << "\n";
    //getchar();

    // for root finding
    double factor = 2.0;
    const boost::uintmax_t maxit = 200; 
    boost::uintmax_t it = maxit; 
    bool is_rising = true;
    boost::math::tools::eps_tolerance<double> tol(15);

    // Fermi Energy 
    int np = 0;
    double gap_tol = bandwidth_ * 1.0E-8;
    fermi_energy_ = ek_list_[idx_[0]];
    for (int i=0; i<ek_list_.size(); ++i) {
      int k = qn_list_[idx_[i]].first;
      //int n = qn_list_[idx_[i]].second;
      int iw = std::nearbyint(blochbasis_.kweight(k));
      double ek = ek_list_[idx_[i]];
      np += iw;
      if (np == num_fill_particles_) {
        // in case all states are filled
        if (np == ek_list_.size()) {
          metallic_ = false;
          fermi_energy_ = ek_list_.back();
          break;
        }
        // all states are not filled
        if (std::abs(ek_list_[idx_[i+1]]-ek) > gap_tol) {
          //std::cout <<  "ek = " << ek << "\n";
          //std::cout << "np, num_fill_particles_ = " << np << "  " << num_fill_particles_ << "\n";
          metallic_ = false;
          fermi_energy_ = 0.5*(ek+ek_list_[idx_[i+1]]);
          break;
        }
        else {
          metallic_ = true;
          fermi_energy_ = ek;
          //std::cout << qn_list_[idx_[i+1]].second << "\n";
          //std::cout << std::abs(ek_list_[idx_[i+1]]-ek); getchar();
          // find fermi energy by solving with smeared levels
          break;
        }
      }
      else if (np > num_fill_particles_) {
        metallic_ = true;
        fermi_energy_ = ek;
        break;
      }
    }

    // in case 'metallic_' fix fermi energy wrt smearing
    if (metallic_) {
    //if (true) {
      std::pair<double,double> r = 
        boost::math::tools::bracket_and_solve_root([this](double x) 
          {return particle_constraint_eqn(x);},
            fermi_energy_, factor, is_rising, tol, it);
      fermi_energy_ = r.first + 0.5*(r.second-r.first);
      if (it == maxit) {
        std::cout << " ** warning: fermi energy solve - iteraction exceeded\n";
      }
    }

    if (iteration_zero_) {
      bandwidth_zero_ = bandwidth_;
      fermi_energy_zero_ = fermi_energy_;
      metallic_zero_ = metallic_;
    }
    iteration_zero_ = false;
    std::cout << "metallic, e_F = "<<metallic_<<"  "<<fermi_energy_<< "\n";
    std::cout << "Bandwidth = " << bandwidth_ << "\n";

    // check
    /*double particle_sum = 0.0;
    double W_INV = 1.0/smear_width_;
    for (int i=0; i<ek_list_.size(); ++i) {
      double x = (ek_list_[idx_[i]]-fermi_energy_)*W_INV;
      double smear_wt = MethfesselPaxton_func(smear_func_order_,x);
      int k = qn_list_[idx_[i]].first;
      double degeneracy = std::round(blochbasis_.kweight(k));
      particle_sum += degeneracy * smear_wt;
    }
    std::cout << "particles = " << num_fill_particles_ 
              << " =? " << particle_sum << "\n";
    getchar();
    */
    

    // ground state
    double gs_energy = 0.0;
    double W_inv = 1.0/smear_width_;
    std::vector<int> shell_nmax(num_symm_kpoints_);
    for (auto& elem : shell_nmax) elem = -1; // invalid default value
    int k, n;
    if (metallic_) {
    //if (true) {
      RealMatrix wt_matrix(ek_list_.size(),kblock_dim_);
      for (int i=0; i<ek_list_.size(); ++i) {
        double x = (ek_list_[idx_[i]]-fermi_energy_)*W_inv;
        double smear_wt = MethfesselPaxton_func(smear_func_order_,x);
        if (std::abs(smear_wt)<1.0E-12) continue;
        std::tie(k,n) = qn_list_[idx_[i]];
        if (shell_nmax[k] < n) shell_nmax[k] = n;
        wt_matrix(k,n) = smear_wt;
        gs_energy += ek_list_[idx_[i]]*smear_wt;
      }
      kshells_up_.clear();
      for (int k=0; k<num_symm_kpoints_; ++k) {
        int nmax = shell_nmax[k];
        if (nmax != -1) {
          realArray1D smear_wts(nmax+1);
          for (int n=0; n<=nmax; ++n) smear_wts(n) = wt_matrix(k,n);
          //std::cout << smear_wts.transpose() << "\n";
          kshells_up_.push_back({k,0,nmax,smear_wts});
        }
      }
    }
    else {
      // insulator state
      int np = 0;
      for (int i=0; i<ek_list_.size(); ++i) {
        std::tie(k,n) = qn_list_[idx_[i]];
        if (shell_nmax[k] < n) shell_nmax[k] = n;
        int iw = std::nearbyint(blochbasis_.kweight(k));
        gs_energy += ek_list_[idx_[i]]*iw;
        np += iw;
        if (np == num_fill_particles_) break;
      }
      // store the filled k-shells
      kshells_up_.clear();
      for (int k=0; k<num_symm_kpoints_; ++k) {
        int nmax = shell_nmax[k];
        if (nmax != -1) {
          realArray1D smear_wts = realArray1D::Ones(nmax+1);
          kshells_up_.push_back({k,0,nmax,smear_wts});
        }
      }
    }
    // total energy per site
    gs_energy /= num_sites_;
    //std::cout << "Ground state energy (per site) = " << gs_energy << "\n";
    /* 
    for (int k=0; k<kshells_up_.size(); ++k) {
      std::cout << kshells_up_[k].k << " " << kshells_up_[k].nmin << "  "
          << kshells_up_[k].nmax << "\n";
    }
    getchar();
    */
//----------------------------------------------------
#ifdef ON

    // k-points to include in the ground state
    int top_filled_level = 0;
    std::vector<int> shell_nmax(num_kpoints_);
    for (auto& elem : shell_nmax) elem = -1; // invalid default value
    std::cout << "W = " << W << "\n";
    double particle_sum = 0.0;
    //double num_particles_d = static_cast<double>(num_fill_particles);
    for (int i=0; i<ek.size(); ++i) {
      double x = (ek[idx[i]]-fermi_energy_)*W_inv;
      double smear_wt = MethfesselPaxton_func(N_order,x);
      std::cout << i << "   " << smear_wt << "\n";
      int k = qn_list[idx[i]].first;
      double degeneracy = std::round(blochbasis_.kweight(k));
      particle_sum += degeneracy * smear_wt;
      // particle number constraint
      //if (std::abs(particle_sum-num_particles_d)<1.0E-8) break;
      // smearing weight cut-off
      //if (i==top_filled_level) break;
      //if (std::abs(smear_wt)<1.0E-15) break;

      int n = qn_list[idx[i]].second;
      if (shell_nmax[k] < n) shell_nmax[k] = n;
    }
    std::cout << "particles = " << num_fill_particles 
              << " =? " << particle_sum << "\n";
    getchar();

    top_filled_level = num_fill_particles-1;
    std::cout << "top level --> " << top_filled_level << "\n";
    for (int i=num_fill_particles-10; i<num_fill_particles+10; ++i) {
      std::cout << i << "  " << idx[i] << "  " << ek[idx[i]] << 
      "   " << qn_list[idx[i]].first << "--" << qn_list[idx[i]].second <<  "\n";
    }
    getchar();


    double degeneracy_tol = 1.0E-12;
    top_filled_level = num_fill_particles-1;
    fermi_energy_ = ek[idx[top_filled_level]];
    int num_valence_states = 1;
    int num_valence_particle = 1;
    // look upward in energy
    for (int i=top_filled_level+1; i<ek.size(); ++i) {
      if (std::abs(fermi_energy_-ek[idx[i]])>degeneracy_tol) break;
      num_valence_states++;
    }
    // look downward in energy
    if (num_valence_states>1) {
      for (int i=top_filled_level-1; i>=0; --i) {
        if (std::abs(fermi_energy_ - ek[idx[i]])>degeneracy_tol) break;
        num_valence_states++;
        num_valence_particle++;
      }
      // warn
      degeneracy_warning_ = true;
      if (degeneracy_warning_) {
        std::cout << ">> warning: Spinon groundstate degeneracy: " << num_valence_particle <<
          " particles in " << num_valence_states << " states." << "\n";
      }
    }
    /* 
      Filled k-shells. A k-shell is a group of energy levels having same 
      value of quantum number k.
    */
    // find 'nmax' values of filled k-shells
    //shell_nmax(num_kpoints_);
    int k, n;
    for (auto& elem : shell_nmax) elem = -1; // invalid default value
    for (int i=0; i<num_fill_particles; ++i) {
      int state = idx[i]; 
      std::tie(k,n) = qn_list[state];
      if (shell_nmax[k] < n) shell_nmax[k] = n;
    }
    // store the filled k-shells
    kshells_up_.clear();
    for (int k=0; k<num_kpoints_; ++k) {
      int nmax = shell_nmax[k];
      if (nmax != -1) kshells_up_.push_back({k,0,nmax});
    }
    /*
    for (unsigned k=0; k<kshells_up_.size(); ++k) {
      std::cout << kshells_up_[k].k << " " << kshells_up_[k].nmin << "  "
          << kshells_up_[k].nmax << "\n";
    }
    */


#endif 
//----------------------------------------------------

  }
  else {
    /* 
      No T.P (Time Reversal * Inversion) symmetry. 
      Need to consider both UP and DN states.
    */
    throw std::range_error("Spinon::construct_groundstate: case not implemented\n");
  }
}

double Spinon::particle_constraint_eqn(const double& mu)
{
  double particle_sum = 0.0;
  double W_inv = 1.0/smear_width_;
  for (int i=0; i<ek_list_.size(); ++i) {
    double x = (ek_list_[idx_[i]]-mu)*W_inv;
    double smear_wt = MethfesselPaxton_func(smear_func_order_,x);
    int k = qn_list_[idx_[i]].first;
    double degeneracy = std::round(blochbasis_.kweight(k));
    particle_sum += degeneracy * smear_wt;
  }
  double f = particle_sum-static_cast<double>(num_fill_particles_);
  //std::cout << "mu,  f = "<< mu << "   " << f << "\n";
  return f;
}

double Spinon::MethfesselPaxton_func(const int& N, const double& x)
{
  assert(N>=0);
  double S_0 = 0.5*(1.0-std::erf(x));
  // O-th order
  if (N == 0) return S_0;
  double A = -0.25;
  double exp_factor = std::exp(-x*x)/std::sqrt(pi());
  std::vector<double> Hermite(2*N);
  Hermite[0] = 1.0;
  Hermite[1] = 2.0*x;
  double S_N = A*Hermite[1];
  // 1-st order
  if (N ==1) return (S_0 + S_N*exp_factor);
  // Higher orders
  for (int n=2; n<2*N; ++n) {
    Hermite[n] = 2.0*(x*Hermite[n-1] - (n-1)*Hermite[n-2]);
  }
  for (int n=2; n<=N; ++n) {
    A *= -0.25/n;
    S_N += A*Hermite[2*n-1];
  }
  return S_0 + S_N*exp_factor;
}

double Spinon::MarzariVenderbilt_smear(const double& x)
{
  double a = std::sqrt(2.0);
  double b = 1.0/a;
  double x2 = (x-b)*(x-b);
  return std::exp(-x2)*(2.0-a*x)/pi();
}

void Spinon::set_particle_num(const input::Parameters& inputs)
{
  int nowarn;
  int particle_per_site = inputs.set_value("particle_per_site",1,nowarn);
  if (nowarn == 0) {
    num_spins_ = particle_per_site*num_sites_;
    //num_dnspins_ = num_spins_/2;
    //num_upspins_ = num_spins_-num_dnspins_;
    band_filling_ = 2.0*static_cast<double>(num_spins_)/num_total_states_;
    hole_doping_ = 1.0 - band_filling_;
    return;
  }
  /*
  if (Model::id()==model::model_id::PYROCHLORE) {
    num_spins_ = 5*num_sites_;
    num_dnspins_ = num_spins_/2;
    num_upspins_ = num_spins_-num_dnspins_;
    band_filling_ = 2.0*static_cast<double>(num_spins_)/num_total_states_;
    hole_doping_ = 1.0 - band_filling_;
    return;
  }
  */
  

  hole_doping_ = inputs.set_value("hole_doping", 0.0);
  if (std::abs(last_hole_doping_-hole_doping_)<1.0E-15) {
    // no change in hole doping, particle number remails same
    return;
  }
  last_hole_doping_ = hole_doping_;
  band_filling_ = 1.0-hole_doping_;
  num_spins_ = 0.5*static_cast<int>(std::round(band_filling_*num_total_states_));
  if (num_spins_<0 || num_spins_>num_total_states_) {
    throw std::range_error("Spinon::set_particle_num:: hole doping out-of-range");
  }
  // take even no of spins
  if (num_spins_%2 !=0 ) num_spins_ += 1;
  //num_dnspins_ = num_spins_/2;
  //num_upspins_ = num_spins_ - num_dnspins_;
  band_filling_ = 2.0*static_cast<double>(num_spins_)/num_total_states_;
  hole_doping_ = 1.0 - band_filling_;
  /*std::cout << "band_filling = " << band_filling_ << "\n";
  std::cout << "N_up = " << num_upspins_ << "\n";
  std::cout << "N_dn = " << num_dnspins_ << "\n";
  */
}

void Spinon::update_terms(void)
{
  update_unitcell_terms();
}

void Spinon::update_site_parameter(const std::string& pname, const double& pvalue)
{
  Model::update_parameter(pname, pvalue);
  for (unsigned i=0; i<usite_terms_.size(); ++i) 
    usite_terms_[i].eval_coupling_constant(Model::parameters(),Model::constants());
}

void Spinon::build_unitcell_terms(const lattice::LatticeGraph& graph)
{
  // take only quadratic & pairing terms
  int num_siteterms = 0;
  int num_bondterms = 0;
  for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
    if (sterm->qn_operator().is_quadratic() || sterm->qn_operator().is_pairing())
      num_siteterms++;
  }
  for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
    if (bterm->qn_operator().is_quadratic() || bterm->qn_operator().is_pairing()) {
      //std::cout << bterm->qn_operator().name() << "\n"; getchar();
      num_bondterms++;
    }
  }
  usite_terms_.resize(num_siteterms);
  ubond_terms_.resize(num_bondterms);
  unsigned i = 0;
  for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
    if (sterm->qn_operator().is_quadratic() || sterm->qn_operator().is_pairing()) {
      usite_terms_[i].build_siteterm(*sterm, graph);
      i++;
    }
  }
  i = 0;
  for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
    if (bterm->qn_operator().is_quadratic() || bterm->qn_operator().is_pairing()) {
      ubond_terms_[i].build_bondterm(*bterm, graph);
      i++;
    }
  }
}

void Spinon::construct_kspace_block(const MF_Params& mf_params, const Vector3d& kvec)
{
  work.setZero(); 
  // bond terms
  for (auto& term : ubond_terms_) {
    if (term.qn_operator().is_quadratic() && term.qn_operator().spin_up()) {
      for (int i=0; i<term.num_out_bonds(); ++i) {
        Vector3d delta = term.bond_vector(i);
        //std::cout << ">>HACK: remormalization of spinon parameters OFF\n";
        work += term.coeff_matrix(i) * std::exp(ii()*kvec.dot(delta));
      }
    }
  }
  // add hermitian conjugate part
  quadratic_block_up_ = work + work.adjoint();

  // site terms 
  for (const auto& term : usite_terms_) {
    if (term.qn_operator().spin_up()) {
      quadratic_block_up_ += term.coeff_matrix();
      //std::cout << term.coeff_matrix() << "\n"; getchar();
    }
  }
  if (!no_spinon_lambda_) {
    for (int i=0; i<mf_params.num_sites(); ++i) {
      int site_dim = mf_params.site(i).dim();
      realArray1D lm_params = mf_params.site(i).lm_params(); 
      //realArray1D lm_params_noint = mf_params.site(i).lm_params_noint(); 
      for (int m=0; m<site_dim; ++m) {
        auto n = mf_params.site(i).state_indices()[m];
        quadratic_block_up_(n,n) += -lm_params(m);
        //std::cout << "lamba["<<n<<"] = " << lm_params[m] << "\n"; 
      }
    }
  }
  //std::cout << "e0 =\n" << quadratic_block_up_.diagonal().transpose() << "\n"; getchar();
  // site terms
  //std::cout << "k = " << kvec.transpose() << "\n";
  //std::cout << "hk =\n" << quadratic_block_up_ << "\n"; getchar();
}

void UnitcellTerm::update_bondterm_cc(const int& term_id, const MF_Params& mf_params)
{
  for (auto& M : coeff_matrices_) M.setZero();
  for (int i=0; i<mf_params.num_bonds(); ++i) {
    int id = mf_params.bond(i).vector_id();
    //std::cout << "id = " << id << "\n";
    int src = mf_params.bond(i).src();
    int tgt = mf_params.bond(i).tgt();
    int rows = mf_params.site(src).dim();
    int cols = mf_params.site(tgt).dim();
    auto renorm_cc = mf_params.bond(i).boson_renormed_cc(term_id);
    //renorm_bond_couplings_[i] = mf_params.bond(i).spinon_renormed_cc(0);
    //std::cout << renorm_cc << "\n" << "\n";
    for (int m=0; m<rows; ++m) {
      auto ii = mf_params.site(src).state_indices()[m];
      for (int n=0; n<cols; ++n) {
        auto jj = mf_params.site(tgt).state_indices()[n];
        coeff_matrices_[id](ii,jj) += renorm_cc(m,n);
      }
    }
  }
}

void UnitcellTerm::update_siteterm_cc(const MF_Params& mf_params)
{
  for (auto& M : coeff_matrices_) M.setZero();
  //coeff_matrices_[0].setZero();
  for (int i=0; i<mf_params.num_sites(); ++i) {
    int rows = mf_params.site(i).dim();
    int cols = mf_params.site(i).dim();
    auto renorm_soc = mf_params.site(i).boson_renormed_soc();
    //renorm_bond_couplings_[i] = mf_params.bond(i).spinon_renormed_cc(0);
    //std::cout << renorm_cc << "\n" << "\n";
    for (int m=0; m<rows; ++m) {
      auto ii = mf_params.site(i).state_indices()[m];
      for (int n=0; n<cols; ++n) {
        auto jj = mf_params.site(i).state_indices()[n];
        coeff_matrices_[0](ii,jj) += renorm_soc(m,n);
        //std::cout << "coeff = " << renorm_soc(m,n) << "\n";
      }
    }
  }
  //std::cout << "coeff =\n" << coeff_matrices_[0] << "\n"; getchar();
}


void Spinon::update_unitcell_terms(void)
{
  for (unsigned i=0; i<ubond_terms_.size(); ++i) 
    ubond_terms_[i].eval_coupling_constant(Model::parameters(),Model::constants());
  for (unsigned i=0; i<usite_terms_.size(); ++i) 
    usite_terms_[i].eval_coupling_constant(Model::parameters(),Model::constants());
}

/* Write a bond term like,
 H = \sum_{Ia,Jb}c^{\dag}_{Ia} t_{Ia,Jb} c_{Jb}
 for lattices with multiple sites per unit cell as
 H = \sum_{I,delta} Psi^{\dag}_{I} M^{\delta} Psi_{I+delta}
 Assumption: 'site's in the Graph are numbered contigously. For
 sites in the unitcell, the 'sl number' is same as 'uid'.
*/
void UnitcellTerm::build_bondterm(const model::HamiltonianTerm& hamterm,
  const lattice::LatticeGraph& graph)
{
  num_basis_sites_ = graph.lattice().num_basis_sites();
  dim_ = graph.lattice().num_basis_orbitals();
  lattice::LatticeGraph::out_edge_iterator ei, ei_end;
  // get number of unique 'cell bond vectors'
  num_out_bonds_ = 0;
  for (unsigned i=0; i<num_basis_sites_; ++i) {
    //std::cout << i << "\n";
    for (std::tie(ei, ei_end)=graph.out_bonds(i); ei!=ei_end; ++ei) {
      unsigned id = graph.vector_id(ei);
      if (id > num_out_bonds_) num_out_bonds_ = id;
    }
  }
  num_out_bonds_++;
  //std::cout << "num_out_bonds_ = " << num_out_bonds_ << "\n";
  bond_vectors_.resize(num_out_bonds_);
  coeff_matrices_.resize(num_out_bonds_);
  for (auto& M : coeff_matrices_) {
    M.resize(dim_, dim_);
    M.setZero();
  }
  expr_matrices_.clear();
  expr_matrices_.resize(num_out_bonds_);
  for (auto& M : expr_matrices_) {
    M.resize(dim_,dim_);
  }
  // initialize expression matrices
  for (int id=0; id<num_out_bonds_; ++id) {
    for (int i=0; i<dim_; ++i) {
      for (int j=0; j<dim_; ++j) {
        expr_matrices_[id](i,j) = "0";
      }
    }
  }

  // operator
  name_ = hamterm.name(); 
  op_ = hamterm.qn_operator();
  // build the matrices (for each 'bond vector')
  for (unsigned i=0; i<num_basis_sites_; ++i) {
    for (std::tie(ei,ei_end)=graph.out_bonds(i); ei!=ei_end; ++ei) {
      unsigned id = graph.vector_id(ei);
      unsigned t = graph.target(ei);
      unsigned j = graph.site_uid(t);
      unsigned btype = graph.bond_type(ei);

      strMatrix expr_mat = hamterm.coupling_expr(btype);
      ComplexMatrix coeff_mat = hamterm.coupling(btype);
      int rows = graph.site_dim(i);
      int cols = graph.site_dim(j);
      if ((expr_mat.rows()!=rows) && (expr_mat.cols()!=cols)) {
        throw std::runtime_error("Dimension mismatch in the 'bondterm'");
      }
      for (int m=0; m<rows; ++m) {
        for (int n=0; n<cols; ++n) {
          int p = graph.lattice().basis_index_number(i,m);
          int q = graph.lattice().basis_index_number(j,n);
          std::string cc_expr(expr_mat(m,n));
          if (cc_expr.size()>0) {
            expr_matrices_[id](p,q) += "+("+cc_expr+")";
          }
          // values
          coeff_matrices_[id](p,q) += coeff_mat(m,n);
          //std::cout << id << " " << coeff_matrices_[id](p,q) << "\n";
        }
      }
      bond_vectors_[id] = graph.vector(ei);
      //std::cout << i << " " << j << "\n"; 
      //std::cout << id << ": " << bond_vectors_[id].transpose() << "\n"; 
      //getchar();
    }
  }
  /*std::cout <<"------"<<op_.name()<<"-------\n"; getchar();
  for (int id=0; id<num_out_bonds_; ++id) {
    std::cout << "delta = " << id << "\n";
    for (int i=0; i<dim_; ++i) {
      for (int j=0; j<dim_; ++j) {
        std::cout<<i<<", "<<j<<"= "<<expr_matrices_[id][i][j]<<"\n";
      }
    }
    std::cout<<"\n";
  }
  getchar();
  */
}

void UnitcellTerm::build_siteterm(const model::HamiltonianTerm& hamterm,
  const lattice::LatticeGraph& graph)
{
  num_basis_sites_ = graph.lattice().num_basis_sites();
  dim_ = graph.lattice().num_basis_orbitals();
  num_out_bonds_ = 1; // dummy, no real contrib
  bond_vectors_.resize(1);
  coeff_matrices_.resize(1);
  coeff_matrices_[0].resize(dim_,dim_);
  coeff_matrices_[0].setZero();
  expr_matrices_.resize(1);
  expr_matrices_[0].resize(dim_,dim_);
  // initialize expression matrix
  for (int i=0; i<dim_; ++i) {
    for (int j=0; j<dim_; ++j) {
      expr_matrices_[0](i,j) = "0";
    }
  }

  // operator
  name_ = hamterm.name(); 
  op_ = hamterm.qn_operator();
  // build the matrix 
  for (unsigned i=0; i<num_basis_sites_; ++i) {
    unsigned stype = graph.site_type(i);
    //coeff_matrices_[0](i,i) = hamterm.coupling(stype);
    // expression
    strMatrix expr_mat = hamterm.coupling_expr(stype);
    ComplexMatrix coeff_mat = hamterm.coupling(stype);
    int rows = graph.site_dim(i);
    int cols = graph.site_dim(i);
    if (expr_mat.cols()!=cols) {
      throw std::runtime_error("Dimension mismatch in the 'siteterm'");
    }
    // diagonal site term
    if (expr_mat.rows()==1) {
      for (int j=0; j<cols; ++j) {
        int n = graph.lattice().basis_index_number(i,j);
        coeff_matrices_[0](n,n) = coeff_mat(0,j);
        std::string cc_expr(expr_mat(0,j));
        boost::trim(cc_expr);
        if (cc_expr.size()>0) {
          expr_matrices_[0](n,n) += "+("+cc_expr+")";
        }
        //n++;
      }
    }
    else {
      // off-diagonal site term
      if (expr_mat.rows()!=rows) {
        throw std::runtime_error("Dimension mismatch in the 'siteterm'");
      }
      for (int m=0; m<rows; ++m) {
        for (int n=0; n<cols; ++n) {
          int p = graph.lattice().basis_index_number(i,m);
          int q = graph.lattice().basis_index_number(i,n);
          std::string cc_expr(expr_mat(m,n));
          if (cc_expr.size()>0) {
            expr_matrices_[0](p,q) += "+("+cc_expr+")";
          }
          // values
          coeff_matrices_[0](p,q) += coeff_mat(m,n);
        }
      }
    }
  }
  bond_vectors_[0] = Vector3d(0.0,0.0,0.0);
}

/*
void UnitcellTerm::eval_coupling_constant(const model::ModelParams& pvals, const model::ModelParams& cvals)
{
  expr::Expression expr;
  expr::Expression::variables vars;
  for (const auto& p : pvals) {
    vars[p.first] = p.second;
    //std::cout << p.first << " = " << p.second << "\n"; getchar();
  }
  for (const auto& c : cvals) vars[c.first] = c.second;
  try { 
    for (unsigned n=0; n<num_out_bonds_; ++n) {
      for (unsigned i=0; i<dim_; ++i) {
        for (unsigned j=0; j<dim_; ++j) {
          std::string cc_expr(expr_matrices_[n][i][j]);
          if (cc_expr.size()>0) {
            coeff_matrices_[n](i,j) = expr.evaluate(cc_expr, vars); 
            //std::cout << "cc = " << coeff_matrices_[n](i,j) << "\n"; getchar();
          }
          else
            coeff_matrices_[n](i,j) = 0.0;
        }
      }
    }
  }
  catch (std::exception& e) 
  { 
    std::string msg = "UnitcellTerm::evaluate_coupling_constant:\n" + std::string(e.what());
    throw std::runtime_error(msg);
  }
}*/

void UnitcellTerm::eval_coupling_constant(const model::ModelParams& pvals, const model::ModelParams& cvals)
{
  expr::ComplexExpr expr;
  for (const auto& p : pvals) expr.add_var(p.first, p.second);
  for (const auto& c : cvals) expr.add_var(c.first, c.second);
  try { 
    for (unsigned n=0; n<num_out_bonds_; ++n) {
      for (unsigned i=0; i<dim_; ++i) {
        for (unsigned j=0; j<dim_; ++j) {
          std::string cc_expr(expr_matrices_[n](i,j));
          if (cc_expr.size()>0) {
            expr.set_expr(cc_expr);
            coeff_matrices_[n](i,j) = expr.evaluate(); 
            //std::cout << "cc = " << cc_expr << "\n"; getchar();
          }
          else coeff_matrices_[n](i,j) = 0.0;
        }
      }
    }
  }
  catch (std::exception& e) 
  { 
    std::string msg = "UnitcellTerm::evaluate_coupling_constant:\n" + std::string(e.what());
    throw std::runtime_error(msg);
  }
}

/*void Spinon::check_xml(void)
{
  std::cout << "Checking XML parser\n";
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("model.xml", pugi::parse_trim_pcdata);
  //std::cout << "Load result: " << result.description() << ", mesh name: " << "\n"; 
  pugi::xml_node model = doc.child("model");
  //std::cout << model.child("parameter").attribute("default").value() << std::endl;
  for (pugi::xml_node p = model.child("parameter"); p; p = p.next_sibling())
  {
    std::cout << "Parameter: ";
    for (pugi::xml_attribute attr = p.first_attribute(); attr; attr = attr.next_attribute())
    {
      std::cout << attr.name() << " = " << attr.as_double() << "\n";
    }
    std::cout << std::endl;
  }
  std::cout << "---------------------------------\n\n";
}*/



} // end namespace ssmf











