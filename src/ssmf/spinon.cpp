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
#include "../expression/complex_expression.h"
//#include "../expression/expression.h"

namespace srmf {

Spinon::Spinon(const input::Parameters& inputs, const model::Hamiltonian& model, 
    const lattice::LatticeGraph& graph, const SB_Params& srparams)
: model::Hamiltonian(model)
, blochbasis_(graph)
{
  // add mean-field terms to the spinon model
  if (graph.lattice().id()==lattice::lattice_id::SQUARE_2BAND) {
  }

  if (graph.lattice().id()==lattice::lattice_id::PYROCHLORE_3D) {
  }
  have_TP_symmetry_ = model.have_TP_symmetry();
  SO_coupling_ = model.is_spinorbit_coupled();
  if (SO_coupling_) spin_multiply_ = 1;
  else spin_multiply_ = 2;
  num_sites_ = graph.num_sites();
  num_bonds_ = graph.num_bonds();
  num_kpoints_ = blochbasis_.num_kpoints();
  kblock_dim_ = blochbasis_.subspace_dimension();
  num_total_states_ = spin_multiply_ * num_kpoints_ * kblock_dim_;
  num_basis_sites_ = graph.lattice().num_basis_sites();
  //dim_ = graph.lattice().num_basis_orbitals();
  quadratic_block_up_.resize(kblock_dim_,kblock_dim_);
  pairing_block_.resize(kblock_dim_,kblock_dim_);
  orbital_en_.resize(spin_multiply_ * kblock_dim_);
  orbital_en_shifted_.resize(kblock_dim_);
  work.resize(kblock_dim_,kblock_dim_);
  build_unitcell_terms(graph);
  set_particle_num(inputs);
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
  return 0;
}

void Spinon::update(const input::Parameters& inputs)
{
  Model::update_parameters(inputs);
  update_terms();
  orbital_en_.setZero();
  for (const auto& term : usite_terms_) {
    //std::cout << " --------- here --------\n";
    if (term.qn_operator().spin_up()) {
      auto e0 = term.coeff_matrix().diagonal().real();
      if (spin_multiply_==2) {
        for (int i=0; i<kblock_dim_; ++i) {
          orbital_en_(i) = e0(i);
          orbital_en_(kblock_dim_+i) = e0(i);
        }
      }
      else orbital_en_ = e0;
      //std::cout << "e0 =" << orbital_en_.transpose() << "\n"; getchar();
    }
  }
  orbital_en_shifted_ = orbital_en_;
}

void Spinon::solve(const lattice::LatticeGraph& graph, SB_Params& srparams)
{
  construct_groundstate(srparams);
  compute_averages(graph, srparams);
}

void Spinon::compute_averages(const lattice::LatticeGraph& graph, SB_Params& srparams)
{
  //Eigen::VectorXcd amplitude_vec1, amplitude_vec2;
  Eigen::Matrix<std::complex<double>,1,Eigen::Dynamic> amplitude_vec1, amplitude_vec2;
  for (int i=0; i<srparams.num_sites(); ++i) {
    srparams.site(i).spinon_density().setZero();
  }
  for (int i=0; i<srparams.bonds().size(); ++i) {
    srparams.bond(i).spinon_ke().setZero();
  }
  if (have_TP_symmetry_) {
    for (int i=0; i<kshells_up_.size(); ++i) {
      unsigned k = kshells_up_[i].k;
      Vector3d kvec = blochbasis_.kvector(k);
      construct_kspace_block(srparams, kvec);
      es_k_up_.compute(quadratic_spinup_block());

      // site density
      unsigned nmin = kshells_up_[i].nmin;
      unsigned nmax = kshells_up_[i].nmax;
      unsigned nbands = nmax-nmin+1;
      Eigen::VectorXcd eigvec_wt(nbands);
      for (int j=0; j<srparams.num_sites(); ++j) {
        unsigned site_dim = srparams.site(j).dim();
        realArray1D n_avg(spin_multiply_*site_dim); 
        for (unsigned m=0; m<site_dim; ++m) {
          auto ii = srparams.site(j).state_indices()[m];
          double norm = 0.0;
          for (unsigned band=nmin; band<=nmax; ++band) {
            norm += std::norm(es_k_up_.eigenvectors().row(ii)[band]);
          }
          n_avg(m) = norm;
        }
        srparams.site(j).spinon_density() += n_avg;
      }

      // bond averages
      amplitude_vec1.resize(nbands);
      amplitude_vec2.resize(nbands);
      for (int j=0; j<srparams.num_bonds(); ++j) {
        Vector3d delta = srparams.bond(j).vector();
        std::complex<double> exp_kdotr = std::exp(ii()*kvec.dot(delta));
        unsigned src = srparams.bond(j).src();
        unsigned tgt = srparams.bond(j).tgt();
        unsigned rows = srparams.site(src).dim();
        unsigned cols = srparams.site(tgt).dim();
        //sr_bond bond(srparams.bonds()[j]);
        cmplArray2D ke_matrix=cmplArray2D::Zero(spin_multiply_*rows, spin_multiply_*cols);
        for (unsigned m=0; m<rows; ++m) {
          auto ii = srparams.site(src).state_indices()[m];
          amplitude_vec1 = es_k_up_.eigenvectors().block(ii,nmin,1,nbands);
          //std::cout << "amplitude_vec1 = " << amplitude_vec1.transpose() << "\n"; getchar();
          for (unsigned n=0; n<cols; ++n) {
            auto jj = srparams.site(tgt).state_indices()[n];
            amplitude_vec2 = es_k_up_.eigenvectors().block(jj,nmin,1,nbands);
            std::complex<double> chi_sum = amplitude_vec1.dot(amplitude_vec2);
            ke_matrix(m,n) = exp_kdotr * chi_sum;
            //std::complex<double> chi_sum(0.0);
            //for (unsigned band=nmin; band<=nmax; ++band) {
            //  chi_sum += exp_kdotr * std::conj(es_k_up_.eigenvectors().row(ii)[band]) * 
            //         es_k_up_.eigenvectors().row(jj)[band];
            //}
            //ke_matrix(m,n) = chi_sum;
          }
        }
        srparams.bond(j).spinon_ke() += ke_matrix;
      }
    }
  }
  else {
    /* 
      No T.P (Time Reversal * Inversion) symmetry. 
      Need to consider both UP and DN states.
    */
    throw std::range_error("Spinon::compute_averages: case not implemented\n");
  }

  // final site density 
  for (int i=0; i<srparams.num_sites(); ++i) {
    //srparams.site(i).spinon_density() /= num_kpoints_;
    realArray1D n_avg = srparams.site(i).spinon_density()/num_kpoints_;
    if (spin_multiply_==2) {
      // for DN spins, same as UP spins
      auto dim = srparams.site(i).dim();
      for (int n=0; n<dim; ++n) n_avg(dim+n) = n_avg(n);
    }
    srparams.site(i).spinon_density() = n_avg;
    // print
    /*
    std::cout<<"site-"<<i<<":"<<"\n";
    for (int m=0; m<srparams.site(i).dim(); ++m) {
      std::cout << "<n_up>["<<m<<"] = " << srparams.site(i).spinon_density()[m] << "\n";
    }
    std::cout << "\n";
    */
  }

  // final KE average 
  for (int i=0; i<srparams.num_bonds(); ++i) {
    cmplArray2D ke_matrix = srparams.bond(i).spinon_ke()/num_kpoints_;
    if (spin_multiply_==2) {
      int m = ke_matrix.rows()/2;
      int n = ke_matrix.cols()/2;
      ke_matrix.block(m,n,m,n) = ke_matrix.block(0,0,m,n);
      ke_matrix.block(0,n,m,n) = cmplArray2D::Zero(m,n);
      ke_matrix.block(m,0,m,n) = cmplArray2D::Zero(m,n);
    }
    srparams.bond(i).spinon_ke() = ke_matrix;
    srparams.bond(i).set_spinon_renormalization();
    //std::cout << "ke_matrix =\n" << ke_matrix << "\n"; getchar();
    // print
    /*std::cout<<"bond-"<<i<<":"<<"\n";
    for (int m=0; m<srparams.bond(i).spinon_ke().rows(); ++m) {
      for (int n=0; n<srparams.bond(i).spinon_ke().cols(); ++n) {
        std::cout << "chi["<<m<<","<<n<<"] = " << srparams.bond(i).spinon_ke()(m,n) << "\n";
      }
    }
    std::cout << "\n";
    */
  }
  //std::cout << "Exiting at Spinon::compute_averages\n"; exit(0);
}

void Spinon::construct_groundstate(const SB_Params& srparams)
{
  if (have_TP_symmetry_) {
    /* Has T.P (Time Reversal * Inversion) symmetry. 
       So we have e_k(up) = e_k(dn).
    */
    std::vector<std::pair<unsigned,unsigned>> qn_list; // list of (k,n)
    std::vector<double> ek;
    for (unsigned k=0; k<num_kpoints_; ++k) {
      Vector3d kvec = blochbasis_.kvector(k);
      construct_kspace_block(srparams, kvec);
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
    double degeneracy_tol = 1.0E-12;
    int top_filled_level = num_upspins_-1;
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
    for (int i=0; i<num_upspins_; ++i) {
      int state = idx[i]; 
      std::tie(k,n) = qn_list[state];
      if (shell_nmax[k] < n) shell_nmax[k] = n;
    }
    // store the filled k-shells
    kshells_up_.clear();
    for (unsigned k=0; k<num_kpoints_; ++k) {
      int nmax = shell_nmax[k];
      if (nmax != -1) kshells_up_.push_back({k,0,static_cast<unsigned>(nmax)});
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

void Spinon::set_particle_num(const input::Parameters& inputs)
{
  hole_doping_ = inputs.set_value("hole_doping", 0.0);
  if (std::abs(last_hole_doping_-hole_doping_)<1.0E-15) {
    // no change in hole doping, particle number remails same
    return;
  }
  last_hole_doping_ = hole_doping_;
  band_filling_ = 1.0-hole_doping_;
  int num_up_states = 0.5*static_cast<int>(num_total_states_); 

  // band_filling = 1 implies 'half-filling' conventionally.
  bool pairing_type_ = false;
  if (pairing_type_) {
    int n = 0.5*static_cast<int>(std::round(band_filling_*num_up_states));
    if (n<0 || n>num_up_states) throw std::range_error("Wavefunction:: hole doping out-of-range");
    num_upspins_ = static_cast<unsigned>(n);
    num_dnspins_ = num_upspins_;
    num_spins_ = num_upspins_ + num_dnspins_;
    band_filling_ = static_cast<double>(2*n)/num_up_states;
  }
  else{
    int n = static_cast<int>(std::round(band_filling_*num_up_states));
    if (n<0 || n>2*num_up_states) throw std::range_error("Wavefunction:: hole doping out-of-range");
    num_spins_ = static_cast<unsigned>(n);
    num_dnspins_ = num_spins_/2;
    num_upspins_ = num_spins_ - num_dnspins_;
    band_filling_ = static_cast<double>(n)/num_up_states;
  }
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

void Spinon::construct_kspace_block(const SB_Params& srparams, const Vector3d& kvec)
{
  work.setZero(); 
  pairing_block_.setZero();
  //work2 = Matrix::Zero(dim_,dim_);
  // bond terms
  //for (const auto& term : uc_bondterms_) {
  for (const auto& term : ubond_terms_) {
    if (term.qn_operator().is_quadratic() && term.qn_operator().spin_up()) {
      //std::cout << term.num_out_bonds() << "\n";
      for (unsigned i=0; i<term.num_out_bonds(); ++i) {
        Vector3d delta = term.bond_vector(i);
        // renormalization by bosonic order parameter
        cmplArray2D boson_bond_avg = srparams.bond(i).boson_ke();
        if (SO_coupling_) {
          cmplArray2D term_matrix = term.coeff_matrix(i).array() * boson_bond_avg  
                                    * std::exp(ii()*kvec.dot(delta));
          work += term_matrix.matrix();
        }
        else {
          int m = boson_bond_avg.rows()/2;
          int n = boson_bond_avg.cols()/2;
          cmplArray2D term_matrix = term.coeff_matrix(i).array()  
                                    * boson_bond_avg.block(0,0,m,n)  
                                    * std::exp(ii()*kvec.dot(delta));
          work += term_matrix.matrix();
          //std::cout << "t = " << term.coeff_matrix(i).array() << "\n";
          //std::cout << "tB = " << term.coeff_matrix(i).array()* boson_bond_avg.block(0,0,m,n) << "\n";
          //getchar();
        }
        //std::cout << "delta = " << delta.transpose() << "\n"; getchar();
      }
    }
    if (term.qn_operator().is_pairing()) {
      for (unsigned i=0; i<term.num_out_bonds(); ++i) {
        Vector3d delta = term.bond_vector(i);
        pairing_block_ += term.coeff_matrix(i) * std::exp(ii()*kvec.dot(delta));
      }
    }
  }
  // add hermitian conjugate part
  quadratic_block_up_ = work + work.adjoint();


  //orbital_en_shifted_(0) = 0.4421;
  //orbital_en_shifted_(1) = -0.4421;

  // site terms 
  for (const auto& term : usite_terms_) {
    //std::cout << " --------- here --------\n";
    if (term.qn_operator().spin_up()) {
      quadratic_block_up_ += term.coeff_matrix();
      //for (int i=0; i<kblock_dim_; ++i) quadratic_block_up_(i,i) += orbital_en_shifted_(i);
      //if (v==1) {std::cout << " sterm =" << term.coeff_matrix().diagonal().transpose() << "\n"; getchar();}
    }
  }
  //std::cout << "e0 =\n" << quadratic_block_up_.diagonal().transpose() << "\n"; 
  // LM parameters/chemical potential
  // site density
  /*
  realArray1D e(2);
  e(0) = 3.73115;
  e(1) = -3.73115;
  */
  for (int i=0; i<srparams.num_sites(); ++i) {
    unsigned site_dim = srparams.site(i).dim();
    realArray1D lm_params = srparams.site(i).lm_params(); 
    realArray1D lm_params_noint = srparams.site(i).lm_params_noint(); 
    for (unsigned m=0; m<site_dim; ++m) {
      auto n = srparams.site(i).state_indices()[m];
      quadratic_block_up_(n,n) += -lm_params(m);
      //std::cout << "lamba["<<n<<"] = " << lm_params[m] << "\n"; 
    }
  }
  //std::cout << "e0 =\n" << quadratic_block_up_.diagonal().transpose() << "\n"; getchar();

  //quadratic_block_up_ += work1.adjoint();
  //pairing_block_ = work2;
  //pairing_block_ += work2.adjoint();
  // site terms
  //std::cout << "k = " << kvec.transpose() << "\n";
  //std::cout << "hk =\n" << quadratic_block_up_ << "\n"; getchar();
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
          int p = graph.lattice().basis_index_number(i, m);
          int q = graph.lattice().basis_index_number(j, n);
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
      //std::cout<<bond_vectors_[id].transpose()<<"\n"; getchar();
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
  op_ = hamterm.qn_operator();
  // build the matrix 
  for (unsigned i=0; i<num_basis_sites_; ++i) {
    unsigned stype = graph.site_type(i);
    //coeff_matrices_[0](i,i) = hamterm.coupling(stype);
    // expression
    strMatrix expr_mat = hamterm.coupling_expr(stype);
    ComplexMatrix coeff_mat = hamterm.coupling(stype);
    if ((expr_mat.rows()!=1) && (expr_mat.cols()!=graph.site_dim(i))) {
      throw std::runtime_error("Dimension mismatch in the 'siteterm'");
    }
    for (int j=0; j<graph.site_dim(i); ++j) {
      unsigned n = graph.lattice().basis_index_number(i, j);
      coeff_matrices_[0](n,n) = coeff_mat(0,j);
      std::string cc_expr(expr_mat(0,j));
      boost::trim(cc_expr);
      if (cc_expr.size()>0) {
        expr_matrices_[0](n,n) += "+("+cc_expr+")";
      }
      n++;
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



} // end namespace srmf











