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
    const lattice::LatticeGraph& graph, const SR_Params& srparams)
: model::Hamiltonian(model)
, blochbasis_(graph)
{
  // add mean-field terms to the spinon model
  if (graph.lattice().id()==lattice::lattice_id::SQUARE_2BAND) {
  }

  if (graph.lattice().id()==lattice::lattice_id::PYROCHLORE_3D) {
  }

  num_sites_ = graph.num_sites();
  num_bonds_ = graph.num_bonds();
  num_kpoints_ = blochbasis_.num_kpoints();
  kblock_dim_ = blochbasis_.subspace_dimension();
  num_total_states_ = 2 * num_kpoints_ * kblock_dim_;
  num_basis_sites_ = graph.lattice().num_basis_sites();
  //dim_ = graph.lattice().num_basis_orbitals();
  quadratic_block_up_.resize(kblock_dim_,kblock_dim_);
  pairing_block_.resize(kblock_dim_,kblock_dim_);
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

void Spinon::solve(const lattice::LatticeGraph& graph, SR_Params& srparams)
{
  construct_groundstate();
  compute_averages(graph, srparams);
}

void Spinon::compute_averages(const lattice::LatticeGraph& graph, SR_Params& srparams)
{
  std::vector<double> nup_avg(kblock_dim_,0.0);
  Eigen::VectorXcd amplitude_vec1, amplitude_vec2;
  ComplexMatrix chi_mat = ComplexMatrix::Zero(2,2);
  if (have_TP_symmetry_) {
    for (int i=0; i<kshells_up_.size(); ++i) {
      unsigned k = kshells_up_[i].k;
      Vector3d kvec = blochbasis_.kvector(k);
      construct_kspace_block(kvec);
      es_k_up_.compute(quadratic_spinup_block());

      // site density
      unsigned nmin = kshells_up_[i].nmin;
      unsigned nmax = kshells_up_[i].nmax;
      Eigen::VectorXcd eigvec_wt(nmax-nmin+1);
      for (unsigned p=0; p<kblock_dim_; ++p) {
        for (unsigned j=nmin; j<=nmax; ++j)
          eigvec_wt[j] = es_k_up_.eigenvectors().row(p)[j];
        nup_avg[p] += eigvec_wt.squaredNorm();
      }

      // bond averages
      amplitude_vec1.resize(nmax-nmin+1);
      amplitude_vec2.resize(nmax-nmin+1);
      //for (int j=0; j<srparams.bonds().size(); ++j) {
      for (int j=0; j<1; ++j) {
        sr_bond bond(srparams.bonds()[j]);
        Vector3d delta = bond.vector();
        std::complex<double> exp_kdotr = std::exp(ii()*kvec.dot(delta));
        for (unsigned m=0; m<bond.src_state_indices().size(); ++m) {
          auto ii = bond.src_state_indices()[m];
          for (unsigned n=0; n<bond.tgt_state_indices().size(); ++n) {
            auto jj = bond.tgt_state_indices()[n];
            std::complex<double> chi_sum(0.0);
            for (unsigned band=nmin; band<=nmax; ++band) {
              chi_sum += exp_kdotr * std::conj(es_k_up_.eigenvectors().row(ii)[band]) * 
                     es_k_up_.eigenvectors().row(jj)[band];
            }
            chi_mat(m,n) += chi_sum;
            //std::cout << "(m,n) = " << m << ", " << n << "\n";
            //std::cout << "chi = " << chi << "\n";
          }
        }
        //getchar();
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
  for (unsigned p=0; p<kblock_dim_; ++p) {
    nup_avg[p] /= num_kpoints_;
    std::cout << "<n_up>["<<p<<"] = " << nup_avg[p] << "\n";
  }

  // final KE average 
  chi_mat /= num_kpoints_;
  for (int m=0; m<chi_mat.rows(); ++m) {
    for (int n=0; n<chi_mat.cols(); ++n) {
      std::cout << "chi["<<m<<","<<n<<"] = " << chi_mat(m,n) << "\n";
    }
  }


  std::cout << "Exiting at Spinon::compute_averages\n"; exit(0);
}

void Spinon::construct_groundstate(void)
{
  if (have_TP_symmetry_) {
    /* Has T.P (Time Reversal * Inversion) symmetry. 
       So we have e_k(up) = e_k(dn).
    */
    std::vector<std::pair<unsigned,unsigned>> qn_list; // list of (k,n)
    std::vector<double> ek;
    for (unsigned k=0; k<num_kpoints_; ++k) {
      Vector3d kvec = blochbasis_.kvector(k);
      construct_kspace_block(kvec);
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
      std::cout << ">> warning: Spinon groundstate degeneracy: " << 2*num_valence_particle <<
        " particles in " << 2*num_degen_states << " states." << "\n";
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

void Spinon::update(const input::Parameters& inputs)
{
  Model::update_parameters(inputs);
  update_terms();
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

void Spinon::construct_kspace_block(const Vector3d& kvec)
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
        work += term.coeff_matrix(i) * std::exp(ii()*kvec.dot(delta));
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
  // site terms 
  //for (const auto& term : uc_siteterms_) {
  for (const auto& term : usite_terms_) {
    //std::cout << " --------- here --------\n";
    if (term.qn_operator().spin_up()) {
      quadratic_block_up_ += term.coeff_matrix();
      //std::cout << " sterm =" << term.coeff_matrix() << "\n"; getchar();
    }
  }
  //quadratic_block_up_ += work1.adjoint();
  //pairing_block_ = work2;
  //pairing_block_ += work2.adjoint();
  // site terms
  //std::cout << "ek = " << quadratic_block_up_ << "\n"; getchar();
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











