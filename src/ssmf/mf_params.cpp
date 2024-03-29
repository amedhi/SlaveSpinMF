/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-03-12 12:20:33
* @Last Modified by:   Amal Medhi
* @Last Modified time: 2022-07-13 12:55:57
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "mf_params.h"

namespace ssmf {

//-------------------------------MF Site-----------------------------
MF_Site::MF_Site() 
  : type_{0}, dim_{0} 
{
  spin_orbitals_.clear(); 
  state_indices_.clear();
  connected_bonds_.clear();
  bond_outgoing_.clear();
} 

MF_Site::MF_Site(const int& type, const int& dim, const idx_list& state_indices)
  : type_{type}, dim_{dim}, state_indices_{state_indices}
{
  connected_bonds_.clear();
  bond_outgoing_.clear();
  /*
    'spin_orbitals' Note:  the spin-orbitals correspond to same
    ordering of quantum as defined in Model hamiltonian. 
    That is terms for DOWN spins should be explicitly added. Then
    only it will be considered here.
  */
  spin_orbitals_.resize(dim_); 
  std::iota(spin_orbitals_.begin(),spin_orbitals_.end(),0);
  spinon_density_.resize(dim_);
  boson_density_.resize(dim_);
  lm_params_.resize(dim_); lm_params_.setZero();
  qp_weights_.resize(dim_); qp_weights_.setOnes();
  spinon_fluct_.resize(dim_,dim_);
  soc_matrix_.resize(dim_,dim_); soc_matrix_.setZero();
  spinon_flip_ampl_.resize(dim_,dim_); spinon_flip_ampl_.setZero();
  boson_flip_ampl_.resize(dim_,dim_); boson_flip_ampl_.setZero();
  spinon_renormed_soc_.resize(dim_,dim_); spinon_renormed_soc_.setZero();
  boson_renormed_soc_.resize(dim_,dim_); boson_renormed_soc_.setZero();
}

void MF_Site::set_soc_matrix(const cmplArray2D& soc_mat)
{
  soc_matrix_ = soc_mat;
  spinon_renormed_soc_ = soc_matrix_;
  boson_renormed_soc_ = soc_matrix_;
  //std::cout << "soc_matrix=" << soc_matrix_ << "\n"; getchar();
}

void MF_Site::set_spinon_renormalization(void)
{
  spinon_renormed_soc_ = soc_matrix_*spinon_flip_ampl_;
}

void MF_Site::set_boson_renormalization(void)
{
  boson_renormed_soc_ = soc_matrix_*boson_flip_ampl_;
  //std::cout << "\n" << soc_matrix_ << "\n";
  //std::cout << "\n" << boson_renormed_soc_ << "\n";
  //getchar();
}


void MF_Site::add_bond(const int& id, const bool& outgoing)
{
  connected_bonds_.push_back(id); 
  bond_outgoing_.push_back(outgoing);  
}
//-------------------------------------------------------------------

//-------------------------------MF Bond-----------------------------
MF_Bond::MF_Bond(const int& type, const bool& is_intracell, const int& src, const int& tgt, 
  const int& vector_id, const Vector3d& vector)
  : type_{type}, is_intracell_{is_intracell}, src_{src}, tgt_{tgt}, vector_id_{vector_id}, 
    vector_{vector}
{
  term_couplings_.clear(); // for all 'bond terms'
  spinon_ke_.clear(); 
  boson_ke_.clear(); 
}

void MF_Bond::add_term_cc(const cmplArray2D& mat) 
{ 
  term_couplings_.push_back(mat); 
  int m = mat.rows();
  int n = mat.cols();
  spinon_renormed_cc_.push_back(cmplArray2D::Zero(m,n)); 
  boson_renormed_cc_.push_back(cmplArray2D::Zero(m,n)); 
  spinon_ke_.push_back(cmplArray2D::Zero(m,n));
  boson_ke_.push_back(cmplArray2D::Zero(m,n));
  //std::cout << spinon_ke_.size() << "\n";
  //std::cout << boson_ke_.size() << "\n\n";
}

void MF_Bond::set_spinon_renormalization(void)
{
  for (int i=0; i<term_couplings_.size(); ++i) {
    spinon_renormed_cc_[i] = term_couplings_[i] * spinon_ke_[i]; 
  }
}

void MF_Bond::set_boson_renormalization(void)
{
  for (int i=0; i<term_couplings_.size(); ++i) {
    boson_renormed_cc_[i] = term_couplings_[i] * boson_ke_[i]; 
  }
}

//-------------------------------MF Params-----------------------------
MF_Params::MF_Params(const input::Parameters& inputs, const lattice::LatticeGraph& graph,
    const model::Hamiltonian& model)
{
  num_basis_sites_ = graph.lattice().num_basis_sites();
  SO_coupling_ = model.is_spinorbit_coupled();
  //num_bondterms_ = model.num_bondterms();

  // add the sites
  sites_.clear();
  for (int i=0; i<num_basis_sites_; ++i) {
    int type = graph.site_type(i);
    int dim = graph.site_dim(i);
  	MF_Site::idx_list state_indices;
    state_indices.resize(dim);
    for (int j=0; j<dim; ++j) 
       state_indices[j] = graph.lattice().basis_index_number(i,j);
    sites_.push_back({type, dim, state_indices});
  }

  // store the bonds in a 'unit cell'
  lattice::LatticeGraph::out_edge_iterator ei, ei_end;
  bonds_.clear();
  bool is_intracell;
  for (int i=0; i<num_basis_sites_; ++i) {
    for (std::tie(ei, ei_end)=graph.out_bonds(i); ei!=ei_end; ++ei) {
      auto type = graph.bond_type(ei);
      auto s = graph.source(ei);
      auto t = graph.target(ei);
      //std::cout << "i="<<i<<", s="<<s<<", t="<<t<<"\n";
      //getchar()

      // src site
      int src = graph.site_uid(s);
      //int src_dim = graph.site_dim(s);

      // tgt site
      int tgt = graph.site_uid(t);
      if (t==tgt) is_intracell = true;
      else is_intracell = false;
      //int tgt_dim = graph.site_dim(t);

      //std::cout << "src="<<src<<", tgt="<<tgt<<"\n"; getchar();
      bonds_.push_back({type,is_intracell,src,tgt,graph.vector_id(ei),graph.vector(ei)});

      // store id of the bond connected to the site
      int id = bonds_.size()-1;
      sites_[src].add_bond(id,true); // outgoing true 
      sites_[tgt].add_bond(id,false); // outgoing false
      //std::cout<<bonds_.back().type()<<": "<< bonds_.back().src()<<" - "<< bonds_.back().tgt()<<"\n";
    }
  }
  num_bonds_ = bonds_.size();

  num_bondterms_ = 0;
  // coupling constants for the Hamiltonian bond term
  for (auto bterm=model.bondterms_begin(); bterm!=model.bondterms_end(); ++bterm) {
    if (bterm->qn_operator().is_quadratic() && bterm->qn_operator().spin_up()) {
    	num_bondterms_++;
    	for (auto& bond : bonds_) {
    		bond.add_term_cc(bterm->coupling(bond.type()));		
    	} 
    }
  }
  //std::cout << "num_bondterms = " << num_bondterms_ << "\n";
  // SOC matrix
  for (auto sterm=model.siteterms_begin(); sterm!=model.siteterms_end(); ++sterm) {
    if (sterm->qn_operator().id()==model::op_id::spin_flip) {
      for (auto& site : sites_) {
        site.set_soc_matrix(sterm->coupling(site.type()));   
      } 
    }
  }
  ke_per_site_ = 0.0;
  pe_per_site_ = 0.0;
}

void MF_Params::update(const model::Hamiltonian& model)
{
  //std::cout << "MF_Params::update\n"; getchar();
  // coupling constants for the Hamiltonian bond term
  for (auto bterm=model.bondterms().begin(); bterm!=model.bondterms().end(); ++bterm) {
    if (bterm->qn_operator().is_quadratic() && bterm->qn_operator().spin_up()) {
      for (auto& bond : bonds_) {
        bond.add_term_cc(bterm->coupling(bond.type()));   
      } 
    }
  }
  //std::cout << "num_bondterms = " << num_bondterms_ << "\n";
  // SOC matrix
  for (auto sterm=model.siteterms().begin(); sterm!=model.siteterms().end(); ++sterm) {
    if (sterm->qn_operator().id()==model::op_id::spin_flip) {
      for (auto& site : sites_) {
        //std::cout << "MF_Params::update::soc_matrix=" << sterm->coupling(0)<<"\n"; getchar();
        site.set_soc_matrix(sterm->coupling(site.type()));   
      } 
    }
  }
  ke_per_site_ = 0.0;
  pe_per_site_ = 0.0;
}

void MF_Params::init_params(void)
{
  for (auto& site : sites_) {
    site.lm_params().setZero();
    site.qp_weights().setOnes();
    site.spinon_flip_ampl().setOnes();
    site.boson_flip_ampl().setOnes();
    site.set_spinon_renormalization();
    site.set_boson_renormalization();
  }
  for (auto& bond : bonds_) {
  	for (int i=0; i<num_bondterms_; ++i) {
    	bond.boson_ke(i).setOnes();
  	  bond.spinon_ke(i).setOnes();
  	}
   	bond.set_spinon_renormalization();
    bond.set_boson_renormalization();
  }
}


} // end namespace ssmf
