/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-29 21:46:50
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-10-05 13:29:19
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "sb_params.h"

namespace srmf {

//-------------------------SR_Bonds------------------------
sb_bond::sb_bond(const unsigned& type, const unsigned& src, const unsigned& src_dim,
  const unsigned& tgt, const unsigned& tgt_dim, const Vector3d& vector,
  const bool& SO_coupled)
  : type_{type}, src_{src}, tgt_{tgt}, vector_{vector}, SO_coupled_{SO_coupled}
{
  if (SO_coupled_) {
    spinon_ke_.resize(src_dim, tgt_dim);
    boson_ke_.resize(src_dim, tgt_dim);
  }
  else {
    spinon_ke_.resize(2*src_dim, 2*tgt_dim);
    boson_ke_.resize(2*src_dim, 2*tgt_dim);
  }
  term_couplings_.clear();
  term_spins_.clear();
  spinon_renormed_couplings_.clear();
  boson_renormed_couplings_.clear();
}

void sb_bond::add_term_cc(const cmplArray2D& mat, const model::spin& s) 
{ 
  term_couplings_.push_back(mat); 
  term_spins_.push_back(s);
  int m = mat.rows();
  int n = mat.cols();
  if (!SO_coupled_) {
    m *= 2;
    n *= 2;
  }
  spinon_renormed_couplings_.push_back(cmplArray2D::Zero(m,n));
  boson_renormed_couplings_.push_back(cmplArray2D::Zero(m,n));
}

void sb_bond::set_spinon_renormalization(void)
{
  if(SO_coupled_) {
    for (int i=0; i<term_couplings_.size(); ++i) {
      spinon_renormed_couplings_[i] = term_couplings_[i] * spinon_ke_;   
    }
  }
  else {
    int m = spinon_ke_.rows();
    int n = spinon_ke_.cols();
    cmplArray2D coupling=cmplArray2D::Zero(m,n);
    for (int i=0; i<term_couplings_.size(); ++i) {
      coupling.block(0,0,m/2,n/2) = term_couplings_[i];
      coupling.block(m/2,n/2,m/2,n/2) = term_couplings_[i];
      spinon_renormed_couplings_[i] = coupling * spinon_ke_;   
    }
  }
}


//-------------------------SB_Params------------------------
SB_Params::SB_Params(const input::Parameters& inputs,const lattice::LatticeGraph& graph,
  const model::Hamiltonian& model)
{
  SO_coupling_ = model.is_spinorbit_coupled();
  //std::cout << "----SB_Params::SB_Params----\n";
  // sites in the unit cell
  num_sites_ = graph.lattice().num_basis_sites();
  sb_site::idx_list state_indices;
  sites_.clear();
  for (unsigned i=0; i<num_sites_; ++i) {
    unsigned type = graph.site_type(i);
    unsigned dim = graph.site_dim(i);
    //sites_.push_back()
    state_indices.resize(dim);
    for (unsigned j=0; j<dim; ++j) 
       state_indices[j] = graph.lattice().basis_index_number(i,j);
    sites_.push_back({type, dim, state_indices, SO_coupling_});
  }

  // store the bonds in a 'unit cell'
  std::vector<unsigned> src_state_indices;
  std::vector<unsigned> tgt_state_indices;
  lattice::LatticeGraph::out_edge_iterator ei, ei_end;
  bonds_.clear();
  //site_links_.clear();
  //site_links_.resize(num_sites_);
  for (auto& site : sites_) site.clear();
  for (unsigned i=0; i<num_sites_; ++i) {
    for (std::tie(ei, ei_end)=graph.out_bonds(i); ei!=ei_end; ++ei) {
      auto type = graph.bond_type(ei);
      auto s = graph.source(ei);
      auto t = graph.target(ei);
      //std::cout << "i="<<i<<", s="<<s<<", t="<<t<<"\n";

      // src site
      unsigned src = graph.site_uid(s);
      unsigned src_dim = graph.site_dim(s);

      //src_state_indices.resize(n);
      //for (unsigned j=0; j<n; ++j) 
      //  src_state_indices[j] = graph.lattice().basis_index_number(m,j);

      // tgt site
      unsigned tgt = graph.site_uid(t);
      unsigned tgt_dim = graph.site_dim(t);
      //tgt_state_indices.resize(q);
      //for (unsigned j=0; j<q; ++j) 
      //  tgt_state_indices[j] = graph.lattice().basis_index_number(p,j);

      //std::cout << "m="<<m<<", n="<<n<<"\n"; getchar();
      bonds_.push_back({type,src,src_dim,tgt,tgt_dim,graph.vector(ei),SO_coupling_});
      // Hamiltonian terms for this bonds
      //ComplexMatrix coeff_mat = hamterm.coupling(btype);
      for (auto bterm=model.bondterms_begin(); bterm!=model.bondterms_end(); ++bterm) {
        if (bterm->qn_operator().is_quadratic()) {
          bonds_.back().add_term_cc(bterm->coupling(type), bterm->qn_operator().sigma());
          //std::cout << "coupling = " << bterm->coupling(type) << "\n";
        }
      }

      // store id of the bond connected to the site
      int id = bonds_.size()-1;
      sites_[src].add_bond(id,true); // outgoing true 
      sites_[tgt].add_bond(id,false); // outgoing false

      //site_links_[m].push_back({id,false}); // outgoing from 'm'
      //site_links_[p].push_back({id,true}); // incoming to 'n'
      //std::cout<<bonds_.back().type()<<": "<< bonds_.back().src()<<" - "<< bonds_.back().tgt()<<"\n";
    }
  }
  num_bonds_ = bonds_.size();
}

void SB_Params::init_mf_params(void)
{
  for (unsigned i=0; i<num_sites_; ++i) {
    sites_[i].lm_params().setZero();
    sites_[i].qp_weights().setOnes();
  }
  for (unsigned i=0; i<num_bonds_; ++i) {
    bonds_[i].boson_ke().setOnes();
    bonds_[i].spinon_ke().setOnes();
    bonds_[i].set_spinon_renormalization();
  }
}


} // end namespace srmf
