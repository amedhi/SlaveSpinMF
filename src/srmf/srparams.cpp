/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-29 21:46:50
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-09-28 10:12:26
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "srparams.h"

namespace srmf {

SR_Params::SR_Params(const input::Parameters& inputs,const lattice::LatticeGraph& graph,
  const model::Hamiltonian& model)
{
  SO_coupling_ = model.is_spinorbit_coupled();
  //std::cout << "----SR_Params::SR_Params----\n";
  // sites in the unit cell
  num_sites_ = graph.lattice().num_basis_sites();
  sr_site::idx_list state_indices;
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
      bonds_.push_back({type,src,src_dim,tgt,tgt_dim,graph.vector(ei)});
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


  /*
  // storages
  sp_bond_ke_.resize(num_bonds_);
  for (int i=0; i<num_bonds_; ++i) {
    int rows = bonds_[i].src_state_indices().size();
    int cols = bonds_[i].tgt_state_indices().size();
    sp_bond_ke_[i].resize(rows,cols);
  }

  sp_site_density_.resize(num_sites_);
  for (int i=0; i<num_sites_; ++i) {
    sp_site_density_[i].resize(graph.site_dim(i));
  }

  //bond_tchi_.resize(num_bonds_);
  rbond_avg_.resize(num_bonds_);
  ssite_avg_.resize(num_sites_);
  rsite_avg_.resize(num_sites_);
  spinon_site_density_.resize(num_sites_);

  for (int i=0; i<num_bonds_; ++i) {
    //bond_tchi_[i] = 1.0;
    rbond_avg_[i] = 1.0;
  }
  for (int i=0; i<num_sites_; ++i) {
    ssite_avg_[i] = 0.0;
    rsite_avg_[i] = 0.0;
    spinon_site_density_[i] = 1.0;
  }

  // spinon density
  */
}



} // end namespace srmf
