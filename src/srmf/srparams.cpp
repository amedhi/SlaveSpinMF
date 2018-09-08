/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-29 21:46:50
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-09-08 16:54:43
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "srparams.h"

namespace srmf {

SR_Params::SR_Params(const input::Parameters& inputs, const lattice::LatticeGraph& graph)
{
  //std::cout << "----SR_Params::SR_Params----\n";
  num_sites_ = graph.lattice().num_basis_sites();
  // store the bonds in a 'unit cell'
  std::vector<unsigned> src_state_indices;
  std::vector<unsigned> tgt_state_indices;
  lattice::LatticeGraph::out_edge_iterator ei, ei_end;
  bonds_.clear();
  site_links_.clear();
  site_links_.resize(num_sites_);
  for (unsigned i=0; i<num_sites_; ++i) {
    for (std::tie(ei, ei_end)=graph.out_bonds(i); ei!=ei_end; ++ei) {
      auto type = graph.bond_type(ei);
      auto s = graph.source(ei);
      auto t = graph.target(ei);
      //std::cout << "i="<<i<<", s="<<s<<", t="<<t<<"\n";

      // src site
      unsigned m = graph.site_uid(s);
      unsigned n = graph.site_dim(s);
      src_state_indices.resize(n);
      for (unsigned j=0; j<n; ++j) 
        src_state_indices[j] = graph.lattice().basis_index_number(m,j);

      // tgt site
      unsigned p = graph.site_uid(t);
      unsigned q = graph.site_dim(t);
      tgt_state_indices.resize(q);
      for (unsigned j=0; j<q; ++j) 
        tgt_state_indices[j] = graph.lattice().basis_index_number(p,j);

      //std::cout << "m="<<m<<", n="<<n<<"\n"; getchar();
      //bonds_.push_back({type,m,n,graph.vector(ei)});
      bonds_.push_back({type,src_state_indices,tgt_state_indices,graph.vector(ei)});
      // store id of the bond connected to the site
      int id = bonds_.size()-1;
      site_links_[m].push_back({id,false}); // outgoing from 'm'
      site_links_[p].push_back({id,true}); // incoming to 'n'
      //std::cout<<bonds_.back().type()<<": "<< bonds_.back().src()<<" - "<< bonds_.back().tgt()<<"\n";
    }
  }
  num_bonds_ = bonds_.size();

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

}



} // end namespace srmf
