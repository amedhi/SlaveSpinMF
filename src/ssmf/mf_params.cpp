/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-03-12 12:20:33
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-04-05 11:02:18
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "mf_params.h"

namespace srmf {

//-------------------------------MF Site-----------------------------
MF_Site::MF_Site() 
  : type_{0}, dim_{0} 
{
  spin_orbitals_.clear(); 
  state_indices_.clear();
  connected_bonds_.clear();
  bond_outgoing_.clear();
} 

MF_Site::MF_Site(const int& type, const int& dim, const idx_list& state_indices, 
    const bool& SO_coupled)
  : type_{type}, dim_{dim}, state_indices_{state_indices}
{
  connected_bonds_.clear();
  bond_outgoing_.clear();
  // spin-oribitals
  if (SO_coupled) spin_orbitals_.resize(dim); 
  else spin_orbitals_.resize(2*dim); // spin + oribitals;
  std::iota(spin_orbitals_.begin(),spin_orbitals_.end(),0);
  spinon_density_.resize(spin_orbitals_.size());
  boson_density_.resize(spin_orbitals_.size());
  lm_params_.resize(spin_orbitals_.size());
  qp_weights_.resize(spin_orbitals_.size());
}

void MF_Site::add_bond(const int& id, const bool& outgoing)
{
  connected_bonds_.push_back(id); 
  bond_outgoing_.push_back(outgoing);  
}
//-------------------------------------------------------------------

//-------------------------------MF Bond-----------------------------
MF_Bond::MF_Bond(const int& type, const int& src, const int& tgt, 
  const int& vector_id, const Vector3d& vector,  const bool& SO_coupled)
  : type_{type}, src_{src}, tgt_{tgt}, vector_id_{vector_id}, 
    vector_{vector}, SO_coupled_{SO_coupled}
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
  if (!SO_coupled_) {
    m *= 2;
    n *= 2;
  }
  spinon_renormed_cc_.push_back(cmplArray2D::Zero(m,n)); 
  boson_renormed_cc_.push_back(cmplArray2D::Zero(m,n)); 
  spinon_ke_.push_back(cmplArray2D::Zero(m,n));
  boson_ke_.push_back(cmplArray2D::Zero(m,n));
  //std::cout << spinon_ke_.size() << "\n";
  //std::cout << boson_ke_.size() << "\n\n";
}

void MF_Bond::set_spinon_renormalization(void)
{
  if(SO_coupled_) {
    for (int i=0; i<term_couplings_.size(); ++i) {
      spinon_renormed_cc_[i] = term_couplings_[i] * spinon_ke_[i]; 
    }
  }
  else {
    for (int i=0; i<term_couplings_.size(); ++i) {
    	int m = spinon_ke_[i].rows();
    	int n = spinon_ke_[i].cols();
    	cmplArray2D coupling=cmplArray2D::Zero(m,n);
      coupling.block(0,0,m/2,n/2) = term_couplings_[i];
      coupling.block(m/2,n/2,m/2,n/2) = term_couplings_[i];
      spinon_renormed_cc_[i] = coupling * spinon_ke_[i];   
    }
  }
}

void MF_Bond::set_boson_renormalization(void)
{
  if(SO_coupled_) {
    for (int i=0; i<term_couplings_.size(); ++i) {
      boson_renormed_cc_[i] = term_couplings_[i] * boson_ke_[i]; 
    }
  }
  else {
    for (int i=0; i<term_couplings_.size(); ++i) {
      int m = boson_ke_[i].rows();
      int n = boson_ke_[i].cols();
      cmplArray2D coupling=cmplArray2D::Zero(m,n);
      coupling.block(0,0,m/2,n/2) = term_couplings_[i];
      coupling.block(m/2,n/2,m/2,n/2) = term_couplings_[i];
      boson_renormed_cc_[i] = coupling * boson_ke_[i];   
      //std::cout << boson_renormed_cc_[i] << "\n"; getchar();
    }
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
    sites_.push_back({type, dim, state_indices, SO_coupling_});
  }

  // store the bonds in a 'unit cell'
  lattice::LatticeGraph::out_edge_iterator ei, ei_end;
  bonds_.clear();
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
      //int tgt_dim = graph.site_dim(t);

      //std::cout << "src="<<src<<", tgt="<<tgt<<"\n"; getchar();
      bonds_.push_back({type,src,tgt,graph.vector_id(ei),graph.vector(ei),
      	SO_coupling_});

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
}

void MF_Params::init_params(void)
{
  for (int i=0; i<num_basis_sites_; ++i) {
    sites_[i].lm_params().setZero();
    sites_[i].qp_weights().setOnes();
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



} // end namespace srmf
