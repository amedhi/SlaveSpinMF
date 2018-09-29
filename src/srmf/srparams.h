/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2018-04-29 21:46:50
* @Last Modified by:   Amal Medhi, amedhi@macbook
* @Last Modified time: 2018-04-29 21:47:49
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef SRPARAMS_H
#define SRPARAMS_H

#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <numeric>
#include "../scheduler/task.h"
#include "../lattice/graph.h"
#include "../lattice/matrix.h"
#include "../model/model.h"

namespace srmf {

using real_siteparms_t = std::vector<realArray1D>;
using cmpl_siteparms_t = std::vector<cmplArray1D>;
using cmpl_bondparms_t = std::vector<cmplArray2D>;

class sr_site
{
public:
  using idx_list = std::vector<unsigned>;
  sr_site() : type_{0}, dim_{0} 
  {
    state_indices_.clear();
    connected_bonds_.clear();
    inout_types_.clear();
  } 
  sr_site(const unsigned& type, const unsigned& dim, const idx_list& state_indices, 
    const bool& SO_coupled)
    : type_{type}, dim_{dim}, state_indices_{state_indices}
  {
    connected_bonds_.clear();
    inout_types_.clear();
    // spin-oribitals
    if (SO_coupled) spin_orbitals_.resize(dim); 
    else spin_orbitals_.resize(2*dim); // spin + oribitals;
    std::iota(spin_orbitals_.begin(), spin_orbitals_.end(), 0);
    spinon_density_.resize(spin_orbitals_.size());
    boson_density_.resize(spin_orbitals_.size());
    lm_params_.resize(spin_orbitals_.size());
  }
  void clear(void) { connected_bonds_.clear(); inout_types_.clear(); }
  const unsigned& type(void) const { return type_; }
  const unsigned& dim(void) const { return dim_; }
  const std::vector<unsigned> spin_orbitals(void) const { return spin_orbitals_; }
  void add_bond(const unsigned& id, const bool& outgoing) 
    { connected_bonds_.push_back(id); inout_types_.push_back(outgoing); } 
  const idx_list& state_indices(void) const { return state_indices_; }
  const idx_list& connected_bonds(void) const { return connected_bonds_; }
  bool is_outgoing_bond(const unsigned& i) const { return inout_types_[i]; }
  realArray1D& spinon_density(void) { return spinon_density_; }
  realArray1D& lm_params(void) { return lm_params_; }
  const realArray1D& lm_params(void) const { return lm_params_; }
  realArray1D& qp_weights(void) { return qp_weights_; }
  const realArray1D& qp_weights(void) const { return qp_weights_; }
private:
  unsigned type_;
  unsigned dim_;
  std::vector<unsigned> spin_orbitals_;
  idx_list state_indices_;
  idx_list connected_bonds_;
  std::vector<bool> inout_types_;
  realArray1D spinon_density_;
  realArray1D boson_density_;
  realArray1D lm_params_;
  realArray1D qp_weights_;
};

class sr_bond
{
public:
  using idx_list = std::vector<unsigned>;
  sr_bond(const unsigned& type, const unsigned& src, const unsigned& src_dim,
    const unsigned& tgt, const unsigned& tgt_dim, const Vector3d& vector,
    const bool& SO_coupled);
  void add_term_cc(const cmplArray2D& mat, const model::spin& s);
  const unsigned& type(void) const { return type_; }
  const unsigned& src(void) const { return src_; }
  const unsigned& tgt(void) const { return tgt_; }
  void set_spinon_renormalization(void);
  void set_boson_ke(void) { boson_ke_.setOnes(); }
  void set_spinon_ke(void) { spinon_ke_.setOnes(); set_spinon_renormalization(); }
  //const idx_list& src_state_indices(void) const { return src_state_indices_; }
  //const idx_list& tgt_state_indices(void) const { return tgt_state_indices_; }
  const Vector3d& vector(void) const { return vector_; }
  cmplArray2D& spinon_ke() { return spinon_ke_; }
  cmplArray2D& boson_ke() { return boson_ke_; }
  const cmplArray2D& spinon_ke() const { return spinon_ke_; }
  const cmplArray2D& boson_ke() const { return boson_ke_; }
  int num_bond_terms(void) const { return term_couplings_.size(); }
  const cmplArray2D& term_cc(const unsigned& i) const { return term_couplings_[i]; }
  const cmplArray2D& spinon_renormed_cc(const unsigned& i) const 
    { return spinon_renormed_couplings_[i]; }
private:
  unsigned type_;
  unsigned src_;
  unsigned tgt_;
  Vector3d vector_;
  bool SO_coupled_;
  //idx_list src_state_indices_;
  //idx_list tgt_state_indices_;
  std::vector<cmplArray2D> term_couplings_; // model term coupling constants 
  std::vector<model::spin> term_spins_; 
  cmplArray2D spinon_ke_;
  std::vector<cmplArray2D> spinon_renormed_couplings_;
  cmplArray2D boson_ke_;
  std::vector<cmplArray2D> boson_renormed_couplings_;
};

class site_link 
{
public:
  site_link(const unsigned& bond_id, const bool& incoming) 
    : id_{bond_id}, incoming_{incoming}  {}
  const unsigned& id(void) const { return id_; }
  const bool& is_incoming(void) const { return incoming_; }
private:
  unsigned id_{0};
  bool incoming_{false};
};

class SR_Params 
{
public:
  //using bond = sr_bond;
  //using site = sr_site;
  using links = std::vector<site_link>;
  SR_Params(const input::Parameters& inputs, const lattice::LatticeGraph& graph,
    const model::Hamiltonian& model);
  ~SR_Params() {}
  //int init(const lattice::Lattice& lattice) override;
  void init_mf_params(void);
  const unsigned& num_sites(void) const { return num_sites_; }
  const unsigned& num_bonds(void) const { return num_bonds_; }
  const sr_site& site(const unsigned& i) const { return sites_[i]; }
  sr_site& site(const unsigned& i) { return sites_[i]; }
  const sr_bond& bond(const unsigned& i) const { return bonds_[i]; }
  sr_bond& bond(const unsigned& i) { return bonds_[i]; }
  const std::vector<sr_site>& sites(void) const { return sites_; }
  const std::vector<sr_bond>& bonds(void) const { return bonds_; }

  //const std::vector<links>& site_links(void) const { return site_links_; }
  //cmplArray2D& sp_bond_ke(const int& i) { return sp_bond_ke_[i]; }
  //realArray1D& sp_site_density(const int& i) { return sp_site_density_[i]; }
  //std::vector<std::complex<double>> sbond_avg(void) { return sbond_avg_; }
  //std::vector<std::complex<double>> rbond_avg(void) { return rbond_avg_; }
  //const double& spinon_density(const int& i) const { return spinon_site_density_[i]; }
  //const cmplArray1D& bond_tchi(void) const { return bond_tchi_; }
  //cmplArray1D& bond_tchi(void) { return bond_tchi_; }
private:
  //using LatticeGraph = lattice::LatticeGraph;
  bool SO_coupling_{false};
  unsigned num_sites_{0}; // no of sites per unitcell
  unsigned num_bonds_{0}; // no of bonds connected to all sites in unitcell
  std::vector<sr_site> sites_; // list of all sites in unitcell
  std::vector<sr_bond> bonds_; // list of all bonds
  //std::vector<links> site_links_; // bonds connecting every site
  //cmplArray1D bond_tchi_;
  //std::vector<cmplArray2D> sp_bond_ke_;
  //std::vector<realArray1D> sp_site_density_;
  //std::vector<std::complex<double>> rbond_avg_;
  //std::vector<double> spinon_site_density_;
  //std::vector<double> ssite_avg_;
  //std::vector<double> rsite_avg_;
};




} // end namespace srmf

#endif
