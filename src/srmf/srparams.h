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
#include "../scheduler/task.h"
#include "../lattice/graph.h"
#include "../lattice/matrix.h"

namespace srmf {

class sr_bond
{
public:
  using idx_list = std::vector<unsigned>;
  sr_bond(const unsigned& type, const unsigned& src, const unsigned& tgt,
    const Vector3d& vector) 
    : type_{type}, src_{src}, tgt_{tgt}, vector_{vector} {}
  sr_bond(const unsigned& type, const idx_list& idx_list1, 
    const idx_list& idx_list2, const Vector3d& vector)
    : type_{type}, src_state_indices_{idx_list1}, tgt_state_indices_{idx_list2}, 
      vector_{vector} {}
  const unsigned& type(void) const { return type_; }
  const unsigned& src(void) const { return src_; }
  const unsigned& tgt(void) const { return tgt_; }
  const idx_list& src_state_indices(void) const { return src_state_indices_; }
  const idx_list& tgt_state_indices(void) const { return tgt_state_indices_; }
  const Vector3d& vector(void) const { return vector_; }
private:
  unsigned type_;
  unsigned src_;
  unsigned tgt_;
  idx_list src_state_indices_;
  idx_list tgt_state_indices_;
  Vector3d vector_;
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
  using bond = sr_bond;
  using links = std::vector<site_link>;
  SR_Params(const input::Parameters& inputs, const lattice::LatticeGraph& graph);
  ~SR_Params() {}
  //int init(const lattice::Lattice& lattice) override;
  const unsigned& num_sites(void) const { return num_sites_; }
  const unsigned& num_bonds(void) const { return num_bonds_; }
  const std::vector<bond>& bonds(void) const { return bonds_; }
  const std::vector<links>& site_links(void) const { return site_links_; }
  ComplexArray& sp_bond_ke(const int& i) { return sp_bond_ke_[i]; }
  ComplexArray1D& sp_site_density(const int& i) { return sp_site_density_[i]; }
  //std::vector<std::complex<double>> sbond_avg(void) { return sbond_avg_; }
  std::vector<std::complex<double>> rbond_avg(void) { return rbond_avg_; }
  const double& spinon_density(const int& i) const { return spinon_site_density_[i]; }
  const ArrayXcd& bond_tchi(void) const { return bond_tchi_; }
  ArrayXcd& bond_tchi(void) { return bond_tchi_; }
private:
  //using LatticeGraph = lattice::LatticeGraph;
  unsigned num_sites_{0}; // no of sites per unitcell
  unsigned num_bonds_{0}; // no of bonds connected to all sites in unitcell
  std::vector<bond> bonds_; // list of all bonds
  std::vector<links> site_links_; // bonds connecting every site
  ArrayXcd bond_tchi_;
  std::vector<ComplexArray> sp_bond_ke_;
  std::vector<ComplexArray1D> sp_site_density_;
  std::vector<std::complex<double>> rbond_avg_;
  std::vector<double> spinon_site_density_;
  std::vector<double> ssite_avg_;
  std::vector<double> rsite_avg_;
};




} // end namespace srmf

#endif
