/*---------------------------------------------------------------------------
* Author: Amal Medhi
* @Date:   2019-03-12 12:20:33
* @Last Modified by:   Amal Medhi
* @Last Modified time: 2019-05-05 19:43:35
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef MF_PARAMS_H
#define MF_PARAMS_H

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

namespace ssmf {

using real_siteparms_t = std::vector<realArray1D>;
using cmpl_siteparms_t = std::vector<cmplArray1D>;
using cmpl_bondparms_t = std::vector<cmplArray2D>;

class MF_Site
{
public:
  using idx_list = std::vector<int>;
  MF_Site();
  MF_Site(const int& type, const int& dim, const idx_list& state_indices);
  void set_soc_matrix(const cmplArray2D& soc_mat); 
  const int& type(void) const { return type_; }
  const int& dim(void) const { return dim_; }
  const std::vector<int> spin_orbitals(void) const { return spin_orbitals_; }
  void add_bond(const int& id, const bool& outgoing);
  const idx_list& state_indices(void) const { return state_indices_; }
  const idx_list& connected_bonds(void) const { return connected_bonds_; }
  bool bond_outgoing(const int& i) const { return bond_outgoing_[i]; }
  realArray1D& spinon_density(void) { return spinon_density_; }
  realArray1D& lm_params(void) { return lm_params_; }
  const realArray1D& lm_params(void) const { return lm_params_; }
  realArray1D& qp_weights(void) { return qp_weights_; }
  const realArray1D& qp_weights(void) const { return qp_weights_; }
  cmplArray2D& spinon_fluct(void) { return spinon_fluct_; }
  const cmplArray2D& spinon_fluct(void) const { return spinon_fluct_; }
  const cmplArray2D& spinon_flip_ampl(void) const { return spinon_flip_ampl_; }
  cmplArray2D& spinon_flip_ampl(void) { return spinon_flip_ampl_; }
  const cmplArray2D& boson_flip_ampl(void) const { return boson_flip_ampl_; }
  cmplArray2D& boson_flip_ampl(void) { return boson_flip_ampl_; }
  void set_spinon_renormalization(void);
  void set_boson_renormalization(void);
  const cmplArray2D& spinon_renormed_soc(void) const 
    { return spinon_renormed_soc_; }
  const cmplArray2D& boson_renormed_soc(void) const 
    { return boson_renormed_soc_; }
private:
  int type_;
  int dim_;
  std::vector<int> spin_orbitals_;
  idx_list state_indices_;
  idx_list connected_bonds_;
  std::vector<bool> bond_outgoing_;
  realArray1D spinon_density_;
  realArray1D boson_density_;
  realArray1D lm_params_;
  realArray1D lm_params_noint_;
  realArray1D qp_weights_;
  cmplArray2D spinon_fluct_;
  cmplArray2D soc_matrix_; 
  cmplArray2D spinon_flip_ampl_; 
  cmplArray2D boson_flip_ampl_; 
  cmplArray2D spinon_renormed_soc_; 
  cmplArray2D boson_renormed_soc_; 
};

class MF_Bond
{
public:
  using idx_list = std::vector<unsigned>;
  MF_Bond(const int& type, const bool& is_intracell, const int& src, 
    const int& tgt, const int& vector_id, const Vector3d& vector);
  ~MF_Bond() {}
	void add_term_cc(const cmplArray2D& mat);
  void set_spinon_renormalization(void);
  void set_boson_renormalization(void);
  const int& type(void) const { return type_; }
  //void make_intracell(void) const { is_intracell_ = true; }
  const bool& is_intracell(void) const { return is_intracell_; }
  const int& src(void) const { return src_; }
  const int& tgt(void) const { return tgt_; }
  const int& vector_id(void) const { return vector_id_; }
  const Vector3d& vector(void) const { return vector_; }
  const cmplArray2D& term_coupling(const int& i) const 
    { return term_couplings_[i]; }
  const cmplArray2D& spinon_renormed_cc(const int& i) const 
    { return spinon_renormed_cc_[i]; }
  const cmplArray2D& boson_renormed_cc(const int& i) const 
    { return boson_renormed_cc_[i]; }
  cmplArray2D& spinon_ke(const int& i) { return spinon_ke_[i]; }
  cmplArray2D& boson_ke(const int& i) { return boson_ke_[i]; }
private:
  int type_;
  bool is_intracell_{false};
  int src_;
  int tgt_;
  int vector_id_;
  Vector3d vector_;
  std::vector<cmplArray2D> term_couplings_; // for all 'bond terms'
  std::vector<cmplArray2D> spinon_renormed_cc_; 
  std::vector<cmplArray2D> boson_renormed_cc_; 
  std::vector<cmplArray2D> spinon_ke_; // for all 'bond terms'
  std::vector<cmplArray2D> boson_ke_; // for all 'bond terms'
};


class MF_Params
{
public:
  MF_Params(const input::Parameters& inputs, const lattice::LatticeGraph& graph,
    const model::Hamiltonian& model);
  ~MF_Params() {}
  bool SO_coupling(void) const { return SO_coupling_; }
  const int& num_sites(void) const { return num_basis_sites_; }
  const int& num_bonds(void) const { return num_bonds_; }
  void update(const model::Hamiltonian& model);
  void init_params(void);
  const MF_Site& site(const int& i) const { return sites_[i]; }
  MF_Site& site(const int& i) { return sites_[i]; }
  const MF_Bond& bond(const int& i) const { return bonds_[i]; }
  MF_Bond& bond(const int& i) { return bonds_[i]; }
  const std::vector<MF_Bond>& bonds(void) const { return bonds_; }
  const double& ke_per_site(void) const { return ke_per_site_; }
  const double& pe_per_site(void) const { return pe_per_site_; }
  double& ke_per_site(void) { return ke_per_site_; }
  double& pe_per_site(void) { return pe_per_site_; }
  double total_energy(void) const { return ke_per_site_+pe_per_site_; }
private:
  int num_basis_sites_;
  int num_bonds_;
  bool SO_coupling_{false};
  int num_bondterms_;
  double ke_per_site_{0.0};
  double pe_per_site_{0.0};
  std::vector<MF_Site> sites_; // list of all sites in unitcell
  std::vector<MF_Bond> bonds_; // list of all bonds
  std::vector<cmplArray2D> spinon_ke_;
  std::vector<cmplArray2D> boson_ke_;
};



} // end namespace ssmf

#endif