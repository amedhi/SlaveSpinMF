/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@macbook
* @Date:   2018-04-21 11:41:01
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-03-07 13:58:38
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <cassert>
#include "boson_basis.h"

namespace srmf {

//--------------------------------SlaveSpinBasis-------------------------------------
void SlaveSpinBasis::construct(const int& num_sites, const int& site_dim)
{
  /*
    basis for single site = |a_1/b_1; a_2/b_2; ....> 
    The subscripts denote spin-orbital indices. Number of spin-orbitals
    per site is 'site_dim'. For each spin-orbital index, the following
    constraint applies: (a^\dag a + b^\dag b) = 1.  

    'a'-particle is represented by bit value: '1' == |+1/2>
    'b'-particle is represented by bit value: '0' == |-1/2>
  */
  //std::cout << "SlaveSpinBasis::construct\n";
  num_sites_ = num_sites;
  site_dim_ = site_dim; 
  num_bits_ = num_sites_ * site_dim_;
  if (num_bits_ > max_bits-1) {
    throw std::runtime_error("SlaveSpinBasis::construct: basis dimension exceeds limits");
  }
  ndim_ = std::pow(2,num_bits_);

  // state indices
  basis_states_.resize(ndim_+1);
  state_indices_.resize(ndim_+1);
  for (unsigned i=0; i<=ndim_; ++i) {
    basis_states_[i] = state_t(i);
    state_indices_[basis_states_[i].to_ullong()] = i;
    //std::cout << "state: " << i << " =" << basis_states_[i].to_string() << "\n";
  }
  null_state_ = state_t(ndim_);
  null_idx_ = ndim_;
}

SlaveSpinBasis::op_result SlaveSpinBasis::apply_adag_a(const size_t& site, 
  const size_t& alpha, const idx_t& idx) const
{
  if (idx == null_idx_) return std::make_pair(0,null_idx_);
  idx_t pos = site * site_dim_ + alpha;
  // a^\dag a 
  state_t state = basis_states_[idx];
  if (state.test(pos)) return std::make_pair(1,idx);
  else return std::make_pair(0,idx);
}

SlaveSpinBasis::op_result SlaveSpinBasis::apply_bdag_b(const size_t& site, 
  const size_t& alpha, const idx_t& idx) const
{
  if (idx == null_idx_) return std::make_pair(0,null_idx_);
  idx_t pos = site * site_dim_ + alpha;
  // b^\dag b 
  state_t state = basis_states_[idx];
  if (state.test(pos)) return std::make_pair(0,idx);
  else return std::make_pair(1,idx);
}

SlaveSpinBasis::op_result SlaveSpinBasis::apply_Sz(const size_t& site, 
  const size_t& alpha, const idx_t& idx) const
{
  if (idx == null_idx_) return std::make_pair(0,null_idx_);
  idx_t pos = site * site_dim_ + alpha;
  // Sz == 1/2(a^\dag a - b^\dag b)
  state_t state = basis_states_[idx];
  if (state.test(pos)) return std::make_pair(+0.5,idx);
  else return std::make_pair(-0.5,idx);
}

SlaveSpinBasis::op_result SlaveSpinBasis::apply_Oplus(const double& c,
  const size_t& site, const size_t& alpha, const idx_t& idx) const
{
  if (idx == null_idx_) return std::make_pair(0,null_idx_);
  idx_t pos = site * site_dim_ + alpha;
  // O+ == (S+ + c*S-)
  state_t state = basis_states_[idx];
  if (state.test(pos)) {
    state.reset(pos);
    return std::make_pair(c,state_indices_[state.to_ullong()]);
  }
  else {
    state.set(pos);
    return std::make_pair(1,state_indices_[state.to_ullong()]);
  }
}

SlaveSpinBasis::op_result SlaveSpinBasis::apply_Ominus(const double& c,
  const size_t& site, const size_t& alpha, const idx_t& idx) const
{
  if (idx == null_idx_) return std::make_pair(0,null_idx_);
  idx_t pos = site * site_dim_ + alpha;
  // O- == (S- + c*S+)
  state_t state = basis_states_[idx];
  if (state.test(pos)) {
    state.reset(pos);
    return std::make_pair(1,state_indices_[state.to_ullong()]);
  }
  else {
    state.set(pos);
    return std::make_pair(c,state_indices_[state.to_ullong()]);
  }
}

SlaveSpinBasis::op_result SlaveSpinBasis::apply_Splus(const size_t& site, 
  const size_t& alpha, const idx_t& idx) const 
{
  if (idx == null_idx_) return std::make_pair(0,null_idx_);
  idx_t pos = site * site_dim_ + alpha;
  // S+ == a^\dag . b
  state_t state = basis_states_[idx];
  if (state.test(pos)) {
    return std::make_pair(0,null_idx_);
  }
  else {
    state.set(pos);
    return std::make_pair(1,state_indices_[state.to_ullong()]);
  }
}


SlaveSpinBasis::op_result SlaveSpinBasis::apply_Sminus(const size_t& site, 
  const size_t& alpha, const idx_t& idx) const 
{
  if (idx == null_idx_) return std::make_pair(0,null_idx_);
  idx_t pos = site * site_dim_ + alpha;
  // S- == b^dag . a
  state_t state = basis_states_[idx];
  if (state.test(pos)) {
    state.reset(pos);
    return std::make_pair(1,state_indices_[state.to_ullong()]);
  }
  else {
    return std::make_pair(0,null_idx_);
  }
}

SlaveSpinBasis::op_result SlaveSpinBasis::apply_Zplus(const double& c,
  const size_t& site, const size_t& alpha, const idx_t& idx) const
{
  if (idx == null_idx_) return std::make_pair(0,null_idx_);
  idx_t pos = site * site_dim_ + alpha;
  // Z+ = <P+> a^\dag b <P-> = c a^\dag b
  state_t state = basis_states_[idx];
  //std::cout << "pos = " << pos << "\n";
  //std::cout << "state = " << std::bitset<6>(state.to_ulong()) << "\n";
  if (state.test(pos)) {
    return std::make_pair(0,null_idx_);
  }
  else {
    // apply b^\dag a
    state.set(pos); 
    return std::make_pair(c,state_indices_[state.to_ullong()]);
  }
}

SlaveSpinBasis::op_result SlaveSpinBasis::apply_Zminus(const double& c, 
  const size_t& site, const size_t& alpha, const idx_t& idx) const
{
  if (idx == null_idx_) return std::make_pair(0,null_idx_);
  idx_t pos = site * site_dim_ + alpha;
  // Z- = <P-> b^\dag a <P+> = c b^\dag a 
  state_t state = basis_states_[idx];
  if (state.test(pos)) {
    // apply a^\dag b
    state.reset(pos); 
    return std::make_pair(c,state_indices_[state.to_ullong()]);
  }
  else return std::make_pair(0,null_idx_);
}



//--------------------------------RotorBasis-------------------------------------
void qbitset::construct(const unsigned& q, const unsigned& num_qbits)
{
  if (q < 1) {
    throw std::invalid_argument("qbitset::construct: invalid value to argument-1");
  }
  q_ = q;
  num_bits_ = num_qbits;
  if (q_ > 1) {
    unsigned n = q_-1;
    block_size_ = 0;
    while (n) {
      n >>= 1;
      block_size_++;
    }
  }
  //std::cout << "block_size_ = " << block_size_ << "\n";
  if (num_bits_*block_size_ > max_bits) {
    throw std::runtime_error("qbitset::construct: qbitset too long");
  }
  // store the bit values
  bitvals_.clear();
  for (unsigned i=0; i<q_; ++i) { 
    bitvals_.push_back({i}); 
    //std::cout << "bit = " << bitvals_[i] << "\n";
  }
  maxbit_ = bitvals_.back();
  bitstring_.reset();
  //std::cout << bitstring_ << "\n";
}

void qbitset::reset()
{
  bitstring_.reset();
}

void qbitset::reset(const qbitset::size_t& pos)
{
  size_t n = pos * block_size_;
  for (unsigned i=0; i<block_size_; ++i) bitstring_.reset(n+i);
}

void qbitset::set()
{
  for (auto i=0; i<num_bits_; ++i) set(i);
}

void qbitset::set(const qbitset::size_t& pos)
{
  size_t n = pos * block_size_;
  for (unsigned i=0; i<block_size_; ++i) bitstring_[n+i] = maxbit_[i];
}

bool qbitset::raise(const qbitset::size_t& pos) 
{
  assert(pos>=0 || pos<num_bits_);
  size_t n = pos * block_size_;
  bitstring_t bit(0); 
  for (unsigned i=0; i<block_size_; ++i) {
    bit[i] = bitstring_[n+i];
  }
  auto p = bit.to_ulong();
  if (p < bitvals_.size()-1) {
    bit = bitvals_[p+1];
    for (unsigned i=0; i<block_size_; ++i)  bitstring_[n+i] = bit[i];
    //std::cout << "bit = " << bitstring_ << " " << p+1 << "\n";
    return true;
  }
  else return false;
}

bool qbitset::lower(const qbitset::size_t& pos) 
{
  assert(pos>=0 || pos<num_bits_);
  size_t n = pos * block_size_;
  bitstring_t bit(0); 
  for (unsigned i=0; i<block_size_; ++i) {
    bit[i] = bitstring_[n+i];
  }
  auto p = bit.to_ulong();
  if (p > 0) {
    bit = bitvals_[p-1];
    for (unsigned i=0; i<block_size_; ++i) bitstring_[n+i] = bit[i];
    return true;
  }
  else return false;
}

qbitset::size_t qbitset::operator[](const unsigned& pos) const
{
  assert(pos>=0 || pos<num_bits_);
  size_t n = pos * block_size_;
  bitstring_t bit(0); 
  for (unsigned i=0; i<block_size_; ++i) {
    bit[i] = bitstring_[n+i];
  }
  return bit.to_ulong();
}

qbitset::size_t qbitset::bitval(const unsigned& pos) const
{
  assert(pos>=0 || pos<num_bits_);
  size_t n = pos * block_size_;
  bitstring_t bit(0); 
  for (unsigned i=0; i<block_size_; ++i) {
    bit[i] = bitstring_[n+i];
  }
  return bit.to_ulong();
}

std::ostream& operator<<(std::ostream& os, const qbitset& qbits_)
{
  for (int i=qbits_.size()-1; i>=0; --i) {
    size_t n = i * qbits_.block_size_;
    for (int j=qbits_.block_size_-1; j>=0; --j) {
      os << qbits_.bitstring_[n+j];
    }
    os << " ";
  }
  return os;
}


//-------------- RotorBasis --------------
void RotorBasis::construct(const unsigned& num_sites, const int& qn_min, const int&qn_max)
{
  assert(qn_min <= qn_max);
  num_sites_ = num_sites;
  theta_min_ = qn_min;
  theta_max_ = qn_max;
  site_dim_ = (theta_max_-theta_min_)+1;
  ndim_ = std::pow(site_dim_, num_sites_);
  basis_states_.resize(ndim_);
  qbitset qbits_(site_dim_, num_sites_);
  //std::cout << num_sites_ << "\n"; getchar();

  qbits_.set();
  state_idx_.resize(qbits_.to_ullong()+1);
  null_idx_ = ndim_;
  for (auto& idx : state_idx_) idx = null_idx_;

  qbits_.reset();
  basis_states_[0] = qbits_;
  state_idx_[qbits_.to_ullong()] = 0;
  //std::cout << "|" << 0 << "> = " << qbits_ << " " << qbits_.to_ullong() << "\n"; 
  auto max_bit = site_dim_-1;
  for (auto n=1; n<ndim_; ++n) {
    for (auto i=0; i<num_sites_; ++i) {
      if (qbits_.bitval(i) == max_bit) {
        qbits_.reset(i);
      }
      else { 
        qbits_.raise(i);
        break;
      }
    }
    basis_states_[n] = qbits_;
    state_idx_[qbits_.to_ullong()] = n;
    //std::cout << "|" << n << "> = " << qbits_ << " " << qbits_.to_ullong() << "\n"; 
  }
}

RotorBasis::op_result RotorBasis::apply_cidag(const RotorBasis::idx_t& idx, 
  const unsigned& site) const
{
  auto new_state = basis_states_[idx];
  if (new_state.raise(site)) {
    //std::cout << new_state << " " << new_state.to_ullong() << "\n";
    return std::make_pair(1, state_idx_[new_state.to_ullong()]);
  }
  else return std::make_pair(1, null_idx_);
}

RotorBasis::op_result RotorBasis::apply_ni(const RotorBasis::idx_t& idx, 
  const unsigned& site) const
{
  int ntheta = theta_min_ + basis_states_[idx].bitval(site);
  return std::make_pair(ntheta, idx);
}


} // end namespace srmf
