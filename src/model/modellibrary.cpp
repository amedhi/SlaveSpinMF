/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-11 13:02:35
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-17 12:02:29
*----------------------------------------------------------------------------*/
#include <cmath>
#include "model.h"
#include <boost/algorithm/string.hpp>

namespace model {

int Hamiltonian::define_model(const input::Parameters& inputs, 
  const lattice::Lattice& lattice)
{
  //int info;
  //unsigned ntypes;
  //std::vector<MatrixElement> matrix_elem(20);
  double defval;
  //unsigned sitetype, change, type, src_type, tgt_type;
  std::string name, path; //, matrixelem, op, qn, site, src, tgt, fact;
  //SiteBasis site_basis;
  //BasisDescriptor basis;
  //QuantumNumber::value_type min, max, step;
  CouplingConstant cc;

  // define the models 
  model_name = inputs.set_value("model", "HUBBARD");
  boost::to_upper(model_name);

  strMatrix expr_mat;
  strMatrix::row_t expr_vec;
  //smat.getfromtxt("./matrix.txt");
  //smat(0,1) = "b";
  //smat(1,0) = "c";
  //smat(1,1) = "d";
  //std::cout << smat << "\n";

  if (model_name == "HUBBARD") {
    mid = model_id::HUBBARD;
    id2_ = model_id2::HUBBARD;
    if (lattice.id()==lattice::lattice_id::TBG) {
      set_TP_symmetry(true);
      add_parameter(name="t1", defval=1.0, inputs);
      add_parameter(name="t2", defval=1.0, inputs);
      add_parameter(name="tp1", defval=1.0, inputs);
      add_parameter(name="tp2", defval=1.0, inputs);
      add_parameter(name="mu", defval=0.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);
      // bond operator terms
      cc.create(5);
      expr_mat.resize(2,2);
      expr_mat(0,0) = "t1-tp1"; expr_mat(0,1) = "0";
      expr_mat(1,0) = "0"; expr_mat(1,1) = "t1+tp1";
      cc.add_type(0, expr_mat);
      expr_mat(0,0) = "t1+0.5*tp1"; expr_mat(0,1) = "(-0.866*tp1)";
      expr_mat(1,0) = "(-0.866*tp1)"; expr_mat(1,1) = "t1-0.5*tp1";
      cc.add_type(1, expr_mat);
      expr_mat(0,0) = "t1+0.5*tp1"; expr_mat(0,1) = "(0.866*tp1)";
      expr_mat(1,0) = "(0.866*tp1)"; expr_mat(1,1) = "t1-0.5*tp1";
      cc.add_type(2, expr_mat);
      expr_mat(0,0) = "t2"; expr_mat(0,1) = "tp2";
      expr_mat(1,0) = "-tp2"; expr_mat(1,1) = "t2";
      cc.add_type(3, expr_mat);
      expr_mat(0,0) = "t2"; expr_mat(0,1) = "-tp2";
      expr_mat(1,0) = "tp2"; expr_mat(1,1) = "t2";
      cc.add_type(4, expr_mat);
      add_bondterm(name="hopping", cc, op::spin_hop());
      // hubbard
      add_siteterm(name="hubbard", cc="U", op::hubbard_int());
    }

    else if (lattice.id()==lattice::lattice_id::CHAIN) {
      set_spinorbit_coupling(true);
      set_TP_symmetry(true);
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);
      int nowarn;
      add_parameter(name="ext_field", defval=0.0, inputs, nowarn);

      if (true) {
        // external site field
        cc.create(2);
        expr_vec.resize(2);
        expr_vec[0] = "-ext_field";
        expr_vec[1] = "ext_field";
        cc.add_type(0, expr_vec);
        expr_vec[0] = "ext_field";
        expr_vec[1] = "-ext_field";
        cc.add_type(1, expr_vec);
        add_siteterm(name="ExtField", cc, op::ni_up());
        // bond operators (Diagonal in either product basis or SOC basis)
        cc.create(2);
        int num_bands = 2;
        expr_mat.resize(num_bands,num_bands);
        for (int i=0; i<num_bands; ++i) {
          for (int j=0; j<num_bands; ++j) {
            if (i==j) expr_mat(i,j) = "-t";
            else expr_mat(i,j) = "0";
          }
        }
        cc.add_type(0, expr_mat);
        cc.add_type(1, expr_mat);
        add_bondterm(name="hopping", cc, op::upspin_hop());
      }
      else {
        // external site field
        cc.create(1);
        expr_vec.resize(2);
        expr_vec[0] = "-ext_field";
        expr_vec[1] = "ext_field";
        cc.add_type(0, expr_vec);
        add_siteterm(name="ExtField", cc, op::ni_up());
        // bond operators (Diagonal in either product basis or SOC basis)
        cc.create(1);
        int num_bands = 2;
        expr_mat.resize(num_bands,num_bands);
        for (int i=0; i<num_bands; ++i) {
          for (int j=0; j<num_bands; ++j) {
            if (i==j) expr_mat(i,j) = "-t";
            else expr_mat(i,j) = "0";
          }
        }
        cc.add_type(0, expr_mat);
        add_bondterm(name="hopping", cc, op::upspin_hop());
      }
    }

    else {
      set_TP_symmetry(true);
      // model parameters
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);
      // bond operator terms
      add_bondterm(name="hopping", cc="-t", op::spin_hop());
      add_siteterm(name="hubbard", cc="U", op::hubbard_int());
    }
  }

  else if (model_name == "TJ") {
    mid = model_id::TJ;
    id2_ = model_id2::TJ;
    int nowarn;
    if (inputs.set_value("no_double_occupancy",true,nowarn))
      set_no_dbloccupancy();
    // model parameters
    add_parameter(name="t", defval=1.0, inputs);
    add_parameter(name="J", defval=0.0, inputs);
    // bond operator terms
    add_bondterm(name="hopping", cc="-t", op::spin_hop());
    add_bondterm(name="exchange", cc="J", op::sisj_plus());
  }

  else if (model_name == "HUBBARD_NBAND") {
    mid = model_id::HUBBARD_NBAND;
    int num_bands;
    switch (lattice.id()) {
      case lattice::lattice_id::SQUARE_NBAND:
        id2_ = model_id2::HUBBARD_NBAND;
        // model parameters
        add_parameter(name="t", defval=1.0, inputs);
        add_parameter(name="U", defval=0.0, inputs);
        add_parameter(name="J", defval=0.0, inputs);
        //add_parameter(name="lambda", defval=0.0, inputs);
        // bond operator terms
        num_bands = inputs.set_value("num_bands", 1);
        cc.create(1);
        expr_mat.resize(num_bands,num_bands);
        for (int i=0; i<num_bands; ++i) {
          for (int j=0; j<num_bands; ++j) {
            if (i==j) expr_mat(i,j) = "-t";
            else expr_mat(i,j) = "0";
          }
        }
        cc.add_type(0, expr_mat);
        add_bondterm(name="hopping", cc, op::spin_hop());
        add_siteterm(name="hubbard", cc="U", op::hubbard_int());
        break;
      default:
        throw std::range_error("*error: modellibrary: model not defined for the lattice"); 
    }
  }

  else if (model_name == "HUBBARD_SQIRIDATE") {
    mid = model_id::HUBBARD_SQIRIDATE;
    if (lattice.id()==lattice::lattice_id::SQUARE_IRIDATE) {
      id2_ = model_id2::HUBBARD_SQIRIDATE;
      set_spinorbit_coupling(true);
      set_TP_symmetry(true);
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="lambda", defval=0.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);
      add_parameter(name="J", defval=0.0, inputs);
      int nowarn;
      add_parameter(name="ext_field", defval=0.0, inputs, nowarn);

      // external site field
      expr_vec.resize(6);
      for (int i=0; i<6; i+=2) expr_vec[i] = "-ext_field";
      for (int i=1; i<6; i+=2) expr_vec[i] = "ext_field";
      cc.create(1);
      cc.add_type(0, expr_vec);
      add_siteterm(name="ExtField", cc, op::ni_up());

      // SOC term in product basis (not diagonal)
      path = "/Users/amedhi/Projects/PhDs/ArunMaurya/PyrochloreIrdidate/ModelParameters/product_basis/";
      cc.create(1);
      expr_mat.resize(6,6);
      expr_mat.getfromtxt(path+"soc_matrix.txt");
      cc.add_type(0,expr_mat);
      add_siteterm(name="spin_flip", cc, op::spin_flip());

      // bond operators (Diagonal in either product basis or SOC basis)
      int num_bands = 6;
      cc.create(2);
      expr_mat.resize(num_bands,num_bands);
      for (int i=0; i<num_bands; ++i) {
        for (int j=0; j<num_bands; ++j) {
          if (i==j) expr_mat(i,j) = "-t";
          else expr_mat(i,j) = "0";
        }
      }
      cc.add_type(0, expr_mat);
      cc.add_type(1, expr_mat);
      add_bondterm(name="hopping", cc, op::upspin_hop());
    }
    else {
      throw std::range_error("*error: modellibrary: model not defined for the lattice"); 
    }
  }

  else if (model_name == "TBI_HUBBARD") {
    mid = model_id::TBI_HUBBARD;
    switch (lattice.id()) {
      case lattice::lattice_id::SQUARE_2BAND:
        id2_ = model_id2::BHZ;
        set_TP_symmetry(true);
        add_parameter(name="e0", defval=1.0, inputs);
        add_parameter(name="t", defval=1.0, inputs);
        add_parameter(name="tsp", defval=1.0, inputs);
        // site operators
        expr_vec.resize(2);
        expr_vec[0] = "-(e0-4*t)";
        expr_vec[1] = "(e0-4*t)";
        cc = expr_vec;
        add_siteterm(name="ni_up", cc, op::ni_up());
        add_siteterm(name="ni_dn", cc, op::ni_dn());

        // bond operators
        cc.create(2);
        expr_mat.resize(2,2);
        expr_mat(0,0) = "-t"; expr_mat(0,1) = "-i*tsp";
        expr_mat(1,0) = "-i*tsp"; expr_mat(1,1) = "t";
        cc.add_type(0, expr_mat);
        expr_mat(0,0) = "-t"; expr_mat(0,1) = "-tsp";
        expr_mat(1,0) = "tsp"; expr_mat(1,1) = "t";
        cc.add_type(1, expr_mat);
        add_bondterm(name="hopping", cc, op::upspin_hop());

        expr_mat(0,0) = "-t"; expr_mat(0,1) = "i*tsp";
        expr_mat(1,0) = "i*tsp"; expr_mat(1,1) = "t";
        cc.add_type(0, expr_mat);
        expr_mat(0,0) = "-t"; expr_mat(0,1) = "-tsp";
        expr_mat(1,0) = "tsp"; expr_mat(1,1) = "t";
        cc.add_type(1, expr_mat);
        add_bondterm(name="hopping", cc, op::dnspin_hop());

        /*
        cc = CouplingConstant({0,"-(e0-4*t)"}, {1,"(e0-4*t)"});
        add_siteterm(name="hopping", cc, op::ni_up());
        add_siteterm(name="hopping", cc, op::ni_dn());
        // upspin hop
        cc = CouplingConstant({0,"-t"}, {1,"-i*tsp"},
          {2,"-i*tsp"}, {3,"t"}, {4,"-tsp"},{5,"tsp"});
        add_bondterm(name="hopping", cc, op::upspin_hop());
        // dnspin hop
        cc = CouplingConstant({0,"-t"}, {1,"i*tsp"},
          {2,"i*tsp"}, {3,"t"}, {4,"-tsp"},{5,"tsp"});
        add_bondterm(name="hopping", cc, op::dnspin_hop());
        */
        // Hubbard interaction
        add_parameter(name="U", defval=0.0, inputs);
        add_siteterm(name="hubbard", cc="U", op::hubbard_int());
        break;

      case lattice::lattice_id::HONEYCOMB:
        id2_ = model_id2::KM;
        /*
        add_parameter(name="t", defval=1.0, inputs);
        add_parameter(name="lambda", defval=1.0, inputs);
        // upspin hop
        cc = CouplingConstant({0,"-t"}, {1,"i*lambda"},
          {2,"-i*lambda"}, {3,"-i*lambda"}, {4,"i*lambda"});
        add_bondterm(name="hopping", cc, op::upspin_hop());
        // dnspin hop
        cc = CouplingConstant({0,"-t"}, {1,"-i*lambda"},
          {2,"i*lambda"}, {3,"i*lambda"}, {4,"-i*lambda"});
        add_bondterm(name="hopping", cc, op::dnspin_hop());
        // Hubbard interaction
        add_parameter(name="U", defval=0.0, inputs);
        add_siteterm(name="hubbard", cc="U", op::hubbard_int());
        */
        break;

      case lattice::lattice_id::KAGOME:
        /*
        add_parameter(name="t", defval=1.0, inputs);
        add_parameter(name="t2", defval=1.0, inputs);
        add_parameter(name="lambda", defval=1.0, inputs);
        add_parameter(name="lambda2", defval=1.0, inputs);
        // upspin hop
        cc = CouplingConstant({0,"-t+i*lambda"},{1,"-t-i*lambda"},
          {2,"-t2+i*lambda2"},{3,"-t2-i*lambda2"});
        add_bondterm(name="hopping", cc, op::upspin_hop());
        // dnspin hop
        cc = CouplingConstant({0,"-t-i*lambda"},{1,"-t+i*lambda"},
          {2,"-t2-i*lambda2"},{3,"-t2+i*lambda2"});
        add_bondterm(name="hopping", cc, op::dnspin_hop());
        // Hubbard interaction
        add_parameter(name="U", defval=0.0, inputs);
        add_siteterm(name="hubbard", cc="U", op::hubbard_int());
        */
        break;

      case lattice::lattice_id::PYROCHLORE_V1:
        /*
        add_parameter(name="t", defval=1.0, inputs);
        add_parameter(name="t2", defval=1.0, inputs);
        add_parameter(name="lambda", defval=1.0, inputs);
        add_parameter(name="lambda2", defval=1.0, inputs);
        add_parameter(name="th", defval=1.0, inputs);
        // upspin hop
        cc = CouplingConstant({0,"-t+i*lambda"},{1,"-t-i*lambda"},
          {2,"-t2+i*lambda2"},{3,"-t2-i*lambda2"},{4,"-th"});
        add_bondterm(name="hopping", cc, op::upspin_hop());
        // dnspin hop
        cc = CouplingConstant({0,"-t-i*lambda"},{1,"-t+i*lambda"},
          {2,"-t2-i*lambda2"},{3,"-t2+i*lambda2"},{4,"-th"});
        add_bondterm(name="hopping", cc, op::dnspin_hop());
        // Hubbard interaction
        add_parameter(name="U", defval=0.0, inputs);
        add_siteterm(name="hubbard", cc="U", op::hubbard_int());
        */
        break;

      case lattice::lattice_id::PYROCHLORE_3D:
        id2_ = model_id2::PYROCHLORE;
        set_spinorbit_coupling(true);
        set_TP_symmetry(true);
        add_parameter(name="t", defval=1.0, inputs);
        add_parameter(name="lambda", defval=0.0, inputs);
        add_parameter(name="U", defval=0.0, inputs);
        add_parameter(name="J", defval=0.0, inputs);

        // site term
        expr_vec.resize(6);
        expr_vec[0] = "lambda"; expr_vec[1] = "lambda";
        expr_vec[2] = "-0.5*lambda"; expr_vec[3] = "-0.5*lambda";
        expr_vec[4] = "-0.5*lambda"; expr_vec[5] = "-0.5*lambda";
        cc.create(4);
        cc.add_type(0, expr_vec);
        cc.add_type(1, expr_vec);
        cc.add_type(2, expr_vec);
        cc.add_type(3, expr_vec);
        add_siteterm(name="onsite", cc, op::ni_up());

        // bond operators
        path = "/Users/amedhi/Projects/PhDs/ArunMaurya/PyrochloreIrdidate/hoppings/";
        cc.create(9);
        expr_mat.resize(6,6);
        expr_mat.getfromtxt(path+"intracell_01.txt");
        cc.add_type(0, expr_mat);
        expr_mat.getfromtxt(path+"intracell_02.txt");
        cc.add_type(1, expr_mat);
        expr_mat.getfromtxt(path+"intracell_03.txt");
        cc.add_type(2, expr_mat);
        expr_mat.getfromtxt(path+"intracell_12.txt");
        cc.add_type(3, expr_mat);
        expr_mat.getfromtxt(path+"intracell_13.txt");
        cc.add_type(4, expr_mat);
        expr_mat.getfromtxt(path+"intracell_23.txt");
        cc.add_type(5, expr_mat);

        expr_mat.getfromtxt(path+"intercell_10.txt");
        cc.add_type(6, expr_mat);
        expr_mat.getfromtxt(path+"intercell_20.txt");
        cc.add_type(7, expr_mat);
        expr_mat.getfromtxt(path+"intercell_30.txt");
        cc.add_type(8, expr_mat);
        add_bondterm(name="hopping", cc, op::upspin_hop());

        // Huubard U
        add_siteterm(name="hubbard", cc="U", op::hubbard_int());
        break;

      default:
        throw std::range_error("*error: modellibrary: model not defined for the lattice"); 
    }
  }

  else if (model_name == "PYROCHLORE") {
    mid = model_id::PYROCHLORE;
    id2_ = model_id2::PYROCHLORE;
    switch (lattice.id()) {
      case lattice::lattice_id::PYROCHLORE_3D:
        set_spinorbit_coupling(true);
        set_TP_symmetry(true);
        add_parameter(name="t", defval=1.0, inputs);
        add_parameter(name="lambda", defval=0.0, inputs);
        add_parameter(name="U", defval=0.0, inputs);
        add_parameter(name="J", defval=0.0, inputs);

        // SOC term in product basis (not diagonal)
        path = "/Users/amedhi/Projects/PhDs/ArunMaurya/PyrochloreIrdidate/ModelParameters/product_basis/";
        cc.create(4);
        expr_mat.resize(6,6);
        expr_mat.getfromtxt(path+"soc_matrix.txt");
        cc.add_type(0,expr_mat);
        cc.add_type(1,expr_mat);
        cc.add_type(2,expr_mat);
        cc.add_type(3,expr_mat);
        add_siteterm(name="spin_flip", cc, op::spin_flip());

        // bond operators
        path = "/Users/amedhi/Projects/PhDs/ArunMaurya/PyrochloreIrdidate/ModelParameters/product_basis/";
        cc.create(9);
        expr_mat.resize(6,6);
        expr_mat.getfromtxt(path+"hopping_intracell_01.txt");
        cc.add_type(0, expr_mat);
        expr_mat.getfromtxt(path+"hopping_intracell_02.txt");
        cc.add_type(1, expr_mat);
        expr_mat.getfromtxt(path+"hopping_intracell_03.txt");
        cc.add_type(2, expr_mat);
        expr_mat.getfromtxt(path+"hopping_intracell_12.txt");
        cc.add_type(3, expr_mat);
        expr_mat.getfromtxt(path+"hopping_intracell_13.txt");
        cc.add_type(4, expr_mat);
        expr_mat.getfromtxt(path+"hopping_intracell_23.txt");
        cc.add_type(5, expr_mat);

        expr_mat.getfromtxt(path+"hopping_intercell_10.txt");
        cc.add_type(6, expr_mat);
        expr_mat.getfromtxt(path+"hopping_intercell_20.txt");
        cc.add_type(7, expr_mat);
        expr_mat.getfromtxt(path+"hopping_intercell_30.txt");
        cc.add_type(8, expr_mat);
        add_bondterm(name="hopping", cc, op::upspin_hop());

        // Huubard U
        add_siteterm(name="hubbard", cc="U", op::hubbard_int());
        break;
      default:
        throw std::range_error("*error: modellibrary: model not defined for the lattice"); 
    }
  }



  /*------------- undefined model--------------*/
  else {
    throw std::range_error("*error: modellibrary: undefined model");
  }

  // if the model has site disorder
  /*
  if (site_disorder) {
    add_disorder_term(name="disorder", op::ni_sigma());
  }*/
  
  return 0;
}

int Hamiltonian::construct(const input::Parameters& inputs, 
  const lattice::Lattice& lattice)
{
  init(lattice);
  define_model(inputs, lattice);
  finalize(lattice);
  return 0;
}


} // end namespace model
