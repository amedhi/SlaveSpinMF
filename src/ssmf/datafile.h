/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-07-02 16:04:32
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-07-02 16:16:29
*----------------------------------------------------------------------------*/
#ifndef DATA_FILE_H
#define DATA_FILE_H

#include <string>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace file {

class DataFile 
{
public:
  DataFile() {}
  DataFile(const std::string& prefix, const std::string& name, 
    const bool& replace_mode=true);
  ~DataFile() {}
  void init(const std::string& prefix, const std::string& name, 
    const bool& replace_mode=true); 
  void init(const std::string& prefix, const std::string& name, 
    const std::string& heading, const bool& replace_mode=true); 
  void print_heading(const std::string& header); 
  void set_file_mode(const bool& replace_mode) { replace_mode_=replace_mode; }
  void open(void); 
  void close(void); 
  bool is_open(void) const { return fs_.is_open(); }
  std::ofstream& fs(void) { return fs_; }  
private:
  std::ofstream fs_;
  std::string dname_{""};
  std::string fname_{""};
  bool heading_printed_{false};
  bool replace_mode_{true};

  void set_fname(const std::string& prefix, const std::string& name);
};


} // end namespace ssmf

#endif
