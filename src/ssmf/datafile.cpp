/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-07-02 16:04:32
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-07-04 10:36:31
*----------------------------------------------------------------------------*/
#include "datafile.h"
#include <locale> 

namespace file {

DataFile::DataFile(const std::string& prefix, const std::string& name, 
  const bool& replace_mode)
{
  init(prefix,name,replace_mode);  
}

void DataFile::init(const std::string& prefix, const std::string& name, 
  const bool& replace_mode)
{
  dname_ = name;
  set_fname(prefix,name);
  replace_mode_ = replace_mode;
}

void DataFile::init(const std::string& prefix, const std::string& name, 
  const std::string& header, const bool& replace_mode)
{
  init(prefix,name,replace_mode);
  print_heading(header);
}

void DataFile::set_fname(const std::string& prefix, const std::string& name)
{
  fname_ = prefix+name;
  std::locale loc;
  for (int i=0; i<fname_.size(); ++i)
    fname_[i] = std::tolower(fname_[i],loc);
  auto pos = fname_.find('^');
  if (pos != std::string::npos) fname_.erase(pos,1);
  //fname_ = "res_"+fname_+".txt";
  fname_ = fname_+".txt";
}

void DataFile::print_heading(const std::string& header)
{
  if (heading_printed_) return;
  if (!replace_mode_) return;
  if (!fs_.is_open()) open();
  fs_ << header;
  fs_ << "# Results:\n";
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << std::flush;
  heading_printed_ = true;
  close();
}

void DataFile::open(void) 
{
  if (fs_.is_open()) return;
  if (replace_mode_) {
    fs_.open(fname_);
    replace_mode_ = false;
  }
  else fs_.open(fname_,std::ios::app);
  if (!fs_.is_open()) 
    throw std::runtime_error("DataFile::open: file open failed");
}

void DataFile::close(void) 
{
  fs_.close();
}

} // end namespace file

