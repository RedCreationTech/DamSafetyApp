#pragma once
// Test-only CSV adapter. Does not stand in for MOOSE framework integration tests.
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
namespace MooseUtils {
class DelimitedFileReader {
 public:
 enum class HeaderFlag {ON}; enum class FormatFlag {COLUMNS};
 explicit DelimitedFileReader(const std::string &p):path(p){}
 void setHeaderFlag(HeaderFlag){} void setFormatFlag(FormatFlag){}
 void read(){
  std::ifstream f(path);if(!f)throw std::runtime_error("missing CSV: "+path);
  std::string line; std::getline(f,line); std::stringstream head(line);std::string cell;
  while(std::getline(head,cell,',')){if(!cell.empty()&&cell.back()=='\r')cell.pop_back();names.push_back(cell);}
  data.resize(names.size());
  while(std::getline(f,line)) {if(line.empty()||line=="\r")continue;std::stringstream row(line);std::size_t i=0;
   while(std::getline(row,cell,',')){if(i>=data.size())throw std::runtime_error("extra column");data[i++].push_back(std::stod(cell));}
   if(i!=data.size())throw std::runtime_error("missing column");
  }
 }
 const std::vector<std::string>& getNames()const{return names;}
 const std::vector<std::vector<double>>& getData()const{return data;}
 const std::vector<double>&getData(const std::string& n)const{auto p=std::find(names.begin(),names.end(),n);if(p==names.end())throw std::runtime_error("missing column");return data[p-names.begin()];}
 private:std::string path;std::vector<std::string>names;std::vector<std::vector<double>>data;
};}
