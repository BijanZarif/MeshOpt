#pragma once
#include "MyParameterList.h"

class CaseParameterList: public MyParameterList{
 public:
  CaseParameterList(){};
  CaseParameterList(const parameter_map_type& pmap){ Initialize(pmap); }

  std::string gmshFileName;
  std::string stepFileName;
  std::string outputFileName = "meshopt.msh";

 private:
  std::set<std::string> GetValidParameters() const;
  void Populate(const parameter_map_type& pmap);

  void Broadcast(){}

  std::string getListName() const{ return "CaseParameterList"; }
};
