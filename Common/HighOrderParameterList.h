#pragma once

#include "MyParameterList.h"

class HighOrderParameterList: public MyParameterList{
 public:
  HighOrderParameterList(){}
  HighOrderParameterList(const parameter_map_type& pmap){ Initialize(pmap); }

  int order = 2;

  double targetMinQuality = 0.9;

  //bool useOptimizedPoints = 0;
    
  bool OptimizedNodeSpacing = false;

 private:
  std::set<std::string> GetValidParameters() const;
  void Populate(const parameter_map_type& pmap);

  void Broadcast(){}

  std::string getListName() const{ return "HighOrderParameterList"; }

};
