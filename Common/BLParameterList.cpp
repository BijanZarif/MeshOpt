#include "BLParameterList.h"
#include "ParameterManip.h"
#include <iostream>

std::set<std::string> BLParameterList::GetValidParameters() const{
  return {"Thickness","Re","DistortionRatio","NLayers","FlowType",
      "MinSpacing","GrowthRatio","GenerateBL","HasBL","BLSurfID",
      "BLTermID","MaxBLGenIters"};
}

void BLParameterList::Populate(const parameter_map_type& pmap){
  ParameterManip::OverrideDouble(thickness,"Thickness",pmap);
  ParameterManip::OverrideDouble(Re,"Re",pmap);
  ParameterManip::OverrideDouble(minSpacing,"MinSpacing",pmap);
  ParameterManip::OverrideDouble(growthRatio,"GrowthRatio",pmap);
  ParameterManip::OverrideInteger(distortionRatio,"DistortionRatio",pmap);
  ParameterManip::OverrideInteger(NLayers,"NLayers",pmap);
  ParameterManip::OverrideString(flowType,"FlowType",pmap);
  ParameterManip::OverrideBool(GenerateBL,"GenerateBL",pmap);
  ParameterManip::OverrideBool(HasBL,"HasBL",pmap);
  ParameterManip::OverrideInteger(BLSurfID,"BLSurfID",pmap);
  ParameterManip::OverrideInteger(BLTermID,"BLTermID",pmap);
  ParameterManip::OverrideInteger(MaxBLGenIters,"MaxBLGenIters",pmap);
  
}
