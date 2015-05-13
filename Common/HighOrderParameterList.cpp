#include "HighOrderParameterList.h"
#include "ParameterManip.h"

std::set<std::string> HighOrderParameterList::GetValidParameters() const{
  return {"Order","TargetMinQuality","OptimizedNodeSpacing"};
}

void HighOrderParameterList::Populate(const parameter_map_type& pmap){
  ParameterManip::OverrideInteger(order,"Order",pmap);
  ParameterManip::OverrideDouble(targetMinQuality,"TargetMinQuality",pmap);
  ParameterManip::OverrideBool(OptimizedNodeSpacing,"OptimizedNodeSpacing",
				  pmap);
}
