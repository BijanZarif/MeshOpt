#include "CaseParameterList.h"
#include "ParameterManip.h"

std::set<std::string> CaseParameterList::GetValidParameters() const{
  return {"GMSHFileName","STEPFileName","OutputFileName"};
}

void CaseParameterList::Populate(const parameter_map_type& pmap){
  ParameterManip::OverrideString(gmshFileName,"GMSHFileName",pmap);
  ParameterManip::OverrideString(stepFileName,"STEPFileName",pmap);
  ParameterManip::OverrideString(outputFileName,"OutputFileName",pmap);
}
