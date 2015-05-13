#include "MyParameterList.h"
//#include <boost/mpi.hpp>
#include <iostream>

void MyParameterList::Initialize(const parameter_map_type& plist){
 
  FindInvalidParameters(plist,GetValidParameters());
  Populate(plist);
  //Broadcast();
  
}

void MyParameterList::
FindInvalidParameters(const parameter_map_type& plist,
		      const std::set<std::string>& valid) const{
  for(auto it = plist.begin(); it != plist.end(); ++it){
    auto fn = valid.find(it->first);
    if(fn == valid.end()){
      std::cout << "valid " << getListName() << " parameters are: "<< std::endl;
      std::cout << "[";
     for(auto it = valid.begin(); it != valid.end(); ++it){
       std::cout << *it;
       if(it != --valid.end()) std::cout << ", ";
      }
     std::cout << "]." << std::endl;
      throw std::runtime_error("Invalid parameter: " + it->first + "!");
    }

  }

}
