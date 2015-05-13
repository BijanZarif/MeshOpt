#include "NewMeshOptimizer.h"

#include "MeshOptimizer.h"
#include "MeshContainer.h"
#include "OptElManager.h"
#include "elementFactory.h"
#include "IdealElement.h"
#include "LBFGSOptimizer.h"
#include "DenseNewtonOptimizer.h"
#include "PlotOptimizationProgress.h"

#include "SubmeshMeritEvaluator2.h"

#include <assert.h>

NewMeshOptimizer::
NewMeshOptimizer(const MeshContainer& mesh, 
		 OptElManager& optel_manager,
		 const std::map<const MEl*,arma::mat>& ideal_map,
		 int el_dim_to_opt):
  mesh(mesh), optel_manager(optel_manager), 
  el_dim_to_opt(el_dim_to_opt){

  //const element_set& elements = mesh.getElementsOfDimConst(el_dim_to_opt);

  int cnt = 0;
  for(auto it = ideal_map.begin(); it != ideal_map.end(); ++it){
    el2ideal[it->first] = cnt++;
    ideal_elements.push_back(CreateIdealElement(it->first,it->second));
  }

}

int NewMeshOptimizer::Optimize(){
  double merit0 = ComputeGlobalMerit();
  std::cout << "Initial merit: " << merit0 << std::endl;
  double merit_old;
  double merit = merit0;
  PlotOptimizationHeader();
  
  //std::cout << "Relative merit \t # active nodes \t diff" << std::endl;
  //std::cout << "--------------------------------------------------------\n";
  
  for(int iter = 0; iter < MaxIts; iter++){
    SetActiveNodes();
    if(active_nodes.size() == 0) break;
    
    OptimizeSubmesh();
    merit_old = merit;
    merit = ComputeGlobalMerit();
    double diff = (merit_old-merit)/merit_old;
    PlotOptimizationProgress(iter,active_nodes.size(),merit/merit0,diff,
			     GetMinQuality(),min_detJ);
    
    //std::cout << merit/merit0 << "\t\t" << active_nodes.size() << "\t\t" <<
    //  diff << std::endl;

    if(diff < 0.1){
      break;
    }

    //std::cout << "merit: " <<  << std::endl;
    //std::cout << "Active nodes sizes: " << active_nodes.size() << std::endl;
  }
}

int NewMeshOptimizer::OptimizeSubmesh(){
  typedef IdealElementBase::idealPair idealPair;
  for(auto nd = active_nodes.begin(); nd != active_nodes.end(); ++nd){


    std::vector<const idealPair*> ideal_pairs;
    for(auto el = node2el[*nd].begin(); el != node2el[*nd].end(); ++el){
      const auto& ideal = ideal_elements[el2ideal[*el]];
      for(int i = 0; i < ideal->getNumSubelements(); i++){
	ideal_pairs.push_back(&ideal->getIdealPair(i));
	if(ideal->getIdealPair(i).first->getDim() == 2){
	  //std::cout << "element dim is 2!" << std::endl;
	}
      }
    }

    //std::cout << "ideal pairs size: " << ideal_pairs.size() << std::endl;
    MNode* node = const_cast<MNode*>(mesh.getNodes().at(*nd).get());
    SubmeshMeritEvaluator2 
      submesh_merit_evaluator(node,ideal_pairs,optel_manager);
   
    //DenseNewtonOptimizer optimizer(submesh_merit_evaluator);
    //optimizer.Optimize();


    LBFGSOptimizer optimizer(submesh_merit_evaluator);
    optimizer.Optimize();
 
  }
}
double NewMeshOptimizer::ComputeGlobalMerit(){
  min_detJ = 1.0/0.0;
  distortions.resize(ideal_elements.size());
  double global_merit = 0.0;
  for(auto it = el2ideal.begin(); it != el2ideal.end(); ++it){
    const auto& pr = ComputeElementDistortion(it->first);
    double dist = pr.first;
    distortions[it->second] = dist;
    min_detJ = std::min(min_detJ,pr.second);
    
    global_merit+= std::pow(dist-1.0,2);
  }
  return global_merit;
}

double NewMeshOptimizer::GetMinQuality(){
  double min_quality = 1.0;
  
  distortions.resize(ideal_elements.size());
  for(auto it = el2ideal.begin(); it != el2ideal.end(); ++it){
    min_quality = std::min(min_quality,1.0/distortions[it->second]);
  }
  return min_quality;
}

int NewMeshOptimizer::InsertIdealElement(const MEl* el){

  el2ideal[el] = ideal_elements.size();
  ideal_elements.push_back(CreateIdealElement
			   (el,optel_manager.CreateIdealMatrix(el)));

  //ideal_elements.push_back(optel_manager.CreateIdealMatrix(el));
}

std::unique_ptr<IdealElementBase> NewMeshOptimizer::
CreateIdealElement(const MEl* el, const arma::mat& ideal){
  int type = el->getElementType();
  if(type == 3){
    return std::unique_ptr<IdealElementBase>(new QuadIdealElement(el,ideal));
  }
  else if(type == 6){
    return std::unique_ptr<IdealElementBase>(new PrismIdealElement(el,ideal));
  }
  else{
    return std::unique_ptr<IdealElementBase>(new IdealElement(el,ideal));
  }
  
}
int NewMeshOptimizer::SetActiveNodes(){
  active_elements.clear();
  active_nodes.clear();
  for(auto it = el2ideal.begin(); it != el2ideal.end(); ++it){
    if(distortions[it->second] > 1.0/quality_threshold){
      active_elements.insert(it->first);
      const gind* elnodes = it->first->getNodes();
      for(int i = 0; i < it->first->NumNodes(); i++){
	if(node_dims_to_opt[mesh.getNodes().at(elnodes[i])->getType()]){
	  active_nodes.insert(elnodes[i]);
	}
      }
    }
  }
  //std::cout << "after instering active nodes" << std::endl;
  
  const element_set& elements = mesh.getElementsOfDimConst(el_dim_to_opt);
  
  std::vector<bool> newnodes(mesh.getNodes().size(),false);
  
  for(auto it = elements.begin(); it != elements.end(); ++it){
    const gind* elnodes = it->get()->getNodes();
    for(int i = 0; i < it->get()->NumNodes(); i++){
      int currnd = elnodes[i];
      //active_nodes.insert(currnd);
      
      auto nodeit = active_nodes.find(currnd);
      if(nodeit != active_nodes.end()){
	
	auto optit = el2ideal.find(it->get());
	if(optit == el2ideal.end()){
	  InsertIdealElement(it->get());
	}
	assert(currnd < newnodes.size());
	
	auto fn = node2el.find(currnd);
	if(fn == node2el.end()){
	  newnodes[currnd] = true;
	  std::vector<const MEl*> ellist = {it->get()};
	  node2el.insert(std::make_pair(currnd,ellist));
	}
	else if(newnodes[currnd]){
	  fn->second.push_back(it->get());
	}
	
      }
    }
  }
}

std::pair<double,double> NewMeshOptimizer::
ComputeElementDistortion(const MEl* el){
  const auto& ideal = ideal_elements[el2ideal[el]];
  int N = ideal->getNumSubelements();
  double distortion = 0.0;
  double min_detJ = 1.0/0.0;
  for(int i = 0; i < N; i++){
    const auto& ideal_pair = ideal->getIdealPair(i);
    const auto& optel = 
      optel_manager.CreateOptEl(ideal_pair.first,ideal_pair.second);
    
    distortion+= optel->computeDistortion();
    min_detJ = std::min(min_detJ,optel->getMinDetJ());
  }
  return std::make_pair(distortion/N,min_detJ);
  //return distortion/N;
  
}
