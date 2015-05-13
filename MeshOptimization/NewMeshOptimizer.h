#pragma once

class MeshContainer;
class OptElManager;
class MEl;

#include <map>
#include <set>
#include <vector>
#include <armadillo>
#include <memory>

class IdealElementBase;

class NewMeshOptimizer{

 public:
  typedef std::unique_ptr<IdealElementBase> idealPtr;
  
  NewMeshOptimizer(const MeshContainer& mesh, 
		   OptElManager& optel_manager,
		   const std::map<const MEl*,arma::mat>& ideal_map,
		   int el_dim_to_opt);


  
  void setNodeDimToOpt(int dim){ node_dims_to_opt[dim] = true; }
  void setMaxIts(int max){ MaxIts = max; }

  int Optimize();

 private:
  int OptimizeSubmesh();
  double ComputeGlobalMerit();
  int SetActiveNodes();
  double GetMinQuality();
  int InsertIdealElement(const MEl* el);
  std::pair<double,double> ComputeElementDistortion(const MEl* el);

  std::unique_ptr<IdealElementBase> CreateIdealElement(const MEl* el, 
						       const arma::mat& ideal);

  const MeshContainer& mesh;
  OptElManager& optel_manager;
  const int el_dim_to_opt;
  int node_dim_to_opt=2;
  int node_dims_to_opt[4] = {false};
  
  std::vector<double> distortions;
  
  //std::map<const MEl*,idealPair> ideal_elements;
  
  std::vector<idealPtr> ideal_elements;
  //std::vector<unsigned int> ideal_index_map;
  std::map<const MEl*, int> el2ideal;
  std::set<const MEl*> active_elements;
  std::set<int> active_nodes;
  std::map<int,std::vector<const MEl*> > node2el;

  int MaxIts = 20;
  double StopTol = 1.0e-2;
  double quality_threshold = 0.9;
  double min_detJ;

  //std::map<const MEl*,unsigned int> ideal_index_map;
  
};
