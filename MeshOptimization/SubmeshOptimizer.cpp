#include "SubmeshOptimizer.h"
#include "SubmeshMeritEvaluator.h"
#include "LocalMeshMeritEvaluator.h"
#include "Mel.h"
#include "MeshContainer.h"
#include "LBFGSOptimizer.h"
#include "DenseNewtonOptimizer.h"
#include "PlotOptimizationProgress.h"

using namespace std;

int SubmeshOptimizer::RandomizeNodes(double factor){
  
  std::unordered_set<const MEl*> watch;
  std::map<gind,gind> activenodes;
  
 
  node_map& nodes = mesh.getNodesNC();
  //const ideal_map& idealElements = mesh.getIdealElements();
  //std::map<MEl*,double>& all_merits = mesh.getAllMeritsNC();

  FindActiveElements(watch,activenodes,0,2,-1);
  for(auto nd = activenodes.begin(); nd != activenodes.end(); ++nd){
    arma::vec3 xyz = nodes.at(nd->first)->getXYZ();
    for(int i=0; i<2; i++){
      xyz(i)+= factor*((double) std::rand() / (RAND_MAX));
    }
    nodes.at(nd->first)->setXYZ(xyz.memptr());
  }

}

int SubmeshOptimizer::OptimizeGlobal(double threshold, int Nnearest, 
				     int type, double factor){

  std::unordered_set<const MEl*> watch;
  std::map<gind,gind> activenodes;
  //const ideal_map& idealElements = mesh.getIdealElements();
  //std::map<MEl*,double>& all_merits = mesh.getAllMeritsNC();

  FindActiveElements(watch,activenodes,threshold,Nnearest,type);

  cout << "watch size: " << watch.size() << endl;
  cout << "active size: " << activenodes.size() << endl;

  std::vector<const MEl*> temp(watch.size());
  int cnt = 0;
  for(auto el = watch.begin(); el != watch.end(); ++el, ++cnt){
    temp[cnt] = *el;
  }
  LocalMeshMeritEvaluator lm_evaluator(activenodes,temp,
				       idealElements,optel_manager,factor);
  
  LBFGSOptimizer global_optimizer(lm_evaluator);
  global_optimizer.plot_exit_flag=true;
  global_optimizer.Optimize();
  {
    int cnt = 0;
    const std::vector<double>& merits = lm_evaluator.getMerits();
    for(auto el = watch.begin(); el != watch.end(); ++el, ++cnt){
      all_merits[*el] = merits[cnt];
    }
  }
}

double SubmeshOptimizer::OptimizeOneByOne(double threshold, int Nnearest,
				       int type, double factor){

  std::unordered_set<const MEl*> watch;

  activenodes.clear();
  
  const node_map& nodes = mesh.getNodes();
  //const ideal_map& idealElements = mesh.getIdealElements();
  //std::map<MEl*,double>& all_merits = mesh.getAllMeritsNC();
  

  arma::wall_clock timer;
  timer.tic();
  FindActiveElements(watch,activenodes,threshold,Nnearest,type);



  if(activenodes.size() == 0){
    return ComputeMerit();
  }

  std::map<gind, std::vector<const MEl*> > Node2ElementList = 
    CreateNode2ElementList(watch,activenodes);
 
  //std::cout << "active nodes size: " << activenodes.size() << std::endl;

  double global_merit = 0.0;
  
  for(auto nd = activenodes.begin(); nd != activenodes.end(); ++nd){

    //std::cout << "h0" << std::endl;
    const std::vector<const MEl*>& elements = Node2ElementList[nd->first];
    const int sz = elements.size();

    //std::cout << "h1" << std::endl;
    arma::vec3 pos_old = nodes.at(nd->first)->xyzvec3();
    const double* pold = nodes.at(nd->first)->getParametricCoords();
    arma::vec3 parametric_old;
    for(int i = 0; i < nodes.at(nd->first)->getType(); i++){
      parametric_old[i] = pold[i];
    } 

    //std::cout << "h2" << std::endl;
    //assert(elements.size() == 2);
    arma::vec debpos = {0.00178445,0.00558035,0};
    bool debug = false;
    bool to_opt = true;

    if(to_opt){
    std::vector<ElIdealPair> temp;
    for(int i=0; i<sz; i++){
      temp.emplace_back(elements[i],idealElements.at(elements[i]));
    }

    double merit0 = 0;
    for(int i = 0; i < elements.size(); i++){
      merit0+= all_merits[elements[i]];
    }

    //std::cout << "nodes: " << std::endl;
    //std::cout << nodes.at(nd->first)->xyzvec3() << std::endl;

    //assert(elements[0] != elements[1]);

    SubmeshMeritEvaluator evaluator(nodes.at(nd->first).get(),
				    temp,
				    optel_manager, all_merits, factor);

    //std::cout << "h4" << std::endl;

    //std::cout << "before optimization" << std::endl;
    //DenseNewtonOptimizer optimizer(evaluator);
    //optimizer.setDebug(debug);
    //optimizer.Optimize();


    //std::cout << "h5" << std::endl;
    LBFGSOptimizer optimizer(evaluator);
    optimizer.Optimize();

    global_merit+= optimizer.getMerit();

    //std::cout << "After optimize" << std::endl;

    //const std::vector<double>& merits = evaluator.getMerits();
    //for(int i=0; i<merits.size(); i++){
    //  all_merits[elements[i]] = merits[i];
    //}


    }
  }


  /*
  double glob_merit = 0.0;
  int cnt=0;
  for(auto el = all_merits.begin(); el != all_merits.end(); ++el){
    if(el->first->getDim() == el_dim_to_opt){
      double fac = 1.0;
      //if(el->first->IsBL()) fac = factor;
      glob_merit+= pow(fac*(el->second-1.0),2);
      //std::cout << el->second << std::endl;
      cnt++;
    }
  }
  */
  return global_merit;
  //return glob_merit/cnt;

  //return activenodes.size();
}
int SubmeshOptimizer::Optimize(double threshold, std::vector<bool> to_opt,
			       std::set<unsigned short int> geo_faces_to_opt_t,
			       std::unordered_set<int> clamped_nodes_t,
			       bool only_interior_t){

  geo_faces_to_opt = geo_faces_to_opt_t;

  //only_opt_node_with_bctag = only_opt_node_with_bctag_t;

  dims_to_opt = to_opt;

  clamped_nodes = clamped_nodes_t;

  only_interior = only_interior_t;




  /*
  const node_map& nodes = mesh.getNodes();
  for(auto nd = nodes.begin(); nd != nodes.end(); ++nd){
    cout << nd->second->getType() << endl;
  }
  */
  //cout << "in  optimize one by one" << endl;

  //RandomizeNodes(0.03);
  //return 1;

  arma::wall_clock timer;
  timer.tic();
  arma::wall_clock inner_timer;

  
  //FindActiveElements(watch,activenodes,threshold,0,-1);

  OptimizeOneByOne(threshold,0,-1,1);

  double merit0 = ComputeMerit();
  std::cout << "Initial Merit: " << merit0 << std::endl;

  PlotOptimizationHeader();
  
  //cout << "initial min Mesh quality: " << minMeshQuality() << endl;
  //std::cout << "initial min detS: " << minMeshDetS() << std::endl;
  //std::cout << "initial min detJ: " << minMeshDetJ() << std::endl;
  
  int max_it = 20;
  //int Nactive = 0;
  //int Nactive_old;
  double Merit = merit0;
  double Merit_old;
  for(int i=0; i<max_it; i++){
    if(activenodes.size() == 0){
      PlotOptimizationProgress(i,activenodes.size(),Merit/merit0,0,
			       minMeshQuality(),minMeshDetJ());
      break;
    }
    
    inner_timer.tic();
    //Nactive_old = Nactive;
    Merit_old = Merit;
    //Merit = OptimizeGlobal(threshold,1,-1,1);
    //activenodes.clear();
    //FindActiveElements(watch,activenodes,threshold,0,-1);
    Merit = OptimizeOneByOne(threshold,0,-1,1);


    UpdateMeshMerits();
    
    //OptimizeGlobal(threshold,2,-1,1.0e0);
    //cout << "Merit" << Merit << endl;
    double diff = std::abs(Merit_old-Merit)/Merit_old;
    double min_quality = minMeshQuality();
    PlotOptimizationProgress(i,activenodes.size(),Merit/merit0,diff,
			     min_quality,minMeshDetJ());
    
    //double diff = 1.0;
    //cout << "diff: " << diff << endl;
    //std::cout << "min quality: " << minMeshQuality() << std::endl;
    //std::cout << "min detS: " << minMeshDetS() << std::endl;
    //std::cout << "min detJ: " << minMeshDetJ() << std::endl;

    if(diff < 1.0e-1 && min_quality > 1.0e-2){
      break;
    }
   
  }
  
  
  UpdateMeshMerits();


  
}

/*
int SubmeshOptimizer::SnapToBLExtent(){

  std::unordered_set<MEl*> watch;
  std::map<gind,gind> activenodes;

  const node_map& nodes = mesh.getNodes();
  std::map<MEl*,double>& all_merits = mesh.getAllMeritsNC();
  const ideal_map& idealElements = mesh.getIdealElements();
  const element_set& elements = mesh.getElements();

  std::map<gind, std::vector<MEl*> > node2elements;
  std::unordered_set<gind> bl_edge_nodes;

  const double min_qual0 = minMeshQuality();
  cout << "min qual 0: " << min_qual0 << endl;


  for(auto nd = newnode_map.begin(); nd != newnode_map.end(); ++nd){
    bl_edge_nodes.insert(nd->second);
  }

 
  watch.clear();
  for(auto el = elements.begin(); el != elements.end(); ++el){
    MEl* currel = el->get();
    if(dims_to_opt[currel->getDim()]){
      const int ncn = currel->numCornerNodes();
      const gind* cn = currel->getCornerNodes();
      for(int i=0; i<ncn; i++){
	auto it = bl_edge_nodes.find(cn[i]);
	if(it != bl_edge_nodes.end()){
	  if(1){
	    //if(dims_to_opt[nodes.at(*it)->getType()]){
	    watch.insert(currel);
	    node2elements[*it].push_back(currel);
	  }
	}
      }
    }
  }

    std::map<gind,arma::vec3> oldpos;
  
    for(auto nd = normal_map.begin(); nd != normal_map.end(); ++nd){
	const gind sn = nd->first;
	const gind bln = newnode_map.at(nd->first);
	MNode* BLNode = nodes.at(bln).get();
	MNode* SurfNode = nodes.at(sn).get();
	if(1){
	  //if(dims_to_opt[BLNode->getType()]){
	  arma::vec3 pos = SurfNode->getXYZ();
	  oldpos[bln] = BLNode->getXYZ();
	  pos-= nd->second.dist*nd->second.normal;
	  BLNode->setXYZ(pos.memptr());
	}
    }

    arma::wall_clock timer;
    timer.tic();
    UpdateMeshMerits(watch);
    cout << "Update merit time: " << timer.toc() << endl;


    while(minMeshQuality() < 0.5*min_qual0){
      //for(auto nd = normal_map.begin(); nd != normal_map.end(); ++nd){
      //const gind bln = newnode_map.at(nd->first);
      for(auto nd = node2elements.begin(); nd != node2elements.end(); ++nd){
	const std::vector<MEl*>& els = nd->second;
	
	bool updated = false;
	for(auto el = els.begin(); el != els.end(); ++el){
	  if(all_merits[*el] > 2.0/min_qual0){
	    MNode* BLNode = nodes.at(nd->first).get();
	    BLNode->setXYZ(oldpos[nd->first].memptr());
	    updated = true;
	    break;

	  }
	}
	if(updated){
	  arma::mat gradMerit;
	  for(auto el = els.begin(); el != els.end(); ++el){
	    std::unique_ptr<OptEl> optel = 
	      optel_manager.CreateOptEl(*el,idealElements.at(*el));
	    double detS;
	    optel->computeGradMerit(gradMerit,all_merits[*el],detS);
	  }

	}

      }
      }
    cout << "Min quality: " << minMeshQuality() << endl;
    
      cout << "after moving back" << endl;
 
}
*/
