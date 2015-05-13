#include "HighOrderMeshGeneratorAndOptimizer.h"
#include "MeshContainer.h"
#include "GeometryContainer.h"
#include "HighOrderMeshGenerator.h"
#include "NodeIndexer.h"
#include "ShapeFunctionMatrices.h"
#include "OptElManager.h"
#include "SubmeshOptimizer.h"
#include "NewMeshOptimizer.h"
#include "IdealElement.h"

HighOrderMeshGeneratorAndOptimizer::
HighOrderMeshGeneratorAndOptimizer(MeshContainer& mesh, 
				   GeometryContainer&geometry):
  mesh(mesh), geometry(geometry){

}

std::shared_ptr<MeshContainer> HighOrderMeshGeneratorAndOptimizer::
GenerateAndOptimize(int order){
  HighOrderMeshGenerator ho_mesh_generator(mesh,geometry);
  
  int curr_order = (*mesh.getElements().begin())->getOrder();


  element_set& elements = mesh.getElementsOfDim(mesh.Dimension());
  for(auto el = elements.begin(); el != elements.end(); ++el){
    const MEl* element = el->get();

  }
  std::shared_ptr<MeshContainer> ho_mesh;
  for(int ord = curr_order+1; ord <= order; ord++){
    std::cout << "\nGenerating mesh of order: " << ord << std::endl;
    
    ho_mesh = ho_mesh_generator.generateHighOrderMeshRecursive(ord);

    mesh = *ho_mesh;

    std::map<const MEl*,arma::mat> ideal_elements;

    NodeIndexFactory index_factory;
    ShapeFunctionMatricesFactory sf_factory;
  
    OptElManager opt_el_manager(mesh.getNodesNC(),geometry,sf_factory,
				index_factory,mesh.getNodeSpacingType(),
				mesh.Dimension());

 
    for(int dim = 2; dim <=mesh.Dimension(); dim++){
      element_set& elements = mesh.getElementsOfDim(dim);
      //element_set& elements = mesh.getElementsOfDim(mesh.Dimension());
      //std::cout <<  "element size for dim 1: " << elements.size() << std::endl;
      std::cout << "Generating elements of dimension: " << dim << std::endl;
      
      for(auto el = elements.begin(); el != elements.end(); ++el){
	const MEl* element = el->get();
	int nc = element->NumChildren();
	bool is_curved = false;
	for(int ch = 0; ch < nc; ch++){
	  const MEl* child = element->getChild(ch);
	  if(child->hasGeoEntity()){
	    is_curved = true;
	    break;
	  }
	}
	//ideal_elements[element] = opt_el_manager.CreateIdealMatrix(element);

	//if(dim == 2){
	if(dim == 2 && element->hasGeoEntity()){
	  if(is_curved){
	    ideal_elements[element] = opt_el_manager.CreateIdealMatrix(element);
	  }
	}
	//else if(1){
	else if(dim == 3 && is_curved){
	  ideal_elements[element] = opt_el_manager.CreateIdealMatrix(element);
	}
	
      }


      std::cout << "Optimizing elements of dimension: " << dim << std::endl;
      SubmeshOptimizer mesh_optimizer(*ho_mesh,opt_el_manager,ideal_elements,
				      dim);
      
      
      std::vector<bool> to_opt(4,false);
      to_opt[dim] = true;
      //to_opt[1] = true;
      mesh_optimizer.Optimize(1.05,to_opt,std::set<unsigned short int>(),
      			      std::unordered_set<int>(),true);
      
      
      /*
      if(dim == mesh.Dimension()){
	NewMeshOptimizer new_mesh_optimizer(*ho_mesh,opt_el_manager,
					    ideal_elements,
					    dim);

	new_mesh_optimizer.setNodeDimToOpt(1);
	new_mesh_optimizer.setNodeDimToOpt(2);
	new_mesh_optimizer.setNodeDimToOpt(3);
	new_mesh_optimizer.setMaxIts(10);
      
	new_mesh_optimizer.Optimize();
      
      }
      */
    }
  }

  //mesh_optimizer.Optimize(1.2,{false,true,true,false});
  return ho_mesh;
}

