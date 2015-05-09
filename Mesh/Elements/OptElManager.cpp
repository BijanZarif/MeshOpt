#include "OptElManager.h"
#include "OptEl.h"
#include "ShapeFunctionMatrices.h"

std::unique_ptr<OptEl> OptElManager::CreateOptEl(const MEl* el, 
						 const arma::mat& ideal){

  int dim = el->getDim();

  int order = el->getOrder();
  int degree_int = dim*(order-1) + dim*order;
  int type = el->getElementType();

  if(el->getElementType() == 3 || el->getElementType() == 4){
    degree_int = dim*(2*order-1)+dim*2*order;
  }
  //degree_int = 5*order-2;
  //degree_int = 3*order;
  //degree_int = 6*order-3;
  degree_int = 6*order-3;

  if(type == 3 || type == 6) degree_int = 6*order;

  //if(el->getElementType() == 3) degree_int = 4*order;

  int interp = interp_type;
  if(type == 2 && dim == 3){
    interp = 2;
  }
  
  const ShapeFunctionMatrices* sf = 
    sf_factory.getShapeFunction(el->getElementType(),
				el->getOrder(),interp,degree_int);

  const ShapeFunctionMatrices* sf_ideal = 
    sf_factory.getShapeFunction(el->getElementType(),
				1,0,degree_int);
  
  if(el->getDim() == 2){
    return std::unique_ptr<OptEl2D>(new OptEl2D(el,node_indexer,
    						nodes,sf,sf_ideal,
						ideal,&geometry));
  }
  else if(el->getDim() == 3){
    return std::unique_ptr<OptEl3D>(new OptEl3D(el,node_indexer,
    						nodes,sf,sf_ideal,
						ideal,&geometry));
  }
}
arma::mat OptElManager::CreateIdealMatrix(const MEl* el){

  ActiveMEl activeEl(el,node_indexer,nodes,NULL);
  return activeEl.getLinearNodes();
  //return activeEl.getNodesMatrix();
  
}
