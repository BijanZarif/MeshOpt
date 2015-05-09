#include "SubmeshMeritEvaluator2.h"
#include "Mel.h"
#include "Mnode.h"
#include "GeometryContainer.h"

SubmeshMeritEvaluator2::
SubmeshMeritEvaluator2(MNode* nd,
		       const std::vector<const idealPair*>& elements,
		       OptElManager& optel_manager):
  nd(nd), elements(elements), optel_manager(optel_manager){

}

void SubmeshMeritEvaluator2::GetCurrentState(double* state){
  const double* temp = nd->getParametricCoords();
  for(int i=0; i<nd->getType(); i++) state[i] = temp[i];  
}

void SubmeshMeritEvaluator2::SetCurrentState(const double* state){
  double temp[] = {0.0,0.0,0.0};
  
  for(int i=0; i<nd->getType(); i++) temp[i] = state[i];


  
  arma::vec3 newpos;
  const short int entity = nd->getGeoEntity();
  if(nd->getType() != 3){
    GeometryContainer& geometry = optel_manager.getGeometry();
    if(nd->getType() == 1){
      std::vector<myEdge>& edges = geometry.getEdgesNC();
      myEdge& edge = edges[entity];
      if(edge.areParamsValid(state)){
	edge.param2xyz(state,newpos.memptr());
	nd->setParametricCoords(temp);
	nd->setXYZ(newpos.memptr());
      }
      else{
	//std::cout << state[0] << std::endl;
	//std::cout << "edge params are not valid" << std::endl;
      }
    }
    else if(nd->getType() == 2){
      std::vector<myFace>& faces = geometry.getFacesNC();
      myFace& face = faces[entity];
      if(face.areParamsValid(state)){
	face.param2xyz(state,newpos.memptr());
	nd->setParametricCoords(temp);
	nd->setXYZ(newpos.memptr());
	//nd->setXYZ(temp);
      }
      else{
	//std::cout << "face params are not valid" << std::endl;
      }
    }
    //nd->setXYZ(newpos.memptr());
  }
  else{
    nd->setParametricCoords(temp);
    nd->setXYZ(temp);
  }
}

double SubmeshMeritEvaluator2::EvaluateMerit(){
  arma::vec grad(nd->getType());
  return EvaluateGradient(grad.memptr());
}

double SubmeshMeritEvaluator2::EvaluateGradient(double* gradient){
  
  double submesh_merit = 0.0;

  const int node_type = nd->getType();
  arma::vec gradVec(gradient,node_type,false);
  gradVec.zeros();

  for(auto it = elements.begin(); it != elements.end(); ++it){
    const MEl* el = (*it)->first;
    const arma::mat& ideal = (*it)->second;

    std::unique_ptr<OptEl> optel = optel_manager.CreateOptEl(el,ideal);

    const gind* elnodes = el->getNodes();
    
    int currnode = nd->getND();

    int elnode = -1;
    for(int i = 0; i < el->NumNodes(); i++){
      if(elnodes[i] == currnode){
	elnode = i;
	break;
      }
    }
    //assert(elnode != -1);
    
    if(elnode >= 0){
      arma::vec grad_temp(node_type);
      grad_temp.zeros();
    
      submesh_merit+= 
	optel->computeGradMeritParamNode(elnode,grad_temp,1);
      gradVec+= grad_temp;
    }

  }

  return submesh_merit;
}

double SubmeshMeritEvaluator2::EvaluateHessian(double* gradient, 
					       double* Hessian){
 double submesh_merit = 0.0;

  const int node_type = nd->getType();
  
  arma::vec gradVec(gradient,node_type,false);
  arma::mat Hess(Hessian,node_type,node_type,false);
  gradVec.zeros();
  Hess.zeros();
  
  for(auto it = elements.begin(); it != elements.end(); ++it){
    const MEl* el = (*it)->first;
    const arma::mat& ideal = (*it)->second;
    
    std::unique_ptr<OptEl> optel = optel_manager.CreateOptEl(el,ideal);

    const gind* elnodes = el->getNodes();
    
    int currnode = nd->getND();

    int elnode = -1;
    for(int i = 0; i < el->NumNodes(); i++){
      if(elnodes[i] == currnode){
	elnode = i;
	break;
      }
    }
    
    if(elnode >= 0){
      arma::vec gradVec_temp = 0*gradVec;
      arma::mat Hess_temp = 0*Hess;
      submesh_merit += 
	optel->computeHessMeritParamNode(elnode,gradVec_temp,Hess_temp,1);

      gradVec+= gradVec_temp;
      Hess+= Hess_temp;
    }
  }

  return submesh_merit;
}
