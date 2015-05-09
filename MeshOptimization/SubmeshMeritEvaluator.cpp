#include "SubmeshMeritEvaluator.h"
#include "Mel.h"
#include "Mnode.h"
#include "GeometryContainer.h"

#include <assert.h>

void SubmeshMeritEvaluator::GetCurrentState(double* state){
  //double* temp = nd->xyzptr();
  const double* temp = nd->getParametricCoords();
  //arma::vec3 temp = nd->xyzvec3();
  
  for(int i=0; i<nd->getType(); i++) state[i] = temp[i];
  //for(int i=0; i<2; i++) temp[i] = state[i];
  
  //state = nd->xyzptr();

}

void SubmeshMeritEvaluator::SetCurrentState(const double* state){
  
  double temp[] = {0.0,0.0,0.0};
  
  for(int i=0; i<nd->getType(); i++) temp[i] = state[i];

  //for(int i=0; i<2; i++) temp[i] = state[i];

  //nd->setXYZ(temp);

  
  arma::vec3 newpos;
  const short int entity = nd->getGeoEntity();
  if(nd->getType() != 3){
    GeometryContainer& geometry = optel_manager.getGeometry();
    if(nd->getType() == 1){
      std::vector<myEdge>& edges = geometry.getEdgesNC();
      myEdge& edge = edges[entity];
      //const double* u = nd->getParametricCoords();
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
      //const double* uv = nd->getParametricCoords();
      //nd->setXYZ(temp); // HACK!
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

double SubmeshMeritEvaluator::EvaluateMerit(){
  arma::vec grad(nd->getType());
  return EvaluateGradient(grad.memptr());
}

double SubmeshMeritEvaluator::EvaluateHessian(double* gradient, 
					      double* Hessian){
  double submesh_merit = 0.0;

  const int node_type = nd->getType();
  
  arma::vec gradVec(gradient,node_type,false);
  arma::mat Hess(Hessian,node_type,node_type,false);
  gradVec.zeros();
  Hess.zeros();
  
  for(auto it = elements.begin(); it != elements.end(); ++it){
    const MEl* el = it->el;

    std::unique_ptr<OptEl> optel = optel_manager.CreateOptEl(el,it->ideal);

    const gind* elnodes = el->getNodes();
    
    int currnode = nd->getND();

    int elnode = -1;
    for(int i = 0; i < el->NumNodes(); i++){
      if(elnodes[i] == currnode){
	elnode = i;
	break;
      }
    }
    assert(elnode != -1);
    
    
    arma::vec gradVec_temp = 0*gradVec;
    arma::mat Hess_temp = 0*Hess;
    submesh_merit += 
      optel->computeHessMeritParamNode(elnode,gradVec_temp,Hess_temp,factor);

    gradVec+= gradVec_temp;
    Hess+= Hess_temp;
  if(nd->getND() == 5562){
    //std::cout << arma::solve(Hess_temp,gradVec_temp) << std::endl;
  }
  }

  if(nd->getND() == 5562){
    //std::cout << gradVec << std::endl;
    //std::cout << arma::solve(Hess,gradVec) << std::endl;
  }

  return submesh_merit;
}

double SubmeshMeritEvaluator::EvaluateGradient(double* gradient){
  
  double submesh_merit = 0.0;

  const int node_type = nd->getType();
  arma::vec gradVec(gradient,node_type,false);
  gradVec.zeros();

  for(auto it = elements.begin(); it != elements.end(); ++it){
    const MEl* el = it->el;

    std::unique_ptr<OptEl> optel = optel_manager.CreateOptEl(el,it->ideal);

    const gind* elnodes = el->getNodes();
    
    int currnode = nd->getND();

    int elnode = -1;
    for(int i = 0; i < el->NumNodes(); i++){
      if(elnodes[i] == currnode){
	elnode = i;
	break;
      }
    }
    assert(elnode != -1);
    

    arma::vec grad_temp(node_type);
    grad_temp.zeros();
    
    submesh_merit += 
      optel->computeGradMeritParamNode(elnode,grad_temp,factor);
    gradVec+= grad_temp;

    if(currnode == 5562){
      //std::cout << grad_temp << std::endl;
    }
  }

  return submesh_merit;

  /*
  using std::cout;
  using std::endl;
  //std::cout << "in eval gradient" << std::endl;

  const node_map& nodes = optel_manager.getNodes();

  //cout << "begenning eval gradient" << endl;
  merits.resize(elements.size());
  double merit = 0.0;
  arma::mat gradMerit;
  //arma::vec3 gradNode = {0.0, 0.0, 0.0};
  arma::vec gradNode(nd->getType(),arma::fill::zeros);

  //arma::vec gradNode(2,arma::fill::zeros);


  if(elements.size() == 2){
    //std::cout << "element size is 2" << std::endl;
  }
  int order;
  int cnt=0;
  for(auto el = elements.begin(); el != elements.end(); ++el, ++cnt){
 
    std::unique_ptr<OptEl> optel = optel_manager.CreateOptEl(el->el,el->ideal);

    order = el->el->getDim();

    double temp, minDetS;
    //double mer = optel->computeGradMerit(gradMerit,merits[cnt], factor,minDetS);

    //double mer = optel->computeMerit();
    //gradMerit = optel->computeGradMeritFD();
    
    //std::cout << "mer: " << mer << std::endl;
    //std::cout << gradMerit << std::endl;

    const gind* elnodes = el->el->getNodes();
    const int nn = el->el->NumNodes();
    int fd_node=-1;
    for(int i=0; i<nn; i++){
      if(nd->getND() == elnodes[i]){
	fd_node = i;
      }
    }
    assert(fd_node != -1);
    double mer = optel->computeGradMeritParam(gradMerit,merits[cnt], factor,
    					      minDetS,fd_node);
    
    //std::cout << mer << std::endl;
    //std::cout << gradMerit << std::endl;
    
    //double mer = optel->computeGradMeritParam(gradMerit,merits[cnt], factor);
  
    //merits[cnt] = (sqrt(mer) + 1.0)/10.0;
    //const gind* cn = el->el->getCornerNodes();
    //const int ncn = el->el->numCornerNodes();



    bool found = false;
    for(int i=0, dofcnt=0; i<nn; i++){
 
      if(elnodes[i] == nd->getND()){
	//std::cout << nd->xyzvec3().t() << std::endl;
	//for(int j = 0; j < 2; j++){
	for(int j=0; j<nd->getType(); j++){
	  gradNode[j]+= gradMerit[dofcnt+j];
	  //std::cout << gradNode[j] << " ";
	}
	//std::cout << std::endl;
	found = true;
      }
      //dofcnt+= 3;
      dofcnt+= nodes.at(elnodes[i])->getType();
    }
    
    assert(found);

    merit+= mer;
  }

  //cout <<  "submesh merit: " << merit << endl;

  //cout << "endign eval gradient" << endl;
  arma::vec3 pos = nd->xyzvec3();
  //pos/=arma::norm(pos,2);
  arma::vec3 gradVec;
  
  //std::cout << nd->xyzvec3().t() << std::endl;
  //std::cout << "grad: " << std::endl;
  for(int i=0; i<nd->getType(); i++){
    gradient[i] = gradNode[i];
    gradVec[i] = gradient[i];
    if(nd->getType() == 3 && order != 1){
      //gradient[i] = pos[i];
    }
  }
  if(nd->getType() == 3 && order != 1){
    //std::cout << arma::dot(gradVec,pos) << std::endl;
  }
  //for(int i = 0; i < 2; i++) gradient[i] = gradNode[i];

  //if(merit != merit) merit = 1.0/0.0;

  return std::abs(merit);
  */
}
