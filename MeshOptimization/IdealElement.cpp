#include "IdealElement.h"

#include "Mel.h"
#include "elementFactory.h"

IdealElement::IdealElement(const MEl* el, const arma::mat& ideal): 
  ideal_pair(std::make_pair(el,ideal)){}

QuadIdealElement::QuadIdealElement(const MEl* el, const arma::mat& ideal){

  //std::cout << "In QuadIdealElement constructor" << std::endl;

  arma::mat subideal(3,3);
  int order = el->getOrder();
  int tridof = (order+1)*(order+2)/2;
  //arma::mat subtri(tridof,3);
  
  std::vector<int> subtriind(tridof);
  
  std::vector<int> subtri(tridof);
  const gind* quadnodes = el->getNodes();

  // First sub triangle
  subideal.col(0) = ideal.col(0);
  subideal.col(1) = ideal.col(1);
  subideal.col(2) = ideal.col(2);
  
  for(int i = 0, cnt = 0; i <= order; i++){
    for(int j = 0; j <= (order-i); j++, cnt++){
      subtriind[cnt] =  i*(order+1)+j;
      subtri[cnt] =  quadnodes[i*(order+1)+j];

    }
  }
  //std::cout << arma::Col<int>(subtriind).t() << std::endl;

  elptrs[0] = elementFactory::Instance()->
    CreateElement(2,subtri,NULL,order,el->getBCTag(),el->getGeoEntity(),
		  el->getGeoType());

  ideal_pairs.emplace_back(std::make_pair(elptrs[0].get(),subideal));
  

  
  // Second sub triangle
  subideal.col(0) = ideal.col(1);
  subideal.col(1) = ideal.col(3);
  subideal.col(2) = ideal.col(0);

  for(int i = 0, cnt = 0; i <=order; i++){
    for(int j = 0; j <= (order-i); j++, cnt++){
      subtriind[cnt] = (order-i) + j*(order+1);
      subtri[cnt] = quadnodes[(order-i) + j*(order+1)];
    }
  }
  //std::cout << arma::Col<int>(subtriind).t() << std::endl;

  elptrs[1] = elementFactory::Instance()->
    CreateElement(2,subtri,NULL,order,el->getBCTag(),el->getGeoEntity(),
		  el->getGeoType());

  ideal_pairs.emplace_back(std::make_pair(elptrs[1].get(),subideal)); 

  // Third sub triangle
  subideal.col(0) = ideal.col(2);
  subideal.col(1) = ideal.col(0);
  subideal.col(2) = ideal.col(3);
  
  for(int i = 0, cnt = 0; i <=order; i++){
    for(int j = 0; j <= (order-i); j++, cnt++){
      subtriind[cnt] = (order-j)*(order+1) + i;
      subtri[cnt] = quadnodes[(order-j)*(order+1) + i];
    }
  }
  //std::cout << arma::Col<int>(subtriind).t() << std::endl;

  elptrs[2] = elementFactory::Instance()->
    CreateElement(2,subtri,NULL,order,el->getBCTag(),el->getGeoEntity(),
		  el->getGeoType());

  ideal_pairs.emplace_back(std::make_pair(elptrs[2].get(),subideal)); 

  // Fourth sub triangle
  subideal.col(0) = ideal.col(3);
  subideal.col(1) = ideal.col(2);
  subideal.col(2) = ideal.col(1);
  
  for(int i = 0, cnt = 0; i <=order; i++){
    for(int j = 0; j <= (order-i); j++, cnt++){
      subtriind[cnt] = (order-j) + (order-i)*(order+1);
      subtri[cnt] = quadnodes[(order-j) + (order-i)*(order+1)];
    }
  }
  //std::cout << arma::Col<int>(subtriind).t() << std::endl;

  elptrs[3] = elementFactory::Instance()->
    CreateElement(2,subtri,NULL,order,el->getBCTag(),el->getGeoEntity(),
		  el->getGeoType());

  ideal_pairs.emplace_back(std::make_pair(elptrs[3].get(),subideal)); 

}

PrismIdealElement::PrismIdealElement(const MEl* el, const arma::mat& ideal){
  arma::mat subideal(3,4);
  int order = el->getOrder();
  int tridof = (order+1)*(order+2)*(order+3)/6;
  //arma::mat subtri(tridof,3);
  
  //std::vector<int> subtriind(tridof);
  
  std::vector<int> subtet(tridof);
  const gind* quadnodes = el->getNodes();

  int subindices[6][4] = {
    {0,1,2,3},
    {1,2,0,4},
    {2,0,1,5},
    {3,5,4,0},
    {4,3,5,1},
    {5,4,3,2}
  };

  for(int se = 0; se < 6; se++){
    for(int nd = 0; nd < 4; nd++){
      subideal.col(nd) = ideal.col(subindices[se][nd]);
      subtet[nd] = quadnodes[subindices[se][nd]];
    }
    elptrs[se] = elementFactory::Instance()->
      CreateElement(4,subtet,NULL,order,el->getBCTag(),el->getGeoEntity(),
		    el->getGeoType());

    ideal_pairs.emplace_back(std::make_pair(elptrs[se].get(),subideal));

  }
  /*
  // First sub triangle
  subideal.col(0) = ideal.col(0);
  subideal.col(1) = ideal.col(1);
  subideal.col(2) = ideal.col(2);
  
  for(int i = 0, cnt = 0; i <= order; i++){
    for(int j = 0; j <= (order-i); j++, cnt++){
      subtriind[cnt] =  i*(order+1)+j;
      subtri[cnt] =  quadnodes[i*(order+1)+j];

    }
  }
  //std::cout << arma::Col<int>(subtriind).t() << std::endl;

  elptrs[0] = elementFactory::Instance()->
    CreateElement(2,subtri,NULL,order,el->getBCTag(),el->getGeoEntity(),
		  el->getGeoType());

  ideal_pairs.emplace_back(std::make_pair(elptrs[0].get(),subideal));
  */
}
