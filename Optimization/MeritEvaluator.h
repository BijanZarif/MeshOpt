#pragma once

//#include "MeshContainer.h"
//#include "NodeIndexer.h"
//#include "ShapeFunctionMatrices.h"
//#include <lbfgs.h>

class MeritEvaluator{

 public:
  virtual int NumDOFs() = 0;
  virtual void GetCurrentState(double* state) = 0;
  virtual void SetCurrentState(const double* state) = 0;
  virtual double EvaluateMerit() = 0;
  virtual double EvaluateGradient(double* gradient) = 0;
  virtual double EvaluateHessian(double* gradient, double* Hessian) = 0;

  //virtual double EvaluateHessian(double* gradient, double* Hessian) = 0;

};


/*
class MeritEvaluator{
private:
  MeshContainer& mesh;
  ShapeFunctionMatricesFactory sf_factory;
  NodeIndexFactory index_factory;

  std::unordered_map<MNode*,int> active_nodes;
  std::unordered_set<MEl*> active_elements;

  arma::mat state;
  arma::mat gradient;
  double computeElementMerit(MEl* el,arma::mat& ideal,node_map& nodes);

public:
  MeritEvaluator(MeshContainer& mesh_t): mesh(mesh_t){}

  void Initialize(int Nlayers, double threshold, int type_to_opt);
  lbfgsfloatval_t EvaluateMerit();
  virtual double EvaluateGradient(lbfgsfloatval_t* grad);
  virtual void GetCurrentState(lbfgsfloatval_t* st);
  virtual void SetCurrentState(const lbfgsfloatval_t* st, double step);
  virtual const int NumDOFS() const{ return 3*active_nodes.size(); }

};
*/
