#pragma once

#include "MeritEvaluator.h"
#include "GlobalDefines.h"
#include "OptElManager.h"
#include <vector>
#include <armadillo>

class MNode;
class MEl;

class SubmeshMeritEvaluator2: public MeritEvaluator{
 public:
  typedef std::pair<const MEl*,const arma::mat> idealPair;
  SubmeshMeritEvaluator2(MNode* nd,
			 const std::vector<const idealPair*>& elements,
			 OptElManager& optel_manager);

  void SetDimToOpt(int dim){ dim_to_opt[dim] = true; }

  int NumDOFs(){ return nd->getType(); }
  void GetCurrentState(double* state);
  void SetCurrentState(const double* state);
  double EvaluateMerit();
  double EvaluateGradient(double* gradient);
  double EvaluateHessian(double* gradient, double* Hessian);

 private:
  MNode* nd;
  const std::vector<const idealPair*>& elements;
  OptElManager& optel_manager;
  bool dim_to_opt[4] = {false};

};
