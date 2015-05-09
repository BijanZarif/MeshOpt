#pragma once
#include "Optimizer.h"

class MeritEvaluator;

class DenseNewtonOptimizer: public Optimizer{
 public:
 DenseNewtonOptimizer(MeritEvaluator& merit_evaluator, double alpha0=1,
		      double rho=0.5, double c = 1.0e-4): 
  Optimizer(merit_evaluator), alpha0(alpha0), rho(rho), c(c){}
  int Optimize();
  int setDebug(bool deb){ debug = deb; }
  
 private:
  double backtrackingLineSearch(double* x, double* direction);

  const double alpha0;
  const double rho;
  const double c;

  bool debug=false;

};
