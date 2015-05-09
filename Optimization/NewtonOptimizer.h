#pragma once
#include "Optimizer.h"
#include "LocalMeshMeritEvaluator.h"

class NewtonOptimizer: public Optimizer{
 private:
  LocalMeshMeritEvaluator& hessian_evaluator;
 protected:

 public:

 NewtonOptimizer(LocalMeshMeritEvaluator& merit_evaluator_t): 
  Optimizer(merit_evaluator_t), hessian_evaluator(merit_evaluator_t){}
  int Optimize();

};
