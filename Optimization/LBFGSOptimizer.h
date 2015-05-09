#pragma once
#include "Optimizer.h"
#include "MeritEvaluator.h"

class LBFGSOptimizer: public Optimizer{
 private:
  double fx;

 protected:
  static double _evaluate(
			  void *instance,
			  const double *x,
			  double *g,
			  const int n,
			  const double step
			  )
  {
    return reinterpret_cast<LBFGSOptimizer*>(instance)->evaluate(x, g, n, step);
  }
  static int _progress(
		       void *instance,
		       const double *x,
		       const double *g,
		       const double fx,
		       const double xnorm,
		       const double gnorm,
		       const double step,
		       int n,
		       int k,
		       int ls
		       )
  {
    return reinterpret_cast<LBFGSOptimizer*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
  }
  double evaluate(
		  const double *x,
		  double *g,
		  const int n,
		  const double step
		  );
  
  int progress(
	       const double *x,
	       const double *g,
	       const double fx,
	       const double xnorm,
	       const double gnorm,
	       const double step,
	       int n,
	       int k,
	       int ls
	       );
 
 public:
  bool plot_exit_flag=false;
 LBFGSOptimizer(MeritEvaluator& merit_evaluator_t): 
  Optimizer(merit_evaluator_t){}
  int Optimize();
  double getMerit() const{ return fx; }

};
