#include "LBFGSOptimizer.h"
#include "MeritEvaluator.h"
#include <lbfgs.h>
#include <iostream>
#include <iomanip>

int PrintLBFGSFlag(int flag){
  using std::cout;
  using std::endl;
  cout << "flag: " << flag << endl;
  switch(flag){
  case LBFGS_SUCCESS:
    cout << "Success!" << endl;
    break;
  case LBFGS_ALREADY_MINIMIZED:
    cout << "Already Minimized" << endl;
    break;
  case LBFGSERR_UNKNOWNERROR:
    cout << "Unknown error" << endl;
    break;
  case LBFGSERR_LOGICERROR:
    cout << "Logic error" << endl;
    break;
  case LBFGSERR_OUTOFMEMORY:
    cout << "Out of memory" << endl;
    break;
  case LBFGSERR_CANCELED:
    cout << "Canceled" << endl;
    break;
  case LBFGSERR_INVALID_N:
    cout << "Invalid N" << endl;
    break;
  case LBFGSERR_INVALID_N_SSE:
    cout << "Invalid_N_SSE" << endl;
    break;
  case LBFGSERR_INVALID_X_SSE:
    cout << "Invalid_X_SSE" << endl;
    break;
  case LBFGSERR_INVALID_EPSILON:
    cout << "Invalid epsilon" << endl;
    break;
  case LBFGSERR_INVALID_TESTPERIOD:
    cout << "Invalid TESTPERIOD" << endl;
    break;
  case LBFGSERR_INVALID_DELTA:
    cout << "Invalid Delta" << endl;
    break;
  case LBFGSERR_INVALID_LINESEARCH:
    cout << "Invalid Linesearch" << endl;
    break;
  case LBFGSERR_INVALID_MINSTEP:
    cout << "Invalid Minstep" << endl;
    break;
  case LBFGSERR_INVALID_MAXSTEP:
    cout << "Invalid Maxstep: " << endl;
    break;
  case LBFGSERR_INVALID_FTOL:
    cout << "Invalid FTOL" << endl;
    break;
  case LBFGSERR_INVALID_WOLFE:
    cout << "Invalid Wolfe" << endl;
    break;
  case LBFGSERR_INVALID_GTOL:
    cout << "Invalid GTOL" << endl;
    break;
  case LBFGSERR_INVALID_XTOL:
    cout << "Invalid XTOL" << endl;
    break;
  case LBFGSERR_INVALID_MAXLINESEARCH:
    cout << "Invalid MaxLineSearch" << endl;
    break;
  case LBFGSERR_INVALID_ORTHANTWISE:
    cout << "Invalid Orthantwise" << endl;
    break;
  case LBFGSERR_INVALID_ORTHANTWISE_START:
    cout << "Invalid Orthantwise start" << endl;
    break;
  case LBFGSERR_INVALID_ORTHANTWISE_END:
    cout << "Invalid Orthantwise end" << endl;
    break;
  case LBFGSERR_OUTOFINTERVAL:
    cout << "Out Of Interval" << endl;
    break;
  case LBFGSERR_INCORRECT_TMINMAX:
    cout << "Incorretct Tmin/max" << endl;
    break;
  case LBFGSERR_ROUNDING_ERROR:
    cout << "Rounding Error" << endl;
    break;
  case LBFGSERR_MINIMUMSTEP:
    cout << "Minimum step" << endl;
    break;
  case LBFGSERR_MAXIMUMSTEP:
    cout << "Maximum step" << endl;
    break;
  case LBFGSERR_MAXIMUMLINESEARCH:
    cout << "Maximum line search" << endl;
    break;
  case LBFGSERR_MAXIMUMITERATION:
    cout << "Maximum Iteration" << endl;
    break;
  case LBFGSERR_WIDTHTOOSMALL:
    cout << "Width too small" << endl;
    break;
  case LBFGSERR_INVALIDPARAMETERS:
    cout << "Invalid Parameters" << endl;
    break;
  case LBFGSERR_INCREASEGRADIENT:
    cout << "Increase gradient" << endl;
    break;
  }



}
double  LBFGSOptimizer::evaluate(
		 const double *x,
		 double *g,
		 const int n,
		 const double step
		 )
{

  using std::cout;
  using std::endl;

  merit_evaluator.SetCurrentState(x);
 
 
  double merit = merit_evaluator.EvaluateGradient(g);

  

  return merit;
}

int LBFGSOptimizer::Optimize(){
  using std::cout;
  using std::endl;

  const int N = merit_evaluator.NumDOFs();
 

  lbfgs_parameter_t param;
  //double fx;
  double *state = lbfgs_malloc(N);

  
  merit_evaluator.GetCurrentState(state);


  lbfgs_parameter_init(&param);
  //param.max_step = 1.0;
  param.min_step = 0.0;
  //param.ftol = 0.1;
  param.past = 2;
  //param.epsilon = 1.0e-3;
  //param.delta = 1.0e-3;
  param.max_iterations = 100;
  //param.xtol = 0.0;

  param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
  //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
  //param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;

  //int ret = lbfgs(N, state, &fx, _evaluate, _progress, this, &param);
  int ret = lbfgs(N, state, &fx, _evaluate, NULL, this, &param);
  //plot_exit_flag = true;
  if(plot_exit_flag){
    PrintLBFGSFlag(ret);
  }
  /* Report the result. */
  //printf("L-BFGS optimization terminated with status code = %d\n", ret);
  //cout << std::setprecision(15) << "fx: " << fx << endl;

  lbfgs_free(state);
 
  return ret;
   
}

int LBFGSOptimizer::progress(
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
    
  std::cout << "Iteration: " << k << " fx: " << fx << 
    " step: " << step << std::endl;

  return 0;
}
