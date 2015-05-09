#include "DenseNewtonOptimizer.h"
#include "MeritEvaluator.h"
#define ARMA_NO_DEBUG

#include <armadillo>

int DenseNewtonOptimizer::Optimize(){
  const int N = merit_evaluator.NumDOFs();
  const int MaxIt = 10;
  
  arma::vec grad(N);
  arma::mat Hess(N,N);
  arma::vec p(N,arma::fill::zeros);
  arma::vec x(N);

  merit_evaluator.GetCurrentState(x.memptr());

  const double merit0 = merit_evaluator.EvaluateMerit();
  

  
  for(int opt_iter = 0; opt_iter < MaxIt; opt_iter++){
    //merit_evaluator.EvaluateGradient(grad.memptr());
    merit_evaluator.EvaluateHessian(grad.memptr(),Hess.memptr());
    //std::cout << Hess << std::endl;
    //if(Hess.is_finite()){
    //std::cout << Hess << std::endl;
      //std::cout << grad << std::endl;
    bool good_solve = arma::solve(p,Hess,grad);
    p = -p;

    if(!good_solve){
      std::cout << Hess << std::endl;
    }
    if(debug){
      std::cout << "p: " << std::endl;
      std::cout << p.t() << std::endl;
    }

    //p = -grad;

    /*
    //bool good_solve = false;
    if(!good_solve || !p.is_finite()){
      p = -grad;
    }
    else{
      //p = p;
    }
    if(arma::norm(p,2) > 1){
      p/=arma::norm(p,2);
    }
    *///}
      //else{
      //p = grad;
      //}
    merit_evaluator.GetCurrentState(x.memptr());
    double alpha = backtrackingLineSearch(x.memptr(),p.memptr());

    if(debug){
      std::cout << merit_evaluator.EvaluateMerit() << " " << alpha << std::endl;
    }
    //x+= alpha*p;
    //merit_evalu
  }
  merit = merit_evaluator.EvaluateMerit();
  //std::cout << merit << " " << merit0 << std::endl;
  return merit < merit0;
  
}

double DenseNewtonOptimizer::backtrackingLineSearch(double* x, 
						    double* direction){

  const int N = merit_evaluator.NumDOFs();
  const double alpha_min = 1.0e-10;

  //std::cout << "N: " << N << std::endl;

  double alpha = alpha0;
  
  arma::vec xvec(x,N);

  
  const arma::vec p(direction,N);
  
  arma::vec gradk(N);
  //const double mer0 = merit_evaluator.EvaluateMerit();

  const double mer0 = merit_evaluator.EvaluateGradient(gradk.memptr());

  if(debug){
    std::cout << "mer0: " << mer0 << std::endl;
  }

  const double gradkp = arma::as_scalar(gradk.t()*p);
  
  arma::vec xnew = xvec + alpha*p;
  merit_evaluator.SetCurrentState(xnew.memptr());
  

  
  while(merit_evaluator.EvaluateMerit() > mer0 && alpha > alpha_min){
    //while(merit_evaluator.EvaluateMerit() > mer0 + c*alpha*gradkp && 
    //	alpha > alpha_min){
    //std::cout << alpha << std::endl;
    alpha*= rho;
    xnew = xvec + alpha*p;
    merit_evaluator.SetCurrentState(xnew.memptr());
    
  }
  
  if(debug){
  std::cout << "merit/mer0: " << merit_evaluator.EvaluateMerit() << " " <<
    mer0 << std::endl;
  }
  if(merit_evaluator.EvaluateMerit() > mer0){
    merit_evaluator.SetCurrentState(xvec.memptr());
  }
  //std::cout << merit_evaluator.EvaluateMerit() << " " << mer0 << 
  //  " " << alpha << std::endl;
  return alpha;
}


