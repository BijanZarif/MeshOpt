#include "NewtonOptimizer.h"
//#include "computeSolutionBelos.h"

//#include <armadillo>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <Teuchos_GlobalMPISession.hpp>

#include <Kokkos_DefaultNode.hpp>

#include <EpetraExt_RowMatrixOut.h>
#include <MatrixMarket_Tpetra.hpp>

typedef double ST;
typedef int LO;
typedef int GO;

typedef Kokkos::DefaultNode::DefaultNodeType node_type;
typedef Tpetra::MultiVector<ST,LO,GO, node_type> multivector_type;
typedef Tpetra::CrsMatrix<ST,LO,GO,node_type> sparse_mat_type;
typedef Tpetra::Map<LO,GO,node_type> map_type;

int NewtonOptimizer::Optimize(){
  using namespace std;

  const int ndof = merit_evaluator.NumDOFs();

  Teuchos::RCP<const Teuchos::Comm<int> > comm =     
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

 Teuchos::RCP<node_type> node = rcp(new node_type());
 //node->init(1);
 
  Teuchos::RCP<const map_type> data_map = 
    Tpetra::createContigMapWithNode<LO,GO>
    (ndof,ndof,comm,node);


  Teuchos::RCP<multivector_type> DU =
    rcp(new multivector_type(data_map,1,true) );

 Teuchos::RCP<multivector_type> GRAD =
    rcp(new multivector_type(data_map,1,true) );

  Teuchos::RCP<sparse_mat_type> HESS =
    rcp(new sparse_mat_type(data_map,0));

  //double mer = hessian_evaluator.EvaluateHessian(HESS,GRAD);
  arma::wall_clock timer;
  
 
  arma::vec state(ndof);
  arma::vec grad(ndof);
  arma::vec dx(ndof);

  const double merit0 = hessian_evaluator.EvaluateGradient(grad.memptr());
  double merit = merit0;
  cout << "MERIT 0: " << merit0 << endl;
  cout << "Ndof: " << ndof << endl;
  double merit_old = merit;
  int cnt=0;
  int good_steps = 0;
  double step = 1.0;
  double diff = 1.0/0.0;
  bool good_step = true;
  while(good_steps < 20 && step > 1.0e-4 & diff/merit > 1.0e-4){

    if(good_step){
      timer.tic();
      hessian_evaluator.EvaluateHessian(HESS,GRAD);
      cout << "Assemble time: " << timer.toc() << endl;
      //EpetraExt::RowMatrixToMatlabFile("testfile",HESS);
      Tpetra::MatrixMarket::Writer<sparse_mat_type>::writeSparseFile("testfile",HESS);
      break;
      timer.tic();
      // Solve the system
      //computeSolutionBelos(HESS,GRAD,DU,1,0,100,100,1.0e-10,false,"RILUK");
      cout << "Solve time: " << timer.toc() << endl;

      Teuchos::ArrayRCP<ST> DU_data = DU->get1dViewNonConst();
      for(int i=0; i<ndof; i++){
	dx[i] = DU_data[i];
      }

    }

    merit_evaluator.GetCurrentState(state.memptr());
    

    state-= step*dx;
    merit_evaluator.SetCurrentState(state.memptr());
    merit_old = merit;
 
    merit = merit_evaluator.EvaluateGradient(grad.memptr());
    cout << "merit: " << merit << " step: " << step << endl;
    if(merit < merit_old){
      diff = merit_old-merit;
      good_steps++;
      step = std::min(2.0*step,1.0);
      good_step = true;
    }
    else{
      state+= step*dx;
      merit_evaluator.SetCurrentState(state.memptr());
      step = 0.5*step;
      merit = merit_old;
      good_step = false;
    }

    cnt++;
  } 
  

  /*
  using std::cout;
  using std::endl;

  const int N = merit_evaluator.NumDOFs();
 
  arma::vec state(N);
  arma::vec grad(N);
  arma::vec grad_temp(N);
  arma::mat Hessian(N,N);
  
  grad.zeros();
  Hessian.zeros();

  const double merit0 = merit_evaluator.EvaluateGradient(grad.memptr());
  double merit = merit0;
  double merit_old = merit;
  int cnt=0;
  int good_steps = 0;
  double step = 1.0;
  bool good_step = true;
  while(good_steps < 5 && step > 1.0e-10){

    if(good_step){
      hessian_evaluator.EvaluateHessianDense(grad.memptr(),Hessian.memptr());
      
      //arma::mat Herm = (Hessian + Hessian.t())/2.0;
      //arma::cx_vec eigval = arma::eig_gen(Herm);
      //cout << eigval << endl;
    }

    merit_evaluator.GetCurrentState(state.memptr());
    
    //cout << "cond hess: " << arma::cond(Hessian) << endl;
    //cout << "Hessian Symmetry: " << norm(Hessian-Hessian.t(),2.0) << endl;
    // Hessian.save("Hessian",arma::raw_ascii);
    
    arma::vec dx = solve(Hessian,1.0*grad);
    
    
    state-= step*dx;
    merit_evaluator.SetCurrentState(state.memptr());
    merit_old = merit;
    merit = merit_evaluator.EvaluateGradient(grad_temp.memptr());
    cout << "merit: " << merit << " step: " << step << endl;
    if(merit < merit_old){
      good_steps++;
      step = std::min(2.0*step,1.0);
      good_step = true;
    }
    else{
      state+= step*dx;
      merit_evaluator.SetCurrentState(state.memptr());
      step = 0.5*step;
      merit = merit_old;
      good_step = false;
    }

    cnt++;
  } 
  */
  
}
