#include "HModel.h"
#include <iostream>
#include <memory>
#include <fstream>
//#include <omp.h>
#include "OCC_Handler.h"

//#include <Tpetra_DefaultPlatform.hpp>
//#include <Teuchos_GlobalMPISession.hpp>

using namespace std;

int main(int argc, char ** argv){

  //Teuchos::GlobalMPISession mpiSession (&argc, &argv);



  //unique_ptr<HModel> hmodel(new HModel);  
  HModel hmodel;

 


  try{
    hmodel.parseCommandLineOptions(argc,argv);    
  }
  catch(std::exception& e){
    cout << "std::runtime error thrown: " << e.what() << endl;
    return 1;
  }

 
  //hmodel->readOCC();
  hmodel.ReadGeometry();

  try{
    //hmodel->readGMSH();
    hmodel.readMesh();
  }
  catch(std::exception& e){
    cout << "std::runtime error thrown: " << e.what() << endl;
    return 1;
  }
  
  hmodel.generateBL();

  hmodel.PrepareMesh();

  hmodel.MeshHighOrder();

    
  hmodel.WriteMesh("GMSH");
  
 


  return 0;
}


