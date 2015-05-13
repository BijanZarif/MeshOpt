#include "PlotOptimizationProgress.h"
#include <iostream>
#include <iomanip>

int PlotOptimizationHeader(){
  std::cout << std::left << std::setw(6) << std::setfill(' ') << "Iter";
  std::cout << std::left << std::setw(12) << std::setfill(' ') << "# nodes";
  std::cout << std::left << std::setw(12) << std::setfill(' ') << "Rel. Merit";
  std::cout<< std::left << std::setw(12) << std::setfill(' ') << "Diff";
  std::cout<< std::left << std::setw(12) << std::setfill(' ') << "Min Q";
  std::cout<< std::left << std::setw(12) << std::setfill(' ') << "Min |J|";
  std::cout << std::endl;
  std::cout << 
    "------------------------------------------------------------" << std::endl;
}

int PlotOptimizationProgress(int iteration, 
			     int num_active,
			     double rel_merit,
			     double diff,
			     double min_quality,
			     double min_detj){

  std::cout << std::left << std::setw(6) << std::setfill(' ') << iteration;
  std::cout << std::left << std::setw(12) << std::setfill(' ') << num_active;
  std::cout << std::left << std::setw(12) << std::setfill(' ') << rel_merit;
  std::cout << std::left << std::setw(12) << std::setfill(' ') << diff;
  std::cout << std::left << std::setw(12) << std::setfill(' ') << min_quality;
  std::cout << std::left << std::setw(12) << std::setfill(' ') << min_detj;
  std::cout << std::endl;

}
