#pragma once

int PlotOptimizationHeader();

int PlotOptimizationProgress(int iteration, 
			     int num_active,
			     double rel_merit,
			     double diff,
			     double min_quality,
			     double min_detj);

