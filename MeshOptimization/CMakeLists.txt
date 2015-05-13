SET(S
  MeshOptimizer.cpp
  NewMeshOptimizer.cpp
#  GlobalMeshOptimizer.cpp
  SubmeshOptimizer.cpp
  LocalMeshMeritEvaluator.cpp
  SubmeshMeritEvaluator.cpp
  SubmeshMeritEvaluator2.cpp
  IdealElement.cpp
  PlotOptimizationProgress.cpp
)

SET(H
  MeshOptimizer.h
  NewMeshOptimizer.h
#  GlobalMeshOptimizer.h
  SubmeshOptimizer.h
  LocalMeshMeritEvaluator.h
  SubmeshMeritEvaluator.h
  SubmeshMeritEvaluator2.h
  IdealElement.h
  PlotOptimizationProgress.h
)

append_src(S H)
