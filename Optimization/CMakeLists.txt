SET(S
  LBFGSOptimizer.cpp
  DenseNewtonOptimizer.cpp
  #NewtonOptimizer.cpp
)

SET(H
  Optimizer.h
  LBFGSOptimizer.h
  DenseNewtonOptimizer.h
  #NewtonOptimizer.h
)

add_library(Optimization ${S} ${H})
target_include_directories(Optimization PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(Optimization_INC ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

install(TARGETS Optimization DESTINATION lib)
install(FILES ${H} DESTINATION include)

#append_src(S H)
