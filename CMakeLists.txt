cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(MeshOpt CXX)


set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;*.la;")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#SET(Trilinos_PREFIX ${Trilinos_Install_DIR} "/lib/cmake/Trilinos")
#SET(CMAKE_PREFIX_PATH ${Trilinos_PREFIX})
# Get Trilinos as one entity
#FIND_PACKAGE(Trilinos REQUIRED)

#set(Boost_USE_STATIC_LIBS   ON)
#FIND_PACKAGE(Boost COMPONENTS program_options serialization mpi chrono iostreams REQUIRED)

include(ExternalProject)

if(NOT BLAS_LAPACK_LIBS)
  ExternalProject_Add(OpenBLAS
    GIT_REPOSITORY https://github.com/xianyi/OpenBLAS.git
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/OpenBLAS
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/OpenBLAS
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/OpenBLAS/src/OpenBLAS
    CONFIGURE_COMMAND  ""
    UPDATE_COMMAND ""
    BUILD_COMMAND make USE_THREAD=0 NO_AFFINITY=1
    INSTALL_COMMAND make PREFIX=${CMAKE_INSTALL_PREFIX}/OpenBLAS install
  )
set(BLAS_LAPACK_LIBS libopenblas.a)
set(BLAS_LAPACK_INC_DIR ${CMAKE_INSTALL_PREFIX}/OpenBLAS/include)
set(BLAS_LAPACK_LIB_DIR ${CMAKE_INSTALL_PREFIX}/OpenBLAS/lib)
endif(NOT BLAS_LAPACK_LIBS)

ExternalProject_Add(Armadillo
  URL  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/armadillo-5.000.1.tar.gz
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/Armadillo
  INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/Armadillo
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}/Armadillo
)

set(ARMADILLO_ROOT ${CMAKE_INSTALL_PREFIX}/Armadillo)

ExternalProject_Add(LBFGS
  GIT_REPOSITORY https://github.com/chokkan/liblbfgs.git
  CONFIGURE_COMMAND ./autogen.sh && ./configure --prefix=${CMAKE_INSTALL_PREFIX}/LBFGS
  INSTALL_COMMAND make install
  BUILD_IN_SOURCE 1
  UPDATE_COMMAND ""
  BUILD_COMMAND make
  INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/LBFGS
)

set(LBFGS_ROOT ${CMAKE_INSTALL_PREFIX}/LBFGS)

SET(CMAKE_CXX_FLAGS "--std=c++11 -O3 -DARMA_DONT_USE_WRAPPER -DARMA_USE_LAPACK -DARMA_USE_BLAS ${CMAKE_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++ -lgfortran")

macro(append_src FILES FILES2)
  foreach(FILE ${${FILES}})
    list(APPEND LIST ${CMAKE_CURRENT_SOURCE_DIR}/${FILE})
  endforeach(FILE)
  set(SOURCES ${SOURCES};${LIST} PARENT_SCOPE)
  foreach(FILE ${${FILES2}})
    list(APPEND LIST2 ${CMAKE_CURRENT_SOURCE_DIR}/${FILE})
  endforeach(FILE)
  set(HEADERS ${HEADERS};${LIST2} PARENT_SCOPE) 
  set(INCDIRS ${INCDIRS};${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
endmacro(append_src)



macro(move_src_up_level)
  set(SOURCES ${SOURCES} PARENT_SCOPE)
  set(HEADERS ${HEADERS} PARENT_SCOPE)
  set(INCDIRS ${INCDIRS} PARENT_SCOPE)
endmacro(move_src_up_level)

INCLUDE_DIRECTORIES(${OCC_ROOT}/inc)
#INCLUDE_DIRECTORIES(${OCC_ROOT})
INCLUDE_DIRECTORIES(${LBFGS_ROOT}/include)
#INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
#INCLUDE_DIRECTORIES(${ARMA_ROOT}/include)
#INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${BLAS_LAPACK_INC_DIR})
INCLUDE_DIRECTORIES(${ARMADILLO_ROOT}/include)

#LINK_DIRECTORIES(${ARMA_ROOT}/lib)
LINK_DIRECTORIES(${OCC_ROOT}/lib)
LINK_DIRECTORIES(${LBFGS_ROOT}/lib)
#LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})
#LINK_DIRECTORIES(${BLAS_LIB_DIRS})
LINK_DIRECTORIES(${BLAS_LAPACK_LIB_DIR})
LINK_DIRECTORIES(${ARMADILLO_ROOT}/lib)

#SET(BLAS_LIBS mkl_intel_lp64 mkl_core mkl_sequential)
#SET(BLAS_LIB_DIRS /home/john/local/intel/mkl/lib/intel64 /home/john/local/intel/lib/intel64)

#message("BLAS_LIBS: !" ${BLAS_LIBS})




FILE(GLOB OCC_LIBS "${OCC_ROOT}/lib/*.a")

#SET(MeshOpt_LIBS armadillo ${Boost_LIBRARIES} Common Mesh FiniteElement Optimization ${LBFGS_ROOT}/lib/liblbfgs.a ${BLAS_LIBS})

#SET(MeshOpt_LIBS armadillo Common Mesh FiniteElement Optimization ${LBFGS_ROOT}/lib/liblbfgs.a ${BLAS_LIBS})

SET(MeshOpt_LIBS Common Mesh Geometry FiniteElement Optimization ${LBFGS_ROOT}/lib/liblbfgs.a ${BLAS_LAPACK_LIBS})

#ADD_SUBDIRECTORY(Common)
#ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/../Solvers 
#  ${CMAKE_CURRENT_BINARY_DIR}/Solvers)

ADD_SUBDIRECTORY(Optimization)
ADD_SUBDIRECTORY(FiniteElement)
ADD_SUBDIRECTORY(Geometry)
ADD_SUBDIRECTORY(Mesh)
ADD_SUBDIRECTORY(ParameterList)
#ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/../Trilinos 
#  ${CMAKE_CURRENT_BINARY_DIR}/Trilinos)

ADD_SUBDIRECTORY(Common)
ADD_SUBDIRECTORY(BL)
ADD_SUBDIRECTORY(MeshOptimization)
ADD_SUBDIRECTORY(HighOrder)


#ADD_SUBDIRECTORY(Debugging)

#move_src_up_level()

message("MeshOpt_LIBS: ${MeshOpt_LIBS}")
message("Mesh_INC: " ${Mesh_INC})
message("sources: " ${SOURCES})
message("Headers: " ${HEADERS})

add_executable(MeshOpt ${SOURCES} ${HEADERS})
#target_include_directories(MeshOpt PUBLIC ${INCDIRS} ${Mesh_INC} ${Solvers_INC} ${Optimization_INC} ${FiniteElement_INC})
target_include_directories(MeshOpt PUBLIC ${INCDIRS} ${Geometry_INC} ${Mesh_INC} ${Optimization_INC} ${FiniteElement_INC})

#target_link_libraries(MeshOpt ${MeshOpt_LIBS})

#target_link_libraries(MeshOpt -Wl,--start-group  ${MeshOpt_LIBS} ${OCC_LIBS} ${OPENMPI_LIBS} ${Trilinos_TPL_LIBRARIES} ${Trilinos_LIBRARIES} -Wl,--end-group -Wl,--whole-archive -Wl,--no-whole-archive pthread m dl)

target_link_libraries(MeshOpt -Wl,--start-group  ${MeshOpt_LIBS} ${OCC_LIBS} -Wl,--end-group -Wl,--whole-archive -Wl,--no-whole-archive pthread m dl gfortran)



install(TARGETS MeshOpt DESTINATION bin)


#install(FILES ${HEADERS} DESTINATION include)


