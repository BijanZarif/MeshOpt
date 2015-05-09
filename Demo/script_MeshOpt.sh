SOURCE_BASE=$HOME/Dropbox/MeshOpt
~/local/CMake/install/bin/cmake \
  -D CMAKE_INSTALL_PREFIX:PATH=$PWD/../install \
  -D CMAKE_BUILD_TYPE=Debug \
  -D OCC_ROOT=/home/john/local/OpenCascade/OpenCASCADE_install \
  -D BLAS_LAPACK_LIBS="mkl_intel_lp64.a;mkl_core.a;mkl_sequential.a" \
  -D BLAS_LAPACK_LIB_DIR=/home/john/local/IntelParallelStudio/install/mkl/lib/intel64 \
  ${SOURCE_BASE}
