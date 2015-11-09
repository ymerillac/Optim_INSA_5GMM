# CMake project file for Apame
cmake_minimum_required (VERSION 2.6)
project (Apame)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# set fortran compiler flags
set (CMAKE_Fortran_FLAGS_RELEASE "-fopenmp -O2")
set (CMAKE_Fortran_FLAGS_DEBUG   "-fopenmp -O0 -g")

# default installation
get_filename_component (default_prefix "." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the Apame directory."
      FORCE)

# build executable
FILE(GLOB APAME_F90 "src/*.f90")
add_executable(apame ${APAME_F90})

target_link_libraries(apame blas)
target_link_libraries(apame lapack)

# install executables and scripts
install (TARGETS apame 
         RUNTIME DESTINATION "bin")
