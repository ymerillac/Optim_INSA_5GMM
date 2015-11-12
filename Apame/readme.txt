Install Apame software, developed by Daniel Filkovic (http://www.3dpanelmethod.com/)

To install:
- by cmake:
   mkdir build
   cd build
   cmake ..
   make
   make install
   
- by "hand":
   gfortran *.f90 -o ../bin/apame -llapack -lblas -fopenmp
