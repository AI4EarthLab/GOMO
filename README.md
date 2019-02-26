# GOMO
Generalized Operator Modelling of the Ocean (GOMO) is a practical ocean model based on OpenArray which is a simple operator library for the decoupling of ocean modelling and parallel computing.

The fundamental equations of GOMO are derived from POM2k. GOMO features a bottom-following, free-surface, staggered Arakawa C grid. To effectively evolve the rapid surface fluctuations, GOMO uses the mode-splitting algorithm to address the fast propagating surface gravity waves and slow propagating internal waves in barotropic (external) and baroclinic (internal) modes, respectively. The details of the continuous governing equations, the corresponding operator expression form and the descriptions of all the variables used in GOMO are listed in the docs folder.

An ideal test case `Seamount` is included and `./bin/data/seamount65_49_21.nc` is the input file.

For more details, please see the paper (https://www.geosci-model-dev-discuss.net/gmd-2019-28/).

# Compile GOMO
Before attempting to compile GOMO, the following software are required:
  1) OpenArray ( It is available at https://github.com/hxmhuang/OpenArray_CXX).
  2) Fortran 90 or Fortran 95 compiler.
  3) GNU make version 3.81 or higher.
  4) NetCDF library.
  5) Message Passing Interface (MPI) library.

First, compile and install OpenArray. Second, change the directory to GOMO, set the right path of OpenArray in makefile and create a new folder named lib. Lastly, compile GOMO. When it finishes, change into ./bin directory. if the executable file GOMO is generated, then the compilation was successful.

# Run GOMO
Within the ./bin directory where GOMO and config.txt files exist, then type:
  `./GOMO` 
or
  `mpirun -np N ./GOMO` 
