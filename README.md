# GOMO
Generalized Operator Modelling of the Ocean (GOMO) is a three-dimensional ocean model based on OpenArray which is a simple operator library for the decoupling of ocean modelling and parallel computing.

The fundamental equations and algorithms of GOMO are derived from POM2k (Blumberg and Mellor, 1987). GOMO features bottom-following, free-surface, staggered Arakawa C-grid. To effectively evolve the rapid surface fluctuations, GOMO uses the mode-splitting algorithm to address the fast propagating surface gravity waves and slow propagating internal waves in barotropic (external) and baroclinic (internal) modes, respectively. The details of the continuous governing equations, the corresponding operator expression form and the descriptions of all the variables used in GOMO are listed in the ./docs folder.

GOMO is composed of 42 Fortran files (.F90), a header file (.h), a single namelist file (.txt), and a makefile. Pre-processing processing package written in Matlab is located in ./pre. You can use the package to produce the input file for the ideal test--seamount. The default input file seamount65_49_21.nc is located at the directory ./bin/data. 

For more details, please see the paper (https://www.geosci-model-dev-discuss.net/gmd-2019-28/) and the simple user mannual of OpenArray located at ./docs.

# Compile GOMO
If the installation of OpenArray is done, it is fairly easy to compile GOMO, the basic steps are as follows:

  1) Download GOMO from GitHub:  
        git clone https://github.com/hxmhuang/GOMO.git   
     or  
        wget https://github.com/hxmhuang/GOMO/archive/master.zip  
  2) Set environment variables to specify path to pnetcdf 
        export PATH=/path/to/pnetcdf/bin:$PATH 
        export CPLUS_INCLUDE_PATH=/path/to/pnetcdf/include:$CPLUS_INCLUDE_PATH
	export C_INCLUDE_PATH=/path/to/pnetcdf/include:$C_INCLUDE_PATH
        export LD_LIBRARY_PATH=/path/to/pnetcdf/lib:$LD_LIBRARY_PATH 
        export LIBRARY_PATH=/path/to/pnetcdf/lib:$LIBRARY_PATH
  3) Edit the Makefile, modify the installation directory of OpenArray;  
        export OPEN_ARRAY=/path/to/OpenArray
  4) Make. Please use the same compiler to complile GOMO and OpenArray.  
        For intel compiler: make -f Makefile.intel;  
        For mpich compiler: make -f Makefile.mpich;  

# Run GOMO
After compiling, the executable file ./bin/GOMO will be generated. Within the directory ./bin where GOMO and config.txt files exist, type:   
  ./GOMO   
or  
  mpirun -np N ./GOMO   
