# GOMO
Generalized Operator Modelling of the Ocean (GOMO) is a three-dimensional ocean model based on OpenArray which is a simple operator library for the decoupling of ocean modelling and parallel computing.

The fundamental equations and algorithms of GOMO are derived from POM2k (Blumberg and Mellor, 1987). The details of the continuous governing equations, the corresponding operator expression form and the descriptions of all the variables used in GOMO are listed in the `./docs` folder.

GOMO is composed of 42 Fortran files (.F90), a header file (.h), and a single namelist file (.txt). Pre-processing processing package written in Matlab is located in `./pre`. You can use the package to produce the input file for the ideal test--seamount. The default input file `seamount65_49_21.nc` is located at the directory `./bin/data`. 

For more details, please refer to the paper (https://www.geosci-model-dev-discuss.net/gmd-2019-28/) and the simple user mannual of OpenArray located at ./docs.

## Requirments
* PnetCDF (version 1.7.0 or higher).  
* OpenArray.  

## Compile 
If the installation of OpenArray is done, it is fairly easy to compile GOMO, the basic steps are as follows:

  1) Download GOMO from GitHub:  

```shell
        wget https://github.com/hxmhuang/GOMO/archive/master.zip  
	unzip master.zip
```

  2) Set environment variables to specify path to pnetcdf   

```shell
        export PATH=/path/to/pnetcdf/bin:$PATH   
        export CPLUS_INCLUDE_PATH=/path/to/pnetcdf/include:$CPLUS_INCLUDE_PATH    
	export C_INCLUDE_PATH=/path/to/pnetcdf/include:$C_INCLUDE_PATH   
        export LD_LIBRARY_PATH=/path/to/pnetcdf/lib:$LD_LIBRARY_PATH    
        export LIBRARY_PATH=/path/to/pnetcdf/lib:$LIBRARY_PATH   
```

  3) Specify the path to OpenArray;  
```shell
        export OPEN_ARRAY=/path/to/OpenArray
```

  4) Check the compilers used to build PnetCDF.

```shell
	which pnetcdf_version
	pnetcdf_version
```

  5) Make. Please use the same compilers to build GOMO, PnetCDF and OpenArray. If PnetCDF and OpenArray are built by mpich, then:  

```shell
	./configure MPICC=mpicc MPICXX=mpicxx MPIF77=mpif77 MPIF90=mpif90
        make  
```

## Run GOMO
After compiling, the executable file `./bin/GOMO` will be generated. Within the directory `./bin` where `GOMO` and `config.txt` files exist, run the ideal case:
   
```shell
  ./GOMO   
```

or
  
```shell
  mpirun -np N ./GOMO   
```
