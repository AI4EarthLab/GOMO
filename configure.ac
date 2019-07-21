AC_PREREQ(2.50)
AC_INIT([GOMO], [0.0])
config_flags="$*"
AM_INIT_AUTOMAKE([subdir-objects no-define no-dependencies foreign])

AS_UNSET([CC])
AS_UNSET([CXX])
AS_UNSET([FC])
AC_PROG_CC([mpicc mpiicc])
AC_PROG_CXX([mpicxx mpiicpc])
AC_PROG_FC([mpif90 mpiifort])

# Append installition path to LDFLAGS and CFLAGS
AC_DEFUN([APPEND_SEARCH_PATH],
        [
            LDFLAGS+=" -L$1/lib "
            CFLAGS+=" -I$1/include "
            CPPFLAGS+=" -I$1/include "
            FCFLAGS+=" -L$1/lib"
            AC_SUBST(FCFLAGS)
            AC_SUBST(LDFLAGS)
            AC_SUBST(CFLAGS)
            AC_SUBST(CPPFLAGS)
        ])

# SEARCH_LIB(DIR, HEADER_FILE, LIB, EXAMPLE_FUNC_IN_LIB, DIR_NOTE, ERROR_NOTE_LIB_NOT_FOUND, ERROR_NOTE_HEADER_NOT_FOUND)
AC_DEFUN([SEARCH_LIB], [
        dnl Check for library
        AC_ARG_VAR($1, [$5])
        old_CFLAGS=${CFLAGS}
        old_LDFLAGS=${LDFLAGS}
        CFLAGS=" -I${$1}/include ${CFLAGS}"
        LDFLAGS=" -L${$1}/lib ${LDFLAGS}"
        echo "trying CFLAGS ${CFLAGS}"
        AC_CHECK_HEADER($2, [AC_CHECK_FUNC([$4], [],
                                            dnl check for function in lib
                                            [AC_SEARCH_LIBS([$4], [$3], [], [AC_MSG_ERROR([$6])])]
                                            )],
                        [AC_MSG_ERROR([$7])])
        CFLAGS=${old_CFLAGS}
        LDFLAGS=${old_LDFLAGS}
        AS_IF([test "x${$1}" != x], [APPEND_SEARCH_PATH([${$1}])], [])
])

SEARCH_LIB(MPI_DIR, mpi.h, [mpi mpich], [MPI_Comm_rank],
            [Specify the mpi installition root directory which contains the include and lib subdirectory, for example /usr/lib/x86_64-linux-gnu/openmpi],
            [
                            -----------------------------------------------------------------------
                                Invalid MPI compiler specified or detected
                                A working MPI C compiler is required. Please specify the location
                                of in the MPI_DIR environment variable. Abort.
                            -----------------------------------------------------------------------],
            [
                            -----------------------------------------------------------------------
                                Invalid MPI compiler specified or detected
                                A working pdfcdf lib is required. Please specify the location 
                                in the PNETCDF_DIR environment variable. Abort.
                            -----------------------------------------------------------------------])

SEARCH_LIB(PNETCDF_DIR, pnetcdf.h, [pnetcdf], [ncmpi_open],
            [Specify the pnetCdf lib installition root directory which contains the include and lib subdirectory, for example /path/to/pnetcdf_dir/],
            [
                            -----------------------------------------------------------------------
                                Invalid pnetcdf lib specified or detected
                                A working pdfcdf lib is required. Please specify the location 
                                in the PNETCDF_DIR environment variable. Abort.
                            -----------------------------------------------------------------------],
            [
                            -----------------------------------------------------------------------
                                Can not find pnetcdf.h in any existing include path,
                                please consider specify the pnetcdf installition folder using PNETCDF_DIR,
                                then call the configure again, for exmaple:
                                $ ./configure PNETCDF_DIR=/path/to/pnetcdf_dir/
                            -----------------------------------------------------------------------])

# check for openarray lib
AC_ARG_VAR(OPEN_ARRAY, [openarray installition directory])
LDFLAGS=" -L$OPEN_ARRAY/lib ${LDFLAGS}"
# echo "$(env)"
AC_SEARCH_LIBS([c_init], [openarray], [], [AC_MSG_ERROR([
                                -----------------------------------------------------------------------
                                Invalid openarray lib specified or detected
                                A working openarray lib is required. Please specify the location 
                                in the OPEN_ARRAY environment variable. Abort.
                            -----------------------------------------------------------------------])],
                            [-lstdc++ -lgfortran -lpnetcdf -lm -ldl -lpthread])
AS_IF([test "x$OPEN_ARRAY" != x], [APPEND_SEARCH_PATH([$OPEN_ARRAY])], [])

AC_SUBST(MPI_DIR)
AC_SUBST(PNETCDF_DIR)
AC_SUBST(OPEN_ARRY)

AC_PROG_RANLIB
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
echo "------------------------------------------------------------------------------"
echo \
"
   ${PACKAGE_NAME} Version ${PACKAGE_VERSION}"

echo "\

   Compilers: MPICC    = ${CC}
              MPICXX   = ${CXX}
              FC       = ${FC}
              CPPFLAGS = ${CPPFLAGS}
              CFLAGS   = ${CFLAGS}
              FCFLAGS  = ${FCFLAGS}
              LDFLAGS  = ${LDFLAGS}"

echo "\

   Now run 'make' to build the library and utility tools.
   Then run 'make @<:@<target>@:>@' for testing and installation, where the
   optional <target> can be:
              GOMO - build a test programs (build only, no run)
              install     - install OpenArray library in ${prefix}

------------------------------------------------------------------------------"
