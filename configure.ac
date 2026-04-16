#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([HydroConstruct], [1.0], [bec.gorton@csiro.au])
AC_CONFIG_SRCDIR([getflatinputdata.c])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC


# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [abs], [],[
           echo "Maths library is required for this program"
           exit -1])
# FIXME: Replace `main' with a function in `-lnetcdf':
AC_CHECK_LIB([netcdf], [ncopen], [],[
           echo "netCDF library is required for this program. Install netcdfg-dev."
           exit -1])
# FIXME: Replace `main' with a function in `-lsjwlib':
#AC_CHECK_LIB([sjwlib], [main])


# Checks for header files.
AC_CHECK_HEADERS([fcntl.h malloc.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor memset pow strrchr strstr])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
