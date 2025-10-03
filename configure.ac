#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])

# Remember to update the version number for a new release
AC_INIT([TCRAD], [0.3], [r.j.hogan@ecmwf.int],[tcrad])

# Most code is C++ but there is some C and some Fortran
AC_LANG([C++])

AC_CONFIG_SRCDIR([src/tcrad_overlap.cpp])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Wno-error foreign])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs
#AC_PROG_AWK
AC_PROG_CC
AC_PROG_CXX
#AC_PROG_F77
#AC_PROG_FC
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_LIBTOOL

# Check for OpenMP
AC_OPENMP
AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")
AC_SUBST(AM_CFLAGS,"$OPENMP_CXXFLAGS")

AX_CHECK_ADEPT([have_adept=yes], [have_adept=no])

# Header directories needed during the compilation
#AC_SUBST(AM_CPPFLAGS,"-I../include")

# Check for NetCDF library: REQUIRED
#AC_SEARCH_LIBS([nc_create], [netcdf], , AC_MSG_ERROR([Unable to find NetCDF library]))

# Checks for standard header files
#AC_CHECK_HEADERS([fenv.h malloc.h memory.h stdint.h stdlib.h string.h strings.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics
AC_TYPE_SIZE_T

# Due to autotools stupidity, this breaks the build:
#AC_FUNC_MALLOC
#AC_CHECK_FUNCS([floor gethostname gettimeofday memset pow sqrt strcasecmp strchr strdup strncasecmp strstr strtol])

# Makefiles to create
AC_CONFIG_FILES([Makefile
                 src/Makefile
		 include/Makefile])

# Provide a summary to the user
AC_MSG_NOTICE([********************* Summary **************************************])
AC_MSG_NOTICE([  CXX = $CXX ])
AC_MSG_NOTICE([  CPPFLAGS = $CPPFLAGS])
AC_MSG_NOTICE([  CXXFLAGS = $CXXFLAGS])
AC_MSG_NOTICE([  LDFLAGS =  $LDFLAGS])
AC_MSG_NOTICE([  LIBS = $LIBS])
AC_MSG_NOTICE([Typing "make; make install" will install the headers and library files in])
AC_MSG_NOTICE([  $prefix])
AC_MSG_NOTICE([********************************************************************])

AC_OUTPUT
