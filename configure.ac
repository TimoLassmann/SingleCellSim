AC_INIT([singleCellsim], [m4_esyscmd_s([git describe --always])])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE
# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL

# Checks for libraries.

# FIXME: Replace `main' with a function in `-lgslcblas':
AC_CHECK_LIB([gslcblas], [cblas_sdot])
# FIXME: Replace `main' with a function in `-lgsl':


AC_CHECK_LIB([gsl],[gsl_integration_qag],[],[AC_MSG_WARN([gsl library missing - some programs will not be compiled!])])

AM_CONDITIONAL([GOT_GSL], [ test "$ac_cv_lib_gsl_gsl_integration_qag" = "yes"])

# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [log])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [pthread_create])


AC_CHECK_HEADERS([float.h stdint.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_INT32_T
AC_TYPE_SIZE_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_LSTAT_FOLLOWS_SLASHED_SYMLINK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([getcwd mkdir pow sqrt])


AC_CONFIG_FILES([Makefile src/Makefile])




tldevel_check_f=$ac_abs_confdir/tldevel/README

AC_CHECK_FILE([$tldevel_check_f], [],[have_tldevel=no])


if test "x${have_tldevel}" = xno; then
   echo "------------------------------------------"
   echo " Unable to find pthreads on this system.  "
   echo " Building a single-threaded version.      "
   echo "------------------------------------------"
   AC_MSG_ERROR([ASFSAF])
fi




LIB_TLDEVELDIR="tldevel"
AC_SUBST(LIB_TLDEVELDIR)



AC_ARG_ENABLE(debugging,[AS_HELP_STRING([--enable-debugging],[include debugging code])
AS_HELP_STRING([--enable-debugging=x],[also set diagnostics verbosity level to <x> (1-3)])],
enable_debugging=$enableval, enable_debugging="no")

case $enable_debugging in
yes)  AC_DEFINE(DEBUGLEVEL, 0,[No debugging. ]);;
1)  AC_DEFINE(DEBUGLEVEL, 1,[Defines debugging level 1.]);;
2)  AC_DEFINE(DEBUGLEVEL, 2,[Defines debugging level 2.]);;
3)  AC_DEFINE(DEBUGLEVEL, 3,[Defines debugging level 3.]);;
no)  AC_DEFINE(DEBUGLEVEL, 0,[No debugging.]);;
*)  AC_MSG_ERROR([Unknown argument to --enable-debugging: $enable_debugging]);;
esac

if test "$enable_debugging" != "no"; then
        AC_DEFINE(DEBUG,1,[Defines debugging .])
        CFLAGS="-ggdb -Wall -m64 -std=c99 $MEM $SIGNAL32"
else
        CFLAGS="-O2 -funroll-loops -Wall -std=c99 $MEM $SIGNAL32"
        DEBUG=0
fi



AC_ARG_ENABLE(valgrind-tests,
    [--enable-valgrind-tests run self tests under valgrind],
    opt_valgrind_tests="$enableval", opt_valgrind_tests="no")
    

OPTS="-q --trace-children=no --error-exitcode=1 --leak-check=yes"
  # Run self-tests under valgrind?
if test "$opt_valgrind_tests" = "yes"; then
        AC_DEFINE(DEBUG,1,[Defines debugging .])   
        CFLAGS="-pg -ggdb -pedantic -Wall -m64 -std=gnu11"
	AC_CHECK_PROGS(VALGRIND, valgrind)
	VALGRIND="$VALGRIND $OPTS"
else
	VALGRIND=
fi
AC_MSG_CHECKING([whether self tests are run under valgrind])
AC_MSG_RESULT($opt_valgrind_tests)



ac_configure_args="$ac_configure_args PACKAGE_VERSION=$PACKAGE_VERSION   PACKAGE_NAME=$PACKAGE_NAME"



AC_CONFIG_SUBDIRS([tldevel])

AC_OUTPUT



if test "$ac_cv_lib_gsl_gsl_integration_qag" = "no";then
echo "     WARNING: GSL library was not found on your system."
echo "     Programs: singleCellSim will not be compiled."
fi
