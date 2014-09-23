# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(mcmc_affine, 0.2.8, alexander.conley@colorado.edu)
AC_LANG_CPLUSPLUS
dnl AC_CONFIG_SRCDIR([src])
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIR([m4])

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_PROG_INSTALL

dnl Initialize CXXFLAGS and CPPFLAGS
dnl CPPFLAGS are preprocessor flags, CXX are compiler flags
dnl Unless --enable-debug is set, then -O3 is used
AC_ARG_WITH(cppflags,
            [  --with-cppflags=VAL     CPPFLAGS value (def: )],
            [CPPFLAGS="$withval"],[CPPFLAGS=""] )
AC_ARG_WITH(cxxflags,[AC_HELP_STRING([--with-cxxflags=VAL],
	[CXXFLAGS value (def: -Wall)])], [ CXXFLAGS=$withval $CXXFLAGS ], 
	[CXXFLAGS="-Wall $CXXFLAGS"])
AC_ARG_ENABLE(debug,
	      AC_HELP_STRING([--enable-debug],
	                     [Add -g to CXXFLAGS]),
              [ CXXFLAGS="-g -O0 $CXXFLAGS"], [ CXXFLAGS="-O3 $CXXFLAGS"])
AC_ARG_ENABLE(profile,
              [  --enable-profile        Add gprof-friendly stuff to CXXFLAGS],
              [ CXXFLAGS="$CXXFLAGS -pg -fprofile-arcs" ],)
AC_ARG_ENABLE(analyze,
              [  --enable-analyze        Run static checker (if using clang)],
              [ CXXFLAGS="$CXXFLAGS --analyze" ],)

dnl Example code
AM_CONDITIONAL(DOEXAMPLES, false)
AC_ARG_ENABLE(examples,    
	[  --enable-examples       Enable examples],
	[ AC_MSG_RESULT(Enabling examples)
          AM_CONDITIONAL(DOEXAMPLES, true) ],)

dnl Tests
AC_LANG_CPLUSPLUS
AM_CONDITIONAL(DOTEST, false)
AC_ARG_ENABLE(test,
	[  --enable-test           Enable testing],
	usetest=$enableval, usetest="")
if test $usetest; then
   AC_MSG_RESULT(Enabling testing, unit and other)
   AM_CONDITIONAL(DOTEST, true)

   dnl Find gtest library
   AC_ARG_WITH(gtestlib,
	       [  --with-gtestlib=DIR       DIR=-LDIR, has libgtest],
               gtestlib="$withval/",gtestlib= )
   if test $gtestlib; then
      AC_MSG_RESULT( Checking User Specified gtest )
      LDFLAGS="$LDFLAGS -L$gtestlib"
      AC_CHECK_LIB( gtest, main, ,
      		    AC_MSG_ERROR("libgtest not in location specified by --with-gtestlib: $gtestlib"),
      		    [ -lgtest -lgtest_main -L$gtestlib ] )
   else
      AC_MSG_RESULT( Checking gtest_main )
      AC_CHECK_LIB( gtest_main, main, ,
      AC_MSG_ERROR("libgtest_main not found -- use --with-gtestlib"),
                   [ -lgtest -lgtest_main ] )
      AC_MSG_RESULT( Checking gtest )
      AC_CHECK_LIB( gtest, main, ,
      AC_MSG_ERROR("libgtest not found -- use --with-gtestlib"),
                   [ -lgtest -lgtest_main ] )
   fi

   dnl Find gtest includes
   AC_ARG_WITH(gtestinc,
        [  --with-gtestinc=VAL      VAL=-IVAL, has gtest include files ],
        [gtestinc="$withval/" CPPFLAGS="$CPPFLAGS -I$gtestinc"],
        [gtestinc=] )
   AC_MSG_RESULT( Checking gtest Headers )
   AC_CHECK_HEADERS([gtest/gtest.h],,
     [ if test $gtestinc; then
         AC_MSG_ERROR("gtest headers not in location specified by --with-gtestinc")
       elif test $gtestlib; then
         AC_CHECK_HEADER($gtestlib/../include/test.h,
            [gtestinc="$gtestlib/../include/" CPPFLAGS="$CPPFLAGS -I$gtestinc"],
             AC_MSG_ERROR("Can't find gtest includes in user specified location; use --with-gtestinc"))
       else
          AC_MSG_ERROR("Can't find gtest includes; use --with-gtestinc")
       fi ] )
fi

dnl HDF5
AC_LANG_C
AC_ARG_WITH(hdfinc,
        [AC_HELP_STRING([--with-hdfinc=VAL],
                        [VAL=-IVAL, has HDF5 include files ]) ], 
        [hdfinc="$withval/" CPPFLAGS="$CPPFLAGS -I$hdfinc"],
        [hdfinc=] )
AC_ARG_WITH(hdflib,
	[AC_HELP_STRING([--with-hdflib=DIR],
			[DIR=-LDIR, has libhdf5])],
	hdflib="$withval/", hdflib= )

dnl HDF5 Header
AC_MSG_RESULT( Checking HDF5 Headers )
AC_CHECK_HEADERS([hdf5.h],,
  [ if test $hdfinc; then
   AC_MSG_ERROR("hdf5 headers not in location specified by --with-hdfinc")
   elif test $hdflib; then
    AC_CHECK_HEADER($hdflib/../include/hdf5.h,
     [hdfinc="$hdflib/../include/" CPPFLAGS="$CPPFLAGS -I$hdfinc"],
      AC_MSG_ERROR("Can't find HDF5 includes in user specified location; use --with-hdfinc"))
    else
      AC_MSG_ERROR("Can't find HDF5 includes; use --with-hdf5inc")
   fi ] )
dnl HDF5 lib
if test $hdflib; then
  AC_MSG_RESULT( Checking User specified HDF5 )
  LDFLAGS="$LDFLAGS -L$hdflib"
  LIBS="$LIBS -lhdf5"
  AC_CHECK_LIB( hdf5, main, ,
      AC_MSG_ERROR("libhdf5 not in location specified by --with-hdflib: $hdflib"),
      [ -lhdf5 -L$hdflib ] )

else
  AC_MSG_RESULT( Checking HDF5 )
  AC_CHECK_LIB( hdf5, main, ,
      AC_MSG_ERROR("libhdf5 not found; use --with-hdflib"), )
fi

# Checks for standard header files.
AC_HEADER_STDC

AC_CONFIG_FILES([Makefile src/Makefile test/Makefile examples/Makefile])
AC_OUTPUT
