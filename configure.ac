AC_PREREQ([2.69])

m4_define([platonyzer_version_major],[2])
m4_define([platonyzer_version_minor],[0])
m4_define([platonyzer_version_micro],[0])
m4_define([platonyzer_version_extra],[])
m4_define([platonyzer_version],[platonyzer_version_major().platonyzer_version_minor().platonyzer_version_micro()platonyzer_version_extra])

AC_INIT([platonyzer], [platonyzer_version], [m.hekkelman@nki.nl])

dnl Switch to a decent C++ compiler, and check if it works.
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_17([noext])

AX_CHECK_COMPILE_FLAG([-fstandalone-debug],
	[
		CXXFLAGS="$CXXFLAGS -fstandalone-debug"
	] , , [-Werror])

AC_CONFIG_SRCDIR([src/platonyzer.cpp])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([config/m4])
AC_CONFIG_HEADERS([src/config.hpp])

AC_PREFIX_DEFAULT(/usr/local)

AC_DEFUN([read_test], [AC_LANG_SOURCE(
	esyscmd(config/tools/m4esc.sh config/test/$1))])

AC_PROG_INSTALL
PKG_PROG_PKG_CONFIG

AX_PTHREAD

AC_ARG_VAR([DEBUG], [Build a debug version of the application])

AX_MRC

dnl revision numbering is something used internally at the NKI
AC_ARG_ENABLE(
	revision,
	[AS_HELP_STRING([--disable-revision], [Create a build number as revision])])

AS_IF([test "x$enable_revision" != "xno" ], [
	UPDATE_REVISION=1
])

AC_SUBST([UPDATE_REVISION], [$UPDATE_REVISION])

AX_FILESYSTEM

AX_BOOST_BASE([1.65.1], [], [AC_MSG_ERROR([Could not find a recent version of boost])])
AX_BOOST_IOSTREAMS
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_DATE_TIME
AX_BOOST_REGEX

AX_LIBCIFPP
AX_LIBZEEP

AX_LIBCLIPPER

AX_LIBPDB_REDO

dnl These are still needed outside the Debian environment
AX_CHECK_LIBRARY([LIBZ], [zlib.h], [z],
	[ LIBS="$LIBS -lz" ],
	[AC_MSG_ERROR([libz not found - compressed files not supported])])
AX_CHECK_LIBRARY([LIBBZ2], [bzlib.h], [bz2],
	[ LIBS="$LIBS -lbz2"],
	[AC_MSG_ERROR([libbz2 not found - compressed files not supported])])

AC_SUBST([LIBS], [$LIBS])

dnl Process Makefile.in to create Makefile
AC_OUTPUT([GNUmakefile])
