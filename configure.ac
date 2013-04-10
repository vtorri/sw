
### Version

m4_define([v_maj], [0])
m4_define([v_min], [0])
m4_define([v_mic], [1])

m4_define([sw_version], [v_maj.v_min.v_mic])

m4_define([lt_cur], m4_eval(v_maj + v_min))
m4_define([lt_rev], v_mic)
m4_define([lt_age], v_min)

AC_INIT([sw], [sw_version], [vincent.torri at univ-evry dot fr])
AC_PREREQ([2.68])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([1.11 dist-bzip2])
AM_SILENT_RULES([yes])

LT_PREREQ([2.4])
LT_INIT([win32-dll disable-static])

version_info="lt_cur:lt_rev:lt_age"
AC_SUBST([version_info])

VMAJ=v_maj
AC_SUBST([VMAJ])


### Default options with respect to host

AC_CANONICAL_HOST


### Additional options to configure


# Unit tests, coverage and benchmarking

SW_CHECK_TESTS([enable_tests="yes"], [enable_tests="no"])
SW_CHECK_COVERAGE([${enable_tests}], [enable_coverage="yes"], [enable_coverage="no"])

if test "x${prefer_assert}" = "xno"; then
   DEBUG_CFLAGS="-DNDEBUG"
else
# use debug symbols and set DEBUG if coverage support is enabled
   DEBUG_CFLAGS="${SW_DEBUG_CFLAGS}"
fi
AC_SUBST([DEBUG_CFLAGS])


### Checks for programs
AC_PROG_CC

# pkg-config
PKG_PROG_PKG_CONFIG

# doxygen program for documentation building

SW_CHECK_DOXYGEN([build_doc="yes"], [build_doc="no"])


### Checks for libraries


### Checks for header files


### Checks for types


### Checks for structures


### Checks for compiler characteristics

AM_PROG_CC_C_O

SW_CHECK_COMPILER_FLAGS([-Wall -Wextra -Wshadow -Wno-long-long -Wno-pointer-sign -Wno-empty-body -Wno-format-zero-length -Wno-nonnull -Wno-overflow -Wno-uninitialized -Wpointer-arith])

### Checks for linker characteristics


### Checks for library functions

AC_CONFIG_FILES([
Makefile
sw.pc
])

AC_OUTPUT

#####################################################################
## Info

echo
echo
echo
echo "------------------------------------------------------------------------"
echo "$PACKAGE_NAME $PACKAGE_VERSION"
echo "------------------------------------------------------------------------"
echo
echo "Configuration Options Summary:"
echo
echo "  Documentation........: ${enable_doc}"
echo "  Tests................: ${enable_tests}"
echo "  Coverage.............: ${enable_coverage}"
echo
echo "  Compilation..........: make"
echo
echo "  Installation.........: make install"
echo
echo "    prefix.............: $prefix"
echo
