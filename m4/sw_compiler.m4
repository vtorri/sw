dnl Copyright (C) 2013 Vincent Torri <vincent dot torri at gmail dot com>
dnl This code is public domain and can be freely used or copied.

dnl Macro that check if compiler of linker flags are available


dnl Macro that checks for a compiler flag availability
dnl
dnl SW_CHECK_COMPILER_FLAG(flag)
dnl AC_SUBST : SW_CFLAGS

AC_DEFUN([SW_CHECK_COMPILER_FLAG],
[
m4_pushdef([UP], m4_translit([[$1]], [-a-z], [_A-Z]))

dnl store in options -Wfoo if -Wno-foo is passed
option=m4_bpatsubst([[$1]], [-Wno-], [-W])

CFLAGS_save="${CFLAGS}"
CFLAGS="${CFLAGS} ${option}"

AC_MSG_CHECKING([whether the compiler supports $1])

AC_COMPILE_IFELSE(
   [AC_LANG_PROGRAM([[]])],
   [have_flag="yes"],
   [have_flag="no"])
AC_MSG_RESULT([${have_flag}])

CFLAGS="${CFLAGS_save}"

if test "x${have_flag}" = "xyes" ; then
   SW_CFLAGS="${SW_CFLAGS} [$1]"
fi

AC_ARG_VAR([SW_CFLAGS], [preprocessor flags for $1])
AC_SUBST([SW_CFLAGS])

m4_popdef([UP])
])

dnl Macro that iterates over a sequence of white separated flags
dnl and that call SW_CHECK_COMPILER_FLAG() for each of these flags
dnl
dnl SW_CHECK_COMPILER_FLAGS(flags)

AC_DEFUN([SW_CHECK_COMPILER_FLAGS],
[
m4_foreach_w([flag], [$1], [SW_CHECK_COMPILER_FLAG(m4_defn([flag]))])
])


dnl Macro that checks for a linker flag availability
dnl
dnl SW_CHECK_LINKER_FLAG(flag)
dnl AC_SUBST : SW_LDFLAGS

AC_DEFUN([SW_CHECK_LINKER_FLAG],
[
m4_pushdef([UP], m4_translit([[$1]], [,-a-z], [__A-Z]))

LDFLAGS_save="${LDFLAGS}"
LDFLAGS="${LDFLAGS} $1"

AC_LANG_PUSH([C])
AC_MSG_CHECKING([whether the linker supports $1])

AC_LINK_IFELSE(
   [AC_LANG_PROGRAM([[]])],
   [have_flag="yes"],
   [have_flag="no"])
AC_MSG_RESULT([${have_flag}])

LDFLAGS="${LDFLAGS_save}"
AC_LANG_POP([C])

if test "x${have_flag}" = "xyes" ; then
   SW_LDFLAGS="${SW_LDFLAGS} [$1]"
fi

AC_SUBST([SW_LDFLAGS])

m4_popdef([UP])
])

dnl Macro that iterates over a sequence of white separated flags
dnl and that call SW_CHECK_LINKER_FLAG() for each of these flags
dnl
dnl SW_CHECK_LINKER_FLAGS(flags)

AC_DEFUN([SW_CHECK_LINKER_FLAGS],
[
m4_foreach_w([flag], [$1], [SW_CHECK_LINKER_FLAG(m4_defn([flag]))])
])
