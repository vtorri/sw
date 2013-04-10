MAINTAINERCLEANFILES = Makefile.in

AM_CPPFLAGS = \
-I$(top_srcdir)/src/lib \
-I$(top_srcdir)/src/lib/algebra \
@CHECK_CFLAGS@

if SW_ENABLE_TESTS

check_PROGRAMS = sw_suite

sw_suite_SOURCES = sw_suite.c
sw_suite_LDADD = @CHECK_LIBS@ $(top_builddir)/src/lib/libsw.la -lm

endif
