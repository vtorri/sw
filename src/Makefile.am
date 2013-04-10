MAINTAINERCLEANFILES = Makefile.in

SUBDIRS = . tests

lib_LTLIBRARIES = lib/libsw.la

install_sw_headersdir = $(pkgincludedir)-@VMAJ@

dist_install_sw_headers_DATA = \
lib/sw.h \
lib/sw_forward.h \
lib/sw_mra.h \
lib/sw_mra.x \
lib/sw_scale_fct.h \
lib/sw_scale_fct_base.x \
lib/sw_spline.h \
lib/sw_sweldens.h \
lib/sw_utils.h \
lib/sw_wavelet.h \
lib/sw_wavelet.x \
lib/sw_weights.h

install_sw_headers_algebradir = $(pkgincludedir)-@VMAJ@/algebra

dist_install_sw_headers_algebra_DATA = \
lib/algebra/sw_rational.h \
lib/algebra/sw_polynomial.h \
lib/algebra/sw_system.h

lib_libsw_la_SOURCES = \
lib/algebra/sw_rational.c \
lib/algebra/sw_polynomial.c \
lib/algebra/sw_system.c

lib_libsw_la_SOURCES += \
lib/sw_mra.c \
lib/sw_scale_fct.c \
lib/sw_spline.c \
lib/sw_sweldens.c \
lib/sw_utils.c \
lib/sw_wavelet.c \
lib/sw_weights2.c

lib_libsw_la_CPPFLAGS = \
-I$(top_srcdir)/src/lib \
-I$(top_srcdir)/src/lib/algebra \
-DSW_BUILD

lib_libsw_la_CFLAGS =

lib_libsw_la_LIBADD = @SW_COVERAGE_LIBS@ -lm

lib_libsw_la_LDFLAGS = -no-undefined -version-info @version_info@

### Binary

bin_PROGRAMS = bin/bench

bin_bench_SOURCES = bin/bench.c
bin_bench_CPPFLAGS = -I$(top_srcdir)/src/lib
bin_bench_CFLAGS =
bin_bench_LDADD = lib/libsw.la


clean-local:
	@rm -rf lib/*.gcno
