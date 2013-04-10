
lib_LTLIBRARIES += src/lib/libsw.la

install_sw_headersdir = $(pkgincludedir)-@VMAJ@

dist_install_sw_headers_DATA = \
src/lib/sw.h \
src/lib/sw_forward.h \
src/lib/sw_mra.h \
src/lib/sw_mra.x \
src/lib/sw_scale_fct.h \
src/lib/sw_scale_fct_base.x \
src/lib/sw_spline.h \
src/lib/sw_sweldens.h \
src/lib/sw_utils.h \
src/lib/sw_wavelet.h \
src/lib/sw_wavelet.x \
src/lib/sw_weights.h

install_sw_headers_algebradir = $(pkgincludedir)-@VMAJ@/algebra

dist_install_sw_headers_algebra_DATA = \
src/lib/algebra/sw_rational.h \
src/lib/algebra/sw_polynomial.h \
src/lib/algebra/sw_system.h

src_lib_libsw_la_SOURCES = \
src/lib/algebra/sw_rational.c \
src/lib/algebra/sw_polynomial.c \
src/lib/algebra/sw_system.c

src_lib_libsw_la_SOURCES += \
src/lib/sw_mra.c \
src/lib/sw_scale_fct.c \
src/lib/sw_spline.c \
src/lib/sw_sweldens.c \
src/lib/sw_utils.c \
src/lib/sw_wavelet.c \
src/lib/sw_weights2.c

src_lib_libsw_la_CPPFLAGS = \
-I$(top_srcdir)/src/lib \
-I$(top_srcdir)/src/lib/algebra \
-DSW_BUILD

src_lib_libsw_la_CFLAGS =
src_lib_libsw_la_LIBADD = @SW_COVERAGE_LIBS@ -lm
src_lib_libsw_la_LDFLAGS = -no-undefined -version-info @version_info@
