
bin_PROGRAMS += src/bin/bench src/bin/tr_test src/bin/sw_advection src/bin/sw_advection_2d

src_bin_bench_SOURCES = src/bin/bench.c
src_bin_bench_CPPFLAGS = -I$(top_srcdir)/src/lib
src_bin_bench_CFLAGS =
src_bin_bench_LDADD = src/lib/libsw.la

src_bin_tr_test_CPPFLAGS = -I$(top_srcdir)/src/lib -I$(top_srcdir)/src/lib/tr
src_bin_tr_test_SOURCES = src/bin/tr_test.c
src_bin_tr_test_LDADD = src/lib/tr/libtransport.la src/lib/libsw.la -lm

src_bin_sw_advection_SOURCES = src/bin/sw_advection.c
src_bin_sw_advection_CPPFLAGS = -I$(top_srcdir)/src/lib
src_bin_sw_advection_CFLAGS =
src_bin_sw_advection_LDADD = src/lib/libsw.la -lm

src_bin_sw_advection_2d_SOURCES = src/bin/sw_advection_2d.c
src_bin_sw_advection_2d_CPPFLAGS = -I$(top_srcdir)/src/lib
src_bin_sw_advection_2d_CFLAGS =
src_bin_sw_advection_2d_LDADD = src/lib/libsw.la -lm
