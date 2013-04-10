
bin_PROGRAMS += src/bin/bench

src_bin_bench_SOURCES = src/bin/bench.c
src_bin_bench_CPPFLAGS = -I$(top_srcdir)/src/lib
src_bin_bench_CFLAGS =
src_bin_bench_LDADD = src/lib/libsw.la
