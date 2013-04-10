
if SW_ENABLE_TESTS

check_PROGRAMS = src/tests/sw_suite

src_tests_sw_suite_SOURCES = src/tests/sw_suite.c
src_tests_sw_suite_CPPFLAGS = \
-I$(top_srcdir)/src/lib \
-I$(top_srcdir)/src/lib/algebra \
@CHECK_CFLAGS@
src_tests_sw_suite_LDADD = @CHECK_LIBS@ src/lib/libsw.la -lm

check-local:
	@./src/tests/sw_suite

else

check-local:
	@echo "reconfigure with --enable-tests"

endif
