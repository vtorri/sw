AUTOMAKE_OPTIONS = subdir-objects

ACLOCAL_AMFLAGS = -I m4

MAINTAINERCLEANFILES = \
Makefile.in \
aclocal.m4 \
compile \
config.guess \
config.h.in \
config.h.in~ \
config.sub  \
configure \
depcomp \
install-sh \
ltmain.sh  \
missing  \
$(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.gz \
$(PACKAGE_TARNAME)-$(PACKAGE_VERSION).tar.bz2 \
$(PACKAGE_TARNAME)-$(PACKAGE_VERSION)-doc.tar.bz2 \
m4/libtool.m4 \
m4/lt~obsolete.m4 \
m4/ltoptions.m4 \
m4/ltsugar.m4 \
m4/ltversion.m4

lib_LTLIBRARIES =
bin_PROGRAMS =
EXTRA_DIST =
SW_CLEANFILES =

include src/lib/Makefile.mk
include src/bin/Makefile.mk
include doc/Makefile.mk

EXTRA_DIST += \
AUTHORS \
COPYING \
README \
TODO \
autogen.sh \
m4/sw_compiler.m4 \
m4/sw_coverage.m4 \
m4/sw_doxygen.m4 \
m4/sw_tests.m4

sed_process = \
$(SED) \
-e 's,@VERSION\@,$(VERSION),g' \
-e 's,@VMAJ\@,$(VMAJ),g' \
-e 's,@prefix\@,$(prefix),g' \
-e 's,@exec_prefix\@,$(exec_prefix),g' \
-e 's,@libdir\@,$(libdir),g' \
-e 's,@includedir\@,$(includedir),g' \
-e 's,@pkgincludedir\@,$(pkgincludedir),g' \
< $< > $@ || rm $@

pc_verbose = $(pc_verbose_@AM_V@)
pc_verbose_ = $(pc_verbose_@AM_DEFAULT_V@)
pc_verbose_0 = @echo "  PC      " $@;

sw.pc: $(top_srcdir)/sw.pc.in Makefile
	$(AM_V_at)rm -f $@
	$(pc_verbose)$(sed_process)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = sw.pc

SW_CLEANFILES += sw.pc

.PHONY: doc coverage

# Unit tests

if SW_ENABLE_TESTS

check-local:
	@./src/tests/sw_suite

else

check-local:
	@echo "reconfigure with --enable-tests"

endif

# Coverage report

if SW_ENABLE_COVERAGE
lcov-reset:
	@rm -rf coverage
	@find . -name "*.gcda" -exec rm {} \;
	@lcov --directory . --zerocounters

lcov-report:
	@mkdir coverage
	@lcov --compat-libtool --directory . --capture --output-file coverage/coverage.info
	@lcov -l coverage/coverage.info | grep -v "`cd $(top_srcdir) && pwd`" | cut -d: -f1 > coverage/remove
	@lcov -r coverage/coverage.info `cat coverage/remove` > coverage/coverage.cleaned.info
	@rm coverage/remove
	@mv coverage/coverage.cleaned.info coverage/coverage.info
	@genhtml -t "$(PACKAGE_STRING)" -o coverage coverage/coverage.info

coverage:
	@make lcov-reset
	@make check
	@make lcov-report

SW_CLEANFILES += coverage

else
lcov-reset:
	@echo "reconfigure with --enable-gcov"

lcov-report:
	@echo "reconfigure with --enable-gcov"

coverage:
	@echo "reconfigure with --enable-tests --enable-gcov"
endif

clean-local:
	rm -rf $(SW_CLEANFILES)
