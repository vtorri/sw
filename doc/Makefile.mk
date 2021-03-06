

.PHONY: doc

PACKAGE_DOCNAME = $(PACKAGE_TARNAME)-$(PACKAGE_VERSION)-doc

if SW_BUILD_DOC

SW_DOC_CLEANFILES = doc/html/ doc/latex/ $(PACKAGE_DOCNAME).tar*
SW_CLEANFILES += $(SW_DOC_CLEANFILES)

doc-clean:
	rm -rf $(SW_DOC_CLEANFILES)

doc: all doc-clean
	$(DOXYGEN)
	rm -rf $(PACKAGE_DOCNAME).tar*
	mkdir -p $(PACKAGE_DOCNAME)/doc
	cp -rf doc/html/ doc/latex/ $(PACKAGE_DOCNAME)/doc
	tar cf $(PACKAGE_DOCNAME).tar $(PACKAGE_DOCNAME)/
	bzip2 -9 $(PACKAGE_DOCNAME).tar
	rm -rf $(PACKAGE_DOCNAME)/
	mv $(PACKAGE_DOCNAME).tar.bz2 $(top_srcdir)

else

doc: all
	@echo "Documentation not built. Run ./configure --help"

endif

EXTRA_DIST += Doxyfile
