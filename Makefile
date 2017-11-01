# Makefile for Sphinx documentation

# Top-level directory of all the ReST source files
GH_SOURCE_DIR     = doc
# Repository branch that contains the source
GH_SOURCE_BRANCH  = develop
# Repository branch that contains the rendered HTML
GH_PUBLISH_BRANCH = gh-pages
# Repository that contains the rendered HTML
GH_UPSTREAM_REPO  = https://github.com/svalinn/DAGMC
# Directory that contains the rendered HTML
BUILDDIR          = gh-build
# Sphinx executable
SPHINXBUILD       = sphinx-build
# Sphinx options
SPHINXOPTS        =

ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) $(GH_SOURCE_DIR)

help:
	@echo "Please use \"make <target\" where <target> is one of"
	@echo "  html     to build HTML in directory \"$(BUILDDIR)\" for testing"
	@echo "  clean    to cleanup HTML build in directory \"$(BUILDDIR)\" after testing"
	@echo "  publish  final build and push from \"$(GH_SOURCE_BRANCH)\" branch to \"$(GH_PUBLISH_BRANCH)\" branch"

clean:
	-rm -rf $(BUILDDIR)

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)."

publish:
	git checkout $(GH_PUBLISH_BRANCH) && \
	git rm -rf * && \
	git checkout $(GH_SOURCE_BRANCH) -- $(GH_SOURCE_DIR) Makefile && \
	git reset HEAD && \
	make clean && \
	make html && \
	rsync -a $(BUILDDIR)/* . && \
	rsync -a $(BUILDDIR)/.* . && \
	rm -rf $(GH_SOURCE_DIR) $(BUILDDIR) && \
	git add . && \
	git commit -m "Generated $(GH_PUBLISH_BRANCH) for `git log $(GH_SOURCE_BRANCH) -1 --pretty=short --abbrev-commit`" && \
	git push $(GH_UPSTREAM_REPO) $(GH_PUBLISH_BRANCH) && \
	git checkout $(GH_SOURCE_BRANCH)
