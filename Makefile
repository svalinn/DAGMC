# Makefile for Sphinx documentation
#

GH_PAGES_SOURCES = doc-src
SOURCE_BRANCH = master
BUILD_BRANCH = gh-pages

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = build

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) doc-src/source
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) source

.PHONY: help clean html gh-revert gh-push

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  gh-revert  to revert changes made to gh-pages branch and switch back to master"
	@echo "  gh-push    to push HTML documentation in gh-pages branch"
	@echo "  html       to make standalone HTML files"


gh-revert:
	git checkout -f --
	rm -rf $(GH_PAGES_SOURCES) build
	git checkout master

gh-push:
	rm -rf $(GH_PAGES_SOURCES) build
	git add -A
	git commit -m "Generated $(BUILD_BRANCH) for `git log master -1 --pretty=short --abbrev-commit`" && git push origin $(BUILD_BRANCH)
	git checkout $(SOURCE_BRANCH)



clean:
	-rm -rf _images _sources people projects papers 
	-rm -rf index.html

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)."

install:
	rsync -a $(BUILDDIR)/* .
	rm -rf $(BUILDDIR)/*



