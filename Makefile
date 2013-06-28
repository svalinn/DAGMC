# Makefile for Sphinx documentation
#

GH_SOURCE_DIRS = doc-src
GH_BUILT_DIRS = _static _sources _downloads neutronics group fuelcycle build
GH_BUILT_FILES = objects.inv genindex.html index.html searchindex.js search.html 

GH_CURRENT_BRANCH = $(shell git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/\1/')
GH_SOURCE_BRANCH = master
GH_BUILD_BRANCH = gh-pages

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
	rm -rf $(GH_SOURCE_DIRS) build
	git checkout $(GH_SOURCE_BRANCH)

gh-push:
	rm -rf $(GH_SOURCE_DIRS) build
	git add -A 
	git commit -m "Generated $(GH_BUILD_BRANCH) for `git log $(GH_SOURCE_BRANCH) -1 --pretty=short --abbrev-commit`" && git push origin $(GH_BUILD_BRANCH)
	git checkout $(GH_SOURCE_BRANCH)


gh-install:
	rsync -a $(BUILDDIR)/html/* .
	rm -rf $(BUILDDIR)/html/*

clean:
	rm -rf $(GH_BUILT_DIRS) $(GH_BUILT_FILES)


html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)."

