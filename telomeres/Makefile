.PHONY: clean clean-test clean-pyc clean-build docs lint help
.DEFAULT_GOAL := help
define BROWSER_PYSCRIPT
import os, webbrowser, sys
try:
	from urllib import pathname2url
except:
	from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT
BROWSER := python -c "$$BROWSER_PYSCRIPT"


##############################
# Global variables
##############################

PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
PROJECT_NAME = telomeres


ifeq (,$(shell which conda))
HAS_CONDA=False
else
HAS_CONDA=True
endif

##############################
# Commands
##############################
clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

lint: ## Lint using flake8
	flake8 src

docs: ## generate Sphinx HTML documentation
	sphinx-apidoc -o docs/
	$(MAKE) -C docs clean
	$(MAKE) -C docs html

showdocs: ## Show the documentation in a html browser
	$(BROWSER) docs/_build/html/index.html

docker: clean ## Make docker image
	docker build -t ${DOCKER_REPO}/${DOCKER_TAG} . --build-arg DOCKER_REPO=${DOCKER_REPO}

dockertag: docker ## Make docker image and tag
	docker tag ${id} ${DOCKER_REPO}/${DOCKER_TAG} ${DOCKER_REPO}/${DOCKER_TAG}:latest

dockerupload: clean dockertag ## Make docker image and upload
	docker push ${DOCKER_REPO}/${DOCKER_TAG}:latest


##############################
# Documentation
##############################
help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)
