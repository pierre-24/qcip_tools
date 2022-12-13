all: help

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  init                        to install python dependencies through pipenv"
	@echo "  sync                        update dependencies of pipenv"
	@echo "  lint                        to lint backend code (flake8)"
	@echo "  test                        to run test suite"
	@echo "  doc                         to build documentation"
	@echo "  help                        to get this help"

init:
	pip3 install -e .

install-dev:
	pip-sync

lint:
	flake8 qcip_tools tests --max-line-length=120 --ignore=N802

test:
	python3 -m unittest discover -s tests

doc:
	cd docs; make html