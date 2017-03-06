all: help

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  install-dependenciess       to install python dependencies through pip"
	@echo "  lint-back                   to lint backend code (flake8)"
	@echo "  test-back                   to run test suite"
	@echo "  doc                         to build documentation"
	@echo "  help                        to get this help"

install-dependencies:
	pip install --upgrade -r requirements.txt -r requirements-dev.txt

lint-back:
	flake8 qcip_tools tests --max-line-length=120 --ignore=N802

test-back:
	python setup.py test

doc:
	cd documentation; make html