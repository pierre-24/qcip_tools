all: help

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  install-dependencies        to install python dependencies through pipenv"
	@echo "  lint                        to lint backend code (flake8)"
	@echo "  tests                       to run test suite"
	@echo "  doc                         to build documentation"
	@echo "  help                        to get this help"

install-dependencies:
	pipenv install --dev

lint:
	pipenv run flake8 qcip_tools tests --max-line-length=120 --ignore=N802

tests:
	pipenv run python -m unittest discover -s tests

doc:
	cd documentation; pipenv run make html