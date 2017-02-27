
install-dependencies:
	pip install --upgrade -r requirements.txt -r requirements-dev.txt

lint-back:
	flake8 qcip_tools tests --max-line-length=120 --ignore=N802

test-back:
	python setup.py test
