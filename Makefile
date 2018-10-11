#TESTS=test_all.py
TESTS=tests

check:
	# No unused imports, no undefined vars
	flake8 pyrmsd.py get_snapshots.py
	pylint pyrmsd.py get_snapshots.py

test:
	py.test -v $(TESTS)

coverage:
	py.test -v --cov=. --cov-report term-missing $(TESTS)

build:
	python setup.py sdist bdist_wheel

docs:
	python setup.py build_sphinx

push:
	git push -u origin master
