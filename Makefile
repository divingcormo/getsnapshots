#TESTS=test_all.py
TESTS=tests

check:
	# No unused imports, no undefined vars
	flake8 --exit-zero src
	pylint src

test:
	py.test -v $(TESTS)

coverage:
	python -m pytest --cov=src --cov-report term-missing $(TESTS)

tox:
	tox -v

build: src/*/*.py setup.py MANIFEST.in
	python setup.py sdist bdist_wheel

push:
	git push -u origin master

clean:
	rm -r build dist 
	rm src/getsnapshots.egg-info/SOURCES.txt
	rm -r docs/build

distclean:
	rm -r build dist
	rm -r src/*.egg-info
	rm -r docs/build
	rm -r .pytest_cache
	rm -r tests/.pytest_cache
	rm -r .tox
	rm .coverage
