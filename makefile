.PHONY: test build upload distribute

test:
	python -m unittest discover tests

build:
	python setup.py sdist bdist_wheel

upload:
	twine upload dist/blazy-$$(python setup.py --version)*

distribute:
	make test
	make build
	make upload