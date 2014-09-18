PYTHON=python
PACKAGE=igakit

.PHONY: build
build:
	$(PYTHON) setup.py build

BACKEND=none
.PHONY: testdemo
testdemo:
	$(PYTHON) demo/plot_crv.py $(BACKEND)
	$(PYTHON) demo/plot_srf.py $(BACKEND)
	$(PYTHON) demo/plot_vol.py $(BACKEND)
	$(PYTHON) demo/pipe.py     $(BACKEND)
	$(PYTHON) demo/bentpipe.py $(BACKEND)
	$(PYTHON) demo/unclamp.py  $(BACKEND)
	$(PYTHON) demo/sweep.py    $(BACKEND)
	$(PYTHON) demo/ruled.py    $(BACKEND)
	$(PYTHON) demo/revolve.py  $(BACKEND)
	$(PYTHON) demo/refine.py   $(BACKEND)
.PHONY: doctest
doctest:
	$(PYTHON) -m doctest src/igakit/transform.py
	$(PYTHON) -m doctest src/igakit/nurbs.py
	$(PYTHON) -m doctest src/igakit/cad.py
	-@$(RM) -r src/igakit/__pycache__
	-@$(RM) -r src/igakit/*.py[co]
.PHONY: unittest
unittest:
	$(PYTHON) test/runtests.py

.PHONY: test testall
test: unittest
testall: test testdemo


.PHONY: install
install:
	$(PYTHON) setup.py install --user

.PHONY: uninstall
uninstall:
	-$(RM) -r $(shell $(PYTHON) -m site --user-site)/$(PACKAGE)
	-$(RM) -r $(shell $(PYTHON) -m site --user-site)/$(PACKAGE)-*-py*.egg-info


.PHONY: clean distclean
clean:
	$(PYTHON) setup.py -q clean
	-$(RM) -r build/src.*
distclean: clean
	-$(RM) -r build dist MANIFEST *.py[co]
	-$(RM) -r __pycache__ */__pycache__
	-find test -name '*.py[co]' -exec rm -f {} ';'
	-find src  -name '*.py[co]' -exec rm -f {} ';'
