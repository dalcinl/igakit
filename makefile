PYTHON=python
PACKAGE=igakit

.PHONY: build
build:
	$(PYTHON) setup.py build

.PHONY: test
test:
	$(PYTHON) test/runtests.py

.PHONY: install install_home install_user
install: install_home
install_home:
	$(PYTHON) setup.py install --home=${HOME}
install_user:
	$(PYTHON) setup.py install --user

.PHONY: uninstall uninstall_home uninstall_user
uninstall: uninstall_home uninstall_user
uninstall_home:
	-$(RM) -r $(HOME)/lib/python/$(PACKAGE)*
	-$(RM) -r $(HOME)/lib64/python/$(PACKAGE)*
uninstall_user:
	-$(RM) -r `$(PYTHON) -m site --user-site`/$(PACKAGE)*


.PHONY: distclean
clean:
	$(PYTHON) setup.py -q clean
	-$(RM) -r build/src.*
distclean: clean
	-$(RM) -r build dist MANIFEST *.py[co]
	-$(RM) -r __pycache__ */__pycache__
	-find test -name '*.py[co]' -exec rm -f {} ';'
	-find src  -name '*.py[co]' -exec rm -f {} ';'
