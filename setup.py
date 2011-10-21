#!/usr/bin/env python
"""
Toolkit for IsoGeometric Analysis (IGA)
"""

NAME    = 'igakit'
VERSION = '0.1'
AUTHOR  = 'Lisandro Dalcin'
EMAIL   = 'dalcinl@gmail.com'
DESCR   = __doc__.strip()

def setup_package():
    from numpy.distutils.core import setup
    from numpy.distutils.core import Extension
    setup(name=NAME,
          version=VERSION,
          author=AUTHOR,
          author_email=EMAIL,
          description=DESCR,
          long_description=DESCR,
          packages = ['igakit'],
          package_dir = {'igakit' : 'src/igakit'},
          ext_modules  = [
            Extension('igakit.igalib',
                      sources = ['src/igakit/igalib.pyf', 
                                 'src/igakit/igalib.f90'],
                      f2py_options = ['--quiet'],
                      define_macros = [#('F2PY_REPORT_ATEXIT', 0),
                                       ('F2PY_REPORT_ON_ARRAY_COPY', 0)]),
            ],
          )
    
if __name__ == "__main__":
    setup_package()
