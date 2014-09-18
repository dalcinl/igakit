#!/usr/bin/env python
"""
Toolkit for IsoGeometric Analysis (IGA)
"""

url = 'https://bitbucket.org/dalcinl/igakit'

setup_args = dict(
    name             = 'igakit',
    version          = '0.1',
    description      = __doc__.strip(),
    long_description = __doc__.strip(),
    author           = 'Lisandro Dalcin',
    author_email     = 'dalcinl@gmail.com',
    license          = 'BSD',
    keywords         = ['FEM', 'IGA'],
    url              = url,
    download_url     = url+'/get/default.tar.gz',
)

def setup_package():
    import sys
    from numpy.distutils.core import setup
    from numpy.distutils.core import Extension
    if 'setuptools' in sys.modules:
        setup_args['install_requires'] = ['numpy']
    setup(packages = ['igakit'],
          ext_modules  = [
              Extension('igakit.igalib',
                        sources = [
                            'igakit/igalib.pyf',
                            'igakit/igalib.f90',
                        ],
                        define_macros = [
                            #('F2PY_REPORT_ATEXIT', 0),
                            #('F2PY_REPORT_ON_ARRAY_COPY', 0),
                        ],
              ),
          ],
          **setup_args)

if __name__ == "__main__":
    setup_package()
