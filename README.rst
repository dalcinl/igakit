igakit: Toolkit for IsoGeometric Analysis (IGA)
===============================================


Overview
--------

This package implements many of the NURBS routines in Piegl's book
using Fortran. It provides Python bindings for these functions using
``f2py``. Finally we provide a NURBS class which uses this
functionality for the simplified, manually creation of geometry for
use in isogeometric analysis.

Installation
------------

Quick installation::

  $ pip install https://bitbucket.org/dalcinl/igakit/get/master.tar.gz

Alternatively, clone the repository hosted at
`Bitbucket <https://bitbucket.org/dalcinl/igakit>`_::

  $ git clone https://bitbucket.org/dalcinl/igakit

Enter the top level directory, next build and install the package
using the standard distutils's ``setup.py`` script::

  $ cd igakit
  $ python setup.py install --user

Acknowledgments
---------------

This project was partially supported by the Center for Numerical
Porous Media, Division of Computer, Electrical, and Mathematical
Sciences & Engineering (`CEMSE <http://cemse.kaust.edu.sa/>`_),
King Abdullah University of Science and Technology (`KAUST
<http://www.kaust.edu.sa/>`_).
