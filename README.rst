igakit: Toolkit for IsoGeometric Analysis (IGA)
===============================================


Overview
--------

XXX To be written ...

Installation
------------

Clone the `Mercurial <http://mercurial.selenic.com/>`_ repository
hosted at `Bitbucket <https://bitbucket.org/dalcinl/igakit>`_ ::

  $ hg clone https://bitbucket.org/dalcinl/igakit

Alternatively, download_ the code as a tarball and unpack it::

  $ curl -O https://bitbucket.org/dalcinl/igakit/get/default.tar.gz
  $ tar -zxf default.tar.gz
  $ mv dalcinl-igakit-XXXXXXXXXXXX igakit

.. _download: https://bitbucket.org/dalcinl/igakit/get/default.tar.gz

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
