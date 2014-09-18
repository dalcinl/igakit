#!/usr/bin/env python
import sys, os
import unittest

def getbuilddir():
    from distutils.util import get_platform
    s = os.path.join("build", "lib.%s-%.3s" % (get_platform(), sys.version))
    if hasattr(sys, 'gettotalrefcount'): s += '-pydebug'
    return s

def bootstrap():
    try:
        import igakit
    except ImportError:
        from os.path import dirname, abspath, join
        top_dir = abspath(join(dirname(__file__), '..'))
        build_dir = join(top_dir, getbuilddir())
        sys.path.insert(0, build_dir)
        import igakit
        sys.stderr.write(
            "Loaded package 'igakit' from build/ directory\n")
    return igakit

bootstrap()
