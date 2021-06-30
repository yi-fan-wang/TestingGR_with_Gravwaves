#!/usr/bin/env python
"""
setup.py file for dipole pycbc waveform plugin package
"""

from setuptools import Extension, setup, Command
from setuptools import find_packages

VERSION = '0.0.dev0'

setup (
    name = 'pycbc-dipole',
    version = VERSION,
    description = 'An waveform plugin for PyCBC',
    author = 'The PyCBC team',
    author_email = 'yifan.wang@aei.mpg.de',
    url = 'http://www.pycbc.org/',
    #download_url = 'https://github.com/gwastro/revchirp/tarball/v%s' % VERSION,
    keywords = ['pycbc', 'signal processing', 'gravitational waves'],
    install_requires = ['pycbc'],
    py_modules = ['dipole'],
    entry_points = {"pycbc.waveform.fd":"IMRPhenomdipole = dipole:genwav",
                    "pycbc.waveform.fd":"SEOBNRv4_ROM_dipole = dipole:genwav_seobnrv4_rom"},
)