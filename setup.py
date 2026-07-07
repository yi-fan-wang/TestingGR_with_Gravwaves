"""
setup.py file for testing general relativity pycbc waveform plugin package
"""

from pathlib import Path

from setuptools import setup, find_packages

VERSION = '1.1.0'
README = Path(__file__).with_name('README.md').read_text(encoding='utf-8')

setup (
    name = 'tgr',
    version = VERSION,
    description = 'A waveform plugin for PyCBC with non-general relativity templates',
    long_description = README,
    long_description_content_type = 'text/markdown',
    author = 'Yifan Wang',
    author_email = 'yifan.wang@aei.mpg.de',
    url = 'https://github.com/yi-fan-wang/TestingGR_with_Gravwaves',
    keywords = ['testing general relativity', 'gravitational waves', 'pycbc'],
    license = 'GPL-3.0-only',
    packages = find_packages(),
    package_data={
        'tgr': ['data/*.npy'],
    },
    include_package_data=True,
    entry_points = {"pycbc.waveform.fd":["birefringence = tgr.birefringence:gen_waveform",
                                         "massive_graviton = tgr.massive_graviton:gen_mg_waveform",
                                         "fta = tgr.fta:gen_waveform",
                                         "ppe = tgr.ppe:gen_ppe_waveform",
                                         "lsa = tgr.lineofsight:gen_waveform"
                                         ],
                    "pycbc.waveform.td":["nrsxs = tgr.nr:gen_sxs_waveform",
                                         "lvcnr = tgr.nr:gen_lvcnr_waveform",
                                         "NRSur7dq4QNM = tgr.nrsurqnm:gen_nrsurqnm",
                                         "NRSur7dq4_remove_qqnm = tgr.nrsurqnm:gen_nrsur_remove_qqnm",
                                         "NRSur7dq4_tdtaper = tgr.nrsurqnm:gen_nrsur7dq4_tdtaper"
                                         ],
                    "pycbc.waveform.length":["birefringence = tgr:length_in_time",
		    			                     "massive_graviton = tgr:length_in_time",
					                         "fta = tgr:length_in_time",
                                             "ppe = tgr:length_in_time",
                                             "lsa = tgr:length_in_time"]},
    python_requires='>=3.11',
    classifiers=[
        'Programming Language :: Python :: 3.11',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    install_requires=[
        "numpy",
        "scipy",
        "pycbc",
    ],
)
