"""
setup.py file for testing birefringence pycbc waveform plugin package
"""

from setuptools import Extension, setup, Command, find_packages

VERSION = '0.2'

setup (
    name = 'PyTGR',
    version = VERSION,
    description = 'A waveform plugin for PyCBC with non-general relativity templates',
    author = 'Yifan Wang',
    author_email = 'yifan.wang@aei.mpg.de',
    url = 'https://github.com/yi-fan-wang/TestingGR_with_Gravwaves',
    keywords = ['testing general relativity', 'gravitational waves', 'pycbc'],
    packages = find_packages(),
    package_data={  
        'pytgr': ['data/*.npy'],
    },
    include_package_data=True,
    #py_modules = ['birefringence'],
    #package_dir = {'':'src'},
    #package_dir={'PyTGR': 'src'},
    entry_points = {"pycbc.waveform.fd":["birefringence = pytgr.birefringence:gen_waveform",
                                         "massivegraviton = pytgr.massivegraviton:gen_waveform",
                                         "fta = pytgr.fta:gen_waveform",
                                         "ppe = pytgr.ppe:gen_waveform",
                                         "lsa = pytgr.lineofsight:gen_waveform"
                                         ],
                    "pycbc.waveform.td":["nrsxs = pytgr.nr:gen_sxs_waveform",
                                         "lvcnr = pytgr.nr:gen_lvcnr_waveform",
                                         "NRSur7dq4QNM = pytgr.nrsurqnm:gen_nrsurqnm",
                                         "NRSur7dq4_quadratic = pytgr.nrsurqnm:gen_nrsur_linearqnm",
                                         "NRSur7dq4_tdtaper = pytgr.nrsurqnm:gen_nrsur7dq4_tdtaper"
                                         ],
                    "pycbc.waveform.length":["birefringence = pytgr:length_in_time",
		    			                     "massivegraviton = pytgr:length_in_time",
					                         "fta = pytgr:length_in_time",
                                             "ppe = pytgr:length_in_time",
                                             "lsa = pytgr:length_in_time"]},
    python_requires='>=3.7',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
    install_requires=[
        "scipy",
        "pycbc",
    ],
)
