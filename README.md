[![Documentation Status](https://readthedocs.org/projects/ansicolortags/badge/?version=latest)](https://yi-fan-wang.github.io/TestingGR_with_Gravwaves)
# PyTGR (PyTesting General Relativity) 

This code collects a variety of gravitational waveform templates in modified gravity, and can be easily called by [PyCBC](http://pycbc.org/) based on the [PyCBC waveform plugin](http://pycbc.org/pycbc/latest/html/waveform_plugin.html). This enables any subsequent gravitational wave data analysis with PyCBC, such as Bayesian parameter estimation with a non-GR waveform template, using matched-filtering to search for non-GR gravitational wave signals , etc. 

Currently, the models supported are:
 - Massive graviton
 - Gravitational wave birefringence `[1,2]`
 - dipole radiation
 - FTA (flexible theory agnostic)
   
Constribution are welcome!

## Usage

First install PyCBC by

`pip install pycbc`

Then install the birefringence waveforms for gravitational waves by

`python setup.py install`

More information should in principle be in the [documentation](https://yi-fan-wang.github.io/TestingGR_with_Gravwaves/html/index.html) and [examples](https://github.com/yi-fan-wang/TestingGR_with_Gravwaves/tree/master/examples), which are under construction.

## References

[1] Tests of Gravitational-Wave Birefringence with the Open Gravitational-Wave Catalog [arXiv](https://arxiv.org/abs/2109.09718). Wang et al.

[2] Gravitational-Wave Implications for the Parity Symmetry of Gravity at GeV Scale [arXiv](https://arxiv.org/abs/2002.05668). Wang et al.
