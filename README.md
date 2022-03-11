# TestingGR_with_Gravwaves

General Relativity is the most beautiful and successful theory to describe gravity in modern physics. Despite its success, there are motivations to go beyond General Relativity both from theoretical and observational sides, such as addressing the singularity problem, quantizing gravity, accounting for the nature of dark matter and dark energy, etc. Since the first detection in 2015, gravitational waves provides a golden oppoturnity to test General Relativity in a brand new regime and look for signs of new physics. One of the effective ways to test General Relativity with gravitational waves is to construct gravitational wave templates predicted by modified gravity theories, then compare templates with data. Within a Bayesian inference framework, the preference of data towards General Relativity or a modified gravity can be established by quantitative indicators, such as Bayes Odds.

Therefore, to test General Relativity, it's crucial to build beyond-General-Relativity waveforms with clear physics meanings first. This repository collects different waveform models, and make use of the [PyCBC waveform plugin](http://pycbc.org/pycbc/latest/html/waveform_plugin.html) to make models avaiable to analyze real gravitational-wave data. Such that it's reproducible and serves as examples for the community to build any new waveform models for data analysis. 

Constribution are welcome!


## Waveform models
Currently, this repository includes the following gravitational wave waveform models for compact binary coalescence:

- Gravitational wave birefreingence `[1,2]`

[1] Tests of Gravitational-Wave Birefringence with the Open Gravitational-Wave Catalog [arXiv](https://arxiv.org/abs/2109.09718). Wang et al.

[2] Gravitational-Wave Implications for the Parity Symmetry of Gravity at GeV Scale [arXiv](https://arxiv.org/abs/2002.05668). Wang et al.

- Gravitational wave dipole radiation (to be added)

## Using the waveform to analyze GW with PyCBC

First install PyCBC by

`pip install pycbc`

Then install the birefringence waveforms for gravitational waves by

`python setup.py install`

## Documentation

 - [yi-fan-wang.github.io/TestingGR_with_Gravwaves](https://yi-fan-wang.github.io/TestingGR_with_Gravwaves) (to be updated)