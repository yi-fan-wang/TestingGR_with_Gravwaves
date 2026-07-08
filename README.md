[![Documentation Status](https://readthedocs.org/projects/ansicolortags/badge/?version=latest)](https://yi-fan-wang.github.io/TestingGR_with_Gravwaves)

# TGR: Testing General Relativity with Gravitational Waves

TGR is a collection of PyCBC waveform plugins for testing general relativity
with gravitational-wave data. It provides waveform approximants for modified
propagation, phenomenological beyond-GR corrections, ringdown models. The package
is intended for analyses that call waveforms through the [PyCBC waveform plugin
interface](http://pycbc.org/pycbc/latest/html/waveform_plugin.html), can be used for
Bayesian parameter estimation, matched-filtering, injection campaigns, etc.

## Installation

TGR is published on PyPI as `tgr`:

```bash
pip install tgr
```

For development from a local checkout:

```bash
git clone https://github.com/yi-fan-wang/TestingGR_with_Gravwaves.git
cd TestingGR_with_Gravwaves
pip install -e .
```

The package requires Python 3.11 or newer and installs PyCBC as a dependency.

## Waveform Plugins

After installation, TGR registers the following PyCBC approximants.

Frequency-domain waveform plugins:

- `birefringence`: gravitational-wave birefringence models
- `massive_graviton`: massive-graviton propagation corrections
- `fta`: flexible theory-agnostic waveform deformations
- `ppe`: parameterized post-Einsteinian corrections
- `lsa`: line-of-sight amplitude corrections

Time-domain waveform plugins:

- `nrsxs`: SXS numerical-relativity waveform interface
- `lvcnr`: LIGO/Virgo/KAGRA numerical-relativity waveform interface
- `NRSur7dq4_remove_qqnm`: NRSur7dq4 waveform with quadratic QNM contributions removed
- `NRSur7dq4_tdtaper`: tapered NRSur7dq4 waveform generation

## Documentation and Examples

Documentation is available at
[yi-fan-wang.github.io/TestingGR_with_Gravwaves](https://yi-fan-wang.github.io/TestingGR_with_Gravwaves).
Example notebooks and scripts are available in the
[examples](https://github.com/yi-fan-wang/TestingGR_with_Gravwaves/tree/master/examples)
directory.

## Usage in Scientific Publications

The `birefringence` is used in Refs.[1,2] for testing gravitational wave parity violation. The `NRSur7dq4_remove_qqnm` approximant is designed for subtracting predicted
quadratic quasi-normal-mode contributions from an `NRSur7dq4` waveform and used in the analysis of nonlinear ringdown evidence in GW250114 [3]. 

## References

[1] Tests of Gravitational-Wave Birefringence with the Open Gravitational-Wave Catalog.
[arXiv:2109.09718](https://arxiv.org/abs/2109.09718).

[2] Gravitational-Wave Implications for the Parity Symmetry of Gravity at GeV Scale.
[arXiv:2002.05668](https://arxiv.org/abs/2002.05668).

[3] A nonlinear voice from GW250114 ringdown.
[arXiv:2601.05734](https://arxiv.org/abs/2601.05734).
