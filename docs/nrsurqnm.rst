.. _nrsurqnm:

Ringdown and quadratic QNM models (``tgr.nrsurqnm``)
======================================================

The :mod:`tgr.nrsurqnm` module builds beyond-GR ringdown waveforms on top of
the numerical-relativity surrogate ``NRSur7dq4``. The full inspiral-merger-ringdown
signal is generated with the surrogate, and the late-time ringdown of selected
spherical-harmonic modes is then *re-modelled* as a superposition of quasi-normal
modes (QNMs). Working at the level of the QNMs allows controlled, physically
meaningful deviations from General Relativity (GR) to be injected mode by mode,
which is exactly what a test of GR requires.

This page describes the physics and the mathematics implemented by the module.
The auto-generated API reference for every function lives under
:mod:`tgr.nrsurqnm` in the :doc:`tgr` package documentation.

Quasi-normal-mode ansatz
------------------------

After merger, each spherical-harmonic mode of the strain rings down as a sum of
damped sinusoids. A single QNM labelled by :math:`(\ell, m, n)` (harmonic
indices and overtone number :math:`n`) is

.. math::
   :label: qnm-ansatz

   \psi_{\ell m n}(t) = A_{\ell m n}\,
       \exp\!\big[-i\,\omega_{\ell m n}\,(t - t_0)\big],
   \qquad t \ge t_0 ,

where :math:`t_0` is the ringdown start time and :math:`A_{\ell m n}\in\mathbb{C}`
is a complex amplitude that fixes both the size and the initial phase of the mode.
The **complex angular frequency** packages the oscillation frequency
:math:`f_{\ell m n}` and the damping time :math:`\tau_{\ell m n}` into a single
quantity,

.. math::
   :label: complex-omega

   \omega_{\ell m n} = 2\pi f_{\ell m n} - \frac{i}{\tau_{\ell m n}} ,

so that the real part drives the oscillation and the imaginary part produces the
exponential decay,

.. math::

   \psi_{\ell m n}(t) = A_{\ell m n}\,
       e^{-i\,2\pi f_{\ell m n}(t-t_0)}\;
       e^{-(t-t_0)/\tau_{\ell m n}} .

QNM spectrum of the remnant: ``get_qnm_freqtau`` and ``QNMTable``
-----------------------------------------------------------------

The QNM frequencies and damping times are uniquely determined by the mass
:math:`M_f` and dimensionless spin :math:`\chi_f` of the remnant black hole. The
helper :func:`tgr.nrsurqnm.get_qnm_freqtau` first predicts :math:`(M_f, \chi_f)`
from the initial binary parameters using the ``NRSur7dq4`` remnant fits, and then
maps them to :math:`(f_{\ell m n}, \tau_{\ell m n})` for each requested mode.

Two kinds of mode labels are supported:

* **Linear (fundamental/overtone) modes** — a 3-character label ``"lmn"`` such as
  ``"220"`` or ``"221"``, evaluated directly as
  :math:`(f_{\ell m n}, \tau_{\ell m n})`.

* **Quadratic modes (QQNMs)** — a 6-character label that concatenates two linear
  modes, e.g. ``"220220"`` for :math:`220\times220`. A quadratic mode arises from
  the second-order coupling of two parent linear modes and rings at the *sum*
  frequency with the *combined* decay rate,

  .. math::
     :label: quadratic-ftau

     f_{(1)(2)} = f_{(1)} + f_{(2)} ,
     \qquad
     \frac{1}{\tau_{(1)(2)}} = \frac{1}{\tau_{(1)}} + \frac{1}{\tau_{(2)}} .

The result is returned as a :class:`tgr.nrsurqnm.QNMTable`, a small frozen
dataclass that carries the remnant parameters (``final_mass``, ``final_spin``)
together with the per-mode dictionaries ``freq`` and ``tau``. It offers three
conveniences used throughout the module:

* ``table.modes`` — the list of tabulated mode labels;
* ``table.omega(mode)`` — the complex angular frequency :eq:`complex-omega`;
* ``table.subset(modes)`` — a new table restricted to the given modes.

Because evaluating the remnant fits has a non-negligible cost and the waveform
generators run inside stochastic samplers, the table is computed **once per
waveform call** and passed down to every routine that needs it.

Least-squares QNM decomposition: ``least_square_qnmfitting``
------------------------------------------------------------

Given a target complex multipole :math:`h(t)` (e.g.
:math:`h_{22} = h^{22}_{\rm real} + i\,h^{22}_{\rm imag}` from the surrogate),
a :class:`~tgr.nrsurqnm.QNMTable`, a start time :math:`t_0`, and an explicit
list of modes to fit,
:func:`tgr.nrsurqnm.least_square_qnmfitting` extracts the complex amplitudes by
a linear least-squares fit over the ringdown window
:math:`t \in [t_0,\, t_0 + T]`, where the window length
:math:`T = \ln(1000)\,\max_k \tau_k` is set by the slowest-decaying fitted mode
(i.e. it extends until the dominant amplitude has dropped by a factor of
:math:`10^3`).

The fitted modes must be given explicitly (``fit_modes``) and every label must
be present in the table — the table may well contain more modes than are
fitted, e.g. the quadratic modes used elsewhere in the same waveform call; a
missing label raises a :class:`ValueError` listing the available modes.

Sampling the templates :eq:`qnm-ansatz` (with unit amplitude and an overall
dynamical-range factor :math:`\eta=10^{-22}` so that the normal equations are
well conditioned for strain-scale data) on the :math:`N` time samples
:math:`t_j` of the window defines the design matrix
:math:`G\in\mathbb{C}^{N\times M}`,

.. math::

   G_{jk} = \eta\,\exp\!\big[-i\,\omega_{k}\,(t_j - t_{\rm grid})\big] ,

where :math:`t_{\rm grid} \le t_0` is the first sample of the window. The
amplitudes :math:`A = (A_1,\dots,A_M)^\top` are the minimiser of the residual
:math:`\lVert G A - h \rVert_2^2`, i.e. the ordinary least-squares solution of
the normal equations

.. math::
   :label: normal-eq

   \big(G^{H} G\big)\, A = G^{H} h ,
   \qquad
   A = \big(G^{H} G\big)^{-1} G^{H} h ,

where :math:`G^{H}` denotes the conjugate transpose. The fitted coefficients
are finally rescaled by :math:`\eta` and propagated from the sample grid to the
requested reference time :math:`t_0` through the phase factor
:math:`e^{-i\omega_k (t_0 - t_{\rm grid})}`.

The routine returns three objects: the dictionary of complex amplitudes
referenced to :math:`t_0`; the best-fit reconstruction
:math:`\sum_k A_k\,e^{-i\omega_k(t - t_{\rm grid})}` over the fit window (a
complex ``TimeSeries``, useful for residual checks); and the slice of the
target waveform that entered the fit.

Replacing the (4,4) ringdown: ``gen_nrsurqnm``
----------------------------------------------

:func:`tgr.nrsurqnm.gen_nrsurqnm` generates the full ``NRSur7dq4`` signal,
keeps every mode except :math:`(4,\pm 4)`, and reconstructs the
:math:`(4,4)` ringdown from its QNM decomposition :eq:`normal-eq`. Deviations from
GR are injected through per-mode fractional amplitude parameters
:math:`\delta_{\ell m n}` (keyword arguments ``delta_<mode>``):

.. math::
   :label: linear-tgr

   A_{\ell m n} \;\longrightarrow\; A_{\ell m n}\,\big(1 + \delta_{\ell m n}\big),
   \qquad \delta_{\ell m n}=0 \ \text{recovers GR.}

The reconstructed mode is mapped back to the two polarizations with the
spin-weighted spherical harmonics of weight :math:`-2` (see
:ref:`nrsurqnm-recombine`).

Removing quadratic modes: ``gen_nrsur_remove_qqnm``
---------------------------------------------------

:func:`tgr.nrsurqnm.gen_nrsur_remove_qqnm` isolates the contribution of the
quadratic modes to the :math:`(4,4)` ringdown so that a GR deviation can be
applied to them. The waveform keywords ``mode22`` and ``mode_quadratic`` are
space-separated label strings selecting, respectively, the parent
:math:`(2,2,n)` overtones entering the least-squares fit and the quadratic
modes to subtract. Amplitude-ratio tables exist for the six quadratic modes in
``QUADRATIC_MODES`` (``220220``, ``220221``, ``221221``, ``220222``,
``221222``, ``222222``); they are interpolated in remnant spin from data files
shipped with the package and loaded once at import time. Both mode lists are
validated up front: an unknown quadratic mode, or a quadratic mode whose parent
is missing from ``mode22``, raises a :class:`ValueError` before any waveform is
generated. Note that ``mode22`` may contain *more* overtones than the quadratic
modes require (e.g. ``"220 221 222 223"``): the extra overtones only stabilise
the fit and do not enter the amplitude products.

**Parent amplitudes.** The parent amplitudes :math:`A_{(1)}, A_{(2)}` are
measured from the surrogate :math:`h_{22}` with
:func:`~tgr.nrsurqnm.least_square_qnmfitting`. Overtone fits are only reliable
for start times :math:`t_0 \ge` ``FIT_TSTART_MIN`` (2 ms after the peak). For
``toffset`` at or beyond this threshold a single fit at :math:`t_0` is used.
For earlier start times the modes are instead fitted on the grid
``FIT_TSTART_GRID`` (six start times between 2.00 and 3.67 ms), each fit is
propagated back to :math:`t_0` with the phase factor
:math:`e^{-i\omega(t_0 - t_{\rm fit})}`, and one amplitude per mode is drawn
from a Gaussian matching the mean and scatter of these fits (real and imaginary
parts independently). This marginalises over the fit-start systematics, at the
price of making the returned waveform *stochastic*; passing the optional
``seed`` keyword makes the draw reproducible.

**Quadratic amplitudes.** The amplitude of a quadratic mode is fixed by the
amplitudes of its two parent linear :math:`(2,2)` modes, a theory-predicted
complex ratio :math:`R_{(1)(2)}(\chi_f)` (interpolated as a function of remnant
spin), and a unit conversion to the surrogate's dimensionless strain,

.. math::
   :label: quadratic-amp

   A_{(1)(2)} = A_{(1)}\,A_{(2)}\,R_{(1)(2)}(\chi_f)\,\mathcal{C},
   \qquad
   \mathcal{C} = \frac{d_L\,(\mathrm{Mpc}\!\to\!\mathrm{m})}
                      {M_f\,\big(M_\odot\!\to\!\mathrm{m}\big)} ,

where :math:`d_L` is the luminosity distance and :math:`\mathcal{C}` rescales the
physical amplitudes back to surrogate units. Each quadratic mode is then built as
a damped sinusoid :eq:`qnm-ansatz` with parameters :eq:`quadratic-ftau` and
subtracted from the :math:`(4,4)` ringdown, from ``toffset`` to the end of the
waveform.

**Amplitude and phase deviation.** The subtracted quadratic contribution is scaled
by a single *complex* deviation factor

.. math::
   :label: qqnm-tgr

   \kappa = a_{\rm TGR}\,e^{\,i\,\varphi_{\rm TGR}} ,

controlled by the keyword arguments ``quadratic_tgr`` (amplitude
:math:`a_{\rm TGR}`) and ``quadratic_tgr_phase`` (phase
:math:`\varphi_{\rm TGR}`, in radians). The GR limit is
:math:`a_{\rm TGR}=1,\ \varphi_{\rm TGR}=0`, i.e. :math:`\kappa = 1`, which removes
exactly the GR-predicted quadratic mode. Tuning :math:`a_{\rm TGR}` rescales how
much of the mode is removed, while :math:`\varphi_{\rm TGR}` rotates it in the
complex plane, probing a phase offset of the quadratic coupling.

**Frequency and damping-time deviation.** Optionally, a quadratic mode with
*non-GR* spectral parameters can be re-injected. If either of the keyword
arguments ``qqnm_deltaf`` and ``qqnm_deltatau`` is given, a quadratic mode with
the full theory amplitude :eq:`quadratic-amp` but fractionally shifted spectrum,

.. math::
   :label: qqnm-ftau-tgr

   f \;\longrightarrow\; f\,(1 + \delta f),
   \qquad
   \tau \;\longrightarrow\; \tau\,(1 + \delta\tau),

is added back after the subtraction, so the residual :math:`(4,4)` ringdown
retains a quadratic mode that oscillates and decays at shifted rates relative to
the GR prediction.

.. _nrsurqnm-recombine:

Recombination into polarizations
--------------------------------

A mode :math:`(\ell, m)` in the co-precessing/source frame is projected onto the
observer frame using the spin-weight :math:`-2` spherical harmonics
:math:`{}_{-2}Y_{\ell m}(\iota, \phi)`, evaluated at the inclination
:math:`\iota` and a reference phase :math:`\phi = \pi/2 - \phi_{\rm coa}`. The
modelled :math:`(4,4)`/:math:`(4,-4)` ringdown is recombined as

.. math::
   :label: recombine

   h_{44}^{\rm obs}(t) =
       h_{44}(t)\, {}_{-2}Y_{4,4}(\iota,\phi)
     + h_{44}^{*}(t)\, {}_{-2}Y_{4,-4}(\iota,\phi) ,

and the complex strain :math:`h = h_+ - i\,h_\times` finally yields the two
polarizations returned by every generator,
:math:`h_+ = \mathrm{Re}\,h` and :math:`h_\times = -\,\mathrm{Im}\,h`.

Tapered surrogate helper
------------------------

For convenience :func:`tgr.nrsurqnm.gen_nrsur7dq4_tdtaper` returns the plain
``NRSur7dq4`` polarizations with a short taper of length ``window``
(default :math:`0.05\,\mathrm{s}`) applied at the start of the time series to
suppress turn-on transients.
