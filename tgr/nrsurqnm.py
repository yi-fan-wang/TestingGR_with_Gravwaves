from dataclasses import dataclass
from importlib.resources import files
import numpy as np
import lal
from scipy.interpolate import interp1d

from pycbc.conversions import get_final_from_initial, get_lm_f0tau
from pycbc.types import (TimeSeries, complex64, zeros)

def gen_nrsurqnm(**kwds):
    from pycbc.waveform import get_td_waveform_modes

    hlm = get_td_waveform_modes(approximant='NRSur7dq4', mass1=kwds['mass1'], mass2=kwds['mass2'],
                            spin1x=kwds['spin1x'], spin1y=kwds['spin1y'], spin1z=kwds['spin1z'],
                            spin2x=kwds['spin2x'], spin2y=kwds['spin2y'], spin2z=kwds['spin2z'],
                            distance = kwds['distance'],
                            delta_t=kwds['delta_t'], f_lower=kwds['f_lower'],
                            mode_array=['22','21','20','33','32','31','30','44','43','42','41','40'])

    h = 0
    for l in range(2,5):
        for m in range(-1*l, l+1):
            if l==4 and abs(m)==4:
                continue # skip 44 mode for ringdown treatment later
            h_modes = hlm[(l,m)][0] + 1j * hlm[(l,m)][1]
            Y_lm = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, l, m)
            h += h_modes * Y_lm
    if 'ringdown_mode' not in kwds or kwds['ringdown_mode'] is None:
        return h.real(), -h.imag()

    h44 = hlm[(4,4)][0] + 1j * hlm[(4,4)][1]
    qnm_par = {}
    qnm_par['final_mass'], qnm_par['final_spin'] = get_final_from_initial(
       kwds['mass1'], kwds['mass2'],
       kwds['spin1x'], kwds['spin1y'], kwds['spin1z'],
       kwds['spin2x'], kwds['spin2y'], kwds['spin2z'],
       approximant='NRSur7dq4',
       f_ref=kwds['f_ref'])

    qnm_par['freq'] = {}
    qnm_par['tau'] = {}
    for mode in kwds['ringdown_mode']:
        if len(mode) == 3:
            l = int(mode[0])
            m = int(mode[1])
            n = int(mode[2])
            qnm_par['freq'][mode], qnm_par['tau'][mode] = get_lm_f0tau(qnm_par['final_mass'],
                                                                         qnm_par['final_spin'],
                                                                         l,m,n)
        elif len(mode) == 6:
            l1 = int(mode[0])
            m1 = int(mode[1])
            n1 = int(mode[2])
            l2 = int(mode[3])
            m2 = int(mode[4])
            n2 = int(mode[5])
            f1, tau1 = get_lm_f0tau(qnm_par['final_mass'], qnm_par['final_spin'], l1,m1,n1)
            f2, tau2 = get_lm_f0tau(qnm_par['final_mass'], qnm_par['final_spin'], l2,m2,n2)
            qnm_par['freq'][mode] = f1 + f2
            qnm_par['tau'][mode] = 1 / (1/tau1 + 1/tau2)
        else:
            raise ValueError("Invalid mode format in rindown_mode")
    #print("QNM parameters:", qnm_par)
    t_final = max(qnm_par['tau'][m] for m in kwds['ringdown_mode']) * np.log(1000)
    ringdown_start_time = kwds['toffset']
    start_idx = int(np.floor(float(ringdown_start_time - h44.start_time) * h44.sample_rate))
    end_idx = int(float(ringdown_start_time + t_final - h44.start_time) * h44.sample_rate)

    start_time = h44.sample_times[start_idx]
    end_time = h44.sample_times[end_idx]
    h44_slice = h44.time_slice(start_time, end_time)
    N = len(h44_slice)

    qnm = {}
    for m in kwds['ringdown_mode']:
        qnm[m] = TimeSeries(zeros(N, dtype=complex64), delta_t=kwds['delta_t'], epoch = start_time)
        sample_times = qnm[m].sample_times.numpy()
        omega = 2 * np.pi * qnm_par['freq'][m] - 1j / qnm_par['tau'][m]
        qnm[m].data = 1e-22 * np.exp(-1j * omega * (sample_times - start_time) )

    M = len(kwds['ringdown_mode'])
    G = np.zeros((N, M), dtype=complex)
    #print(N,M)
    for k, m in enumerate(kwds['ringdown_mode']):
        G[:, k] = qnm[m].data

    G_H = G.conj().T
    A = np.linalg.solve(G_H @ G, G_H @ h44_slice.data)

    allqnm = 0
    for i, m in enumerate(kwds['ringdown_mode']):
        delta_param = 'delta_' + m
        if delta_param in kwds and kwds[delta_param] is not None:
            #print("Applying deviation", delta_param, kwds[delta_param])
            A[i] *= (1 + kwds[delta_param])
        qnm[m].data *= A[i]
        allqnm += qnm[m]
    #print(allqnm.data)
    Y_44 = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, 4, 4)
    Y_4m4 = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, 4, -4)
    qnm44 = allqnm *  Y_44 + np.conj(allqnm) * Y_4m4

    h.data[start_idx:start_idx+len(qnm44)] += qnm44
    #print(qnm44.data)
    return h.real(), -h.imag()


def gen_nrsur7dq4_tdtaper(**kwds):
    from pycbc.waveform import get_td_waveform
    import pycbc.waveform.utils

    if 'approximant' in kwds:
        kwds.pop("approximant")
    hp, hc = get_td_waveform(approximant='NRSur7dq4', **kwds)

    if 'window' not in kwds or kwds['window'] is None:
        window = 0.05
    else:
        window = kwds['window']

    hp_taper = pycbc.waveform.utils.td_taper(hp, hp.start_time, hp.start_time + window)
    hc_taper = pycbc.waveform.utils.td_taper(hc, hc.start_time, hc.start_time + window)

    return hp_taper, hc_taper

@dataclass(frozen=True)
class QNMTable:
    '''QNM frequencies and damping times for a set of modes, and the remnant's
    mass and spin'''
    final_mass: float
    final_spin: float
    freq: dict
    tau: dict

    @property
    def modes(self):
        return list(self.freq)

    def omega(self, mode):
        '''Complex frequency of a QNM mode.'''
        return 2 * np.pi * self.freq[mode] - 1j / self.tau[mode]

    def subset(self, modes):
        return QNMTable(self.final_mass, self.final_spin,
                        {m: self.freq[m] for m in modes},
                        {m: self.tau[m] for m in modes})

def get_qnm_freqtau(qnm_modes, **kwds):
    '''Get the QNM parameters (frequency and damping time) for given modes from
    the initial binary parameters using NRSur7dq4 fits.

    Parameters
    ----------
    qnm_modes : list of str
        List of mode labels to get the QNM parameters for, e.g. ["220","221"] for
        linear modes, or ["220220", "220221"] for quadratic modes.
    kwds : dict
        Dictionary containing the initial binary parameters and f_ref for NRSur7dq4.

    Returns
    -------
    QNMTable
    '''
    final_mass, final_spin = \
        get_final_from_initial(kwds['mass1'], kwds['mass2'],
                               kwds['spin1x'], kwds['spin1y'], kwds['spin1z'],
                               kwds['spin2x'], kwds['spin2y'], kwds['spin2z'],
                               approximant='NRSur7dq4',
                               f_ref=kwds['f_ref'])
    freq = {}
    tau = {}
    for mode in qnm_modes:
        if len(mode) == 3:
            l = int(mode[0])
            m = int(mode[1])
            n = int(mode[2])
            freq[mode], tau[mode] = get_lm_f0tau(final_mass, final_spin, l,m,n)
        elif len(mode) == 6:
            l1 = int(mode[0])
            m1 = int(mode[1])
            n1 = int(mode[2])
            l2 = int(mode[3])
            m2 = int(mode[4])
            n2 = int(mode[5])
            f1, tau1 = get_lm_f0tau(final_mass, final_spin, l1,m1,n1)
            f2, tau2 = get_lm_f0tau(final_mass, final_spin, l2,m2,n2)
            freq[mode] = f1 + f2
            tau[mode] = 1 / (1/tau1 + 1/tau2)
        else:
            raise ValueError("Invalid mode format in rindown_mode")
    return QNMTable(final_mass, final_spin, freq, tau)

def least_square_qnmfitting(fit_modes, qnm_par, ringdown_start_time, h_target):
    '''Least-squares fit of a sum of QNMs to a target waveform::

        h_target(t) ~ sum_m A_m exp[-i omega_m (t - t0)]

    The complex QNM frequencies omega_m = 2 pi f_m - i / tau_m are fixed by
    qnm_par, and the complex amplitudes A_m solved by linear least square
    fitting. The fit window starts at t0 and lasts max_m(tau_m) * ln(1000),
    i.e. until the slowest-decaying mode has damped by a factor of 1000.

    Parameters
    ----------
    fit_modes : list of str
        Mode labels to fit, e.g. ['220', '221', '222']
    qnm_par : QNMTable
        Only freq and tau of the requested modes are used.
    ringdown_start_time : float
        Start time t0 of the ringdown model in unit of seconds. The fit
        itself starts on the sample grid at the last sample <= t0; the
        returned amplitudes are then propagated to exactly t0 via
        exp[-i omega_m (t0 - t_grid)].
    h_target : pycbc TimeSeries
        Multipole waveform to decompose, e.g. h22 = hlm_real + 1j * hlm_cross.
        Internally the templates are rescaled with a factor of 1e-22 to
        keep the normal equations well conditioned for GW strain scale
        data; this factor is undone in the returned amplitudes.

    Returns
    -------
    A_modes : dict
        mode label -> complex amplitude A_m referenced to t0
    sum_fit_qnm : pycbc TimeSeries
        Best-fit reconstruction sum_m A_m exp[-i omega_m (t - t_grid)] over
        the fit window (epoch = first fitted sample).
    h_target_slice : pycbc TimeSeries
        The slice of h_target used in the fit
    '''

    missing = [m for m in fit_modes if m not in qnm_par.modes]
    if missing:
        raise ValueError(f"modes {missing} not found in the QNMTable "
                         f"(available: {qnm_par.modes})")

    t_duration = max(qnm_par.tau[m] for m in fit_modes) * np.log(1000)
    start_idx = int(np.floor(float(ringdown_start_time - h_target.start_time) * h_target.sample_rate))
    end_idx = int(float(ringdown_start_time + t_duration - h_target.start_time) * h_target.sample_rate)

    start_time = h_target.sample_times[start_idx]
    end_time = h_target.sample_times[end_idx]
    h_target_slice = h_target.time_slice(start_time, end_time)
    N = len(h_target_slice) # number of samples in the fit window

    qnm = {}
    dynamical_range = 1e-22
    for m in fit_modes:
        qnm[m] = TimeSeries(zeros(N, dtype=complex64), delta_t=h_target.delta_t, epoch = start_time)
        sample_times = qnm[m].sample_times.numpy()
        qnm[m].data = dynamical_range * np.exp(-1j * qnm_par.omega(m) * (sample_times - start_time) )

    M = len(fit_modes) # number of modes to fit
    G = np.zeros((N, M), dtype=complex) # N by M matrix
    for k, m in enumerate(fit_modes):
        G[:, k] = qnm[m].data
    G_H = G.conj().T
    A = np.linalg.solve(G_H @ G, G_H @ h_target_slice.data)

    A_modes = {}
    sum_fit_qnm = 0
    for i, m in enumerate(fit_modes):
        # templates were rescaled at dynamical_range, so A[i] * template is physical strain
        qnm[m].data *= A[i]
        sum_fit_qnm += qnm[m]
        # reported amplitude: shift to t0 (not necessarily on the sample grid)
        A_modes[m] = A[i] * dynamical_range \
            * np.exp(-1j * qnm_par.omega(m) * (ringdown_start_time - start_time))

    return A_modes, sum_fit_qnm, h_target_slice

# earliest start time with a reliable overtone fit
FIT_TSTART_MIN = 0.002
# fit start time grid used when toffset < FIT_TSTART_MIN
FIT_TSTART_GRID = np.linspace(0.002, 0.0036665, 6)

QUADRATIC_MODES = ['220220', '220221','221221','220222','221222','222222']
_interpolation_cache = {}

for mode in QUADRATIC_MODES:
    filename = files(__package__).joinpath('data', f'spin_{mode}_interpolation.npy')
    data = np.load(filename, allow_pickle=True).item()
    spin_loaded = data['spin']
    complex_ratio_loaded = data[mode+'_complex_ratio']
    interp_loaded = interp1d(spin_loaded, complex_ratio_loaded, kind='cubic', fill_value="extrapolate")
    _interpolation_cache[mode] = interp_loaded

def load_interpolation_function(label):
    """
    load interpolation function for quadratic modes from cache
    """
    return _interpolation_cache[label]

def gen_nrsur_remove_qqnm(**kwds):
    '''Generate an NRSur7dq4 waveform with quadratic QNMs (QQNMs) subtracted.

    The QQNM amplitudes are not free parameters. The parent (2,2,n) overtone
    amplitudes A_22n are fitted from the NRSur7dq4 h22 by least-squares fitting,
    and each quadratic mode amplitude follows the theory prediction::

        A_QQNM = R(chi_f) * A_22n * A_22m * D / (G M_f / c^2)

    where R(chi_f) is the tabulated complex quadratic-to-linear amplitude
    ratio interpolated in final spin. Each QQNM is a damped sinusoid at the
    Kerr frequency f_n + f_m and damping time (1/tau_n + 1/tau_m)^-1. It is
    subtracted from the (4,4) multipole from toffset onwards, and the
    modified (4,4) mode is recombined with all other l <= 4 modes into the
    plus/cross polarizations.

    Parameters
    ----------
    mass1, mass2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z : float
        Component masses and spins passed to NRSur7dq4.
    distance, inclination, coa_phase, delta_t, f_lower, f_ref : float
        Extrinsic and generation parameters passed to NRSur7dq4.
    mode22 : str
        Space-separated parent (2,2,n) overtone labels used in the least-squares
        fit, may contain more overtones than the quadratic modes need
        (e.g. '220 221 222 223'), the extra overtones only stabilize the fit and
        do not enter the amplitude products.
    mode_quadratic : str
        Space-separated quadratic modes to subtract, e.g. "220220 220221".
        Must be a subset of QUADRATIC_MODES.
    toffset : float
        Ringdown start time t0 in seconds, e.g. 0.002. If toffset >= FIT_TSTART_MIN,
        the parent amplitudes come from a single fit at t0. Otherwise they are
        fitted on a grid of start times (2.0 - 3.67 ms), propagated back to
        t0, and each amplitude is a Gaussian random draw over the scatter
        of these fits -- the returned waveform is then stochastic.
    quadratic_tgr : float, optional
        Amplitude factor of the subtracted GR QQNM: 1 (default) subtracts
        the full GR prediction, 0 subtracts nothing.
    quadratic_tgr_phase : float, optional
        Extra phase (in radians) applied to the subtracted QQNM. Default 0.
    qqnm_deltaf, qqnm_deltatau : float, optional
        If either is given, a QQNM with frequency f*(1+qqnm_deltaf), damping
        time tau*(1+qqnm_deltatau) and the full theory amplitude is added
        back after the subtraction, so the waveform then contains a QQNM
        with shifted frequency/damping time. Default 0.
    seed : int, optional
        Seed for the Gaussian amplitude draw used when toffset < FIT_TSTART_MIN,
        making that branch reproducible.

    Returns
    -------
    hp, hc : pycbc TimeSeries (real)
        Plus and cross polarizations.
    '''
    # requested modes
    fit_mode22 = kwds['mode22'].split()
    quadratic_modes = kwds['mode_quadratic'].split()
    unknown = [m for m in quadratic_modes if m not in QUADRATIC_MODES]
    if unknown:
        raise ValueError(f"no amplitude-ratio table for {unknown}; "
                         f"available quadratic modes: {QUADRATIC_MODES}")
    needed_parents = {m[:3] for m in quadratic_modes} | {m[3:] for m in quadratic_modes}
    missing_parents = sorted(needed_parents - set(fit_mode22))
    if missing_parents:
        raise ValueError(f"quadratic modes require parent modes {missing_parents}, "
                         f"but mode22 = '{kwds['mode22']}'")
    needed_parents = sorted(needed_parents)

    # Generate the NRSur7dq4 waveform multipole modes
    from pycbc.waveform import get_td_waveform_modes
    hlm = get_td_waveform_modes(approximant='NRSur7dq4',
                                mass1=kwds['mass1'], mass2=kwds['mass2'],
                                spin1x=kwds['spin1x'], spin1y=kwds['spin1y'], spin1z=kwds['spin1z'],
                                spin2x=kwds['spin2x'], spin2y=kwds['spin2y'], spin2z=kwds['spin2z'],
                                distance = kwds['distance'],
                                delta_t=kwds['delta_t'],
                                f_lower=kwds['f_lower'],
                                f_ref=kwds['f_ref'],
                                mode_array=['22','21','20','33','32','31','30','44','43','42','41','40'])

    qnm_par = get_qnm_freqtau(fit_mode22 + quadratic_modes, **kwds)

    # parent amplitudes fitted from h22
    t0 = kwds['toffset']
    h22 = hlm[(2,2)][0] + 1j * hlm[(2,2)][1]
    rng = np.random.default_rng(kwds.get('seed'))
    if t0 >= FIT_TSTART_MIN:
        A_modes_22, _, _ = least_square_qnmfitting(fit_mode22, qnm_par, t0, h22)
    else:
        fit_A_modes_22 = {m: [] for m in needed_parents}
        for t_fit in FIT_TSTART_GRID:
            this_A, _, _ = least_square_qnmfitting(fit_mode22, qnm_par, t_fit, h22)
            # shift to t0
            for m in needed_parents:
                this_A[m] *= np.exp(-1j * qnm_par.omega(m) * (t0 - t_fit))
                fit_A_modes_22[m].append(this_A[m])

        A_modes_22 = {}
        for m in needed_parents:
            re = np.real(fit_A_modes_22[m])
            im = np.imag(fit_A_modes_22[m])
            A_modes_22[m] = rng.normal(re.mean(), re.std()) + 1j * rng.normal(im.mean(), im.std())

    # getting quadratic modes amplitude
    A_modes_quadratic = {}
    conversion = kwds['distance'] * 1e6 * lal.PC_SI / (qnm_par.final_mass * lal.MRSUN_SI)
    for mode in quadratic_modes:
        quad_linear_ratio_func = load_interpolation_function(mode)
        ratio = quad_linear_ratio_func(qnm_par.final_spin)
        mode1 = mode[:3]
        mode2 = mode[3:]
        A_modes_quadratic[mode] = A_modes_22[mode1] * A_modes_22[mode2] * ratio * conversion

    def par(key, default):
        # this function is in case a key exist but is None,
        # kwds.get(key, default) would return None in that case
        # but we want to return the default
        v = kwds.get(key)
        return default if v is None else v

    # TGR deviation parameters (GR values: amplitude 1, phase 0, deltaf = deltatau = 0)
    tgr_factor = par('quadratic_tgr', 1.0) * np.exp(-1j * par('quadratic_tgr_phase', 0.0))
    shift_freqtau = (kwds.get('qqnm_deltaf') is not None
                     or kwds.get('qqnm_deltatau') is not None)
    deltaf = par('qqnm_deltaf', 0.0)
    deltatau = par('qqnm_deltatau', 0.0)

    # subtract quadratic modes from (4,4) mode
    h44 = hlm[(4,4)][0] + 1j * hlm[(4,4)][1]
    start_idx = int(np.floor(float(t0 - h44.start_time) * h44.sample_rate))
    start_time = h44.sample_times[start_idx]
    end_time = h44.sample_times[-1]
    h44_slice = h44.time_slice(start_time, end_time)
    t = h44_slice.sample_times.numpy()

    for m in quadratic_modes:
        # QQNM waveform; the leading minus sign is the 'm/2 pi + pi' phase convention of NRSur7dq4
        qqnm_gr = -A_modes_quadratic[m] * np.exp(-1j * qnm_par.omega(m) * (t - t0) )
        # Subtract the quadratic mode contribution from (4,4) mode
        h44_slice.data -= tgr_factor * qqnm_gr

        if shift_freqtau:
            omega_shifted = 2 * np.pi * qnm_par.freq[m] * (1 + deltaf) - 1j / (qnm_par.tau[m] * (1 + deltatau))
            h44_slice.data += -A_modes_quadratic[m] * np.exp(-1j * omega_shifted * (t - t0) )

    # recombine the modified (4,4) mode into the polarizations
    h_pol = 0
    for l in range(2,5):
        for m in range(-1*l, l+1):
            if l==4 and abs(m)==4:
                continue # skip 44 mode for ringdown treatment later
            h_modes = hlm[(l,m)][0] + 1j * hlm[(l,m)][1]
            Y_lm = lal.SpinWeightedSphericalHarmonic(
                kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, l, m)
            h_pol += h_modes * Y_lm

    Y_44 = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, 4, 4)
    Y_4m4 = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, 4, -4)
    pol44 = h44_slice *  Y_44 + np.conj(h44_slice) * Y_4m4
    h_pol.data[start_idx:start_idx+len(pol44)] += pol44
    return h_pol.real(), -h_pol.imag()
