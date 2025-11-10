import pkg_resources
import numpy as np
import lal
from scipy.interpolate import interp1d

from pycbc.conversions import get_final_from_initial, get_lm_f0tau
from pycbc.waveform import get_td_waveform_modes
from pycbc.types import (TimeSeries, complex64, zeros)

def gen_nrsurqnm(**kwds):
    hlm = get_td_waveform_modes(approximant='NRSur7dq4', mass1=kwds['mass1'], mass2=kwds['mass2'],
                            spin1x=kwds['spin1x'], spin1y=kwds['spin1y'], spin1z=kwds['spin1z'],
                            spin2x=kwds['spin2x'], spin2y=kwds['spin2y'], spin2z=kwds['spin2z'],
                            distance = kwds['distance'],
                            delta_t=kwds['delta_t'], f_lower=kwds['f_lower'],mode_array=['22','21','20','33','32','31','30','44'])
    
    h = 0
    for l in range(2,4):
        for m in range(-1*l, l+1):
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
    start_idx = int(float(ringdown_start_time - h44.start_time) * h44.sample_rate)
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


def qnm_decomposition(qnm_modes, qnm_par, ringdown_start_time, h_target, **kwds):
    '''
    Decompose the target waveform h_target into the given QNM modes starting from ringdown_start_time
    '''
    t_duration = max(qnm_par['tau'][m] for m in qnm_modes) * np.log(1000)
    start_idx = int(float(ringdown_start_time - h_target.start_time) * h_target.sample_rate)
    end_idx = int(float(ringdown_start_time + t_duration - h_target.start_time) * h_target.sample_rate)

    start_time = h_target.sample_times[start_idx]
    end_time = h_target.sample_times[end_idx]
    h_target_slice = h_target.time_slice(start_time, end_time)
    N = len(h_target_slice)

    qnm = {}
    dynamical_range = 1e-22
    for m in qnm_modes:
        qnm[m] = TimeSeries(zeros(N, dtype=complex64), delta_t=h_target.delta_t, epoch = start_time)
        sample_times = qnm[m].sample_times.numpy()
        omega = 2 * np.pi * qnm_par['freq'][m] - 1j / qnm_par['tau'][m]
        qnm[m].data = dynamical_range * np.exp(-1j * omega * (sample_times - start_time) )

    M = len(qnm_modes)
    G = np.zeros((N, M), dtype=complex)
    for k, m in enumerate(qnm_modes):
        G[:, k] = qnm[m].data
    G_H = G.conj().T
    A = np.linalg.solve(G_H @ G, G_H @ h_target_slice.data)

    A_modes = {}
    for i, m in enumerate(qnm_modes):
        A_modes[m] = A[i]

    allqnm = 0
    for m in qnm_modes:
        delta_param = 'delta_' + m
        if delta_param in kwds and kwds[delta_param] is not None:
            #print("Applying deviation", delta_param, kwds[delta_param])
            A_modes[m] *= (1 + kwds[delta_param])
        qnm[m].data *= A_modes[m]
        allqnm += qnm[m]

    for m in qnm_modes:
        A_modes[m] = A_modes[m] * dynamical_range

    return A_modes, allqnm, h_target_slice

def get_qnmpar(qnm_modes, **kwds):
    qnm_par = {}
    qnm_par['final_mass'], qnm_par['final_spin'] = \
        get_final_from_initial(kwds['mass1'], kwds['mass2'],
                               kwds['spin1x'], kwds['spin1y'], kwds['spin1z'],
                               kwds['spin2x'], kwds['spin2y'], kwds['spin2z'],
                               approximant='NRSur7dq4',
                               f_ref=kwds['f_ref'])
    qnm_par['freq'] = {}
    qnm_par['tau'] = {}
    for mode in qnm_modes:
        if len(mode) == 3:
            l = int(mode[0])
            m = int(mode[1])
            n = int(mode[2])
            qnm_par['freq'][mode], qnm_par['tau'][mode] = \
                get_lm_f0tau(qnm_par['final_mass'], qnm_par['final_spin'],
                            l,m,n)
        elif len(mode) == 6:
            l1 = int(mode[0])
            m1 = int(mode[1])
            n1 = int(mode[2])
            l2 = int(mode[3])
            m2 = int(mode[4])
            n2 = int(mode[5])
            f1, tau1 = get_lm_f0tau(qnm_par['final_mass'],qnm_par['final_spin'],l1,m1,n1)
            f2, tau2 = get_lm_f0tau(qnm_par['final_mass'],qnm_par['final_spin'],l2,m2,n2)
            qnm_par['freq'][mode] = f1 + f2
            qnm_par['tau'][mode] = 1 / (1/tau1 + 1/tau2)
        else:
            raise ValueError("Invalid mode format in rindown_mode")
    return qnm_par

QUADRATIC_MODES = ['220220', '220221','221221']
_interpolation_cache = {}

for mode in QUADRATIC_MODES:
    filename = pkg_resources.resource_filename(__name__, f'data/spin_{mode}_interpolation.npy')
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

def gen_nrsur_linearqnm(**kwds):
    '''
    Generate a NRSur7dq4 waveform with linear QNM in 44 mode by removing the quadratic modes
    '''
    hlm = get_td_waveform_modes(approximant='NRSur7dq4', 
                                mass1=kwds['mass1'], mass2=kwds['mass2'],
                                spin1x=kwds['spin1x'], spin1y=kwds['spin1y'], spin1z=kwds['spin1z'],
                                spin2x=kwds['spin2x'], spin2y=kwds['spin2y'], spin2z=kwds['spin2z'],
                                distance = kwds['distance'],
                                delta_t=kwds['delta_t'],
                                f_lower=kwds['f_lower'],
                                f_ref=kwds['f_ref'],
                                mode_array=['22','21','20','33','32','31','30','44'])
    
    h = 0
    for l in range(2,4):
        for m in range(-1*l, l+1):
            h_modes = hlm[(l,m)][0] + 1j * hlm[(l,m)][1] 
            Y_lm = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, l, m)
            h += h_modes * Y_lm

    h22 = hlm[(2,2)][0] + 1j * hlm[(2,2)][1]
    h44 = hlm[(4,4)][0] + 1j * hlm[(4,4)][1]
    qnm_par = get_qnmpar(['220','221','222'] + QUADRATIC_MODES, **kwds)
    
    # construct QNM from (2,2) mode
    qnm_start_time = kwds['toffset']
    if qnm_start_time >= 0.002:
        A_modes_22, _, _ = qnm_decomposition(['220','221','222'], qnm_par, kwds['toffset'], h22, **kwds)
    else:
        A_modes_22, _, _ = qnm_decomposition(['220','221','222'], qnm_par, 0.002, h22, **kwds)
        for m in A_modes_22:
            omega = 2 * np.pi * qnm_par['freq'][m] - 1j / qnm_par['tau'][m]
            A_modes_22[m] *= np.exp(-1j * omega * (qnm_start_time-0.002))

    # removing quadratic modes
    A_modes_quadratic = {}
    for mode in QUADRATIC_MODES:
        load_interp = load_interpolation_function(mode)
        chi_f = qnm_par['final_spin']
        ratio = load_interp(chi_f)
        mode1 = mode[:3]
        mode2 = mode[3:]
        conversion_factor = kwds['distance'] * 1e6 * lal.PC_SI / qnm_par['final_mass'] / lal.MRSUN_SI
        A_modes_quadratic[mode] = A_modes_22[mode1] * A_modes_22[mode2] * ratio * conversion_factor

    # subtract quadratic modes from (4,4) mode
    start_idx = int(float(kwds['toffset'] - h44.start_time) * h44.sample_rate)
    start_time = h44.sample_times[start_idx]
    end_time = h44.sample_times[-1]
    h44_slice = h44.time_slice(start_time, end_time)
    N = len(h44_slice)

    qnm = {}
    if 'quadratic_tgr' in kwds and kwds['quadratic_tgr'] is not None:
        tgr = kwds['quadratic_tgr']
    else:
        tgr = 1.0
    for m in QUADRATIC_MODES:
        qnm[m] = TimeSeries(zeros(N, dtype=complex64), delta_t=h44.delta_t, epoch = start_time)
        sample_times = qnm[m].sample_times.numpy()
        omega = 2 * np.pi * qnm_par['freq'][m] - 1j / qnm_par['tau'][m]
        qnm[m].data = A_modes_quadratic[m] * np.exp(-1j * omega * (sample_times - start_time) )
        h44_slice -= tgr * qnm[m]

    Y_44 = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, 4, 4)
    Y_4m4 = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, 4, -4)
    qnm44 = h44_slice *  Y_44 + np.conj(h44_slice) * Y_4m4
    h.data[start_idx:start_idx+len(qnm44)] += qnm44
    return h.real(), -h.imag()
