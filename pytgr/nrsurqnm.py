import numpy as np
import lal

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
        qnm[m].data *= A[i]
        allqnm += qnm[m]
    #print(allqnm.data)
    Y_44 = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, 4, 4)
    Y_4m4 = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, 4, -4)
    qnm44 = allqnm *  Y_44 + np.conj(allqnm) * Y_4m4
	
    h.data[start_idx:start_idx+len(qnm44)] += qnm44
    #print(qnm44.data)
    return h.real(), -h.imag()