from pycbc.waveform import get_td_waveform_modes

import numpy as np
import postmerger
fit = postmerger.load_fit('3dq8_20M')

import lal

from pycbc.types import (TimeSeries, float64, zeros)

def gen_nrsurqnm(**kwds):
    hlm = get_td_waveform_modes(approximant='NRSur7dq4', mass1=kwds['mass1'], mass2=kwds['mass2'],
                                spin1x=kwds['spin1x'], spin1y=kwds['spin1y'], spin1z=kwds['spin1z'],
                                spin2x=kwds['spin2x'], spin2y=kwds['spin2y'], spin2z=kwds['spin2z'],
                                delta_t=kwds['delta_t'], f_lower=kwds['f_lower'],mode_array=['22','33','44'])
    h44 = hlm[(4,4)][0] + 1j * hlm[(4,4)][1]
    qnm_par = {}
    qnm_par['final_mass'] = postmerger.final_mass(kwds['mass1'], kwds['mass2'], kwds['spin1z'], kwds['spin2z'], aligned_spins=True)
    qnm_par['final_spin'] = postmerger.final_spin(kwds['mass1'], kwds['mass2'], kwds['spin1z'], kwds['spin2z'], aligned_spins=True)
    qnm_par['freq'] = {}
    qnm_par['tau'] = {}
    qnm_par['freq']['440'], qnm_par['tau']['440'] = postmerger.qnm_Kerr(qnm_par['final_mass'], qnm_par['final_spin'],(4,4,0), prograde=1, SI_units=True, qnm_method='interp')

    mode = ['440']
    t_final = max(qnm_par['tau'][m] for m in mode) * np.log(1000) # to decay by 1e-3 for the slowest decaying mode
    kmax = int(t_final / kwds['delta_t']) + 1

    trans_time = kwds['trans_time']
    start_idx = int(float(trans_time - h44.start_time) * h44.sample_rate)
    out = TimeSeries(zeros(kmax, dtype=float64), delta_t=kwds['delta_t'], epoch = h44.sample_times[start_idx])
    sample_times = out.sample_times.numpy()

	
    for m in mode:
        amp = 1e-22
        phi = 0
        omega = 2 * np.pi * qnm_par['freq'][m] - 1j / qnm_par['tau'][m]
        ringdown = amp * np.exp(-1j * phi) * np.exp(-1j * omega * (sample_times - trans_time) )
        out += ringdown
    
    h44_slice = h44.time_slice(float(out.start_time), float(out.end_time))
    optimal_phase = -np.angle(h44_slice.inner(out))
    optimal_amp = np.abs(h44_slice.inner(out)) / np.real(out.inner(out))
    out *= optimal_amp * np.exp(1j * optimal_phase)

    h = 0
    for l in range(2,4):
        for m in range(-1*l, l+1):
            h_modes = hlm[(l,m)][0] + 1j * hlm[(l,m)][1] 
            Y_lm = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, l, m)
            h += h_modes * Y_lm
	
    Y_44 = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, 4, 4)
    Y_4m4 = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], np.pi/2 - kwds['coa_phase'], -2, 4, -4)
    qnm44 = out *  Y_44 + np.conj(out) * Y_4m4
	
    h.data[start_idx:start_idx+len(qnm44)] += qnm44
	
    return h.real(), -h.imag()