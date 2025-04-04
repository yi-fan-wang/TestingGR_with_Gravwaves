import numpy
import lal
import sxs
import scipy.interpolate
import pycbc.types
import pycbc.waveform.utils

def gen_sxs_waveform(sxs_id, extrapolation_order=2, download=False, **kwds):
	wf = sxs.load(sxs_id + "/Lev/rhOverM",extrapolation_order=extrapolation_order,download=download)
	metadata = sxs.load(sxs_id + "/Lev/metadata.json", download=download)
	
	reference_time = metadata['reference_time']
	reference_index = wf.index_closest_to(reference_time)
	
	wf_sliced = wf[reference_index:]	
	time = wf_sliced.t - wf_sliced.max_norm_time()

	h = 0
	for l in range(2,9):
		for m in range(-1*l, l+1):
			h_modes = wf_sliced[:, wf_sliced.index(l, m)]
			Y_lm = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], kwds['coa_phase'], -2, l, m)
			h += h_modes * Y_lm
	
	# convert to physical units
	mtotal = kwds['mass1'] + kwds['mass2']
	distance = kwds['distance']

	time *= mtotal * lal.MTSUN_SI
	h *= mtotal * lal.MRSUN_SI / distance / 1e6 / lal.PC_SI

	interp_cache = scipy.interpolate.interp1d(time, h, kind='cubic')
	linear_time= numpy.arange(time[0], time[-1], kwds['delta_t'])
	h_interp = interp_cache(linear_time)

	hp = pycbc.types.TimeSeries(numpy.real(h_interp), delta_t = kwds['delta_t'], epoch = linear_time[0])
	hc = pycbc.types.TimeSeries(-numpy.imag(h_interp), delta_t = kwds['delta_t'], epoch = linear_time[0])

	window = 0.1
	hp_taper = pycbc.waveform.utils.td_taper(hp, hp.start_time, hp.start_time + window)
	hc_taper = pycbc.waveform.utils.td_taper(hc, hc.start_time, hc.start_time + window)
	
	return hp_taper, hc_taper