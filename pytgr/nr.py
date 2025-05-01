import numpy
import lal
import lalsimulation as lalsim
import sxs
import scipy.interpolate
import pycbc.types
import pycbc.waveform.utils
import h5py

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

def gen_lvcnr_waveform(data_file, **kwds):
	'''
	Generate a LVC NR waveform from a file path
	'''
	with h5py.File(data_file, "r") as f:
		m1 = f.attrs["mass1"] / (f.attrs["mass1"] + f.attrs["mass2"])
		m2 = f.attrs["mass2"] / (f.attrs["mass1"] + f.attrs["mass2"])

	fStart = kwds['f_lower']
	fRef = kwds['f_ref']
	deltaT = kwds['delta_t']

	mtotal = kwds['mass1'] + kwds['mass2']
	m1 *= (mtotal * lal.MSUN_SI)
	m2 *= (mtotal * lal.MSUN_SI)
	s1x, s1y, s1z, s2x, s2y, s2z = lalsim.SimInspiralNRWaveformGetSpinsFromHDF5File(
					fStart, mtotal, data_file)
	inclination = kwds['inclination']
	distance = kwds['distance'] * 1e6 * lal.PC_SI
	phiRef = kwds['coa_phase']
	
	params = lal.CreateDict()
	lalsim.SimInspiralWaveformParamsInsertNumRelData(params, data_file)
	
	hp, hc = lalsim.SimInspiralChooseTDWaveform(m1, m2,
                                            	s1x, s1y, s1z, s2x, s2y, s2z,
                                            	distance, inclination,
                                            	phiRef, 0,
                                            	0, 0, deltaT,
                                            	fStart, fRef,
                                            	params, approximant=lalsim.NR_hdf5)

	hp = pycbc.types.TimeSeries(hp.data.data, delta_t=deltaT, epoch=hp.epoch)
	hc = pycbc.types.TimeSeries(hc.data.data, delta_t=deltaT, epoch=hc.epoch)
	
	return hp, hc