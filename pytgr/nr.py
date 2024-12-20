import sxs
import h5py
import lal
import numpy
import glob
import scipy.interpolate
import pycbc.types
import pycbc.waveform.utils

def gen_sxs_waveform(sxs_tag, sxs_dir, **kwds):

	#sxs_tag = kwds['sxs_tag']
	sxs.load(sxs_tag+"/Lev/rhOverM")

	filedir = glob.glob(sxs_dir+'/'+sxs_tag+'v*/*/rhOverM_Asymptotic_GeometricUnits_CoM.h5')
	nr_file = h5py.File(filedir[0],'r')['OutermostExtraction.dir']
		
	h=0
	time = None
	for l in range(2,9):
	    for m in range(-1*l, l+1):
	        h_modes = nr_file["Y_l"+str(l)+"_m"+str(m)+".dat"]
	        hreal = h_modes[:,1]
	        himag = h_modes[:,2]
	        Y_lm = lal.SpinWeightedSphericalHarmonic(kwds['inclination'], kwds['coa_phase'], -2, l, m)
	        h += (hreal + 1j*himag) * Y_lm
	        if time is None:
	            time = h_modes[:,0]
	        
	interp_cache = scipy.interpolate.interp1d(time, h, kind='cubic')
	linear_time= numpy.arange(time[0], time[-1], kwds['delta_t'])
	h_interp = interp_cache(linear_time)

	tmerge = linear_time[numpy.argmax(numpy.abs(h_interp))]
	align_time = linear_time - tmerge 

	# convert NR unit to physical unit
	mtotal = kwds['mass1'] + kwds['mass2']
	distance = kwds['distance']

	phy_time = align_time * mtotal * lal.MTSUN_SI
	phy_h = h_interp * lal.MRSUN_SI / distance / 1e6 / lal.PC_SI

	hp = pycbc.types.TimeSeries(numpy.real(phy_h), delta_t = phy_time[1]-phy_time[0], epoch = phy_time[0])
	hc = pycbc.types.TimeSeries(numpy.imag(phy_h), delta_t = phy_time[1]-phy_time[0], epoch = phy_time[0])

	window = 0.1
	hp_win = pycbc.waveform.utils.td_taper(hp, hp.start_time, hp.start_time + window)
	hc_win = pycbc.waveform.utils.td_taper(hc, hc.start_time, hc.start_time + window)
	
	return hp_win, hc_win