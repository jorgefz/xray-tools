import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import argparse
import os
import sys
import math
from numpy.lib.recfunctions import append_fields

from .utils import Utils
from .resources import Resources

"""
	----- X-ray Tools -----

	For analysing stellar data in the X-ray regime.
	
	Author: Jorge Fernandez, University of Warwick, Coventry, UK
	Date: 17th December 2020

	Functions:
		rotlaw
			Calculates the star's Rossby number based on its rotation period,
			and plots it along with a luminosity-ratio/rossby-number model by Wright et al 2011.

		agerel
			Fits the star to a luminosity-ratio/colour-index/age model by Jackson et al 2012.

		xraydiv
			Calculates the energy boundary between soft and hard X-rays.

		lcratio
			Plots the ratio between soft and hard light curves.

		cumlc
			Plots a cumulative light curve

		euv
			Estimates the EUV flux based on the X-ray flux in a given range.

"""


PyDir = os.path.dirname(os.path.realpath(sys.argv[0]))

	


def specplot(spectrum, rmf=None):
	"""
	Plots input spectrum.
	Inputs:
	- spectrum:	(string) Filename of the spectrum.
	- rmf:		(string) Filename of the rmf file.
	"""
	pass

def lcplot(lc):
	"""
	Plots input light-curve
	"""
	pass


def rotation_tracks(age_myr=700, lum_x=1.8E28, colors='red'):
	"""
	Stellar rotation tracks by Tu et al (2015).
	Input:
		age_myr:	stellar ages in Myr
		lum_x:		stellar x-ray luminosities in erg/s
		colors:		colors in which to plot the data points (optional)
	"""
	# Making sure input data is a list
	if not Utils.indexable(age_myr):
		age_myr = [age_myr]

	if not Utils.indexable(lum_x):
		lum_x = [lum_x]

	if isinstance(colors, str):
		colors = [colors]*len(age_myr)


	ages = np.linspace(start=1, stop=4600, num=1000) # ages in Myr

	lbol_sun = 6.58E32/0.171
	lx_sat = 10**(-3.13) * lbol_sun
	rot_rate_sun = 2.9E-6 # Rotation rate of Sun in rad/s
	t_sat = [5.7, 23, 226]

	lx_10th = [ 2.0E31 * t**(-1.12) if t>t_sat[0] else lx_sat for t in ages]
	lx_50th = [ 2.6E32 * t**(-1.42) if t>t_sat[1] else lx_sat for t in ages]
	lx_90th = [ 2.3E36 * t**(-2.5) if t>t_sat[2] else lx_sat for t in ages]

	# Plotting Lx vs Age
	plt.plot(ages, lx_10th, color='red')
	plt.plot(ages, lx_50th, color='green')
	plt.plot(ages, lx_90th, color='blue')
	plt.xscale("log")
	plt.yscale("log")
	plt.xlabel('$\log_{10}$ Age (myr)')
	plt.ylabel('$\log_{10}\,L_X$ (erg/s)')
	plt.title("Tu et al (2015) stellar rotation tracks")
	for i in range(len(age_myr)):
			plt.plot(age_myr[i], lum_x[i], marker='.', color=colors[i], markersize=3)
	plt.show()



def rotlaw(lumratio=0, rossby=0, prot=0, vkcolor=None, plot_dataset=True):
	"""
	This function plots a star's rotation period and X-ray luminosity
	versus bolumetric luminosity ratio along with a model by Wright et al 2011
	(https://ui.adsabs.harvard.edu/abs/2011ApJ...743...48W/abstract)

	Args:
		lumratio:	(float) X-ray / bolumetric luminositiy ratio.
		prot:   	(float) Star's rotation period in days.
		rossby: 	(float) Star's Rossby number (rotation period / turnover time)
		vkcolor:     	(float) Star's V-Ks colour index.
		plot_dataset:	(bool) Plot dataset from Wight et al 2011 alongside model fit.

	The luminosity ratio is always required.
	If the Rossby number is defined, the plot can be generated straight away.
	Otherwise, you will need to input both the star's rotation period in days as well
	as its V - Ks colour index.
	Note: the X-ray luminosity is expectd to be in the ROSAT energy range (0.1 to 2.4 keV).
	They used the PIMMS tool to convert all their fluxes to ROSAT.
	"""
	# Input checking
	includes_rossby = True
	# The luminosity ratio is always required
	if(lumratio <= 0):
		print("[ROTLAW] Error: Lx/Lbol is invalid or undefined.")
		return
	# Decide to use either Rossby number or estimate it.
	if(rossby <= 0):
		includes_rossby = False
		if(prot <= 0 or vkcolor is None):
			print( ("[ROTLAW] Error: undefined parameters.\n"
				"Either Rossby number or both rotation period and V-Ks color required.") )
			return
		else:
			print("[ROTLAW] No Rossby number defined. Value will be estimated.")

	# Estimating turnover time and Rossby number, if undefined
	if(includes_rossby == False):
		rossby, turnover_time = rossby_number(vkcolor, prot)
		print(f" Estimated turnover time: {turnover_time:.3f} days")
		print(f" Estimated Rossby number: {rossby:.3f}")

	# Model parameters
	powerlaw = -2.7
	pwl_err = 0.16
	ro_sat = 0.13	# Saturation Rossby number
	ro_sat_err = 0.02
	rx_sat = 10**(-3.13)
	rx_sat_err = 0.08	# Saturation luminosity ratio (log10)

	# Rossby number points for plotting
	ross = np.linspace(0.001, 10, 1000)

	# Building Rossby number model
	const = rx_sat / (ro_sat**powerlaw)
	model = np.array([const*r**powerlaw if r>ro_sat else rx_sat for r in ross])

	rx_sat_uerr = rx_sat*10**rx_sat_err
	rx_sat_lerr = rx_sat/10**rx_sat_err
	ro_sat_uerr = ro_sat+ro_sat_err
	ro_sat_lerr = ro_sat-ro_sat_err

	const_uerr = (rx_sat_uerr) / (ro_sat_uerr**(powerlaw+pwl_err))
	const_lerr = (rx_sat_lerr) / (ro_sat_lerr**(powerlaw-pwl_err))
	model_uerr = np.array([const_uerr*r**(powerlaw+pwl_err) if r>ro_sat_uerr else rx_sat_uerr for r in ross])
	model_lerr = np.array([const_lerr*r**(powerlaw-pwl_err) if r>ro_sat_lerr else rx_sat_lerr for r in ross])

	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.set_xscale('log')
	# Format the X axis to be in log scale but keeping the original values (e.g. 0.01, 0.1, 1, 10, etc)
	ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
	
	# Plotting Wright dataset
	if(plot_dataset):
		dataset = Resources.wright_dset()
		dataset_rossby = [rossby_number(dataset['V-K'][i], dataset['Prot'][i])[0] for i in range(len(dataset['Prot']))]
		dataset_rx = [10**d for d in dataset['Lx/bol']]
		dataset_cluster = dataset['Cluster']
		print(len(dataset_cluster))
		
		for i in range(len(dataset_rx)):
			if(str(dataset_cluster[i]).strip() == 'Hyades'):
				plt.plot(dataset_rossby[i], dataset_rx[i], color='k', marker='*', markersize=2.5, linestyle='None')
			else:
				plt.plot(dataset_rossby[i], dataset_rx[i], color='0.5', marker='.', markersize=1.5, linestyle='None')

		#plt.plot(dataset_rossby[0], dataset_rx[0], color='black', marker='+', markersize=7.5)

	# Plotting model
	plt.plot(ross, model, 'b-', linewidth=1.2)
	plt.plot(ross, model_uerr, 'b--', linewidth=1.2)
	plt.plot(ross, model_lerr, 'b--', linewidth=1.2)
	plt.xlabel("Rossby number")
	plt.title("Rotation Law by Wright et al (2011)\n$R_x = (%.2E) Ro^{%.2f\\pm %.2f}$" % (const, powerlaw, pwl_err) )

	# Plot our star
	plt.plot(rossby, lumratio, 'r*')

	# Plot decorations
	plt.ylabel("$L_x / L_{bol}$")
	plt.yscale("log")
	plt.show()



def agerel(lumratio=None, bvcolor=None, age_myr=None, plot_dataset=True):
	"""
	This function plots a star's colour index (B-V) and age
	along with an age relation by Jackson et al (2012):
	https://arxiv.org/abs/1111.0031

	Args:
		lumratio:    	(float) X-ray / bolometric luminosity ratio of the star.
		bvcolor:     	(float) Star's B-V colour index.
		age_myr:      	(float) Star's age in Myr.
		plot_dataser:	(bool) Plots dataset of stars from the paper alongside your star.

	Both parameters are required, as well as the predefined model parameters
	in the file 'jackson_model.txt'.
	Note: the X-ray luminosity is expectd to be in the ROSAT energy range (0.1 to 2.4 keV).
	They used the PIMMS tool to convert all their fluxes to ROSAT.
	"""

	if(lumratio is None or bvcolor is None or age_myr is None):
		print("[AGEREL] Error: Star colour index or age undefined.")
		return

	model_data = Resources.jackson_model()
	"""
	Data fields:
		icolor: initial colou index, model boundary 
		fcolor: final colour idnex, model boundary
		satlum: saturation luminosity ratio. After this value, the ratio
				goes from being constant to following a powerlaw.
		satage: saturation age, after which the luminosity ratio follows a power law.
		powerlaw: power law index
		errpowerlaw: error in the power law index 
	"""
	
	# Find which colour index boundary the star applies to.
	ind = -1
	for i,m in enumerate(model_data):
		if(m['icolor'] < bvcolor < m['fcolor']):
			ind = i
			break
	if(ind < 0):
		min_color = min(model_data['icolor'])
		max_color = max(model_data['fcolor'])
		print(f"[AGEREL] Error: star's colour index outside model bounds ({min_color:.2f} to {max_color:.2f})")
		return
	
	# Fit and plot
	satlum = 10**model_data['satlum'][ind]
	satage = 10**model_data['satage'][ind]
	powerlaw = model_data['powerlaw'][ind]

	ages = np.linspace(10**6.5, 10**10, 1000)
	const = satlum / (satage**(-powerlaw))
	model = np.array([const*a**(-powerlaw) if a>satage else satlum for a in ages])

	if(plot_dataset):
		dataset = Resources.jackson_dset()
		dataset_age = [d for i,d in enumerate(dataset['log(Age)']) if model_data['icolor'][ind] <= dataset['(B-V)0'][i] <= model_data['fcolor'][ind]]
		dataset_rx = [d for i,d in enumerate(dataset['log(Lx/Lb)']) if model_data['icolor'][ind] <= dataset['(B-V)0'][i] <= model_data['fcolor'][ind]]
		plt.plot(dataset_age, dataset_rx, color='0.5', marker='.', markersize=2.5, linestyle='None')

	plt.plot(np.log10(ages), np.log10(model), 'b-')
	plt.plot(math.log10(age_myr)+6, np.log10(lumratio), 'r*', markersize = 7.5)

	plt.xlabel("Age log10(yr)")
	plt.ylabel("$log10(L_x / L_{bol})$")
	plt.title("Age relation by Jackson et al (2012)\n%.2f < (B-V) < %.2f" % (model_data['icolor'][ind], model_data['fcolor'][ind]))
	plt.show()
	plt.clf()




def xraydiv(file=None, rmf=None, erange=(0.15,2.4), bkg = None, verbose=True):
	"""
	This function extracts spectrum data from a star and calculates the energy
	boundary between soft and hard X-rays (energy at which 50% of the counts are
	on the left, and the remaining 50% on the right).

	Args:
		file:	(string) file where the spectrum data is stored.
		rmf:	(string) RMF response file of the star.
		erange:	(float,float) (optional) Range of energies to include in plot (keV)
		bkg:	(string) file with background counts. If defined, it will be subtarcted from the main spectrum.
		verbose: (bool) (optional) Prints information and plots spectrum if true.

	Returns:
		[0] 	(float) Energy boundary between soft and hard X-rays in keV.
		[0] 	(None) If necessary input arguments are missing.
	"""
	if(not file or not rmf):
		print("[XRAYDIV] Error: spectrum or rmf filenames undefined")
		return None

	if not os.path.exists(file):
		print(f"[XRAYDIV] Error: file '{file}' not found.")
		return None
	if not os.path.exists(rmf):
		print(f"[XRAYDIV] Error: file '{rmf}' not found.")
		return None
	if(not Utils.indexable(erange) or len(erange)!=2 or erange[0]>=erange[1]):
		print("[XRAYDIV] Error: energy range must be an array of two floats (min, max)")
		return None

	data = None 	# Spectrum data (counts on each bin, as well as grouping, quality, etc)
	scaling = 1		# Factor to which to scale counts based on the area of the aperture
	channel_data = None # Energy range for each bin in keV
	exptime = 0		# Exposure time in seconds of the spectrum

	with fits.open(file) as image:
		data = np.copy(image['SPECTRUM'].data)
		scaling = image['SPECTRUM'].header['BACKSCAL']
		# Converting counts to float lets us divide by the scaling factor and get decimal number
		data = data.astype(dtype=[('CHANNEL', '>i2'), ('COUNTS', '<f8'), ('GROUPING', '>i2'), ('QUALITY', '>i2')])
		# Scaling data by aperture area
		#data['COUNTS'] = data['COUNTS']/float(scaling)
		
		# Calculating total exposure time in seconds
		import re
		"""This splits a date string from the header (e.g. "2018-12-5T16:28:33")
		into a list of ints (2018, 12, 5, 16, 28, 33) using a regular expression
		on the separators: dash, space, T, and semicolon.
		"""
		expstart = filter(None, re.split("[T \-:]+", image[0].header['EXPSTART']))
		expstop = filter(None, re.split("[T \-:]+", image[0].header['EXPSTOP']))
		from datetime import datetime
		exptime = (datetime(*map(int,expstop)) - datetime(*map(int,expstart))).total_seconds()

	with fits.open(rmf) as rmf_data:
		channel_data = np.copy(rmf_data['EBOUNDS'].data)

	# Subtract background counts (if file provided)
	if(bkg):
		bkg_data = None
		bkg_scaling = 1
		with fits.open(bkg) as image:
			bkg_data = np.copy(image['SPECTRUM'].data)
			bkg_scaling = image['SPECTRUM'].header['BACKSCAL']
		# Scaling background counts
		bkg_data = bkg_data.astype(dtype=[('CHANNEL', '>i2'), ('COUNTS', '<f8')])
		bkg_data['COUNTS'] = bkg_data['COUNTS']*float(scaling)/float(bkg_scaling)
		# Subtracting background
		data['COUNTS'] = data['COUNTS'] - bkg_data['COUNTS']
	
	data['COUNTS'][data['COUNTS'] < 0] = 0

	# Get average energy of each bin
	energies = (channel_data['E_MAX'] + channel_data['E_MIN'])/2.0

	# Crop spectra to input energy range
	count_start = [i for i,e in enumerate(energies) if e > erange[0]][0]
	count_end = [item[0] for item in enumerate(energies) if item[1] < erange[1]][-1]
	energies = energies[count_start:count_end]
	data = data[count_start:count_end]

	total_counts = data['COUNTS'].sum()
	if(verbose):
		print(f"[XRAYDIV] Total scaled counts = {total_counts:.5g}")
		print(f"[XRAYDIV] Exposure time: {exptime} s")

	# Calculating middle point between soft and hard X-rays.
	cmcounts = 0
	middle_energy = 0
	for i,c in enumerate(data['COUNTS']):
		cmcounts += c
		if(cmcounts > total_counts/2.0):
			middle_energy = energies[i]
			break
	if verbose:
		print(f"[XRAYDIV] Soft/hard X-rays boundary is at {middle_energy:.3g} keV")

	if(verbose == False):
		return middle_energy

	counts_binned = []
	energies_binned = []
	energy_per_bin = []

	# Get indices of where bins begin, ignoring bad channels
	bin_ind = [ i for i,g in enumerate(data['GROUPING']) if(data['QUALITY'][i] != 1 and g == 1) ]

	for i in range(len(bin_ind)-1):
		total = data['COUNTS'][bin_ind[i]:bin_ind[i+1]].sum()
		counts_binned.append( total )
		energies_binned.append( energies[bin_ind[i]:bin_ind[i+1]].mean() )
		energy_per_bin.append( np.abs(energies[bin_ind[i]] - energies[bin_ind[i+1]]) )

	# Normalizing counts
	counts_binned = [d/exptime/energy_per_bin[i] for i,d in enumerate(counts_binned)]

	# Plot resulting spectrum
	plt.step(energies_binned, counts_binned, color='k', linewidth=1)
	plt.title(f"X-Ray Spectrum\n {file}")
	plt.xlabel("Energy (keV)")
	plt.ylabel("Normalized counts s$^{-1}$ keV$^{-1}$")
	plt.axvline(x=middle_energy, color='red', linestyle='-', linewidth=1)
	plt.text(middle_energy*1.01,max(counts_binned)*0.9,f'{middle_energy:.2f} keV',rotation=0)
	plt.show()

	return middle_energy



def lc_hardness(fname, emin, ebound, emax, rebin=0, legend=True):
	"""
	Plots the ratio between soft and hard X-ray lightcurves.
	Input:
		fname:	(str) Name used to generate filenames for hard and soft source and background lightcurves.
		emin:	(float) Min energy of soft X-rays.
		ebound:	(float) Energy boundaryb between soft and hard X-rays.
		emax:	(float) Max energy of hard X-rays.
		rebin:	(float) Seconds per bin to which to rebin the data
		legend: (bool) Plots a legend on the cumulative lightcurve
	"""
	
	# Input checking
	def file_exists(fname):
		if not os.path.exists(fname):
			print(f"[lc-hardness] Error: {fname} not found")
			return None
		return True
	
	fsoft = f"lightcurves/source_{fname}_softxray.lc"
	fhard = f"lightcurves/source_{fname}_hardxray.lc"
	bsoft = f"lightcurves/bkg_{fname}_softxray.lc"
	bhard = f"lightcurves/bkg_{fname}_hardxray.lc"

	if(not file_exists(fsoft)): return	
	if(not file_exists(fhard)): return
	if(not file_exists(bhard)): return
	if(not file_exists(bsoft)): return
	
	# Generating legend tags
	erange = [emin, ebound, emax]
	softname = f"Soft X-rays ({erange[0]:.2g}-{erange[1]:.2g} keV)"
	hardname = f"Hard X-rays ({erange[1]:.2g}-{erange[2]:.2g} keV)"
	
	# Extracting data from fits files
	hdata = sdata = hbkg = sbkg = None
	scaling = bscaling = 0 # Scaling factors for source and background lightcurves
	
	with fits.open(fhard) as image:
		hdata = np.copy(image['RATE'].data)
		scaling = _header_extract_aperture(image['RATE'].header)
	with fits.open(bhard) as bkg:
		hbkg = np.copy(bkg['RATE'].data) 
		bscaling = _header_extract_aperture(bkg['RATE'].header)	
	with fits.open(fsoft) as image:
		sdata = np.copy(image['RATE'].data)
	with fits.open(bsoft) as bkg:
		sbkg = np.copy(bkg['RATE'].data) 

	# Removing NANs
	old_size = len(hdata)
	hdata,sdata,hbkg,sbkg = remove_all_nans(hdata, sdata, hbkg, sbkg)
	print(f"[lc-hardness] {old_size-len(hdata)} NAN values found and removed")
	
	# Scaling apertures
	hdata['RATE'] = hdata['RATE'] / scaling
	hdata['ERROR'] = hdata['ERROR'] / scaling
	sdata['RATE'] = sdata['RATE'] / scaling
	sdata['ERROR'] = sdata['ERROR'] / scaling

	hbkg['RATE'] = hbkg['RATE'] / bscaling
	sbkg['RATE'] = sbkg['RATE'] / bscaling	
	hbkg['ERROR'] = hbkg['ERROR'] / bscaling
	sbkg['ERROR'] = sbkg['ERROR'] / bscaling

	print(f"[lc-hardness] Scaled background count rates by a factor of {scaling/bscaling:.2g}")
	
	# Normalize times and convert to kiloseconds
	hdata['TIME'] = (hdata['TIME'] - hdata['TIME'][0])/1000.0
	sdata['TIME'] = hdata['TIME']

	# Save source and background light-curves
	totdata = np.copy(hdata) # Total uncorrected lightcurve
	totdata['RATE'] = hdata['RATE'] + sdata['RATE']
	totdata['ERROR'] = np.sqrt(hdata['ERROR']**2 + sdata['ERROR']**2)

	bkgdata = np.copy(hbkg) # Background lightcurve
	bkgdata['RATE'] = hbkg['RATE'] + sbkg['RATE']
	bkgdata['ERROR'] = np.sqrt(hbkg['ERROR']**2 + sbkg['ERROR']**2)
	bkgdata['TIME'] = hdata['TIME']

	# Removing background
	hdata['RATE'] = hdata['RATE'] - hbkg['RATE']
	hdata['ERROR'] = np.sqrt( hdata['ERROR']**2 + hbkg['ERROR']**2 )
	sdata['RATE'] = sdata['RATE'] - sbkg['RATE']
	sdata['ERROR'] = np.sqrt( sdata['ERROR']**2 + sbkg['ERROR']**2 )	
	
	# Generating time errors between datapoints
	terrors = timebin_errors(hdata['TIME'])
	hdata = append_fields(hdata, ['TERROR_L', 'TERROR_R'], terrors)
	sdata = append_fields(sdata, ['TERROR_L','TERROR_R'], terrors)

	# Generate cumulative lightcurves
	hdata = append_fields(hdata, 'CUMUL', np.cumsum(hdata['RATE']))	
	sdata = append_fields(sdata, 'CUMUL', np.cumsum(sdata['RATE']))

	# Saving originals before binning
	orig_hdata = np.copy(hdata)
	orig_hbkg = np.copy(hbkg)

	orig_sdata = np.copy(sdata)
	orig_sbkg = np.copy(sbkg)

	# Binning
	binning = np.average(hdata['TERROR_L']+hdata['TERROR_R'])*1000
	print(f"[lc-hardness] Data is binned to {binning} seconds per bin")
	
	factor = int(rebin//binning)
	if(rebin == 0):
		pass
	elif(factor < 1 or (rebin%binning)!=0):
		print(f"[lc-hardness] Warning: you can only rebin to multiples of the original binning and to lower resolutions")
	else:
		print(f"[lc-hardness] Rebinning to {rebin} second bins...")
		sz = len(hdata)
		# Functions to calculate the new combined bins
		errorsum = lambda x: np.sqrt(sum([d**2 for d in x ]))/len(x) 
		elsesum = lambda x: sum(x)/len(x)
		# Applying error propagation in rebinning
		err_hdata = rebin_array(hdata['ERROR'], factor, func=errorsum)
		err_sdata = rebin_array(sdata['ERROR'], factor, func=errorsum)
		# Rebinning all other fields
		hdata = rebin_data(hdata, factor, func=elsesum)
		sdata = rebin_data(sdata, factor, func=elsesum)
		hdata['ERROR'] = err_hdata
		sdata['ERROR'] = err_sdata
		# Recalculating time errors
		hdata['TERROR_L'], hdata['TERROR_R'] = timebin_errors(hdata['TIME'])
		sdata['TERROR_L'] = hdata['TERROR_L']
		sdata['TERROR_R'] = hdata['TERROR_R']
		# Binning total source and background lightcurves
		"""
		totdata = rebin_data(totdata, factor, func=elsesum)
		totdata['TIME'] = hdata['TIME']
		bkgdata = rebin_data(bkgdata, factor, func=elsesum)
		bkgdata['TIME'] = hdata['TIME']
		"""

	# Combined background-subtracted lightcurve
	cdata = np.copy(hdata)
	cdata['RATE'] = hdata['RATE'] + sdata['RATE']
	cdata['ERROR'] = np.sqrt(hdata['ERROR']**2 + sdata['ERROR']**2)
	cdata_terrors = timebin_errors(hdata['TIME'])
	cdata['TERROR_L'] = cdata_terrors[0]
	cdata['TERROR_R'] = cdata_terrors[1]

	# Hardness ratio lightcurve
	hrdata = np.copy(hdata)
	hrdata['RATE'] = (hdata['RATE']-sdata['RATE'])/(hdata['RATE']+sdata['RATE'])
	hrdata['ERROR'] = error_hardness_ratio(hdata['RATE'], sdata['RATE'], hdata['ERROR'], sdata['ERROR'])
	
	# ======================= PLOTS =======================
	
	# Plot 1: total line lightcurve
	fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
	# --> Top: total background-subtracted lightcurve
	axs[0].plot(totdata['TIME'], totdata['RATE']-bkgdata['RATE'], "b-")
	axs[0].set(ylabel = "Counts / s")
	# --> Bottom: source and background lightcurves.
	axs[1].plot(totdata['TIME'], totdata['RATE'], "b-")
	axs[1].plot(bkgdata['TIME'], bkgdata['RATE'], "r-")
	axs[1].set(ylabel = "Counts / s", xlabel = "Time (ks)")
	fig.subplots_adjust(hspace=0)
	plt.show()

	# Plot 2: dot-errorbar lightcurve
	fig, axs = plt.subplots(nrows=1, ncols=1)
	axs.errorbar(cdata['TIME'], cdata['RATE'],
			xerr=[cdata['TERROR_L'],cdata['TERROR_R']],
			yerr=cdata['ERROR'], fmt='b.')
	axs.set(ylabel = "Counts / s", xlabel = "Time (ks)")
	axs.set_title("Total background-subtracted light-curve")
	fig.subplots_adjust(bottom=0.45)
	plt.show()

	# Plot 3: simple hardness
	fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
	# 		NOTE: HARD X-RAYS ARE BLUE, SOFT X-RAYS ARE RED.
	# --> Top: Soft & Hard lightcurves
	axs[0].errorbar(hdata['TIME'], hdata['RATE'],
			xerr=[hdata['TERROR_L'],hdata['TERROR_R']],
			yerr=hdata['ERROR'],
			fmt='b.', linewidth=1, markersize=2)
	axs[0].errorbar(sdata['TIME'], sdata['RATE'],
			xerr=[sdata['TERROR_L'],sdata['TERROR_R']],
			yerr=sdata['ERROR'],
			fmt='r.', linewidth=1, markersize=2)
	axs[0].set(ylabel = "Counts / s")
	# --> Bottom: Hardness ratio (H/S)
	axs[1].plot(hdata['TIME'], hdata['RATE']/sdata['RATE'], "g-")
	axs[1].set(ylabel = "H/S", xlabel = "Time (ks)")
	axs[1].axhline(y=1, color='k', linestyle='--', linewidth=0.5)
	fig.subplots_adjust(hspace=0)
	plt.show()

	# Plot 4: reduced hardness
	fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
	# 		NOTE: HARD X-RAYS ARE BLUE, SOFT X-RAYS ARE RED.
	# --> Top: Soft & Hard lightcurves
	axs[0].errorbar(hdata['TIME'], hdata['RATE'],
			xerr=[hdata['TERROR_L'],hdata['TERROR_R']],
			yerr=hdata['ERROR'],
			fmt='b.', linewidth=1, markersize=2)
	axs[0].errorbar(sdata['TIME'], sdata['RATE'],
			xerr=[sdata['TERROR_L'],sdata['TERROR_R']],
			yerr=sdata['ERROR'],
			fmt='r.', linewidth=1, markersize=2)
	axs[0].set(ylabel = "Counts / s")
	# --> Bottom: Hardness ratio (H/S)
	hardness = (hdata['RATE'] - sdata['RATE']) / (hdata['RATE'] + sdata['RATE'])
	axs[1].plot(hdata['TIME'], hardness, "g-")
	axs[1].set(ylabel = "(H-S)/(H+S)", xlabel = "Time (ks)")
	axs[1].axhline(y=0, color='k', linestyle='--', linewidth=0.5)
	fig.subplots_adjust(hspace=0)
	plt.show()

	# Plot 5: cumulative lightcurves
	hardplot, = plt.plot(hdata['TIME'], hdata['CUMUL'], 'b-')
	softplot, = plt.plot(sdata['TIME'], sdata['CUMUL'], 'r-')
	plt.xlabel(xlabel="Time (ks)")
	plt.ylabel(ylabel="Cumulative rate")
	plt.title("Cumulative light-curves")
	if legend: plt.legend([hardplot,softplot],[hardname, softname])
	plt.show()

	# ------------------------------------------------------
	# Custom plot. Total lightcurve + hardness

	# 1 - Source and background: 1500 s
	# 2 - total bkg-subtracted: 1500 s
	# 3 - hardness ratio: 2500 s

	# ===== Custom rebin totdata and bkgdata =====
	errorsum = lambda x: np.sqrt(sum([d**2 for d in x ]))/len(x) 
	avgsum = lambda x: sum(x)/len(x)
	# Applying error propagation in rebinning
	factor = 2000//500
	err_totdata = rebin_array(totdata['ERROR'], factor, func=errorsum)
	err_bkgdata = rebin_array(bkgdata['ERROR'], factor, func=errorsum)
	# Rebinning all other fields
	totdata = rebin_data(totdata, factor, func=avgsum)
	bkgdata = rebin_data(bkgdata, factor, func=avgsum)
	totdata['ERROR'] = err_totdata
	bkgdata['ERROR'] = err_bkgdata
	tbin_err = timebin_errors(totdata['TIME'])

	# ===== Custom rebin soft and hard xray data =====
	
	factor = 1500//500
	sz = len(orig_hdata)
	# Applying error propagation in rebinning
	err_orig_hdata = rebin_array(orig_hdata['ERROR'], factor, func=errorsum)
	err_orig_sdata = rebin_array(orig_sdata['ERROR'], factor, func=errorsum)
	# Rebinning all other fields
	orig_hdata = rebin_data(orig_hdata, factor, func=avgsum)
	orig_sdata = rebin_data(orig_sdata, factor, func=avgsum)
	orig_hdata['ERROR'] = err_orig_hdata
	orig_sdata['ERROR'] = err_orig_sdata
	# Recalculating time errors
	orig_hdata['TERROR_L'], orig_hdata['TERROR_R'] = timebin_errors(orig_hdata['TIME'])
	orig_sdata['TERROR_L'] = orig_hdata['TERROR_L']
	orig_sdata['TERROR_R'] = orig_hdata['TERROR_R']

	# ============== Plotting ===============
	fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True)

	orig_hdata['RATE'] *= 1e5
	orig_hdata['ERROR'] *= 1e5
	orig_sdata['RATE'] *= 1e5
	orig_sdata['ERROR'] *= 1e5

	totdata['RATE'] *= 1e5
	totdata['ERROR'] *= 1e5
	bkgdata['RATE'] *= 1e5
	bkgdata['ERROR'] *= 1e5

	# --> Top: Source and background lightcurves
	axs[0].errorbar(totdata['TIME'], totdata['RATE'],
			xerr=tbin_err,
			yerr=totdata['ERROR'],
			fmt='k.', linewidth=1, markersize=2)
	axs[0].errorbar(bkgdata['TIME'], bkgdata['RATE'],
			xerr=tbin_err,
			yerr=bkgdata['ERROR'],
			fmt='g.', linewidth=1, markersize=2)
	axs[0].set(ylabel = "Count rate\n($10^{-5}$ s$^{-1}$)")
	axs[0].get_yaxis().set_label_coords(-0.075,0.5)

	# --> Middle: Soft & hard lightcurves
	axs[1].errorbar(orig_hdata['TIME'], orig_hdata['RATE'],
			xerr=[orig_hdata['TERROR_L'], orig_hdata['TERROR_R']],
			yerr=orig_hdata['ERROR'],
			fmt='b.', linewidth=1, markersize=2)
	axs[1].errorbar(orig_sdata['TIME'], orig_sdata['RATE'],
			xerr=[orig_sdata['TERROR_L'], orig_sdata['TERROR_R']],
			yerr=orig_sdata['ERROR'],
			fmt='r.', linewidth=1, markersize=2)
	axs[1].set(ylabel = "Count rate\n($10^{-5}$ s$^{-1}$)")
	axs[1].get_yaxis().set_label_coords(-0.075,0.5)

	# --> Bottom: Hardness ratio (H/S)
	hardness = hdata['RATE'] / sdata['RATE']
	hardness_err = np.sqrt( (hdata['ERROR']/sdata['RATE'])**2 + (sdata['ERROR']*hdata['ERROR']/sdata['RATE']**2)**2 )
	axs[2].errorbar(hdata['TIME'], hardness, xerr=[hdata['TERROR_L'],hdata['TERROR_R']], yerr=hardness_err, fmt='g.' )
	axs[2].set(ylabel = "Hardness ratio\n(H/S)", xlabel = "Time (ks)")
	axs[2].axhline(y=1, color='k', linestyle='--', linewidth=0.5)
	axs[2].set_ylim([-0.5,2.5])
	axs[2].get_yaxis().set_label_coords(-0.075,0.5)
	fig.subplots_adjust(hspace=0)
	plt.show()



def lc_bkg(lcfile, bkg):
	"""
	"""
	pass

# def lcratio( lcfile1, lcfile2, lcname1, lcrange1, lcbkg1, lcname2, lcrange2, lcbkg2)
def lcratio(fsoftx=None, fhardx=None, softrange=None, hardrange=None, rlog=False, extra=True):
	"""
	Generates and plots a ratio between the soft and hard X-ray
	light curves.

	Args:
		fsoftx    	(string) Light curve data file of soft xrays.
		fhardx   	(string) Light curve data file of hard xrays.
		softrange	(optional) (float, float) Soft xrays energy range in keV. For plotting purposes.
		hardrange	(optional) (float, float) Hard xrays energy range in keV. For plotting purposes.
		rlog		(optional) (bool) If true, plots the soft/hard x-ray ratio lightcurve in a log10 scale.
		extra		(optional) (bool) Plots cumulative lightcurves as well as hard-soft difference.

	"""
	# Input checking
	if(not fhardx or not fsoftx):
		print("[LCRATIO] Error: light curve filenames undefined")
		return None

	if not os.path.exists(fhardx):
		print(f"[LCRATIO] Error: file '{fhardx}' not found.")
		return None
	if not os.path.exists(fsoftx):
		print(f"[LCRATIO] Error: file '{fsoftx}' not found.")
		return None

	plot_srange = False
	plot_hrange = False
	if(Utils.indexable(softrange) and len(softrange)==2):
		plot_srange = True
	if(Utils.indexable(hardrange) and len(hardrange)==2):
		plot_hrange = True

	softlc_data = None
	hardlc_data = None
	scaling_soft = 1
	scaling_hard = 1
	with fits.open(fhardx) as image:
		hardlc_data = np.copy(image['RATE'].data)
		scaling_hard = _header_extract_aperture(image['RATE'].header)

	with fits.open(fsoftx) as image:
		softlc_data = np.copy(image['RATE'].data)
		scaling_soft = _header_extract_aperture(image['RATE'].header)

	# Ensure both lightcuves have the same number of data points
	if( len(hardlc_data['TIME']) != len(softlc_data['TIME']) ):
		print("[LCRATIO] Error: the input lightcurves don't have the same number of data points or time bins")
		return

	# Time bin errors
	tbinerr_raw = timebin_errors(hardlc_data['TIME'])

	# Gather indices where value is NAN or <= 0 from both light curves
	rm_ind = []
	for i in range(len(hardlc_data['TIME'])):
		if(math.isnan(hardlc_data['TIME'][i]) or math.isnan(hardlc_data['RATE'][i])
			or math.isnan(softlc_data['TIME'][i]) or math.isnan(softlc_data['RATE'][i])):
			rm_ind.append(i)
		#if( hardlc_data['RATE'][i] <= 0 or softlc_data['RATE'][i] <= 0 ):
		#	rm_ind.append(i)

	# Remove bad data from collected indices
	h_times_raw = np.delete(hardlc_data['TIME'], rm_ind)
	h_energies = np.delete(hardlc_data['RATE'], rm_ind)
	h_errors = np.delete(hardlc_data['ERROR'], rm_ind)
	s_times_raw = np.delete(softlc_data['TIME'], rm_ind)
	s_energies = np.delete(softlc_data['RATE'], rm_ind)
	s_errors = np.delete(softlc_data['ERROR'], rm_ind)
	tbinerr = [np.delete(tbinerr_raw[0], rm_ind), np.delete(tbinerr_raw[1], rm_ind)]
	
	print(f"[LCRATIO] Warning: {len(rm_ind)} invalid values removed")

	# Scaling lightcurves due to having different apertures
	if(not scaling_soft or not scaling_hard):
		print("[LCRATIO] Warning: unable to find aperture area in headers. No scaling will be applied.")
	else:
		print("[LCRATIO] Scaling factor: ", scaling_soft/scaling_hard)
		h_energies = h_energies * scaling_soft/scaling_hard
		h_errors = h_errors * scaling_soft/scaling_hard
	
	# Calibrating time origin
	s_times = s_times_raw - s_times_raw[0]

	#count_ratio = h_energies/s_energies
	count_ratio = (h_energies-s_energies)/(h_energies+s_energies)

	# Plotting lightcurves and their ratio

	if(softrange and hardrange):
		legend_text = ["Soft X-rays" + plot_srange * f"({softrange[0]:.2f}-{softrange[1]:.2f} keV)",
					   "Hard X-rays" + plot_hrange * f"({hardrange[0]:.2f}-{hardrange[1]:.2f} keV)"]
	else:
		legend_text = ["Soft X-rays", "Hard X-rays"]

	fig, axs = plt.subplots(nrows=2, ncols=1)

	softx_plot = axs[0].errorbar(np.asarray(s_times)/1000, s_energies, fmt='b.', xerr=np.asarray(tbinerr)/1000, yerr=s_errors)
	hardx_plot = axs[0].errorbar(np.asarray(s_times)/1000, h_energies, fmt='r.', xerr=np.asarray(tbinerr)/1000, yerr=h_errors)
	axs[0].set_title(" Soft and hard X-ray light-curves")
	axs[0].set(ylabel = "Counts", xlabel = "Time (ks)")
	axs[0].legend([softx_plot,hardx_plot], legend_text)
	
	#ratio_err = [np.sqrt((h_errors[i]/s_energies[i])**2+(s_errors[i]*h_energies[i]/s_energies[i]**2)**2) for i in range(len(s_errors))]
	ratio_err = [0] * len(s_errors)
	axs[1].errorbar(np.asarray(s_times)/1000, count_ratio, fmt='k.', xerr=np.asarray(tbinerr)/1000, yerr=ratio_err, linewidth=0.5, markersize=1)
	axs[1].plot(np.asarray(s_times)/1000, count_ratio, 'g-')
	axs[1].axhline(y=0, color='k', linestyle='--', linewidth=0.8)
	axs[1].set(ylabel = "$(H-S)/(H+S)$")
	axs[1].set_title(" Hardness lightcurve")
	axs[1].set(xlabel = "Time (ks)")
	ylim = max( np.abs(np.max(count_ratio)), np.abs(np.min(count_ratio)) )
	if(rlog): axs[1].set_yscale('log')
	#else: axs[1].set_ylim(ymax= 1.1*ylim + 1, ymin= -1.1*ylim + 1)

	# Increases spacing between the two subplots so that title and axis label don't overlap.
	fig.subplots_adjust(hspace=0.45)
	plt.show()

	if(not extra): return
	
	# Plotting cumulative lightcurves
	s_cumul = np.cumsum(s_energies)
	h_cumul = np.cumsum(h_energies)

	softx_cplot, = plt.plot(np.asarray(s_times)/1000, s_cumul, 'b-')
	hardx_cplot, = plt.plot(np.asarray(s_times)/1000, h_cumul, 'r-')
	plt.xlabel(xlabel="Time (s)")
	plt.ylabel(ylabel="Total counts")
	plt.title("Cumulative hard and soft light-curves")
	plt.legend([softx_cplot,hardx_cplot], legend_text)
	plt.show()

	

def cumlc(lcfile=None, plotlc=True):
	"""
	Generates and plots a cumulative light curve.

	Args:
		lcfile: 	(string) fits file with lightcurve data
		plotlc:		(optional) (bool) Plots original lightcurve

	"""
	if(not lcfile):
		print("[CUMLC] Error: no lightcurve filename provided")
		return

	if not os.path.exists(lcfile):
		print(f"[CUMLC] Error: file '{lcfile}' not found.")
		return

	lc_data = None
	with fits.open(lcfile) as image:
		lc_data = np.copy(image['RATE'].data)

	times_raw = np.nan_to_num(lc_data['TIME'], 0)
	times = times_raw - times_raw[0]
	energies = np.nan_to_num(lc_data['RATE'], 0)
	energy_err = np.nan_to_num(lc_data['ERROR'], 0)

	# Plotting normal lightcurve
	if(plotlc == True):
		plt.plot(times/1000, energies)
		plt.title(f" Counts plot across time \n {lcfile}")
		plt.xlabel("Time (ks)")
		plt.ylabel("Counts / s")
		plt.show()
		plt.clf()
	
	# Gets the cumulative sum of the energies across time
	cumulative_energies = np.cumsum(energies)

	plt.errorbar(times/1000, cumulative_energies, yerr=energy_err)
	plt.title(f" Cumulative counts plot across time \n {lcfile}")
	plt.xlabel("Time (ks)")
	plt.ylabel("Cumulative rate (counts / s)")
	plt.show()





def euv(lx=0, radius=0, rmin=0, rmax=0, lx_err = 0, radius_err =0, verbose = True, calc_rosat=True):
	"""
	Applies a relation by King (2018) to convert X-ray range fluxes
	to EUV fluxes.

	Args:
		lx: 	(float) X-ray luminosity of the star (erg/s)
		radius: (float) Radius of the star in solar radii.
		rmin:	(float) Lower bound of energy range in keV
		rmax:	(float) Higher bound of energy range in keV
		lx_err: (float) Uncertainty of the X-ray luminosity
		radius_err (float) Uncertainty of the stellar radius.
		verbose:(bool) Prints extra information

	Returns:
		euv_lum:	(float) Estimated EUV luminosity from 0.0136 keV to the minimum input energy range. 
		(-1):		(int) if input energy range is not covered by the model, or other error occurs.
	"""
	euv_data = Resources.king_euv()

	if(lx <= 0 or radius <= 0):
		print("[EUV] Error: X-ray luminosity and radius must be greater than zero.")
		return -1
	fx = lx / (4 * math.pi * (6.957e10 * radius )**2 ) #radius converted to cm2
	fx_err = math.sqrt( ( lx_err/(4*math.pi*(6.957e10*radius)**2 ))**2 + ( radius_err*lx/(4*math.pi*(6.957e10**2)*(radius)**3) )**2 )
	
	if( fx < 10):
		print("[EUV] Warning: the flux is quite low. Are your sure you're using the surface X-ray flux in erg/s/cm^2 ?")

	which_range = -1
	for i in range(len(euv_data['xrange_i'])):
		if rmin == euv_data['xrange_i'][i] and rmax == euv_data['xrange_f'][i]:
			which_range = i
			break

	if(which_range == -1):
		print("[EUV] Error: input X-ray energy range not covered by model.")
		return -1
	
	euv_ratio = euv_data['const'][i] * (fx ** euv_data['pwlaw'][i])
	euv_ratio_err = fx_err * euv_data['const'][i] * euv_data['pwlaw'][i] * fx ** (euv_data['pwlaw'][i] - 1)
	
	euv_flux = euv_ratio * fx
	euv_flux_err = math.sqrt( (euv_ratio_err * fx)**2 + (fx_err * euv_ratio)**2 )
	
	euv_lum = euv_flux * 4 * math.pi * (6.957e10 * radius )**2
	euv_lum_err = math.sqrt( (euv_flux_err*4*math.pi*(6.957e10*radius)**2)**2 + (radius_err*euv_flux*8*math.pi*radius*(6.957e10)**2)**2 )

	if verbose:
		print(f"[EUV] Using energy ranges: EUV 0.0136 to {rmin:.3f}, and X-ray {rmin:.3f} to {rmax:.3f} keV")
		print(f"[EUV] F_EUV / F_x = {euv_ratio:.3e} +/- {euv_ratio_err:.3e}")
		print(f"[EUV] L_EUV = {euv_lum:.3e} +/- {euv_lum_err:.3e} erg/s")
		print(f"[EUV] L_XUV = {(euv_lum+lx):.3e} +/- {(math.sqrt((euv_lum_err)**2+(lx_err)**2)):.3e} erg/s")

	# Calculating ROSAT X-ray flux
	if(calc_rosat):
		from scipy.optimize import fsolve
		r_ind = 0
		func = lambda fx_rosat : fx_rosat + euv_data['const'][r_ind]*fx_rosat**(euv_data['pwlaw'][r_ind]+1) - (euv_flux+fx)
		guess = 2*fx
		solution, = fsolve(func, guess)
		rosat_lum = solution * 4 * math.pi * (6.957e10 * radius )**2
		print(f"[EUV] Estimated ROSAT X-ray flux: {rosat_lum:.3e} erg/s")
		rosat_xuv = euv(rosat_lum, radius, 0.1, 2.4, verbose=False, calc_rosat=False) + rosat_lum
		#print(f"[EUV] Residual {(abs(rosat_xuv - (euv_lum+lx))):.3e} to expected XUV flux.")

	return euv_lum





