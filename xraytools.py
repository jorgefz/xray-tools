import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import argparse
import os
import math


"""
	----- X-ray Tools -----

	For analysing stellar data in the X-ray regime.
	
	Author: Jorge Fernandez, University of Warwick, Coventry, UK
	Date: 6th December 2020

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

"""


def rotlaw(star_lratio=0, star_rossby=0, star_rot=0, vkcolor=None):
	"""
	This function plots a star's rotation period and X-ray luminosity
	versus bolumetric luminosity ratio along with a model by Wright et al 2011
	(https://ui.adsabs.harvard.edu/abs/2011ApJ...743...48W/abstract)

	Args:
		star_lratio:	(float) X-ray / bolumetric luminositiy ratio.
		star_rot:   	(float) Star's rotation period in days.
		star_rossby: 	(float) Star's Rossby number (rotation period / turnover time)
		vkcolor:     	(float) Star's V-Ks colour index.

	The luminosity ratio is always required.
	If the Rossby number is defined, the plot can be generated straight away.
	Otherwise, you will need to input both the star's rotation period in days as well
	as its V - Ks colour index.
	"""

	# Input checking
	includes_rossby = True
	# The luminosity ratio is always required
	if(star_lratio <= 0):
		print("[ROTLAW] Error: Lx/Lbol is invalid or undefined.")
		return
	# Decide to use either Rossby number or estimate it.
	if(star_rossby <= 0):
		includes_rossby = False
		if(star_rot <= 0 or vkcolor is None):
			print( ("[ROTLAW] Error: undefined parameters.\n"
				"Either Rossby number or both rotation period and V-Ks color required.") )
			return
		else:
			print("[ROTLAW] No Rossby number defined. Value will be estimated.")

	# Estimating turnover time and Rossby number, if undefined
	if(includes_rossby == False):
		turnover_time = 0
		# Check if colour is within valid range
		if(vkcolor < 1.1 or vkcolor > 6.6):
			print("[ROTLAW] Warning: V-Ks colour index outside valid range (1.1, 6.6)")
		if(vkcolor < 3.5):
			turnover_time = 10**(0.73 + 0.22*vkcolor) #days
			star_rossby = star_rot / turnover_time
		else:
			turnover_time = 10**(-2.16 + 1.5*vkcolor - 0.13*vkcolor**2)
		print(f" Estimated turnover time: {turnover_time:.3f} days")
		print(f" Estimated Rossby number: {star_rossby:.3f}")

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

	#Plotting model & star
	plt.plot(ross, model, 'b-')
	plt.plot(ross, model_uerr, 'b--', ross, model_lerr, 'b--')
	plt.plot(star_rossby, star_lratio, 'r*')
	plt.xlabel("Rossby number")
	plt.title("Rotation Law by Wright et al (2011)\n$R_x = (%.2E) Ro^{%.2f\\pm %.2f}$" % (const, powerlaw, pwl_err) )

	# Plot decorations
	plt.ylabel("$L_x / L_{bol}$")
	plt.yscale("log")
	plt.xscale("log")
	plt.show()
	plt.clf()





def agerel(star_lum=None, star_color=None, star_age=None):
	"""
	This function plots a star's colour index (B-V) and age
	along with an age relation by Jackson et al (2012):
	https://arxiv.org/abs/1111.0031

	Args:
		star_lum:	(float) X-ray / bolometric luminosity ratio of the star.
		star_color:	(float) Star's B-V colour index.
		star_age:	(float) Star's age in log10(Myr).

	Both parameters are required, as well as the predefined model parameters
	in the file 'jackson_model.txt'.
	"""

	if(star_lum is None or star_color is None or star_age is None):
		print("[AGEREL] Error: Star colour index or age undefined.")
		return
	if(not os.path.isfile("jackson_model.txt")):
		print("[AGEREL] Error: model parameter file 'jackson_model.txt' not found.")
		return

	model_data = np.genfromtxt("jackson_model.txt",
									delimiter=';', names=True,
									dtype=float, skip_header=1)
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
		if(star_color >= m['icolor'] and star_color < m['fcolor']):
			ind = i
			break
	if(ind < 0):
		min_color = min(model_data['icolor'])
		max_color = max(model_data['fcolor'])
		print(f"[AGEREL] Error: star's colour index outside model bounds ({min_color:.2f} to {max_color:.2f})")
		return
	
	# Fit and plot
	satlum = model_data['satlum'][ind]
	satage = model_data['satage'][ind]
	# It only looks like the plots in the paper by doubling the pw index! Why??
	powerlaw = model_data['powerlaw'][ind]*2

	ages = np.linspace(6.5, 10, 1000)
	const = satlum / (satage**powerlaw)
	model = np.array([const*a**(powerlaw) if a>satage else satlum for a in ages])

	plt.plot(ages, model, 'b-')
	plt.plot(star_age+6, np.log10(star_lum), 'r*')

	plt.xlabel("Age log10(yr)")
	plt.ylabel("$log10(L_x / L_{bol})$")
	plt.title("Age relation by Jackson et al (2012)\n%.2f < (B-V) < %.2f" % (model_data['icolor'][ind], model_data['fcolor'][ind]))
	plt.show()
	plt.clf()

	# Plots all models together
	"""
	for m in model_data:
		satlum = m['satlum']
		satage = m['satage']
		powerlaw = m['powerlaw']*2
		ages = np.linspace(6.5, 10, 1000)
		const = satlum / (satage**powerlaw)
		model = np.array([const*a**(powerlaw) if a>satage else satlum for a in ages])
		plt.plot(ages, model, 'b-')
	plt.show()
	"""


	



def xraydiv(file=None, rmf=None, erange=None, verbose=False):
	"""
	This function extracts spectrum data from a star and calculates the energy
	boundary between soft and hard X-rays (energy at which 50% of the counts are
	on the left, and the remaining 50% on the right).

	Args:
		file:	(string) file where the spectrum data is stored.
		rmf:	(string) RMF response file of the star.
		erange:	(float,float) (optional) Range of energies to include in plot (keV)
		verbose: (bool) (optional) Plots the spectrum with the energy divide.

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

	data = None
	channel_data = None
	exptime = 0	#Exposure time in seconds
	with fits.open(file) as image:
		data = np.copy(image['SPECTRUM'].data)
		# Calculating total exposure time in seconds
		import re
		"""This splits a date string from the header (e.g. "2018-12-5T16:28:33")
		into a list of ints (2018, 12, 5, 16, 28, 33) using a regular expression
		on the separators: dash, space, T, and semicolon.
		Then, the delta time is calculated using the datetime module.
		"""
		expstart = filter(None, re.split("[T \-:]+", image[0].header['EXPSTART']))
		expstop = filter(None, re.split("[T \-:]+", image[0].header['EXPSTOP']))
		from datetime import datetime
		exptime = (datetime(*map(int,expstop)) - datetime(*map(int,expstart))).total_seconds()

	with fits.open(rmf) as rmf_data:
		channel_data = np.copy(rmf_data['EBOUNDS'].data)

	# Calculating middle point between soft and hard X-rays.
	total_counts = data['COUNTS'].sum()
	if(verbose):
		print(f"[XRAYDIV] Total counts = {total_counts}")
		print(f"[XRAYDIV] Exposure time: {exptime} s")

	energies = (channel_data['E_MAX'] + channel_data['E_MIN'])/2.0

	cmcounts = 0
	middle_energy = 0
	for i,c in enumerate(data['COUNTS']):
		cmcounts += c
		if(cmcounts > total_counts/2):
			middle_energy = energies[i]
			break

	print(f"[XRAYDIV] Soft/hard X-rays boundary is at {middle_energy:.3f} keV")

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
	counts_binned = [d/total_counts/exptime/energy_per_bin[i] for i,d in enumerate(counts_binned)]

	# Plot only energy in defined range
	start_ind = 0
	stop_ind = len(counts_binned)
	if(erange is not None and isinstance(erange,(tuple, list)) ):
		start_ind_lst = [i for i,d in enumerate(counts_binned) if energies_binned[i] >= erange[0] ]
		stop_ind_lst = [i for i,d in enumerate(counts_binned) if energies_binned[i] >= erange[1] ]
		if(len(start_ind_lst) > 0):
			start_ind = start_ind_lst[0]
		if(len(stop_ind_lst) > 0):
			stop_ind = stop_ind_lst[0]


	plt.step(energies_binned[start_ind:stop_ind], counts_binned[start_ind:stop_ind], color='k', linewidth=1)
	plt.title(f" Spectrum\n {file}")
	plt.xlabel("Energy (keV)")
	plt.ylabel("Normalized counts s$^{-1}$ keV$^{-1}$")
	plt.axvline(x=middle_energy, color='red', linestyle='-', linewidth=1)
	plt.text(middle_energy*1.01,max(counts_binned)*0.9,f'{middle_energy:.2f} keV',rotation=0)
	plt.show()
	plt.clf()

	return middle_energy





def lcratio(fsoftx=None, fhardx=None, softrange=None, hardrange=None, rlog=False, extra=False):
	"""
	Generates and plots a ratio between the soft and hard X-ray
	light curves.

	Args:
		fhardx   	(string) Light curve data file of hard xrays.
		fsoftx    	(string) Light curve data file of soft xrays.
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
	if(isinstance(softrange, (list,tuple) ) and len(softrange)==2):
		plot_srange = True
	if(isinstance(hardrange, (list,tuple) ) and len(hardrange)==2):
		plot_hrange = True


	softlc_data = None
	hardlc_data = None
	with fits.open(fhardx) as image:
		hardlc_data = np.copy(image['RATE'].data)
	with fits.open(fsoftx) as image:
		softlc_data = np.copy(image['RATE'].data)


	# Ensure both lightcuves have the same number of data points
	if( len(hardlc_data['TIME']) != len(softlc_data['TIME']) ):
		print("[LCRATIO] Error: the input lightcurves don't have the same number of data points or time bins")
		return

	# Remove all NANs and values <= 0 from both light curves
	rm_ind = []
	for i in range(len(hardlc_data['TIME'])):
		if(math.isnan(hardlc_data['TIME'][i]) or math.isnan(hardlc_data['RATE'][i])
			or math.isnan(softlc_data['TIME'][i]) or math.isnan(softlc_data['RATE'][i])):
			rm_ind.append(i)
		if( hardlc_data['RATE'][i] <= 0 or softlc_data['RATE'][i] <= 0 ):
			rm_ind.append(i)

	print(f"[LCRATIO] Warning: {len(rm_ind)} invalid values removed")

	h_times_raw = np.delete(hardlc_data['TIME'], rm_ind)
	h_energies = np.delete(hardlc_data['RATE'], rm_ind)
	s_times_raw = np.delete(softlc_data['TIME'], rm_ind)
	s_energies = np.delete(softlc_data['RATE'], rm_ind)

	s_times = s_times_raw - s_times_raw[0]
	h_times = h_times_raw - h_times_raw[0]

	s_energy_err = softlc_data['ERROR']
	h_energy_err = hardlc_data['ERROR']

	count_ratio = h_energies/s_energies

	# Plotting lightcurves and their ratio
	fig, axs = plt.subplots(2)

	axs[0].set_title(" Soft and hard X-ray light-curves")
	softx_plot, = axs[0].plot(s_times, s_energies, 'b-')
	hardx_plot, = axs[0].plot(s_times, h_energies, 'r-')
	axs[0].set(ylabel = "Counts / s")

	axs[1].set_title(" Soft / Hard ratio light curve")
	axs[1].set(xlabel = "Time (s)")
	if(rlog):
		ratio_plot, = axs[1].plot(s_times, np.log10(count_ratio), 'g-')
		axs[1].axhline(y=0, color='k', linestyle='--', linewidth=0.8)
		axs[1].set(ylabel = "log10 $R_{hard}/R_{soft}$")
	else:
		ratio_plot, = axs[1].plot(s_times, count_ratio, 'g-')
		axs[1].axhline(y=1, color='k', linestyle='--', linewidth=0.8)
		axs[1].set(ylabel = "$R_{hard}/R_{soft}$")

	if(plot_srange and plot_hrange):
		fig.legend([softx_plot,hardx_plot] ,
			[f"Soft X-rays ({softrange[0]}-{softrange[1]} keV)",
			 f"Hard X-rays ({hardrange[0]}-{hardrange[1]} keV)"])
	else:
		fig.legend([softx_plot,hardx_plot],['Soft X-rays','Hard X-rays'])
	
	plt.show()
	plt.clf()

	if(not extra):
		return
	# Plotting cumulative lightcurves
	scumul = np.cumsum(s_energies)
	hcumul = np.cumsum(h_energies)
	plt.plot(s_times, scumul, 'b-', s_times, hcumul, 'r-', s_times, scumul/hcumul, 'g-')
	plt.title(f" Cumulative counts plot across time")
	plt.xlabel("Time (s)")
	plt.ylabel("Cumulative rate (counts / s)")
	plt.show()
	plt.clf()

	# Plotting subtracted lightcurve Hard - Soft
	plt.plot(s_times, h_energies-s_energies, 'b-')
	plt.title(f" Hard - Soft X-ray lightcurve")
	plt.xlabel("Time (s)")
	plt.ylabel("Counts / s")
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
	if(isinstance(plotlc, bool) and plotlc):
		plt.plot(times, energies)
		plt.title(f" Counts plot across time \n {lcfile}")
		plt.xlabel("Time (s)")
		plt.ylabel("Counts / s")
		plt.show()
		plt.clf()
	
	# Gets the cumulative sum of the energies across time
	cumulative_energies = np.cumsum(energies)

	plt.errorbar(times, cumulative_energies, yerr=energy_err)
	plt.title(f" Cumulative counts plot across time \n {lcfile}")
	plt.xlabel("Time (s)")
	plt.ylabel("Cumulative rate (counts / s)")
	plt.show()







