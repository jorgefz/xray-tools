import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import argparse
import os
import math


helpme = """
X-Ray Tools python package

A compilation of python scripts for facilitating X-ray analysis of stars.
Author: Jorge Fernandez, University of Warwick, Coventry, UK

Functions:
rotlaw 		:plots Wright et al 2011 model for Rossby number vs X-ray/bolometric ratio.
agerel 		:plots Jackson et al 2012 model for colour index vs X-ray/bolometric ratio vs age.
xraydiv 	:finds the soft/hard x-ray divide of a star's X-ray spectrum.
lcratio 	:plots the ratio between the soft and hard X-ray light curves.
cumlc 		:plots a cumulative light curve.

"""


def rotlaw(star_lratio=0, star_rot=0, star_rossby=0):
	"""
	This function plots a star's rotation period and X-ray luminosity
	versus bolumetric luminosity ratio along with a model by Wright et al 2011
	(https://ui.adsabs.harvard.edu/abs/2011ApJ...743...48W/abstract)

	Args:
		star_lratio:	(float) X-ray luminosity - bolumetric luminosity ratio.
		star_rot:   	(float) Star's rotation period in days.
		star_rossby: 	(float) Star's Rossby number (rotation period/turnover time)

	The luminosity ratio is always required.
	Then, either the star rotation or the Rossby number need to be defined.
	If the Rossby number is defined, the star is fit to the Rossby number-based model.
	On the other hand, if only the rotation period is defined, the star is plotted
	on the less accurate rotation period model, and it's expected Rossby number and
	convective turnover time based on the model are printed.
	Finally, if both are defined, the Rossby number takes preference.
	"""

	# Input checking
	includes_rossby = True

	if(star_lratio == 0):
		print("[ROTLAW] Error: Lx/Lbol undefined.")
		return
	if(star_rot == 0 and star_rossby == 0):
		print("[ROTLAW] Error: rotation period and Rossby number undefined.")
		return
	# Cannot retrieve or calculate Rossby number. Estimate from model.
	if(star_rot != 0 and star_rossby == 0):
		print("[ROTLAW] Warning: no Rossby number provided.")
		includes_rossby = False

	"""
	Defining model parameters:
			 { rx_sat 				(for Rx < ro_sat) 
		Rx = {
			 { C * Ro^(powerlaw) 	(for Rx >= ro_sat) 
	"""
	powerlaw = -2.7
	pwl_err = 0.16
	ro_sat = 0.13	# Saturation Rossby number
	ro_sat_err = 0.02
	rx_sat = 10**(-3.13)
	rx_sat_err = 0.08	# Saturation luminosity ratio (log10)

	# Rossby number v Rx model
	ro_init = 0.001
	ro_end = 10
	ro_points = 1000
	ro = np.linspace(ro_init, ro_end, ro_points)

	# Rotaton period v Rx model
	rot_init = 0.1
	rot_end = 100
	rot_points = 1000
	rot = np.linspace(rot_init, rot_end, rot_points)
	rot_sat = 5 #days


	# Building Rossby number model
	const = rx_sat / (ro_sat**powerlaw)
	model = np.array([const*r**powerlaw if r>ro_sat else rx_sat for r in ro])

	rx_sat_uerr = rx_sat*10**rx_sat_err
	rx_sat_lerr = rx_sat/10**rx_sat_err
	ro_sat_uerr = ro_sat+ro_sat_err
	ro_sat_lerr = ro_sat-ro_sat_err

	const_uerr = (rx_sat_uerr) / (ro_sat_uerr**(powerlaw+pwl_err))
	const_lerr = (rx_sat_lerr) / (ro_sat_lerr**(powerlaw-pwl_err))
	model_uerr = np.array([const_uerr*r**(powerlaw+pwl_err) if r>ro_sat_uerr else rx_sat_uerr for r in ro])
	model_lerr = np.array([const_lerr*r**(powerlaw-pwl_err) if r>ro_sat_lerr else rx_sat_lerr for r in ro])

	# Building Rotation period model
	const_rot = rx_sat / rot_sat**powerlaw
	model_rot = np.array([const_rot*r**powerlaw if r>rot_sat else rx_sat for r in rot])
	model_rot_uerr = np.array([const_rot*r**(powerlaw+pwl_err) if r>rot_sat else rx_sat for r in rot])
	model_rot_lerr = np.array([const_rot*r**(powerlaw-pwl_err) if r>rot_sat else rx_sat for r in rot])

	# Expected Rossby number and turnover time for the star based on the model.
	if(includes_rossby == False):
		expected_rossby = (star_lratio/const) ** (1/powerlaw)
		exp_ross_err = np.sqrt( (np.abs(const-const_uerr))**2*(expected_rossby/powerlaw/const)**2 + (pwl_err)**2*(expected_rossby*np.log(star_lratio/const)/powerlaw**2)**2 )
		expected_conv = star_rot/expected_rossby
		exp_conv_err = exp_ross_err * star_rot / expected_rossby**2
		print(f" Model fits:\n   Rossby number: {expected_rossby:.3f}+/-{exp_ross_err:.3f}\n   Turnover time: {expected_conv:.3f}+/-{exp_conv_err:.3f} days")

	#Plotting model & star
	if(includes_rossby == True):
		plt.plot(ro, model, 'b-')
		plt.plot(ro, model_uerr, 'b--', ro, model_lerr, 'b--')
		plt.plot(star_rossby, star_lratio, 'r*')
		plt.xlabel("Rossby number")
		plt.title("Rotation Law by Wright et al (2011)\n$R_x = (%.2E) Ro^{%.2f\\pm %.2f}$" % (const, powerlaw, pwl_err) )
	else:
		plt.plot(rot, model_rot, 'b--')
		plt.plot(star_rot, star_lratio, 'r*')
		plt.xlabel("Rotation period (days)")
		plt.title("Rotation Law by Wright et al (2011)\n$R_x = (%.2E) P_{rot}^{%.2f\\pm %.2f}$" % (const_rot, powerlaw, pwl_err) )

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
	powerlaw = model_data['powerlaw'][ind]

	ages = np.linspace(6.5, 10, 1000)
	const = satlum / (satage**powerlaw)
	model = np.array([const*(a**powerlaw) if a>satage else satlum for a in ages])

	plt.plot(ages, model, 'b.-')
	plt.plot(star_age+6, np.log10(star_lum), 'r*')

	plt.xlabel("Age log10(yr)")
	plt.ylabel("$log10(L_x / L_{bol})$")
	plt.title("Age relation by Jackson et al (2012)\n%.2f < (B-V) < %.2f" % (model_data['icolor'][ind], model_data['fcolor'][ind]))
	plt.show()
	



def xraydiv(file=None, rmf=None, verbose=False):
	"""
	This function extracts spectrum data from a star and calculates the energy
	boundary between soft and hard X-rays (energy at which 50% of the counts are
	on the left, and the remaining 50% on the right).

	Args:
		file:	fits filename where the spectrum data is stored.
		rmf:	response fits file of the star.

	Returns:
		div:	Energy boundary between soft and hard X-rays in keV.
		None: 	If necessary input arguments are missing.
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
	with fits.open(file) as image:
		data = np.copy(image['SPECTRUM'].data)

	with fits.open(rmf) as rmf_data:
		channel_data = np.copy(rmf_data['EBOUNDS'].data)

	# Calculating middle point between soft and hard X-rays.
	total_counts = data['COUNTS'].sum()
	print(f"Total counts = {total_counts}")

	energies = (channel_data['E_MAX'] + channel_data['E_MIN'])/2.0
	energy_err = channel_data['E_MAX'] - channel_data['E_MIN']

	cmcounts = 0
	middle_energy = 0
	for i,c in enumerate(data['COUNTS']):
		cmcounts += c
		if(cmcounts > total_counts/2):
			middle_energy = energies[i]
			break

	print(f"Soft/hard X-rays boundary is at {middle_energy:.3f} keV")

	if(verbose == False):
		return middle_energy

	counts_binned = []
	energies_binned = []

	# Get indices of where bins begin, ignoring bad channels
	bin_ind = [ i for i,g in enumerate(data['GROUPING']) if(data['QUALITY'][i] != 1 and g == 1) ]

	for i in range(len(bin_ind)-1):
		total = data['COUNTS'][bin_ind[i]:bin_ind[i+1]].sum()
		counts_binned.append( total )
		energies_binned.append( energies[bin_ind[i]:bin_ind[i+1]].mean() )

	# Normalizing counts
	counts_binned = [d/total_counts/energies_binned[i] for i,d in enumerate(counts_binned)]

	# Index at which energies are greater than 2 keV
	stop_ind = [i for i,d in enumerate(counts_binned) if energies_binned[i] >= 2 ][0]


	plt.plot(energies_binned[0:stop_ind], counts_binned[0:stop_ind], 'g.--', markersize=1)
	plt.title(f" Spectrum\n {file}")
	plt.xlabel("Energy (keV)")
	plt.ylabel("Normalized counts  keV$^{-1}$")
	plt.axvline(x=middle_energy, color='red', linestyle='-', linewidth=1)
	plt.text(middle_energy*1.01,max(counts_binned)*0.9,f'{middle_energy:.2f} keV',rotation=0)
	plt.show()
	plt.clf()





def lcratio(fsoftx=None, fhardx=None, softrange=None, hardrange=None, rlog=False, cumlc=False):
	"""
	Generates and plots a ratio between the soft and hard X-ray
	light curves.

	Args:
		fhardx   	Light curve data file of hard xrays.
		fsoftx    	Light curve data file of soft xrays.
		softrange	(optional) List of two floats, soft xrays energy range in keV.
					For plotting purposes.
		hardrange	(optional) List of two floats, hard xrays energy range in keV.
					For plotting purposes.
		rlog		(optional) If true, plots the soft/hard x-ray ratio lightcurve
					in a log10 scale.
		cumlc 		(optional) Plots cumulative lightcurves of soft and hard xray lightcurves

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
		fig.legend([softx_plot,hardx_plot] ,
			['Soft X-rays','Hard X-rays'])
	
	plt.show()
	plt.clf()

	if(not cumlc):
		return
	# Plotting cumulative lightcurves
	plt.clf()
	scumul = np.cumsum(s_energies)
	hcumul = np.cumsum(h_energies)
	plt.plot(s_times, scumul, 'b-', s_times, hcumul, 'r-', s_times, scumul/hcumul, 'g-')
	plt.title(f" Cumulative counts plot across time")
	plt.xlabel("Time (s)")
	plt.ylabel("Cumulative rate (counts / s)")
	plt.show()
	
	







def cumlc(lcfile=None, plotlc=True):
	"""
	Generates and plots a cumulative light curve.

	Args:
		lcfile: 	fits file with lightcurve data
		plotlc:		(optional) Plots original lightcurve

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
	if(plotlc):
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







