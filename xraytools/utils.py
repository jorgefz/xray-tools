
import numpy as np
import argparse
import os
import math

"""
	----- X-ray Tools -----

	For analysing stellar data in the X-ray regime.

	--- Utils library ---
	
	Author: Jorge Fernandez, University of Warwick, Coventry, UK
	Date: 17th December 2020

	Functions:
		indexable
		rossby_number
		timebin_errors
		header_extract_aperture
		remove_all_nans
		rebin_array
		rebin_data
		error_division
		error_hardness_ratio

"""

from xraytools.resources import Resources

class Utils:

	def indexable(obj):
		"""
		Returns true if an input object is subscriptable or idnexable (e.g. if obj[0] or len(obj) works)
		"""
		try:
			_ = obj[:]
		except TypeError:
			return False
		else:
			return True


	def rossby_number(vkcolor = 0, prot = 0):
		"""
		Estimates the Rossby number of a star given its V-K color index and its rotation period in days.
		It returns the Rossby number and the turnover time in days.
		The empirical relation comes from the Wright et al 2011 paper:
			(https://ui.adsabs.harvard.edu/abs/2011ApJ...743...48W/abstract)

		Input:
			vkcolor:	V-K color index
			prot:		Rotation period in days
		Returns:
			[0]: rossby number
			[1]: turnover time (days)
		"""
		rossby = 0
		turnover_time = 0
		if(3.5 >= vkcolor >= 1.1):
			turnover_time = 10**(0.73 + 0.22*vkcolor)
			rossby = prot / turnover_time
		elif(6.6 >= vkcolor > 3.5):
			turnover_time = 10**(-2.16 + 1.5*vkcolor - 0.13*vkcolor**2)
			rossby = prot / turnover_time
		return rossby, turnover_time


	def timebin_errors(data):
		"""
		Generates error bars based on input array of data
		as half the spacing between each two data points
		"""
		right_err = [ np.abs(data[i]-data[i+1])/2 for i in range(len(data)-1)]
		left_err = [ np.abs(data[i]-data[i-1])/2 for i in range(1,len(data))]
		right_err.insert(-1, left_err[-1])
		left_err.insert(0, right_err[0])
		return [left_err, right_err]


	def header_extract_aperture(headers):
		"""
		Given a list of headers, it seeks the circle aperture within
		the selection expression, exctracts the radius, and returns
		the aperture area in pixels squared.
		"""
		selexpr = headers['SLCTEXPR']
		if("circle" in selexpr): circle_ind = selexpr.find("circle")
		elif("CIRCLE" in selexpr): circle_ind = selexpr.find("CIRCLE")
		else: return None
		start = circle_ind + len("circle(")
		end = selexpr[start:].find(")") + start - 1
		data = [float(i) for i in selexpr[start:end].split(',')]
		area = np.pi * data[2]**2
		return area


	def remove_all_nans(*datasets):
		"""
		Removes all NAN values from a list of numpy arrays,
		or numpy structured arrays.
		"""
		nans = []
		for dset in datasets:
			if not dset.dtype.names:
				nans += list(np.nonzero(np.isnan(dset[f]))[0])
				continue
			for f in dset.dtype.names:
				nans += list(np.nonzero(np.isnan(dset[f]))[0])
		datasets = list(datasets)
		for i,dset in enumerate(datasets):
			datasets[i] = np.delete(dset, nans)
		return datasets


		
	def rebin_array(input_data, factor, func=None):
		"""
		Rebins an array to a lower resolution by a factor,
		and returns the new array.
		Leftover values are added up on the last bin.
		"""
		if not func: func = sum
		factor = int(factor)
		if(factor < 1):
			raise TypeError("binning factor must be an integer greater than 1")
		data = np.copy(input_data)
		leftover = len(data) % factor #Extra values that don't fit reshape
		leftover_ind = len(data) - leftover
		ndata = data[:leftover_ind].reshape((len(data)//factor, factor))
		ndata = np.array([func(newbin) for newbin in ndata])
		# Append leftover
		if (leftover > 0):
			ndata = np.append(ndata, func(data[leftover_ind:]))
		return ndata

		
	def rebin_data(data, factor, fields=None, func=None):
		"""
		Rebins an structured array to a lower resolution by a factor
		and returns the new array.
		"""
		if not data.dtype.names:
			return rebin_array(input_data, factor, func)
		ndata = data[0:len(data)//factor+(len(data)%factor > 0)]
		colnames = data.dtype.names
		if fields: colnames = fields
		for n in colnames:
			ndata[n] = rebin_array(data[n], factor, func)
		return ndata
			

	def error_division(x, y, xerr, yerr):
		"""
		Given the equation z = x / y, where x and y have uncertainties,
		it calculates and returns the uncertainty of z.
		"""
		dx = 1/y
		dy = x/y**2
		return np.sqrt( (xerr*dx)**2 + (yerr*dy)**2 )

	def error_hardness_ratio(x, y, xerr, yerr):
		"""
		Given the equation z = (x-y)/(x+y), where x and y have uncertainties,
		it calculates and returns the uncertainty of z.
		"""
		return np.abs(2*y/(x+y)**2) * np.sqrt(xerr**2+yerr**2)


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
			rosat_xuv = Utils.euv(rosat_lum, radius, 0.1, 2.4, verbose=False, calc_rosat=False) + rosat_lum
			#print(f"[EUV] Residual {(abs(rosat_xuv - (euv_lum+lx))):.3e} to expected XUV flux.")

		return euv_lum
