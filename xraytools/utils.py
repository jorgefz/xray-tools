
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