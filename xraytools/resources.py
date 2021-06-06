
import os
import sys
import numpy as np


class Resources:

	PyDir = os.path.dirname(os.path.realpath(__file__))

	path_boltable      = PyDir+"\\resources\\"+"bolometric_corrections.csv"
	path_wright_dset   = PyDir+"\\resources\\"+"wright_dataset.txt"
	path_jackson_dset  = PyDir+"\\resources\\"+"jackson_dataset.txt"
	path_jackson_model = PyDir+"\\resources\\"+"jackson_model.txt"
	path_king_euv      = PyDir+"\\resources\\"+"king2018_euv.txt"
	path_fip_laming    = PyDir+"\\resources\\"+"fip_laming.txt"
	path_hyades_dset   = PyDir+"\\resources\\"+"hyades_dataset.npy"
	path_hyades_freund = PyDir+"\\resources\\"+"hyades_freund_full.txt"
	path_hyades_tic    = PyDir+"\\resources\\"+"hyades_freund_tic.txt"

	def print_paths():
		print(path_boltable)
		print(path_boltable)
		print(path_wright_dset)
		print(path_jackson_dset)
		print(path_jackson_model)
		print(path_king_euv)
		print(path_fip_laming)
		print(path_hyades_dset)
		print(path_hyades_freund)
		print(path_hyades_tic)

	def boltable():
		"""
		Boltable: mamajek et al?
		"""
		dataset = np.genfromtxt(path_boltable, delimiter=',',
			filling_values=None, autostrip=True, names=True,
			deletechars='', dtype=None, encoding=None)
		return dataset

	def wright_dset():
		"""
		Paper: Wright et al (2011)
		"""
		dataset = np.genfromtxt(path_wright_dset, delimiter='|',
			filling_values=None, autostrip=True, names=True,
			deletechars='', dtype=None, encoding=None)
		return dataset

	def jackson_dset():
		"""
		Paper: Jackson et al (2012)
		"""
		dataset = np.genfromtxt(path_jackson_dset, delimiter='|',
			filling_values=None, autostrip=True, names=True,
			deletechars='', dtype=None, encoding=None)
		return dataset

	def jackson_model():
		"""
		Paper: Jackson et al (2012)
		Data fields:
			icolor: initial colou index, model boundary 
			fcolor: final colour idnex, model boundary
			satlum: saturation luminosity ratio. After this value, the ratio
					goes from being constant to following a powerlaw.
			satage: saturation age, after which the luminosity ratio follows a power law.
			powerlaw: power law index
			errpowerlaw: error in the power law index 
		"""
		model_data = np.genfromtxt(path_jackson_model,
					delimiter=';', names=True,
					dtype=float, skip_header=1)
		return model_data

	def king_euv():
		euv_data = np.genfromtxt(path_king_euv,
				delimiter=';', names=True,
				dtype=float, skip_header=2)
		return euv_data

	def fip_laming():
		dataset = np.genfromtxt(path_fip_laming,
				delimiter=',', names=True,
				dtype=float, skip_header=1)
		return dataset

	def hyades_dset():
		return np.load(path_hyades_dset, allow_pickle=True)

	def hyades_freund():
		pass

	def hyades_tic():
		pass



