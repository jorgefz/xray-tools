import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import argparse
import os
import math
import subprocess as sp

from .utils import Utils


class Sas:
	"""
	Calls SAS functions in a pythonian way.
	"""

	# values for K2-136
	star_x = 24602
	star_y = 24379
	star_r = 300

	bkg_x = 23227
	bkg_y = 25419
	bkg_r = 1480

	def version():
		process = sp.run("sasversion", shell=True)
		if(process.returncode != 0):
			print('[sas] Error: failed to get SAS version. Are you sure it is initialised?')

	def set_star_position(x, y, r):
		sas.star_x = int(x)
		sas.star_y = int(y)
		sas.star_r = int(r)

	def set_bkg_position(bx, by, br):
		sas.bkg_x = int(bx)
		sas.bkg_y = int(by)
		sas.bkg_r = int(br)

	def lightcurve(filename, event_list, bins, erange=(0.15,2.4), trange=None, x=None, y=None, r=None, bx=None, by=None, br=None, dump=None, path='lightcurves/'):
		"""
		Generates a light-curve using SAS commands.
		Input:
			filename:	(string) Name tag for the light-curve. Will serve as a basis to generate other filenames.
			event_list:	(string) Path and filename of EPIC event list file.
			bins:		(int) Size of each time bin in seconds to which to bin light-curves.
			erange:		(float, float) Energy range (initial, final] to select in keV.
			trange:		(float, float, ...) Time ranges to collect from the start of the observation, in kiloseconds, as pairs of (start,stop,start,stop,...) etc.
			x,y,z:		(int, optional) Position and radius of circle from which to extract source.
			bx,by,bz:	(int, optional) Position and radius of circle from which to extract background.
			dump:		(float) Filename where to dump stdout and stderr if SAS process fails.
			path:       (string) Path where to generate lightcurve files
		"""
		if not x: x=sas.star_x
		else: x = int(x)
		if not y: y=sas.star_y
		else: y = int(y)
		if not r: r=sas.star_r
		else: r = int(r)
		if not bx: bx=sas.bkg_x
		else: bx = int(bx)
		if not by: by=sas.bkg_y
		else: by = int(by)
		if not br: br=sas.bkg_r
		else: br = int(br)

		if(Utils.indexable(erange) and len(erange)==2):
			pass
		else:
			print("[sas-lightcurve] Error: Input energy range must be an array of two floats")
			return None
		
		fdump = sp.PIPE
		if(dump):
			fdump = open(dump, "w")
				
		if (path and path[-1] != '/'): path += '/'

		source_filename = path + f"source_{filename}.lc"
		bkg_filename = path + f"bkg_{filename}.lc"
		final_filename = path + f"final_{filename}.lc"
		
		timeselect = ""
		if(Utils.indexable(trange) and len(trange)%2 == 0):
			# Obtain start and end dates of observation from event list file.
			evlist = fits.open(event_list)
			obs_start = evlist['EVENTS'].data['TIME'][0]
			#obs_end = evlist['EVENTS'].data['TIME'][-1]
			evlist.close()
			timeselect = " && (TIME in "
			for i in range(len(trange)//2):
				t_start = int(1000*trange[i*2] + obs_start)
				t_end = int(1000*trange[i*2+1] + obs_start)
				timeselect += f"[{t_start}:{t_end}],"
			timeselect = timeselect[:-1] # Remove last comma
			timeselect += ") "
		
		print("[sas-lightcurve] Time selection expression:")
		print("   ", timeselect)

		source_cmd = (""
			f"evselect table={event_list} "
			f"energycolumn=PI expression='#XMMEA_EP&&(PATTERN<=4) {timeselect}  && "
			f"((X,Y) IN circle({x},{y},{r})) && (PI in [{int(1000*erange[0])}:{int(1000*erange[1])}])' "
			f"withrateset=yes rateset={source_filename} "
			f"timebinsize={bins} maketimecolumn=yes makeratecolumn=yes ")

		bkg_cmd = (""
			f"evselect table={event_list} "
			f"energycolumn=PI expression='#XMMEA_EP&&(PATTERN<=4) {timeselect}  && "
			f"((X,Y) IN circle({bx},{by},{br})) && (PI in [{int(1000*erange[0])}:{int(1000*erange[1])}])' "
			f"withrateset=yes rateset={bkg_filename} "
			f"timebinsize={bins} maketimecolumn=yes makeratecolumn=yes ")


		final_cmd = (""
			f"epiclccorr srctslist={source_filename} "
			f"eventlist={event_list} outset={final_filename} "
			f"bkgtslist={bkg_filename} "
			"withbkgset=yes applyabsolutecorrections=yes")

		print(f"[sas-lightcurve] Generating source lightcurve '{source_filename}'...")
		process = sp.run(source_cmd, shell=True, stdout=fdump, stderr=fdump)
		if(process.returncode != 0):
			print("[sas-lightcurve] Error: process terminated with code ", process.returncode)
			return None
		
		print(f"[sas-lightcurve] Generating background lightcurve '{bkg_filename}'...")
		process = sp.run(bkg_cmd, shell=True, stdout=fdump, stderr=fdump)
		if(process.returncode != 0):
			print("[sas-lightcurve] Error: process terminated with code ", process.returncode)
			return None

		print(f"[sas-lightcurve] Generating final lightcurve '{final_filename}'...")
		process = sp.run(final_cmd, shell=True, stdout=fdump, stderr=fdump)
		if(process.returncode != 0):
			print("[sas-lightcurve] Error: process terminated with code ", process.returncode)
			return None

		print("[sas-lightcurve] Done!")


	def final_spectrum(filename, countbin, path='spectra/'):
		"""
		Generates a final spectrum using SAS commands.
		Input:
			filename:	(string) Name tag for the spectrum. Used to obtain the names for the source, bkg, rmf, and arf.
			countbin:	(int) Number of counts per bin to which to bin the spectrum.
		"""
		if (path and path[-1] != '/'): path += '/'

		source = path + f"source_{filename}.fits"
		bkg = path + f"bkg_{filename}.fits"
		rmf = path + f"{filename}.rmf"
		arf = path + f"{filename}.arf"
		final = path + f"final_{filename}_maxc{countbin}.fits"

		final_cmd = (f""
			f"specgroup spectrumset={source} "
			f"mincounts={countbin} oversample=10000 "
			f"rmfset={rmf} arfset={arf} "
			f"backgndset={bkg} groupedset={final}")
		
		print(f"[sas-spectrum] Generating final spectrum '{final}'...")
		process = sp.run(final_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
		if(process.returncode != 0):
			print("[sas-spectrum] Error: process terminated with code ", process.returncode)
			return None


	def spectrum(filename, event_list, countbin, trange=None, x=None, y=None, r=None, bx=None, by=None, br=None, dump=None, path='spectra/'):
		"""
		Generates a spectrum using SAS commands.
		Input:
			filename:	(string) Name tag for the spectrum. Will serve as a basis to generate spectra filenames.
			event_list:	(string) Path and filename of EPIC event list file.
			countbin:	(int) Number of counts per bin to which to bin the spectrum.	
			trange:		(float, float, ...) Time ranges to collect from the start of the observation, in kiloseconds, as pairs of (start,stop,start,stop,...) etc.
			x,y,z:		(int, optional) Position and radius of circle from which to extract source.
			x,y,z:		(int, optional) Position and radius of circle from which to extract source.
			bx,by,bz:	(int, optional) Position and radius of circle from which to extract background.
			dump:		(float) Filename where to dump stdout and stderr if SAS process fails.	
		"""
		if not x: x=sas.star_x
		else: x = int(x)
		if not y: y=sas.star_y
		else: y = int(y)
		if not r: r=sas.star_r
		else: r = int(r)
		if not bx: bx=sas.bkg_x
		else: bx = int(bx)
		if not by: by=sas.bkg_y
		else: by = int(by)
		if not br: br=sas.bkg_r
		else: br = int(br)

		# Hypothetical filename generation
		#fname = "spectra/" + {source|bk|final} + "_maxc" + countbin + "_" + {fpbc|no_fpbc} + ".fits"			

		fdump = sp.PIPE
		if(dump):
			fdump = open(dump, "w")
		
		timeselect = ""
		if(Utils.indexable(trange) and len(trange)%2 == 0):
			# Obtain start and end dates of observation from event list file.
			evlist = fits.open(event_list)
			obs_start = evlist['EVENTS'].data['TIME'][0]
			#obs_end = evlist['EVENTS'].data['TIME'][-1]
			evlist.close()
			timeselect = " && (TIME in "
			for i in range(len(trange)//2):
				t_start = int(1000*trange[i*2] + obs_start)
				t_end = int(1000*trange[i*2+1] + obs_start)
				timeselect += f"[{t_start}:{t_end}],"
			timeselect = timeselect[:-1] # Remove last comma
			timeselect += ") "
		
		print("[sas-spectrum] Time selection expression:")
		print("   ", timeselect)

		source_filename = f"spectra/source_{filename}.fits"
		bkg_filename = f"spectra/bkg_{filename}.fits"
		rmf_filename = f"spectra/{filename}.rmf"
		arf_filename = f"spectra/{filename}.arf"
		final_filename = f"spectra/final_{filename}.fits"

		source_cmd = (f""
			f"evselect table={event_list} "
			f"withspectrumset=yes spectrumset={source_filename} "
			f"energycolumn=PI spectralbinsize=5 "
			f"withspecranges=yes specchannelmin=0 specchannelmax=20479 "
			f"expression='(FLAG==0) && (PATTERN<=4) && ((X,Y) IN circle({x},{y},{r})) {timeselect}'")
		
		bkg_cmd = (f""
			f"evselect table={event_list} "
			f"withspectrumset=yes spectrumset={bkg_filename} "
			f"energycolumn=PI spectralbinsize=5 "
			f"withspecranges=yes specchannelmin=0 specchannelmax=20479 "
			f"expression='(FLAG==0) && (PATTERN<=4) && ((X,Y) IN circle({bx},{by},{br})) {timeselect}'")
				
		arf_cmd = (""
			f"arfgen spectrumset={source_filename} arfset={arf_filename} "
			f"withrmfset=yes rmfset={rmf_filename} "
			f"badpixlocation={event_list} detmaptype=psf")
				

		final_cmd = (f""
			f"specgroup spectrumset={source_filename} "
			f"mincounts={countbin} oversample=10000 "
			f"rmfset={rmf_filename} arfset={arf_filename} "
			f"backgndset={bkg_filename} groupedset={final_filename}")
		
		print(f"[sas-spectrum] Generating source spectrum '{source_filename}'...")
		process = sp.run(source_cmd, shell=True, stdout=fdump, stderr=fdump)
		if(process.returncode != 0):
			print("[sas-spectrum] Error: process terminated with code ", process.returncode)
			return None
		
		print(f"[sas-spectrum] Generating background spectrum '{bkg_filename}'...")
		process = sp.run(bkg_cmd, shell=True, stdout=fdump, stderr=fdump)
		if(process.returncode != 0):
			print("[sas-spectrum] Error: process terminated with code ", process.returncode)
			return None
		
		print(f"[sas-spectrum] Backscaling source and background...")
		process = sp.run(f"backscale spectrumset={source_filename} badpixlocation={event_list}", shell=True, stdout=fdump, stderr=fdump)
		if(process.returncode != 0):
			print("[sas-spectrum] Error: process terminated with code ", process.returncode)
			return None
		process = sp.run(f"backscale spectrumset={bkg_filename} badpixlocation={event_list}", shell=True, stdout=fdump, stderr=fdump)
		if(process.returncode != 0):
			print("[sas-spectrum] Error: process terminated with code ", process.returncode)
			return None
		
		print(f"[sas-spectrum] Generating rersponse file '{rmf_filename}'...")
		process = sp.run(f"rmfgen spectrumset={source_filename} rmfset={rmf_filename}", shell=True, stdout=fdump, stderr=fdump)
		if(process.returncode != 0):
			print("[sas-spectrum] Error: process terminated with code ", process.returncode)
			return None
		
		print(f"[sas-spectrum] Generating ancillary file '{arf_filename}'...")
		process = sp.run(arf_cmd, shell=True, stdout=fdump, stderr=fdump)
		if(process.returncode != 0):
			print("[sas-spectrum] Error: process terminated with code ", process.returncode)
			return None
		
		if(not sas.final_spectrum(filename, countbin)):
			return None

		print(f"[sas-spectrum] Done!")
		
		







