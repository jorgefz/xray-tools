# xray-tools

Compilation of useful python functions for analysing stellar spectra and light-curves in the X-ray regime.

## Requirements

* scipy
* matplotlib
* astropy

## Installation

First, clone this repository somewhere in your machine
```bash
git clone https://github.com/jorgefz/xraytools xraytools
```

Then install the module
```bash
pip install -e xraytools
```

It should now be available to be imported from both the interpreter and python scripts.
```python
>>> import xraytools as xrt
```

## Features

This module contains five main submodules:
* plots: 
* utils: varios functions to calculate physical quantities (Rossby number, EUV flux, etc)
* resources: access relevant datasets and relations (Wright et al (2011) relation, Hyades dataset by Freund et al (2020), etc)
* sas: runs XMM-SAS tasks (more info at https://www.cosmos.esa.int/web/xmm-newton/sas)
* boltable: interpolates stellar quantities from the table by Eric Mamajek (https://www.pas.rochester.edu/~emamajek/EEM\_dwarf\_UBVIJHK\_colors\_Teff.txt)

### Utils

```python
from xraytools.utils import Utils
```

List of utilities (WIP):
* rossby\_number
* euv
* ...

### Boltable

```python
>>> from xraytools.boltable import Boltable
```

Obtain the approximate spectral type of a star from its mass:
```python
>>> Boltable.spt('Msun', 0.7) # Returns the spectral type of a 0.7 solar mass star
'K4.5V'
```

Obtain a dictionary with all interpolated quantities:
```python
>>> stardata = Boltable.interpolate('Msun', 0.71)
>>> stardata
{'Msun': 0.71, 'SpT': 'K4.5V', 'Teff': 4540.0, 'logT': 3.657, 'logL': -0.68, 'Mbol': 6.44, 'BCv': -0.6, 'Mv': 7.04, 'B-V': 1.116, 'Bt-Vt': 1.328, 'G-V': -0.425, 'Bp-Rp': 1.38, 'G-Rp': 0.72, 'M_G': 6.62, 'b-y': nan, 'U-B': 1.028, 'V-Rc': 0.654, 'V-Ic': 1.216, 'V-Ks': 2.781, 'J-H': 0.552, 'H-Ks': 0.127, 'Ks-W1': 0.04, 'W1-W2': -0.073, 'W1-W3': -0.03, 'W1-W4': 0.017, 'M_Ks': 4.26, 'i-z': nan, 'z-Y': nan, 'R_Rsun': 0.737}
```

If you want to interpolate using a known quantity other than mass, use:
```python
>>> Boltable.fields()
```
to get a list of available fields.


### Plots

```python
import xraytools.plots as Plots
```
* rotation\_tracks
* lcratio
* lc\_hardness
* ...


### Resources

```python
from xraytools.resources import Resources
```

### Sas

```python
from xraytools.sas import Sas
```



## Sources (WIP)
Wright et al (2011)
Jackson et al (2012)
King et al (2017)
Freund et al (2020)
Mamajek (2020)
