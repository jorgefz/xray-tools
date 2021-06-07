# xray-tools

Compilation of useful python functions for analysing stellar spectra and light-curves in the X-ray regime.

## Installation

First, clone this repository somewhere in your machine
```bash
git clone https://github.com/jorgefz/xraytools xraytools
```

Then install the python module
```bash
pip install -e xraytools
```

It should now be available to import from anywhere
```python
>>> from xraytools.boltable import Boltable
>>> Boltable.spt('Msun', 0.7) # Returns the spectral type of a 0.7 solar mass star
'O3V'
```

## Sources
Wright et al (2011)
Jackson et al (2012)
King et al (2017)
