
#from astropy.io import fits
import numpy as np
import os
from scipy.interpolate import interp1d
from xraytools.resources import Resources

class Boltable:

        _dset = Resources.boltable()
        _fields = _dset.dtype.names
        interp_funcs = dict()

        def fields():
                return Boltable._fields

        def spt(field, value):
                """
                Retrieves nearest spectral types of a star
                """
                if field not in Boltable._fields:
                    print(f" Error: unknown field '{field}'")
                    print(" Use 'Boltable.fields()' to get a list of available fields")
                    return None
                data_field = Boltable._dset[field]
                try: idx = np.nanargmin( np.abs(data_field - value) )
                except ValueError: return None
                return Boltable._dset['SpT'][idx]

        def interpolate(using_field, using_value):
                """
                From a field and a value, returns a row with all data in boltable
                It doesn't work with spectral types
                """
                if using_field not in Boltable._fields:
                    print(f" Error: unknown field '{using_field}'")
                    print(" Use 'Boltable.fields()' to get a list of available fields")
                    return None
                data_field = Boltable._dset[using_field]
                row = dict()
                row[using_field] = using_value
                row['SpT'] = Boltable.spt(using_field, using_value)
                
                # Reusing past interpolating functions
                if using_field in Boltable.interp_funcs.keys():
                        for i,f in enumerate(Boltable._fields[1:-1]):
                                if f == using_field: continue
                                func = Boltable.interp_funcs[using_field][f]
                                if not callable(func): continue
                                row[f] = float(func(using_value))
                        return row


                ifuncs = dict()
                for i,f in enumerate(Boltable._fields[1:-1]):
                        if f == using_field: continue
                        try:
                                interp = interp1d(data_field, Boltable._dset[f])
                        except ValueError:
                                #print(f"[xraytools.boltable.get_row] Warning: unable to interpolate '{f}' with '{using_field}'")
                                row[f] = np.nan
                                ifuncs[f] = 0.0
                        else:
                                row[f] = float(interp(using_value))
                                ifuncs[f] = interp
                Boltable.interp_funcs[using_field] = ifuncs

                return row
