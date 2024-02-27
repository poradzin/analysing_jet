#/usr/bin/env python
# -*-Python-*-
# Created by buchanj at 28 Nov 2016  15:06

# This script allows the user to create a PPF containing the
# input signals required for a JETTO run from the output n# file of a TRANSP run.

# The generated PPF has DDA name TRAU and is created under the
# users id.


import ppf
from netCDF4 import Dataset
import os.path
import sys
import getpass
import numpy as np
# ================================================================================================================

DDAname = 'TRAU'


# Signals to be integrated and saved

signals_to_int = {
    'EHEAT': 'EPIN',
    'IHEAT': 'IPIN',
    'PRAD'  : 'PRIN',
    'EPTR_OBS': 'EPFL',
    'EETR_OBS': 'EEFL',
    'IETR_OBS': 'IEFL'
}

units_to_int = {
    'EHEAT': 'W',
    'IHEAT' : 'W',
    'PRAD':'W',
    'EPTR_OBS': 'W',                                                                             
    'EETR_OBS': 'W',                                                                             
    'IETR_OBS': 'W' 
}
# List of additional signals requested
signals = {
    'BZXR' : 'BZXR',
    'RMNMP': 'RMNM',
    'RMJMP': 'RMJM',
    'CONDE': 'CHIE',
    'CONDI': 'CHII',
    'NEUTT' : 'PFUS',
    'OMEGDATA':'OMD',
    'TAUEA': 'TEA',
    'TAUEE': 'TEE',
    'TAUES': 'TES',
    'TAUPHI':'TAM',
    'CHPHI': 'CHM',
    'TAUA1': 'TE1',
    'EHEAT': 'EPOW',
    'IHEAT': 'IPOW',
    'PRAD' : 'PRAD',
    'UFASTPP': 'UFPP',
    'UFASTPA': 'UFPA',
    'UTHRM' : 'UTH',
    'V': 'VLP',
    'ND':'ND',
    'NH':'NH',
    'NT':'NT',
}

conversions = {
    'BZXR': 1.0e-3,
    'RMNMP': 1.0e-2,
    'RMJMP': 1.0e-2,
    'CONDE': 1.0e-4,
    'CONDI': 1.0e-4,
    'CONDEF': 1.0e-4,
    'NEUTT' : 17.59e6*1.602176e-19,
    'OMEGDATA':1.0,
    'TAUEA': 1.0,
    'TAUEE': 1.0,
    'TAUES': 1.0,
    'TAUA1': 1.0,
    'TAUPHI':1.0,
    'CHPHI':1.0e-4,
    'EHEAT': 1.0e6,
    'IHEAT': 1.0e6,
    'PRAD' : 1.0e6,
    'UFASTPP': 1.0e6,
    'UFASTPA': 1.0e6,
    'UTHRM' :  1.0e6,
    'V':1.0,
    'ND':1.0e6,
    'NT':1.0e6, 
    'NH':1.0e6, 
}

newunits = {
    'BZXR': 'TESLA*m',
    'RMNMP': 'm',
    'RMJMP': 'm',
    'CONDE': 'm**2/s',
    'CONDI': 'm**2/s',
    'CONDEF': 'm**2/s',
    'CHPHI': 'm**2/s',
    'NEUTT': 'W',
    'OMEGDATA':'rad/sec',
    'TAUEA': 's',
    'TAUEE': 's',
    'TAUES': 's',
    'TAUA1': 's',
    'TAUPHI': 's',
    'EHEAT': 'W*m**-3',
    'IHEAT': 'W*m**-3',
    'PRAD' : 'W*m**-3',
    'UFASTPP': 'W*m**-3',
    'UFASTPA': 'W*m**-3',
    'UTHRM' : 'W*m**-3',
    'V': 'V',
    'ND': 'm**-3',
    'NH': 'm**-3',
    'NT': 'm**-3',
}

# JETTO default list 2D signals from netCDF file

signals_JETTO = {
    'BDENS': 'BDN',
    'BDENS_H': 'BDNH',
    'BDENS_D': 'BDND',
    'BDENS_T': 'BDNT',
    'BDEP_D': 'BDPD',
    'CURB': 'CURB',
    'ECCUR': 'ECC',
    'LHCUR': 'LHC',
    'NE': 'NE',
    'NI': 'NI',
    'NIMP': 'NIMP',
    'NMINI': 'NMIN',
    'OMEGA': 'OME',
    'PBE': 'PBE',
    'PBI': 'PBI',
    'PBTH': 'PBTH',
    'PEECH': 'PEC',
    'PEICH': 'QRFE',
    'PIICH': 'QRFI',
    'PRAD': 'PRAD',
    'Q': 'Q',
    'SBE': 'SBE',
    'SBE_D': 'SBED',
    'SBE_T': 'SBET',
    'SBTH_H': 'SBH',
    'SBTH_D': 'SBD',
    'SBTH_T': 'SBT',
    'SCEV': 'SCEV',
    'SDEP_D': 'SDPD',
    'SVH': 'SVH',
    'SVD': 'SVD',
    'SVT': 'SVT',
    'TE': 'TE',
    'TI': 'TI',
    'TQIN': 'TQIN',
    'UBPAR_H': 'UPAH',
    'UBPRP_H': 'UPPH',
    'UBPAR_D': 'UPAD',
    'UBPRP_D': 'UPPD',
    'UBPAR_T': 'UPAT',
    'UBPRP_T': 'UPPT',
    'UMINPA': 'UMPA',
    'UMINPP': 'UMPP',
    'ZEFFP': 'ZEF',
}


conversions_JETTO = {
    'BDENS': 1.0e6,
    'BDENS_H': 1.0e6,
    'BDENS_D': 1.0e6,
    'BDENS_T': 1.0e6,
    'BDEP_D': 1.0e6,
    'CURB': 1.0e4,
    'ECCUR': 1.0e4,
    'LHCUR': 1.0e4,
    'NE': 1.0e6,
    'NI': 1.0e6,
    'NIMP': 1.0e6,
    'NMINI': 1.0e6,
    'OMEGA': 1.0,
    'PBE': 1.0e6,
    'PBI': 1.0e6,
    'PBTH': 1.0e6,
    'PEECH': 1.0e6,
    'PEICH': 1.0e6,
    'PIICH': 1.0e6,
    'PRAD': 1.0e6,
    'Q': 1.0,
    'SBE': 1.0e6,
    'SBE_D': 1.0e6,
    'SBE_T': 1.0e6,
    'SBTH_H': 1.0e6,
    'SBTH_D': 1.0e6,
    'SBTH_T': 1.0e6,
    'SCEV': 1.0e6,
    'SDEP_D': 1.0e6,
    'SVH': 1.0e6,
    'SVD': 1.0e6,
    'SVT': 1.0e6,
    'TE': 1.0,
    'TI': 1.0,
    'TQIN': 1.0e6,
    'UBPAR_H': 1.0e6,
    'UBPRP_H': 1.0e6,
    'UBPAR_D': 1.0e6,
    'UBPRP_D': 1.0e6,
    'UBPAR_T': 1.0e6,
    'UBPRP_T': 1.0e6,
    'UMINPA': 1.0e6,
    'UMINPP': 1.0e6,
    'ZEFFP': 1.0,
}

newunits_JETTO = {
    'BDENS': 'm**-3',
    'BDENS_H': 'm**-3',
    'BDENS_D': 'm**-3',
    'BDENS_T': 'm**-3',
    'BDEP_D': 'm**-3/s',
    'CURB': 'A*m**-2',
    'ECCUR': 'A*m**-2',
    'LHCUR': 'A*m**-2',
    'NE': 'm**-3',
    'NI': 'm**-3',
    'NIMP': 'm**-3',
    'NMINI': 'm**-3',
    'OMEGA': 'rad/s',
    'TE': 'ev',
    'TI': 'ev',
    'Q': '',
    'PBE': 'W*m**-3',
    'PBI': 'W*m**-3',
    'PBTH': 'W*m**-3',
    'PEECH': 'W*m**-3',
    'PEICH': 'W*m**-3',
    'PIICH': 'W*m**-3',
    'PRAD': 'W*m**-3',
    'SBE': 'm**-3/s',
    'SBE_D': 'm**-3/s',
    'SBE_T': 'm**-3/s',
    'SBTH_H': 'm**-3/s',
    'SBTH_D': 'm**-3/s',
    'SBTH_T': 'm**-3/s',
    'SCEV': 'm**-3/s',
    'SDEP_D': 'm**-3/s',
    'SVH': 'm**-3/s',
    'SVD': 'm**-3/s',
    'SVT': 'm**-3/s',
    'TQIN': 'Nm/m3',
    'UBPAR_H': 'J*m**-3',
    'UBPRP_H': 'J*m**-3',
    'UBPAR_D': 'J*m**-3',
    'UBPRP_D': 'J*m**-3',
    'UBPAR_T': 'J*m**-3',
    'UBPRP_T': 'J*m**-3',
    'UMINPA': 'J*m**-3',
    'UMINPP': 'J*m**-3',
    'ZEFFP': '',
}

# merging dictionaries
# since JETTO output has priority, any repetitions in dictionaries will be
# overidden by JETTO output dictionaries

signals = {**signals, **signals_JETTO}
conversions = {**conversions, **conversions_JETTO}
newunits = {**newunits, **newunits_JETTO}

# Routine for summing two signals and writing them to a PPF =======================================================

def convert_unit(data, unit_str):
    # Define a dictionary of known unit conversions
    conversion_dict = {
        'N/CM**2': ['Pa', 10000],
        '/SEC': ['1/s', 1],
        'Pascals': ['Pa', 1],
        'rad/sec': ['rad/s', 1],
        'Tesla**2': ['T**2', 1],
        'V/CM': ['V/m', 100],
        'V':['V',1],
        'N/CM3/SEC': ['N/m**3/s', 1e6],
        'WEBERS': ['Wb', 1],
        'V*s': ['Wb', 1],
        'Tesla*cm': ['T*m', 1e-1],
        'RAD/SEC': ['rad/s', 1],
        'CM': ['m', 0.01],
        'CM**2': ['m**2', 1e-4],
        'CM**3': ['m**3',1e-6],
        'AMPS/CM2': ['A/m**2', 1e4],
        'TESLA/CM2': ['T', 1e4],
        'VOLT*TESLA/CM': ['V*T/m', 1e2],
        'GRAMS/CM3': ['kg/m**3', 1e3],
        'Tesla': ['T', 1],
        'T*(cm/sec)': ['T*m/s', 0.01],
        'N/SEC': ['N/s', 1],
        '#ptcls': ['', 1],
        'AMPS': ['A', 1],
        'N': ['N', 1],
        'Nt-M/CM3/(RAD/S)': ['Nt-M/m**3/(rad/s)', 1e6],
        'WATTS': ['W', 1],
        'eV': ['eV', 1.0],
        'CM**2/SEC': ['m**2/s', 1e-4],
        '1/SEC': ['1/s', 1],
        'ATOMS/SEC': ['ATOMS/SEC', 1],
        'CM**-1': ['m**-1', 1e2],
        'CM**-4': ['m**-4', 1e-8],
        'SECONDS': ['s', 1],
        '#/CM**3': ['#/m**3', 1e6],
        'W/cm**3': ['W/m**3', 1e6],
        'WATTS/CM3': ['W/m**3', 1e6],
        'NT-M': ['NT-M', 1],
        'TESLA*CM': ['T*m', 1e-2],
        'JLES/CM3': ['J/m**3', 1e6],
        'WATTS/CM3/EV': ['W/m**3/eV)', 1e6],
        'EV': ['eV', 1.0],
        'AMP*TESLA/CM2': ['A*T/m**2',1e4],
        '': ['',1]
    }
    return (data*conversion_dict[unit_str][1], conversion_dict[unit_str][0])


def linear_composite(rootgrp, new_signal, signals):
    """
    Creates a linear combination of signals and returns numpy array .

    Parameters:
    -----------
    rootgrp: netCDF4.Dataset object
        The netcdf4 file object.
    new_signal: str
        The name of the new signal to be created.
    signals: dict
        A dictionary of signals to be combined in the form {signal_name: coefficient}.

    Returns:
    --------
    bool:
        Returns True if the signal was added successfully, and False otherwise.
    """
    # Check that all signals exist in the rootgrp
    signals_copy = signals.copy()
    for signal in signals_copy.keys():
        if signal not in rootgrp.variables:
            print(f"{signal} not found in NETCDF file")
            signals.pop(signal)
            #del_signals[signal]

    # Check if signals is empty
    if not signals:
        print("No valid signals were found in the NetCDF file.")
        return False

    # Check if the new signal key already exists in the NetCDF file
    if new_signal in rootgrp.variables:
        print(f"{new_signal} already exists in the NetCDF file.")
        return False
    
    # Check that all signals have the same shape
    shapes = [rootgrp.variables[signal].shape for signal in signals.keys()]
    #print(f'Shapes: {shapes}')
    if len(set(shapes)) > 1:
        print("Signals have different shapes!")
        return
    else:
        shape = shapes[0]
        #print(f'Shape: {shape}')
    
    # Check that all signals have the same units
    units = [rootgrp.variables[signal].units.strip() for signal in signals.keys()]
    if len(set(units)) > 1:
        print("Signals have different units!{signals.keys()}:{units}")
        return
    else:
        # Get the units for a new signal
        units = units[0]

    # check dimensions 
    dims = [len(rootgrp.variables[signal].dimensions) for signal in signals.keys()]
    if len(set(dims)) > 1:
        print("Signals have different dimensions! {signals.keys()}:{dims}")
        return 
    else:
        dims = dims[0]
        #print(f'Dimensions: {dims}')

    #check dimension keys if the same for all signals
    dim_keys = [None]*dims
    dim_keys_all = [None]*dims
    for ind,x  in enumerate(dim_keys_all):
        dim_keys_all[ind] = [rootgrp.variables[signal].dimensions[ind] for signal in signals.keys()]
        if len(set(dim_keys_all[ind])) > 1:
            print(f"Signals have different {ind} dimension! {dim_keys_all[ind]}")
            return
        else:
            dim_keys[ind] = dim_keys_all[ind][0]
        
    #print(f'Dimension keys: {dim_keys}')


    # define tunits and xunits 
    tunits = rootgrp.variables[dim_keys[0]].units.strip() 
    if len(dim_keys)==2:
        xunits = rootgrp.variables[dim_keys[1]].units.strip()  # X dimension units
    else:
        xunits = None

    # Create the new signal as a linear combination of the given signals
    data = sum([rootgrp.variables[signal][:] * coeff for signal, coeff in signals.items()])
    
           
    #print(f'units: {units}')
    #print(f'shape(data): {np.shape(data)}')

    # Create along name
    long_name = ""
    for signal, coeff in signals.items():
        if coeff < 0:
            long_name += f"-{abs(coeff):.2f}*{signal}"
        else:
            long_name += f"+{coeff:.2f}*{signal}"

    # Remove the leading plus sign from the long_name
    if long_name[0] == "+":
        long_name = long_name[1:]
    #print(f'{new_signal} long_name: {long_name}')
    desc = long_name
    time_key = dim_keys[0]
    if len(dim_keys)==2:
        x_key = dim_keys[1]
    #print(f'time_key: {time_key}, x_key: {x_key}')
    return (data, shape,dim_keys, units, xunits, tunits, desc)


import numpy as np

def check_signal(rootgrp, signal_name):
    """
    Check if a signal exists in a NetCDF file and verifies its properties.
    
    Args:
    - rootgrp: netCDF4.Dataset object representing the NetCDF file
    - signal_name: string representing the name of the signal to check
    
    Returns:
    - A dictionary containing the following keys and values:
      - "exists": a boolean indicating whether the signal exists in the NetCDF file
      - "units": a string representing the units of the signal (or None if the signal does not exist)
      - "long_name": a string representing the long_name of the signal (or None if the signal does not exist)
      - "shape": a tuple representing the shape of the signal (or None if the signal does not exist)
      - "min_value": a float representing the minimum value of the signal (or None if the signal does not exist)
      - "max_value": a float representing the maximum value of the signal (or None if the signal does not exist)
      - "avg_value": a float representing the average value of the signal (or None if the signal does not exist)
    """
    signal_exists = signal_name in rootgrp.variables
    signal_properties = {
        "exists": signal_exists,
        "units": None,
        "long_name": None,
        "shape": None,
        "min_value": None,
        "max_value": None,
        "avg_value": None
    }
    
    if signal_exists:
        signal = rootgrp.variables[signal_name]
        signal_properties["units"] = signal.units.strip()
        signal_properties["long_name"] = signal.long_name.strip()
        signal_properties["shape"] = signal.shape
        signal_values = signal[:]
        signal_properties["min_value"] = np.min(signal_values)
        signal_properties["max_value"] = np.max(signal_values)
        signal_properties["avg_value"] = np.mean(signal_values)
    
    return signal_properties

def addtwo(pulse, signal1, signal2, signal3, rootgrp, DDAname, signals, conversions, newunits):

    if signal1 not in rootgrp.variables or signal2 not in rootgrp.variables:
        return

    nx1 = rootgrp.variables[signal1].shape[1]  # Size of X dimension - signal 1
    nt1 = rootgrp.variables[signal1].shape[0]  # Size of T dimension - signal 1

    nx2 = rootgrp.variables[signal2].shape[1]  # Size of X dimension - signal 2
    nt2 = rootgrp.variables[signal2].shape[0]  # Size of T dimension - signal 2

    if nx1 != nx2 or nt1 != nt2:

        print("Signals incompatible!")
        return

    xunits = rootgrp.variables[rootgrp.variables[signal1].dimensions[1]].units.strip()  # X dimension units
    tunits = rootgrp.variables[rootgrp.variables[signal1].dimensions[0]].units.strip()  # T dimension units
    desc = signals[signal1] + ' + ' + signals[signal2]  # Description
    # New units
    dunits = newunits[signal1]

    # Write signal
    irdat = ppf.ppfwri_irdat(nx1, nt1, refx=-1, reft=-1, user=0, system=0)
    ihdat = ppf.ppfwri_ihdat(dunits, xunits, tunits, "f", "f", "f", desc)

    data = rootgrp.variables[signal1][:] + rootgrp.variables[signal2][:]

    # scale data by conversion factor
    data = [datum * conversions[signal1] for datum in data]

    # Change time vector by 40s to convery to 'real time'
    times = [t + 40.0 for t in rootgrp.variables[rootgrp.variables[signal1].dimensions[0]][:]]

    # Write Signal
    iwdat, ier = ppf.ppfwri(
        pulse, DDAname, signal3, irdat, ihdat, data, rootgrp.variables[rootgrp.variables[signal1].dimensions[1]][:], times
    )

def add_composite_signal(pulse, rootgrp,new_signal, composite, integrate=False):


    (data,shape, dim_keys, units, xunits, tunits, desc) = linear_composite(rootgrp,new_signal, composite)

    data, units = convert_unit(data, units)    
    # for 1D time traces nx=1
    nt, nx = 1, 1
    size = [nt, nx]
    # first position in shape is time dimension, second position in shape is x dimension
    for i, var in enumerate(shape):
        size[i] = var

    [nt, nx] = size  
    if integrate and nx>1:
        print(f'{new_signal} integration along X')  
        dvol = rootgrp.variables['DVOL'][:]
        dvol_units = rootgrp.variables['DVOL'].units.strip() 
        dvol, dvol_units = convert_unit(dvol,dvol_units)
        # upper triangle array with ones
        upper = np.triu(np.ones((nx,nx)))
        # cumulative integration
        data_cum = np.dot(dvol*data,upper)
        data = data_cum[:,-1]
        units = units.replace('/m**3', '', 1).replace('m**-3', '', 1)

        nx=1
    elif integrate and nx==1:
        print(f'{new_signal}: 1D signal - NOT integrating in X dimension')
    # Write a signal
    # refx =-1  a new vector is to be written
    # refx = 0  NO vector is to be written, x_values will be ignored in ppfwri
    if nx > 1:
        refx = -1
        x_values = rootgrp.variables[dim_keys[1]][:]
    else:
        refx = 0
        x_values = None

    # Change time vector by 40s to convery to 'real time'
    times = [t + 40.0 for t in rootgrp.variables[dim_keys[0]][:]]

    # New units
    #dunits = newunits[signal1]
    dunits = units
    # Write signal
    
    irdat = ppf.ppfwri_irdat(nx, nt, refx=refx, reft=-1, user=0, system=0)
    ihdat = ppf.ppfwri_ihdat(dunits, xunits, tunits, "f", "f", "f", desc)

    #data = rootgrp.variables[signal1][:] + rootgrp.variables[signal2][:]

    # scale data by conversion factor
    #data = [datum * conversions[signal1] for datum in data]

    ##Write Signal
    print(f'Writing {new_signal}')
    iwdat, ier = ppf.ppfwri(
        pulse, DDAname, new_signal, irdat, ihdat, data, x_values, times
    )
# subroutine for space cumulative integral integral

def cum_int(pulse, rootgrp, signals_2D):
    '''Cumulated integral of a 2D signal'''

    dvol = rootgrp.variables['DVOL']
    dvol_units = rootgrp.variables['DVOL'].units.strip()

    for key in signals_2D:
        if key not in rootgrp.variables:
            print(key + ' not found in NETCDF file')
            continue
        nt, nx = 1, 1
        size = [nt, nx]
        # first position in shape is time dimension, second position in shape is x dimension
        for i, var in enumerate(rootgrp.variables[key].shape):
            size[i] = var

        [nt, nx] = size
        dunits = rootgrp.variables[key].units.strip()  # Data units
        if '/CM3' in dunits:
            dunits = dunits.replace('/CM3','')
           
        tunits, xunits = '', ''
        txunits = [tunits, xunits]
        for i, var in enumerate(rootgrp.variables[key].dimensions):
            txunits[i] = rootgrp.variables[var].units.strip()
            [tunits, xunits] = txunits

        data = rootgrp.variables[key][:]
        # upper triangle array with ones
        upper = np.triu(np.ones((nx,nx)))
        # cumulative integration
        data_cum = np.dot(dvol*data,upper)

        # Description - Original TRANSP signal name + TRANSP description
        desc = f'{key}:cumulative integral'

        if nx > 1:
            refx = -1
            x_values = rootgrp.variables[rootgrp.variables[key].dimensions[1]][:]
        else:
            refx = 0
            x_values = None

        irdat = ppf.ppfwri_irdat(nx, nt, refx=refx, reft=-1, user=0, system=0)
        ihdat = ppf.ppfwri_ihdat(dunits, xunits, tunits, "f", "f", "f", desc)

        # Change time vector by 40s to conversion to 'real time'
        times = [t + 40.0 for t in rootgrp.variables[rootgrp.variables[key].dimensions[0]][:]]
        # Write Signal
        print(f'writing {signals_2D[key]}')
        iwdat, ier = ppf.ppfwri(pulse, DDAname, signals_2D[key], irdat, ihdat, data_cum, x_values, times)

# Open new PPF
def openPPF(pulse, date, time, comment, status=0):
    ier = ppf.ppfopn(pulse, date, time, comment, status=status)
    if ier != 0:
        print("Unable to open new PPF for writing! Exiting.\n")
        OMFITx.End()
    return None


# Add signals - Need X as first dimension, then T. This is opposite to how they seem to be stored in the NetCDF file.
def save_signal(pulse, DDAname, rootgrp, signals):
    for key in signals:

        if key not in rootgrp.variables:

            print(key + ' not found in NETCDF file')
            continue

        # for 1D time traces nx=1
        nt, nx = 1, 1
        size = [nt, nx]
        # first position in shape is time dimension, second position in shape is x dimension
        for i, var in enumerate(rootgrp.variables[key].shape):
            size[i] = var

        [nt, nx] = size

        dunits = rootgrp.variables[key].units.strip()  # Data units

        tunits, xunits = '', ''
        txunits = [tunits, xunits]
        for i, var in enumerate(rootgrp.variables[key].dimensions):
            txunits[i] = rootgrp.variables[var].units.strip()
        [tunits, xunits] = txunits

        # Description - Original TRANSP signal name + TRANSP description
        desc = f'{key}:{rootgrp.variables[key].long_name}'
        # New units
        dunits = newunits[key]

        # Write a signal
        # refx =-1  a new vector is to be written
        # refx = 0  NO vector is to be written, x_values will be ignored in ppfwri
        if nx > 1:
            refx = -1
            x_values = rootgrp.variables[rootgrp.variables[key].dimensions[1]][:]
        else:
            refx = 0
            x_values = None

        irdat = ppf.ppfwri_irdat(nx, nt, refx=refx, reft=-1, user=0, system=0)
        ihdat = ppf.ppfwri_ihdat(dunits, xunits, tunits, "f", "f", "f", desc)
        data = rootgrp.variables[key][:]
        # scale data by conversion factor
        data = data * conversions[key]

        # Change time vector by 40s to conversion to 'real time'
        times = [t + 40.0 for t in rootgrp.variables[rootgrp.variables[key].dimensions[0]][:]]
        # Write Signal
        print(f'writing {signals[key]}')
        iwdat, ier = ppf.ppfwri(pulse, DDAname, signals[key], irdat, ihdat, data, x_values, times)
    return None


def create_composite(pulse, rootgrp, DDAname, signals, conversions, newunits):
    addtwo(pulse, 'PBI', 'PBTH', 'QNBI', rootgrp, DDAname, signals, conversions, newunits)
    addtwo(pulse, 'SBE', 'SCEV', 'NBS', rootgrp, DDAname, signals, conversions, newunits)
    addtwo(pulse, 'UBPAR_D', 'UBPERP_D', 'WPNB', rootgrp, DDAname, signals, conversions, newunits)
    return None


# Close PPF
def close_PPF(pulse, rootgrp):
    seq, ier = ppf.ppfclo(pulse, "OMFIT", 1)
    if ier != 0:
        print("Unable to close PPF. Exiting.\n")
        OMFITx.End()
    else:
        user=ppf.ppfgid(rw='w')
        mystr = f"{transpcid} was saved to {pulse}/{DDAname}/{user}/{seq}"
        #mystr = f"Sequence number of new PPF is {seq}"
        print(mystr)
        print()

    rootgrp.close()


###############################################################################
def main(pulse,runid,DDAname, signals, conversions, newunits):
    
    
    ier = 0

    #pulse = root['SETTINGS']['EXPERIMENT']['shot']
    #runid = str(root['SETTINGS']['PHYSICS']['TRANSPCID'])
    tr_seq = runid.replace(str(pulse), '')

    # Construct path to results directory for this run
    path = f"/common/transp_shared/Data/result/JET/{pulse:d}/{tr_seq.upper()}/"

    # Get information from NetCDF file
    # Path to netCDF file
    cdfpath = f'{path}{pulse:d}{tr_seq.upper()}.CDF'

    # Check file exists in results directory
    if not os.path.isfile(cdfpath):
        print(f'TRANSP CDF not found in {cdfpath}')
        # Is it a problem with the case
        cdfpath = f'{path}{pulse:d}{tr_seq.lower()}.cdf'
        if not os.path.isfile(cdfpath):

            # Check warehouse directory
            whshot = f'{pulse:07d}'
            whousedir = f'/common/transp_shared/Data/whouse/JET/{whshot[0:4]}___/{whshot}/{tr_seq.lower()}'
            cdfpath = f'{whousedir}/{whshot}{tr_seq.lower()}.cdf'

            if not os.path.isfile(cdfpath):

                # Check if netcdf file in OMFIT
                if 'TRANSP_OUTPUT' in root['OUTPUTS']:
                    print('Using netCDF file stored in OMFIT tree\n')
                    cdfpath = root['OUTPUTS']['TRANSP_OUTPUT']
                else:
                    print('No netCDF file found. Quitting.\n')
                    OMFITx.End()

            else:
                print('Using netCDF file found in warehouse directory.\n')
        else:
            print('Using netCDF file found in results directory.\n')
    else:
        print('Using netCDF file found in results directory.\n')

    # Open existing NetCDF file and read signals from it.
    rootgrp = Dataset(cdfpath, "r", format="NETCDF4")

    # Get user ID to write new PPF under
    #user = MainSettings['SERVER']['CCFE_username']
    user = getpass.getuser()
    user_len = len(user)

    #Set user ID


    ier = ppf.ppfuid(user, "w")

    comment = f"TRANSP inputs to JETTO from run {runid}"

    time = '000000'
    date = '000000'

    openPPF(pulse, date, time, comment, status=0)
    # save signals here
    save_signal(pulse, DDAname, rootgrp, signals)

    create_composite(pulse, rootgrp, DDAname, signals, conversions, newunits)
    cum_int(pulse, rootgrp, signals_to_int)
    
    # WDIA = 1.5*UFASTPP+ UTHRM
    composite = {'UFASTPP':1.5, 'UTHRM':1.0}
    new_signal = 'WDIA'
    add_composite_signal(pulse, rootgrp,new_signal, composite, integrate=True)

    # WMHD = 0.75*UFASTPP + 1.5*(UFASTPA+UPHI) + UTHRM
    composite = {'UFASTPP':0.75,'UFASTPA':1.5,'UPHI':1.5, 'UTHRM':1.0}
    new_signal = 'WMHD'
    add_composite_signal(pulse, rootgrp,new_signal, composite, integrate=True)
    # end of saving signals
    close_PPF(pulse, rootgrp)

    #composite = {'UFASTPP':1.5, 'UTHRM':-1.0,'PLUPLU':1.0}
    #linear_composite(rootgrp,'WDIA', composite)

if __name__ =='__main__':
    
    transpcid=''.join(sys.argv[1:])
    pulse=int(transpcid[:-3])
    runid = transpcid
    print(f'Pulse: {pulse}, transpcid: {transpcid}')
    print(f'user: {getpass.getuser()}')
    main(pulse,runid,DDAname, signals, conversions, newunits)
