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


# Integrated signals to be saved

signals_to_int = {
    'EHEAT': 'EPIN',
    'IHEAT': 'IPIN',
    'PRAD'  : 'PRIN',
}

units_to_int = {
    'EHEAT': 'W',
    'IHEAT' : 'W',
    'PRAD':'W',
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
    'TAUA1': 'TE1',
    'EHEAT': 'EPOW',
    'IHEAT': 'IPOW',
    'PRAD' : 'PRAD',
    'UFASTPP': 'UFPP',
    'UFASTPA': 'UFPA',
    'UTHRM' : 'UTH',
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
    'EHEAT': 1.0e6,
    'IHEAT': 1.0e6,
    'PRAD' : 1.0e6,
    'UFASTPP': 1.0e6,
    'UFASTPA': 1.0e6,
    'UTHRM' :  1.0e6,
}

newunits = {
    'BZXR': 'TESLA*m',
    'RMNMP': 'm',
    'RMJMP': 'm',
    'CONDE': 'm**2/s',
    'CONDI': 'm**2/s',
    'CONDEF': 'm**2/s',
    'NEUTT': 'W',
    'OMEGDATA':'rad/sec',
    'TAUEA': 's',
    'TAUEE': 's',
    'TAUES': 's',
    'TAUA1': 's',
    'EHEAT': 'W*m**-3',
    'IHEAT': 'W*m**-3',
    'PRAD' : 'W*m**-3',
    'UFASTPP': 'W*m**-3',
    'UFASTPA': 'W*m**-3',
    'UTHRM' : 'W*m**-3',
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
            dunits = dunits[:-4]

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
    #if user_len < 1:
    #    OMFITx.End()

    #cum_int(pulse, rootgrp, signals_to_int)
    #Set user ID
    ier = ppf.ppfuid(user, "w")

    comment = f"TRANSP inputs to JETTO from run {runid}"

    time = '000000'
    date = '000000'

    openPPF(pulse, date, time, comment, status=0)

    save_signal(pulse, DDAname, rootgrp, signals)

    create_composite(pulse, rootgrp, DDAname, signals, conversions, newunits)
    cum_int(pulse, rootgrp, signals_to_int)
    close_PPF(pulse, rootgrp)

if __name__ =='__main__':
    
    #for idx, arg in enumerate(sys.argv):
    #    print(f'Argument #{idx} is {arg}. type: {type(arg)}')
    transpcid=''.join(sys.argv[1:])
    pulse=int(transpcid[:-3])
    runid = transpcid
    print(f'Pulse: {pulse}, transpcid: {transpcid}')
    print(f'user: {getpass.getuser()}')
    main(pulse,runid,DDAname, signals, conversions, newunits)
