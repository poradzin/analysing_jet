from netCDF4 import Dataset
import os.path
import sys
import getpass

transpcid=''.join(sys.argv[1:])
pulse=int(transpcid[:-3])
runid = transpcid
print(f'Pulse: {pulse}, transpcid: {transpcid}')
print(f'user: {getpass.getuser()}')

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
#rootgrp = Dataset(cdfpath, "r", format="NETCDF4")


# open the netCDF4 file
with Dataset(cdfpath, "r", format="NETCDF4") as rootgrp:
    # get a list of all variables in the file
    signals = list(rootgrp.variables.keys())

    # get a list of all unique units for each variable
    units = set([rootgrp.variables[signal].units.strip() for signal in signals])

    # print out the list of unique units
    print("Unique Units:")
    for unit in units:
        print(unit)