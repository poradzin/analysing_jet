import os
import sys
import re
import netCDF4 as nc

def search_transp(pulse_number, run_id, search_string):
    # Directory structure
    base_dir = '/common/transp_shared/Data/result/JET'
    pulse_dir = os.path.join(base_dir, str(pulse_number))
    run_dir = os.path.join(pulse_dir, run_id)

    # NetCDF file name
    nc_file_1 = os.path.join(run_dir, str(pulse_number) + run_id + '.CDF')
    nc_file_2 = os.path.join(run_dir, run_id + '.CDF')

    if os.path.exists(nc_file_1):
        nc_file = nc_file_1
    elif os.path.exists(nc_file_2):
        print(f"NetCDF file '{nc_file_1}' not found.")
        nc_file = nc_file_2
    else:
        print(f"NetCDF file '{nc_file_1}' not found.")
        print(f"NetCDF file {nc_file_2} not found.")
        return
    try:
        # Open NetCDF file
        dataset = nc.Dataset(nc_file)

        # Search for the string in keys and field descriptions
        found_items = []
        for key in dataset.variables.keys():
            variable = dataset.variables[key]
            if search_string.lower() in key.lower() or search_string.lower() in variable.long_name.lower():
                units = variable.units if hasattr(variable, 'units') else "N/A"
                found_items.append((key, units, variable.long_name))

        # Print the matching keys, descriptions, and units
        if found_items:
            print(f"Matching items for '{search_string}':")
            for item in found_items:
                print(f"Key: {item[0]}\nUnits: {item[1]}\nDescription: {item[2]}\n")
        else:
            print(f"No matching items found for '{search_string}'.")
    except Exception as e:
        print(f"An error occurred while reading the NetCDF file: {str(e)}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python search_transp.py pulse_number run_id")
        sys.exit(1)

    pulse_number = int(sys.argv[1])
    run_id = sys.argv[2]
    print(f'Pulse no.: {pulse_number}')
    print(f'run_id:    {run_id}')
    while True:
        search_string = input("Enter the string to search (press Enter to exit): ")
        if search_string:  # proceed only if the search string is not empty
            search_transp(pulse_number, run_id,search_string)
        else:
            break    


