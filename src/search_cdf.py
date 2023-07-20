import os
import sys
import re
import netCDF4 as nc
import argparse

def search_transp(search_string, file_path=None):
    if not os.path.exists(file_path):
        print(f"NetCDF file '{file_path}' not found.")
        return

    try:
        # Open NetCDF file
        dataset = nc.Dataset(file_path)

        # Search for the string in keys and field descriptions
        found_items = []
        for key in dataset.variables.keys():
            variable = dataset.variables[key]
            long_name = variable.long_name if hasattr(variable, 'long_name') else ""
            units = variable.units if hasattr(variable, 'units') else "N/A"
            dimensions = variable.dimensions
            if search_string.lower() in key.lower() or search_string.lower() in long_name.lower():
                found_items.append((key, units, long_name, dimensions))

        # Print the matching keys, descriptions, units and dimensions
        if found_items:
            print(f"Matching items for '{search_string}':\n")
            for item in found_items:
                print(f"Key: {item[0]}\nUnits: {item[1]}\nDescription: {item[2]}\nDimensions: {item[3]}\n")
        else:
            print(f"No matching items found for '{search_string}'.")
    except Exception as e:
        print(f"An error occurred while reading the NetCDF file: {str(e)}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Search within NetCDF files')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--transp', metavar='runID', help="TRANSP run ID to search in the default path. Example: '--transp 99972V01'")
    group.add_argument('--file', metavar='filepath', help='The file path to search')

    args = parser.parse_args()

    if args.transp:
        base_dir = '/common/transp_shared/Data/result/JET'
        pulse_dir = args.transp[:-3]
        run_dir = args.transp[-3:]
        file_path = os.path.join(base_dir, pulse_dir, run_dir, f"{args.transp}.CDF")
    else:
        file_path = args.file

    while True:
        search_string = input("Enter the string to search (press Enter to exit): ")
        if search_string:  # proceed only if the search string is not empty
            search_transp(search_string, file_path)
        else:
            break

