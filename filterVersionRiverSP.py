import json
import os

import pandas as pd

def filterVersionRiverSP(directories):
    """
    Reads in all filenames, sorts them, and retrieves the best version/
    counter for each file (e.g. PGC0 over PIC0, PIC0_03 over PIC0_02).
    Writes out a json with the filenames for filtering upon read in
    to the /nas/cee-water/cjgleason/data/SWOT/ directory.

    Parameters:
    directories (list): List of directories to search for best files.
    
    Returns:
    None
    """
    
    # Outpath for json
    outpath = '/nas/cee-water/cjgleason/data/SWOT/'
    
    # List to store all .shp file paths
    files = []
    
    # Get all .shp file names from directories
    for directory in directories:
        for file in os.listdir(directory):
            if file.endswith(".shp"):
                files.append(file)
    
    print(f"There are {str(len(files))} original .shp files in directory.")
    
    # Make DataFrame of filenames
    granules = pd.DataFrame({'files': files})
    granules['cycle'] = granules['files'].str.slice(25, 28)
    granules['pass'] = granules['files'].str.slice(29, 32)  
    granules['version'] = granules['files'].str.slice(-11, -7)
    granules['counter'] = granules['files'].str.slice(-6, -4)     

    # Sort the files
    granules = granules.sort_values(by=['cycle', 'pass', 'version', 'counter'],
                                    ascending=[True, True, True, False])    

    # Keep only the best version of each granule
    granules = granules.drop_duplicates(subset=['cycle', 'pass'],
                                        keep='first')    

    # Extract the file names of files passing the test
    best_files = list(granules['files'])
    
    print(f"There are {str(len(best_files))} best .shp files in directory.")
    
    # Extract base names (file name without extensions) from the list
    base_names = set(os.path.splitext(os.path.basename(file))[0] for file in best_files)
    
    all_best_files = []
    # Loop over the directories to find all files with matching base names
    for directory in directories:
        for file in os.listdir(directory):
            # Get the file's base name
            base_name, extension = os.path.splitext(os.path.basename(file))

            # If the base name is in the list of best files,
            # append that file to the list of files to keep
            if base_name in base_names:
                all_best_files.append(file)
    
        # Split filepath for naming json
        pieces = directory.split('/')

        # Write out best files as json
        with open(os.path.join(outpath, pieces[6] + '_' + pieces[7] +
                               '_filtered.json'), 'w', encoding='utf-8') as f:
            json.dump(all_best_files, f)

        print(f"Wrote out the unique and most recently processed {str(len(all_best_files))} files to {outpath}{pieces[6]}_filtered.json")

dirs = ['/nas/cee-water/cjgleason/data/SWOT/RiverSP_v2_0_20240624/GR/Reach']
filterVersionRiverSP(directories=dirs)