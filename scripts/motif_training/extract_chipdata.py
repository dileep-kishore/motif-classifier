"""Script to extract ChIP-seq data"""

import pandas as pd

def read_chipdata(fname, sheetnum=1):
    """Reading binding data from excel sheet"""
    bind_data = pd.read_excel(fname, sheetname=sheetnum, skiprows=[0])
    return bind_data

def parse_chipdata(chip_data):
    """Parse chip_data to extract position information"""
    data_cols = ['Symbol', 'EcoCyc Locus', 'Coverage', 'Type', 'Position']
    return chip_data[data_cols]

def sort_chipdata(chip_data):
    """Sorts the chip_data table based on coverage"""
    sorted_chipdata = chip_data.sort_values('Coverage', ascending=False, inplace=False)
    return sorted_chipdata

def extract_chipdata(chip_file, n, outfile):
    """Write csv file with sorted chip data and return table"""
    data = read_chipdata(chip_file)
    parsed_data = parse_chipdata(data)
    sorted_data = sort_chipdata(parsed_data)
    # Delete rows with same positions
    sorted_data = sorted_data.groupby((sorted_data["Position"] != sorted_data["Position"].shift()).cumsum().values).first()
    if n == 'all':
        sorted_data.to_csv(outfile, index=False)
        return sorted_data
    selected_data = sorted_data.head(n)
    selected_data.to_csv(outfile, index=False)
    return selected_data
