import sys
import argparse
import json
import pathlib
import time

# This script creates a dictionary that holds data.
# The dictionary is then printed via json to the next script
# Example use:
# $ zcat ./PA_chr01.depth.gz | head -n 1000000 | python createDict.py --samples sample.tex --reference reference.tex
# reference.tex contains the index of references in PA_chr01.depth.gz, if the file looks like
# <chromosome name> <position> <reference> <sample1> <sample2> ...
# the reference have index 2.
# In reference.tex and sample.tex each index must be listed on seperatly on the rows

# The dictionary have 4 keys:
# data_dict["chromosome"] -> list of chromosome name for each position
# data_dict["position"] -> list of the postion of each base pair
# data_dict["reference"] -> 2D matrix of size r*n where r is the number of references and n is the number of base pairs
# data_dict["samples"] -> 2D matrix of size s*n where s is the number of samples and n is the number of base pairs

# index_reference is a list that have index of references
# index_sample is a list that have index of samples
def create_dictionary(index_reference, index_sample):
    index_chromosome = 0
    index_position = 1
    data_dict = {"chromosome": [], "position": [], "reference": [], "sample": []}
    for line in sys.stdin:
        data_line = line.split()
        data_dict["chromosome"].append(data_line[index_chromosome])
        data_dict["position"].append(int(data_line[index_position]))
        data_dict["reference"].append([int(data_line[i]) for i in index_reference])
        data_dict["sample"].append([int(data_line[i]) for i in index_sample])

    # Transpose the 2D list
    data_dict["reference"] = list(map(list, zip(*data_dict["reference"])))
    data_dict["sample"] = list(map(list, zip(*data_dict["sample"])))
    return data_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse some arguments")

    parser.add_argument('--sample', metavar='S', type=pathlib.Path, nargs=1, help='index of samples', required=True)
    parser.add_argument('--reference', metavar='R', type=pathlib.Path, nargs=1, help='index of references', required=True)

    # Handle for parsed arguments
    args = parser.parse_args()

    # Create list with index of reference
    with args.reference[0].open() as file:
        # Read the file and then convert string to int
        index_reference = list(map(int, file.readlines()))
    # Create list with index of samples
    with args.sample[0].open() as file:
        # Read the file and then convert string to int
        index_sample = list(map(int, file.readlines()))

    # Create dictionary from the piped data
    data_dict = create_dictionary(index_reference, index_sample)

    print(json.dumps(data_dict))
