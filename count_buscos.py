#!/usr/bin/env python

# 2023 Jamie McGowan <jamie.mcgowan@earlham.ac.uk>
# 
# Reports how many species a BUSCO family was found to be single-copy in
# Only looks at single-copies, ignores if present as multi copy or fragmented

# Usage: python count_buscos.py -d BUSCO_runs

import argparse, os, sys
from os import chdir, listdir
from os.path import abspath, basename, isdir, join
from time import gmtime, strftime

def main():
    parser = argparse.ArgumentParser(description="Reports how many species a BUSCO sequences was found to be single-copy in")
    parser.add_argument("-i", "--input", type=str, help="Input directory containing completed BUSCO runs", required=True)
    parser.add_argument("--busco_version_3", action="store_true", help="Flag to indicate that BUSCO version 3 was used (which has slighly different output structure)")
    args = parser.parse_args()

    working_directory = os.path.abspath(args.input)

    print_message("Looking for BUSCO runs in " + working_directory)

    chdir(working_directory)

    busco_samples = []
    busco_sample_names = []

    # output directories of BUSCO version 3 are structured different to those of version 4 and 5
    if args.busco_version_3:      
        for i in listdir("."):
            i = abspath(i)
            if isdir(i):
                for j in listdir(i):
                    k = join(i, j)
                    if isdir(k) and "single_copy_busco_sequences" in k:
                        busco_samples.append(i)
                        busco_sample_names.append(basename(i))
    else:
        for i in listdir("."):
            i = abspath(i)

            if isdir(i):
                for j in listdir(i):
                    j = join(i, j)
                    if isdir(j) and "run_" in j:
                        busco_samples.append(j)
                        busco_sample_names.append(basename(i)) 

    print_message("Found", str(len(busco_samples)), "BUSCO directories:")

    for busco_sample in busco_samples:
        print(busco_sample)

    print()

    if len(busco_samples) == 0:
        print_message("Exiting as didn't find any BUSCO directories")
        sys.exit()
    elif len(busco_samples) == 1:
        print_message("Exiting as only found one BUSCO directory")
        sys.exit()
    else:
        print_message("Counting BUSCO sequences:")

    buscos_per_species = {}
    all_buscos = set()

    for busco_sample_name, busco_sample in zip(busco_sample_names, busco_samples):
        buscos = []
        
        if args.busco_version_3:
            chdir(join(busco_sample, "single_copy_busco_sequences"))
        else:
            chdir(join(busco_sample, "busco_sequences", "single_copy_busco_sequences"))

        for f in listdir("."):
            if f.endswith(".faa"):
                busco_name = f.rstrip(".faa")
                buscos.append(busco_name)
                all_buscos.add(busco_name)
        
        buscos_per_species[busco_sample_name] = buscos

        chdir(working_directory)    

    all_buscos_count = {}
    for busco in all_buscos:
        all_buscos_count[busco] = 0

        for busco_sample_name in busco_sample_names:
            if busco in buscos_per_species[busco_sample_name]:
                all_buscos_count[busco] += 1

    print("Name\tComplete and single-copy BUSCOs")

    for busco_sample_name in busco_sample_names:
        print(busco_sample_name, len(buscos_per_species[busco_sample_name]))

    print()
    print("BUSCO\t#Species\t%Species")

    present_in_all_species = []

    busco_counts = []

    for busco in all_buscos:
        count = all_buscos_count[busco]

        if count == len(busco_sample_names):
            present_in_all_species.append(busco)

        percent = ((count / float(len(busco_sample_names))) * 100)
        percent = "{:.2f}".format(percent)

        busco_counts.append([busco, count, percent])

    busco_counts.sort(key = lambda x:x[1], reverse = True)

    for i in busco_counts:
        print("\t".join(map(str, i)))


    print()

    if len(present_in_all_species) == 0:
        print_message("No BUSCO sequences were found to be complete and single-copy in all species")
    else:
        print_message(len(present_in_all_species), "BUSCO sequences were found to be complete and single-copy in all species")
        print(present_in_all_species)

    print()
    print_message("Printing presence/absence matrix of all complete single-copy BUSCOs (!NOTE THIS EXCLUDES DUPLICATED/FRAGMENTED BUSCOS!)")
    
    for busco in all_buscos:
        line = [busco]

        for busco_sample_name in busco_sample_names:
            if busco in buscos_per_species[busco_sample_name]:
                line.append("Y")
            else:
                line.append("N")

        print("\t".join(line))

    print()
    print_message("Complete")


def print_message(*message):
    print(strftime("%d-%m-%Y %H:%M:%S", gmtime()) + "\t" + " ".join(map(str, message)))


if __name__ == "__main__":
    main()
