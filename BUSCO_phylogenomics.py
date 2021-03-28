#!/usr/bin/env python

# BUSCO_phylogenomics.py
# 2019 Jamie McGowan <jamie.mcgowan@mu.ie>
#
# Utility script to construct species phylogenies using BUSCO results.
# Can perform ML supermatrix or generate datasets for supertree methods.
# Works directly from BUSCO output, as long as the same BUSCO dataset
# has been used for each genome
#
# Dependencies:
#   - BioPython
#   - MUSCLE
#   - trimAL
#   - IQ-TREE
#

import argparse
import multiprocessing as mp
import os
import sys
from time import gmtime, strftime

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# If these programs aren't in $PATH, replace the string below with full
# paths to the programs, including the program name
muscle = "muscle"
iqtree = "iqtree"
trimal = "trimal"


# astral = "astral.jar"

# TODO Add FastTree support

def main():
    parser = argparse.ArgumentParser(description="Perform phylogenomic reconstruction using BUSCOs")
    parser.add_argument("--supermatrix",
                        help="Concatenate alignments of ubuquitious single copy BUSCOs and perform supermatrix "
                             "species phylogeny reconstruction using IQTREE/ML",
                        action="store_true")
    parser.add_argument("--supertree",
                        help="Generate individual ML phylogenies of each BUSCO persent in at least 4 genomes for "
                             "supertree species phylogeny reconstruction with ASTRAL",
                        action="store_true")
    parser.add_argument("-t", "--threads", type=int, help="Number of threads to use", required=True)
    parser.add_argument("-d", "--directory", type=str, help="Directory containing completed BUSCO runs", required=True)
    parser.add_argument("-o", "--output", type=str, help="Output directory to store results", required=True)
    parser.add_argument("-l", "--lineage", type=str, help="Name of lineage used to run BUSCO", required=False)
    parser.add_argument("-psc", "--percent_single_copy", type=float, action="store", dest="psc",
                        help="BUSCOs that are present and single copy in N percent of species will be included in the "
                             "concatenated alignment")
    parser.add_argument("--stop_early",
                        help="Stop pipeline early after generating datasets (before phylogeny inference)",
                        action="store_true")
    args = parser.parse_args()

    start_directory = os.path.abspath(args.directory)
    working_directory = os.path.abspath(args.output)
    threads = int(args.threads)
    supermatrix = args.supermatrix
    supertree = args.supertree
    stop_early = args.stop_early
    lineage = args.lineage

    if args.psc is None:
        percent_single_copy = 100
        print(percent_single_copy)
    else:
        percent_single_copy = float(args.psc)
        print(percent_single_copy)

    if not supermatrix and not supertree:
        print("Error! Please select at least one of '--supermatrix' or '--supertree'")
        sys.exit(1)

    # Check input directory exists
    if os.path.isdir(start_directory):
        os.chdir(start_directory)
    else:
        print("Error! " + start_directory + " is not a directory!")

    # Check if output directory already exists
    if os.path.isdir(working_directory):
        print("Error! " + working_directory + " already exists")
        sys.exit(1)
    else:
        os.mkdir(working_directory)

    if lineage == None:
        lineage = ""

    # TODO check dependencies are installed

    print_message("Starting BUSCO Phylogenomics Pipeline")

    # Scan start directory to identify BUSCO runs (begin with 'run_')
    busco_dirs = []

    for item in os.listdir("."):
        if item[0:4] == "run_":
            if os.path.isdir(item):
                busco_dirs.append(item)

    print("Found " + str(len(busco_dirs)) + " BUSCO runs:")

    for directory in busco_dirs:
        print("\t" + directory)

    print("")

    buscos = {}
    all_species = []

    for directory in busco_dirs:
        os.chdir(start_directory)

        species = directory.split("run_")[1]
        all_species.append(species)

        os.chdir(directory)
        # os.chdir("run_" + lineage) # Issue with BUSCO version >= 4?
        os.chdir("busco_sequences")
        os.chdir("single_copy_busco_sequences")
        
        print(species)

        for busco in os.listdir("."):
            if busco.endswith(".faa"):
                #print(busco)
                busco_name = busco[0:len(busco) - 4]
                record = SeqIO.read(busco, "fasta")
                new_record = SeqRecord(Seq(str(record.seq)), id=species, description="")

                if busco_name not in buscos:
                    buscos[busco_name] = []

                buscos[busco_name].append(new_record)

    print("BUSCO\t # Species Single Copy")
    for busco in buscos:
        print(busco + " " + str(len(buscos[busco])))

    print_message((str(len(buscos))) + " BUSCOs were found")
    print("")

    if supertree:
        print_message("Beginning SUPERTREE Analysis")
        print("")

        # Identify BUSCOs that are present (single copy) in at least 4 species
        four_single_copy = []

        for busco in buscos:
            if len(buscos[busco]) >= 4:
                four_single_copy.append(busco)

        if len(four_single_copy) == 0:
            print_message("0 BUSCOs are present and single copy in at least 4 species")
            # Should break out or quit here
        else:
            print_message(str(len(four_single_copy)) + " BUSCOs are single copy and present in at least 4 species")

            os.chdir(working_directory)
            os.mkdir("proteins_4")
            os.mkdir("alignments_4")
            os.mkdir("trimmed_alignments_4")
            os.mkdir("trees_4")
            os.mkdir("trees_4/iqtree_files")

            print("")

            print_message("Writing protein sequences to: " + os.path.join(working_directory, "proteins_4"))

            for busco in four_single_copy:
                busco_seqs = buscos[busco]

                SeqIO.write(busco_seqs, os.path.join("proteins_4", busco + ".faa"), "fasta")

            print("")
            print_message("Aligning protein sequences using MUSCLE with", threads, "threads to:",
                          os.path.join(working_directory, "alignments_4"))

            mp_commands = []
            for busco in four_single_copy:
                mp_commands.append(
                    [os.path.join("proteins_4", busco + ".faa"), os.path.join("alignments_4", busco + ".aln")])

            pool = mp.Pool(processes=threads)
            results = pool.map(run_muscle, mp_commands)

            print("")
            print_message("Trimming alignments using trimAl (-automated1) with", threads, "threads to: ",
                          os.path.join(working_directory, "trimmed_alignments_4"))

            mp_commands = []

            for busco in four_single_copy:
                mp_commands.append([os.path.join("alignments_4", busco + ".aln"),
                                    os.path.join("trimmed_alignments_4", busco + ".trimmed.aln")])

            pool = mp.Pool(processes=threads)
            results = pool.map(run_trimal, mp_commands)

            print("")
            print_message("Generating phylogenies using IQ-TREE (with model testing) for each BUSCO family with",
                          threads, "threads to:", os.path.join(working_directory, "trees_4"))

            mp_commands = []

            for busco in four_single_copy:
                mp_commands.append([os.path.join("trimmed_alignments_4", busco + ".trimmed.aln")])

            pool = mp.Pool(processes=threads)
            results = pool.map(run_iqtree, mp_commands)

            # Move all IQ-TREE generated files to trees_4 folder
            os.system("mv trimmed_alignments_4/*.treefile trees_4")
            os.system("mv trimmed_alignments_4/*.trimmed.aln.* trees_4/iqtree_files")

            print("")
            print_message("Concatenating all TREEs to: ", os.path.join(working_directory, "ALL.trees"))

            os.chdir(working_directory)
            os.system("cat trees_4/*.treefile > ALL.trees")

            print("")

            print_message("Finished generating dataset for supertree analysis. Use programs such as Astral or CLANN "
                          "to infer species tree from trees_4/ALL.trees")

            print("")

    if supermatrix:
        single_copy_buscos = []
        if args.psc is None:
            print_message("Identifying BUSCOs that are single copy in all " + str(len(all_species)) + " species")

            for busco in buscos:
                if len(buscos[busco]) == len(all_species):
                    single_copy_buscos.append(busco)

            if len(single_copy_buscos) == 0:
                print_message("0 BUSCO families were present and single copy in all species")
                print_message("Exiting")
                sys.exit(0)
            else:
                print(str(len(single_copy_buscos)) + " BUSCOs are single copy in all " + str(len(all_species)) + " species")
        else:
            psc = args.psc
            # Identify BUSCOs that are single copy and present in psc% of species

            for busco in buscos:
                percent_species_with_single_copy = (len(buscos[busco]) / (len(all_species) * 1.0)) * 100

                if percent_species_with_single_copy >= psc:
                    single_copy_buscos.append(busco)

            print(str(len(single_copy_buscos)) + " BUSCOs are single copy in >= " + str(psc) + " of species")

        os.chdir(working_directory)
        os.mkdir("proteins")
        os.mkdir("alignments")
        os.mkdir("trimmed_alignments")

        print("")

        print_message("Writing protein sequences to: " + os.path.join(working_directory, "proteins"))

        for busco in single_copy_buscos:
            busco_seqs = buscos[busco]

            SeqIO.write(busco_seqs, os.path.join(working_directory, "proteins", busco + ".faa"), "fasta")

        print("")
        print_message("Aligning protein sequences using MUSCLE with", threads, "threads to: ",
                      os.path.join(working_directory))

        mp_commands = []

        for busco in single_copy_buscos:
            mp_commands.append([os.path.join(working_directory, "proteins", busco + ".faa"),
                                os.path.join(working_directory, "alignments", busco + ".aln")])

        pool = mp.Pool(processes=threads)
        results = pool.map(run_muscle, mp_commands)

        print("")
        print_message("Trimming alignments using trimAl (-automated1) with", threads, "threads to: ",
                      os.path.join(working_directory, "trimmed_alignments"))

        mp_commands = []

        for busco in single_copy_buscos:
            mp_commands.append([os.path.join(working_directory, "alignments", busco + ".aln"),
                                os.path.join(working_directory, "trimmed_alignments", busco + ".trimmed.aln")])

        pool = mp.Pool(processes=threads)
        results = pool.map(run_trimal, mp_commands)

        print("")
        print_message("Concatenating all trimmed alignments for SUPERMATRIX analysis")

        os.chdir(os.path.join(working_directory, "trimmed_alignments"))
        alignments = {}

        for species in all_species:
            alignments[species] = ""

        # if psc isn't set, or is == 100, we can simple just concatenate alignments
        if args.psc is None:
            for alignment in os.listdir("."):
                for record in SeqIO.parse(alignment, "fasta"):
                    alignments[str(record.id)] += str(record.seq)
        else:
            # We need to check if a species is missing from a family, if so append with "-" to represent missing data
            for alignment in os.listdir("."):
                # Keep track of which species are present and missing
                check_species = all_species[:]

                for record in SeqIO.parse(alignment, "fasta"):
                    alignments[str(record.id)] += str(record.seq)
                    check_species.remove(str(record.id))

                if len(check_species) > 0:
                    # There are missing species, fill with N * "?"
                    seq_len = len(str(record.seq))
                    for species in check_species:
                        alignments[species] += ("?" * seq_len)

        os.chdir(working_directory)
        fo = open("SUPERMATRIX.aln", "w")

        for species in alignments:
            fo.write(">" + species + "\n")
            fo.write(alignments[species] + "\n")

        fo.close()

        print_message("Supermatrix alignment is " + str(len(alignments[species])) + " amino acids in length")

        if stop_early:
            print_message("Stopping early")
            sys.exit(0)

        print_message("Reconstructing species phylogeny using IQ-TREE with model selection from ModelFinder, "
                      "1000 ultrafast bootstrap approximations and 1000 SH-aLRTs: SUPERMATRIX.aln.treefile")
        print("")

        os.system("iqtree -s SUPERMATRIX.aln -bb 1000 -alrt 1000 -nt AUTO -ntmax " + str(threads) + " > /dev/null")

        print("")
        print_message("SUPERMATRIX phylogeny construction complete! See treefile: SUPERMATRIX.aln.treefile")


def run_muscle(io):
    os.system("muscle -in " + io[0] + " -out " + io[1] + " > /dev/null 2>&1")

def run_trimal(io):
    os.system("trimal -in " + io[0] + " -out " + io[1] + " -automated1 ")

def run_iqtree(io):
    os.system("iqtree -s " + io[0] + " > /dev/null 2>&1")

def print_message(*message):
    print(strftime("%d-%m-%Y %H:%M:%S", gmtime()) + "\t" + " ".join(map(str, message)))

if __name__ == "__main__":
    main()
