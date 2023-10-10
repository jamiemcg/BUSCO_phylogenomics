#!/usr/bin/env python
# 
# BUSCO_phylogenomics
# Utility script to construct species phylogenies using BUSCO results.
# Assumes the same BUSCO dataset has been used on each genome

# 2023 Jamie McGowan <jamie.mcgowan@earlham.ac.uk>
# https://github.com/jamiemcg/BUSCO_phylogenomics


import argparse
import multiprocessing as mp
from os import listdir, chdir, mkdir, system
from os.path import abspath, basename, isdir, join
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from time import gmtime, strftime

def main():
    parser = argparse.ArgumentParser(description="Perform phylogenomic reconstruction using BUSCO sequences")

    parser.add_argument("-i", "--input", type=str, help="Input directory containing completed BUSCO runs", required=True)
    parser.add_argument("-o", "--output", type=str, help="Output directory to store results", required=True)
    parser.add_argument("-t", "--threads", type=int, help="Number of threads to use", required=True)
    parser.add_argument("--supermatrix_only", help="Don't generate gene trees", action="store_true")
    parser.add_argument("--gene_trees_only", help="Don't perform supermatrix analysis", action="store_true")
    parser.add_argument("--nt", help="Align nucleotide sequences instead of amino acid sequences", action="store_true")  
    parser.add_argument("-psc", "--percent_single_copy", type=float, action="store", dest="psc", default=100.0,
                        help="BUSCO presence cut-off. BUSCOs that are complete and single-copy in at least [-psc] percent of species will be included in the contatenated alignment [default=100.0]")
    parser.add_argument("--trimal_strategy", type=str, action="store", dest="trimal_strategy", default="automated1",
                        help="trimal trimming strategy (automated1, gappyout, strict, strictplus) [default=automated1]")
    parser.add_argument("--missing_character", type=str, action="store", dest="missing_character", help="Character to represent missing data [default='?']", default="?")
    parser.add_argument("--gene_tree_program", type=str, action="store", dest="gene_tree_program", default="fasttree", help="Program to use to generate gene trees (fasttree or iqtree) [default=fasttree]")
    parser.add_argument("--busco_version_3", action="store_true", help="Flag to indicate that BUSCO version 3 was used (which has slighly different output structure)")
    
    args = parser.parse_args()

    print_message("Starting BUSCO Phylogenomics Pipeline")
    print_message("User provided arguments:", sys.argv)
    print_message("Parsed arguments:", vars(args))

    input_directory = abspath(args.input)
    working_directory = abspath(args.output)
    threads = args.threads

    # Check if the input directory exists
    if not isdir(input_directory):
        print_message("ERROR. Input BUSCO directory", input_directory, "not found")
        sys.exit()

    # Check if the output directory already exists
    if isdir(working_directory):
        print_message("ERROR. Output directory", working_directory, "already exists")
        sys.exit()

    # Check trimal parameter
    trimal_strategy = args.trimal_strategy.lower()

    if trimal_strategy not in ["automated1", "gappyout", "strict", "strictplus"]:
        print_message("ERROR. Didn't understand --trimal_strategy parameter. Accepted options = automated1, gappyout, strict, strictplus")
        sys.exit()

    trimal_strategy = "-" + trimal_strategy

    gene_tree_program = args.gene_tree_program.lower()

    print(gene_tree_program)
    if gene_tree_program not in ["fasttree", "iqtree"]:
        print_message("ERROR. Didn't understand --gene_tree_program parameter. Accepted options = fasttree, iqtree")
        sys.exit()

    if args.supermatrix_only and args.gene_trees_only:
        print_message("ERROR. You cannot select both --supermatrix_only and --gene_trees_only")
        sys.exit()

    if args.nt:
        sequence_file_extension = ".fna"
        sequence_type = "nucleotide"
    else:
        sequence_file_extension = ".faa"
        sequence_type = "protein"

    print_message("Looking for BUSCO runs in", input_directory)

    busco_samples = []
    busco_sample_names = []

    chdir(input_directory)

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


    if len(busco_samples) == 0:
        print_message("ERROR. Didn't find any BUSCO directories")
        sys.exit()

    print_message("Found", len(busco_samples), "BUSCO directories:")
    for busco_sample_name, busco_sample in zip(busco_sample_names, busco_samples):
        print(busco_sample_name, busco_sample)

    print()
    print_message("Identifying complete and single-copy BUSCO sequences")

    all_buscos = set()
    buscos = {}
    buscos_per_species = {}

    # Identify single-copy BUSCOs per species

    for busco_sample_name, busco_sample in zip(busco_sample_names, busco_samples):
        buscos_per_species[busco_sample_name] = []

        if args.busco_version_3:
            chdir(join(busco_sample, "single_copy_busco_sequences"))
        else:
            chdir(join(busco_sample, "busco_sequences", "single_copy_busco_sequences"))

        for f in listdir("."):
            if f.endswith(sequence_file_extension):
                busco_name = f.rstrip(sequence_file_extension)
                all_buscos.add(busco_name)
                if busco_name not in buscos:
                    buscos[busco_name] = []

                # TODO This slows things down a bit as it reads all sequences even if they aren't going to be used later
                record = SeqIO.read(f, "fasta")
                new_record = SeqRecord(Seq(str(record.seq)), id=busco_sample_name, description="")
                
                buscos_per_species[busco_sample_name].append(busco_name)
                buscos[busco_name].append(new_record)

    print("Name\tComplete and single-copy BUSCOs:")
    for busco_sample_name in busco_sample_names:
        print(busco_sample_name, len(buscos_per_species[busco_sample_name]))
                
    print()
    print_message(len(all_buscos), "unique BUSCO sequences considered")
    print()

    print("BUSCO\tNumber of species (complete and single-copy only)")
    for busco in buscos:
        print(busco + "\t" + str(len(buscos[busco])))

    print()

    mkdir(working_directory)
    chdir(working_directory)

    if not args.gene_trees_only:

        if args.psc == 100.0:
            print_message("Identifying BUSCOs that are complete and single-copy in all species")
        else:
            print_message("Identifying BUSCOs that are complete and single-copy in at least", args.psc, "percent of species")

        single_copy_buscos = []

        for busco in buscos:
            percent_present = (len(buscos[busco]) / len(busco_sample_names) * 100)

            if percent_present >= args.psc:
                single_copy_buscos.append(busco)

        if len(single_copy_buscos) == 0:
            if args.psc == 100.0:
                print_message("ERROR. Didn't identify any BUSCO sequences that are complete and single-copy in all species")
            else:
                print_message("ERROR. Didn't identify any BUSCO sequences that are complete and single-copy in at least", args.psc, "percent of species")
            print_message("You may want to adjust the --percent_single_copy parameter to allow for greater amounts of missing data if your dataset is patchy")
            sys.exit()

        if args.psc == 100.0:
            print_message("Identified", len(single_copy_buscos), "BUSCO sequences that are complete and single-copy in all species:")
        else:
            print_message("Identified", len(single_copy_buscos), "BUSCO sequences that are complete and single-copy in at least", args.psc, "percent of species:")
        print(",".join(single_copy_buscos))
        print()

        mkdir("supermatrix")
        chdir("supermatrix")
        
        mkdir("sequences")
        print_message("Writing " + sequence_type + " sequences to", join(working_directory, "supermatrix", "sequences"))

        for busco in single_copy_buscos:
            busco_records = buscos[busco]
            SeqIO.write(busco_records, join(working_directory, "supermatrix", "sequences", busco + sequence_file_extension), "fasta")

        mkdir("alignments")
        print_message("Aligning " + sequence_type + " sequences using MUSCLE with", threads, "parallel jobs to:", join(working_directory, "supermatrix", "alignments"))

        mp_commands = []
        for busco in single_copy_buscos:
            mp_commands.append([join(working_directory, "supermatrix", "sequences", busco + sequence_file_extension),
                                join(working_directory, "supermatrix", "alignments", busco + ".aln")])

        pool = mp.Pool(processes=threads)
        results = pool.map(run_muscle, mp_commands)

        pool.close()
        pool.join()


        mkdir("trimmed_alignments")
        print_message("Trimming alignments using trimAL (" + trimal_strategy + ") with", threads, "parallel jobs to:", join(working_directory, "supermatrix", "trimmed_alignments"))

        mp_commands = []
        for busco in single_copy_buscos:
            mp_commands.append([join(working_directory, "supermatrix", "alignments", busco + ".aln"),
                                join(working_directory, "supermatrix", "trimmed_alignments", busco + ".trimmed.aln")])

        pool = mp.Pool(processes=threads)
        results = pool.map(run_trimal, mp_commands)

        pool.close()
        pool.join()

        print_message("Concatenating all trimmed alignments")

        chdir("trimmed_alignments")

        alignments = {}
        partitions = []
        for busco_sample_name in busco_sample_names:
            alignments[busco_sample_name] = ""

        start = 1 # tracking lengths for partitions file

        # If psc == 100; we can just concatenate alignments (no missing data)
        if args.psc == 100:
            for alignment in listdir("."):
                if alignment.endswith(".trimmed.aln"):
                    for record in SeqIO.parse(alignment, "fasta"):
                        alignments[str(record.id)] += str(record.seq)
                    partitions.append([alignment.replace(".trimmed.aln", ""), start, start + len(str(record.seq)) - 1])
                    start += len(record.seq)
        else:
            # we need to handle missing data here
            for alignment in listdir("."):
                if alignment.endswith(".trimmed.aln"):
                    check_samples = busco_sample_names[:]

                    for record in SeqIO.parse(alignment, "fasta"):
                        alignments[str(record.id)] += str(record.seq)
                        check_samples.remove(str(record.id))

                    partitions.append([alignment.replace(".trimmed.aln", ""), start, start + len(str(record.seq)) - 1])
                    start += len(record.seq)
                    
                    if len(check_samples) > 0:
                        # This means some species were missing this busco, fill alignment with missing character ("?" is default)
                        seq_len = len(str(record.seq))
                        for sample in check_samples:
                            alignments[sample] += (args.missing_character * seq_len)
        
        chdir(working_directory)
        chdir("supermatrix")

        alignment_records = []
        for sample in alignments:
            record = SeqRecord(Seq(alignments[sample]), id = sample, description = "")
            alignment_records.append(record)

        print_message("Writing concatenated supermatrix alignment to fasta format:", join(working_directory, "supermatrix", "SUPERMATRIX.fasta"))
        SeqIO.write(alignment_records, "SUPERMATRIX.fasta", "fasta")

        print_message("Writing concatenated supermatrix alignment to phylip format:", join(working_directory, "supermatrix", "SUPERMATRIX.phylip"))
        SeqIO.write(alignment_records, "SUPERMATRIX.phylip", "phylip-relaxed")

        print_message("Writing NEXUS partitions file to:", join(working_directory,  "supermatrix", "SUPERMATRIX.partitions.nex"))
        fo = open(join(working_directory,  "supermatrix", "SUPERMATRIX.partitions.nex"), "w")

        fo.write("#nexus\n")
        fo.write("begin sets;\n")

        for p in partitions:
            fo.write("\tcharset " + p[0] + " = SUPERMATRIX.phylip: " + str(p[1]) + "-" + str(p[2]) + ";\n")
        
        fo.write("end;\n")

        fo.close()

        print_message("Supermatrix alignment is", str(len(alignments[sample])), "amino acids in length")

        print()

    if not args.supermatrix_only:
        print_message("Identifying BUSCOs that are complete and single-copy in at least 4 species")

        single_copy_buscos_4_species = []

        chdir(working_directory)

        for busco in buscos:
            if len(buscos[busco]) >= 4:
                single_copy_buscos_4_species.append(busco)
                # print(busco, len(buscos[busco]))

        print_message("Identified", len(single_copy_buscos_4_species), "BUSCO sequences that are complete and single-copy in at least 4 species:")
        print(",".join(single_copy_buscos_4_species))

        mkdir("gene_trees_single_copy")
        chdir("gene_trees_single_copy")

        print()

        mkdir("sequences_4")
        print_message("Writing " + sequence_type + " sequences to", join(working_directory, "gene_trees", "sequences_4"))

        for busco in single_copy_buscos_4_species:
            busco_records = buscos[busco]
            SeqIO.write(busco_records, join(working_directory, "gene_trees_single_copy", "sequences_4", busco + sequence_file_extension), "fasta")

        mkdir("alignments_4")
        print_message("Aligning " + sequence_type + " sequences using MUSCLE with", threads, "parallel jobs to:", join(working_directory, "gene_trees_single_copy", "alignments_4"))

        mp_commands = []
        for busco in single_copy_buscos_4_species:
            mp_commands.append([join(working_directory, "gene_trees_single_copy", "sequences_4", busco + sequence_file_extension),
                                join(working_directory, "gene_trees_single_copy", "alignments_4", busco + ".aln")])
        
        pool = mp.Pool(processes=threads)
        results = pool.map(run_muscle, mp_commands)

        mkdir("trimmed_alignments_4")
        print_message("Trimming alignments using trimAL (" + trimal_strategy + ") with", threads, "parallel jobs to:", join(working_directory, "gene_trees_single_copy", "trimmed_alignments_4"))

        mp_commands = []
        for busco in single_copy_buscos_4_species:
            mp_commands.append([join(working_directory, "gene_trees_single_copy", "alignments_4", busco + ".aln"),
                                join(working_directory, "gene_trees_single_copy", "trimmed_alignments_4", busco + ".trimmed.aln")])

        pool = mp.Pool(processes=threads)
        results = pool.map(run_trimal, mp_commands)

        mkdir("trees_4")

        if gene_tree_program == "fasttree":
            
            mp_commands = []
            for busco in single_copy_buscos_4_species:
                mp_commands.append([join(working_directory, "gene_trees_single_copy", "trimmed_alignments_4", busco + ".trimmed.aln"),
                                    join(working_directory, "gene_trees_single_copy", "trees_4", busco + ".tree")])

            print_message("Generating gene trees using fasttree with", threads, "parallel jobs to:", join(working_directory, "gene_trees_single_copy", "trees_4"))

            pool = mp.Pool(processes=threads)
            results = pool.map(run_fasttree, mp_commands)
            
            concatenate_commant = "cat " + join(working_directory, "gene_trees_single_copy", "trees_4", "*.tree") + " > " + join(working_directory, "gene_trees_single_copy", "ALL.tree")
        elif gene_tree_program == "iqtree":
            
            mp_commands = []
            for busco in single_copy_buscos_4_species:
                mp_commands.append([join(working_directory, "gene_trees_single_copy", "trimmed_alignments_4", busco + ".trimmed.aln"),
                                    join(working_directory, "gene_trees_single_copy", "trees_4", busco)])

            print_message("Generating gene trees using iqtree with", threads, "parallel jobs to:", join(working_directory, "gene_trees_single_copy", "trees_4"))

            pool = mp.Pool(processes=threads)
            results = pool.map(run_iqtree, mp_commands)

            concatenate_commant = "cat " + join(working_directory, "gene_trees_single_copy", "trees_4", "*.treefile") + " > " + join(working_directory, "gene_trees_single_copy", "ALL.tree")


        print_message("Concatenating all", len(single_copy_buscos_4_species), "gene trees to:", join(working_directory, "gene_trees_single_copy", "ALL.tree"))

        system(concatenate_commant)

    print_message("Done")

def run_muscle(io):
    system("muscle -threads 1 -align " + io[0] + " -output " + io[1] + " > /dev/null 2>&1")

def run_trimal(io):
    system("trimal -in " + io[0] + " -out " + io[1] + " -automated1 ")

def run_iqtree(io):
    system("iqtree --quiet -T 1 -s " + io[0] + " --prefix " + io[1])

def run_fasttree(io):
    system("fasttree " + io[0] + " > " + io[1] + " 2> /dev/null")

def print_message(*message):
    print(strftime("%d-%m-%Y %H:%M:%S", gmtime()) + "\t" + " ".join(map(str, message)))

if __name__ == "__main__":
    main()
