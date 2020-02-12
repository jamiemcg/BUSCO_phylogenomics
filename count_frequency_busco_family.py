#!/usr/bin/env python

# count_frequence_busco_family.py
# 2020 Jamie McGowan <jamie.mcgowan@mu.ie>
# 
# Reports how many species a BUSCO family was found to be single copy in
# ONLY LOOKS AT SINGLE COPIES, IGNORES IF PRESENT AS MULTI COPY
# Usage: python count_frequency_busco_family.py busco_working_directory

import os, sys
	
if len(sys.argv) < 2:
	print("Usage: python count_frequency_busco_family.py busco_working_directory")
	sys.exit(0)

working_directory = os.path.abspath(sys.argv[1])

busco_dirs = []

for item in os.listdir("."):
	if item[0:4] == "run_":
		if os.path.isdir(item):
			busco_dirs.append(item)

n_busco_runs = len(busco_dirs)

print("Found " + str(n_busco_runs) + " BUSCO runs")
print("")
print("BUSCO Run\t#Single Copy BUSCOs")

all_species = []
all_buscos = set()
busco_per_species = {}

for directory in busco_dirs:
	os.chdir(working_directory)
	species = directory.split("run_")[1]
	all_species.append(species)
	buscos = []

	os.chdir(directory)
	os.chdir("single_copy_busco_sequences")


	for busco in os.listdir("."):
		if busco.endswith(".faa"):
			busco_name = busco[0:len(busco) - 4]
			buscos.append(busco_name)
			all_buscos.add(busco_name)

	print(species + "\t" + str(len(buscos)))
	busco_per_species[species] = buscos

print("\n\n")

all_buscos_count = {}
for busco in all_buscos:
	all_buscos_count[busco] = 0
	
	for species in all_species:
		if busco in busco_per_species[species]:
			all_buscos_count[busco] += 1

print("BUSCO\t#Species\t%Species")
for busco in all_buscos:
	percent = ((all_buscos_count[busco] / float(n_busco_runs)) * 100)
	percent = "{:.2f}".format(percent)
	print(busco + "\t" + str(all_buscos_count[busco]) + '\t' + percent)


print("\n\n")

# Print presence/absence matrix of all found buscos per species
line = ["BUSCO"] + all_species
print ("\t".join(line))

for busco in all_buscos:
	line = [busco]
	for species in all_species:
		if busco in busco_per_species[species]:
			line.append("Y")
		else:
			line.append("N")

	print ("\t".join(line))

