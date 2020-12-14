# BUSCO_Phylogenomics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4320788.svg)](https://doi.org/10.5281/zenodo.4320788)


<a href="https://jamiemcgowan.ie" target="_blank">Jamie McGowan - 2020</a>

Utility script to construct species phylogenies using BUSCOs. Works directly from BUSCO output and can be used for supermatrix or supertree methods.

This pipeline runs directly on the output from BUSCO. Move results directories from each BUSCO run (begins with "run_" by default) into the same directory. Example structure, where `INPUT_DIRECTORY` is passed to the `-d` parameter of the pipeline:


* INPUT_DIRECTORY
	* run_species1
	* run_species2
	* run_species3
	* run_species4
	* run_species5
	* run_species6
	* ........


The majority of steps are parallelizable (e.g. family alignments) so running the pipeline with multiple threads leads to a dramatic decrease in runtime.

### Usage
	python BUSCO_Phylogenomics.py -d INPUT_DIRECTORY -o OUTPUT_DIRECTORY --supermatrix --threads 20
	
	
	
### Required parameters:
* `-d --directory` input directory containing BUSCO runs
* `-o --output` output directory
* `-t --threads` number of threads to use
* `--supermatrix` and/or `--supertree` choose to run supermatrix and/or supertree methods


### Optional parameters
* `-psc` BUSCO families that are present and single-copy in N% of species will be included in supermatrix analysis [default = 100%]
* `--stop_early` stop pipeline early before phylogenetic inference (i.e. for the supermatrix approach this will stop after generating the concatenated alignment and for the supertree approach this will stop after generating the trimmed alignments for each BUSCO family)


### Requirements
* [Python](https://www.python.org/)
* [BioPython](https://biopython.org/)
* [MUSCLE](https://www.drive5.com/muscle/)
* [trimAl](http://trimal.cgenomics.org/)
* [IQ-TREE](http://www.iqtree.org/)


`muscle`, `trimal` and `iqtree` should be in `$PATH`


### Pipeline
![BUSCO Phylogenomics pipeline](./pipeline.png)


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4320788.svg)](https://doi.org/10.5281/zenodo.4320788)
