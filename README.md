# BUSCO_Phylogenomics


<a href="https://jamiemcgowan.ie" target="_blank">Jamie McGowan - 2020</a>

Utility script to construct species phylogenies using BUSCOs. Works directly from BUSCO output and can be used for supermatrix or supertree/coalescent methods.

![BUSCO Phylogenomics pipeline](./pipeline.png)

This pipeline runs directly on the output from BUSCO. Move results directories from each BUSCO run into the same directory. Example structure, where `INPUT_DIRECTORY` is passed to the `-d` parameter of the pipeline:

```
* INPUT_DIRECTORY
	* run_species1
	* run_species2
	* run_species3
	* run_species4
	* run_species5
	* run_species6
	* ........

```


The majority of steps are parallelizable (e.g. family alignments) so running the pipeline with multiple threads should dramatically decrease runtime.

### Usage
	python BUSCO_Phylogenomics.py -d INPUT_DIRECTORY -o OUTPUT_DIRECTORY --supermatrix --threads 20

I reccommend using the `--stop_early` flag which stops the pipeline just after generating the concatenated alignment. Then you can take a look at the concatenated alignment first and manually choose parameters for phylogenetic inference. Similarly, you may want to change the alignment trimming strategy which currently uses the `-automated1` parameter from `trimal`
	
	
### Required parameters
* `-d --directory` input directory containing BUSCO runs
* `-o --output` output directory
* `-t --threads` number of threads to use
* `--supermatrix` and/or `--supertree` choose to run supermatrix and/or supertree methods


### Optional parameters
* `-psc` BUSCO families that are present and single-copy in N% of species will be included in supermatrix analysis (default = 100%). Families that are missing for a species will be replaced with missing characters ("?").
* `--stop_early` stop pipeline early before phylogenetic inference (i.e., for the supermatrix approach this will stop after generating the concatenated alignment). This is **recommended** so you can manually choose your own parameters (e.g., bootstrapping/model selection methods) or manually processing/filtering the alignments further when running IQ-Tree, etc..



### Requirements
* [Python](https://www.python.org/)
* [BioPython](https://biopython.org/)
* [MUSCLE (v5)](https://www.drive5.com/muscle/)
* [trimAl](http://trimal.cgenomics.org/)
* [IQ-TREE](http://www.iqtree.org/)


`muscle`, `trimal` and `iqtree` should be in `$PATH`


### Citation

These scripts were initially written to generate species phylogenies for the following publications:

- [McGowan, J., & Fitzpatrick, D. A. (2020). Recent advances in oomycete genomics. Advances in Genetics. **DOI: 10.1016/bs.adgen.2020.03.001**](https://www.sciencedirect.com/science/article/abs/pii/S0065266020300043)
- [McGowan, J., O’Hanlon, R., Owens, R. A., & Fitzpatrick, D. A. (2020). Comparative Genomic and Proteomic Analyses of Three Widespread Phytophthora Species: Phytophthora chlamydospora, Phytophthora gonapodyides and Phytophthora pseudosyringae. Microorganisms, 8(5), 653. **DOI: 10.3390/microorganisms8050653**](https://www.mdpi.com/2076-2607/8/5/653)


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4320788.svg)](https://doi.org/10.5281/zenodo.4320788)
