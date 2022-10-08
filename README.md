# AncestralPaths

AncestralPaths is a local ancestry painting technique based on genealogies produced by [RELATE](https://myersgroup.github.io/relate/). The method aims to assign 'path' ancestries to segments of chromosomes the decribe the path back through a population history that segment has taken. This code is based on a model of the population history of Europeans described in this ![figure](Model_schematic.pdf), but the code can be adapted to different populaiton histories. Steps of the method as as follows:

1. Simulate chromosomes from the model of European population structure.
2. Run Relate on simulated data.
3. Train a neural network to predict the path a haplotype has taken from a feature of the genealogy relating sample haplotypes.
4. Test the neural network for classification accuracy by population.
5. Apply the classifer to RELATE genealogies constructed from real genomes to paint chromosomes with path ancestries. 

The code present here is for creating a classifier trained on a model of European population structure but the number of samples and their ages simulated can be specified. Once painted chromosomes of real genomes are created, there are add on modules to estimate population size and selection signals along paths and to estimate the time since admixture of individual samples. 

## Requirements

The code requires RELATE to by installed. Follow the installation intructions can be found in the [documentation](https://myersgroup.github.io/relate/). RELATE requires that genomes are phased.

The admixture times analysis require R version 4 or later.  

## Input

Population ages file: To specify the number and ages of each sample to simulate and train over a sample ages file is required. This should be a two column file, tab separated, with the population label and age in units of generations, one line per sample. The population labels should be one of *present_day, Bronze_Age, Neo, Yam, WHG, EHG, CHG, Ana, BAA* (caps-sensitive). These are populations from the model of European population structure. This ![figure](Model_schematic.pdf) describes the relationships of these ancient groups and present day Europeans. Ancient genomes from all of these groups are publically available.

## Usage

To run the simulation, training and testing of a classifier over a speicified number of samples:

```
./pipeline.sh <output_prefix> <population_ages_file> <path_to_relate>
```
where:

- <output_prefix> specifies the prefix for all the output files and output directory
- <population_ages_file> is the file specifying the population and ages of each sample described above
- <path_to_relate> is the relative path the the RELATE module

### Ouput 

Two diretories, *simulated_files* and *relate_files* contain the output tree sequences, VCFs and associated files for simulating and constructing RELATE genealogies.

The simulation produces a *output*.yaml file which describes the simulated demographic structure. 

The model_*output*.h5 is the trained classifier that is used to paint real chromosomes in the next step. 

The results from testing the classifier in the form of confusion matrices are saved. One for each population as *output*_*population*_confusion.txt, and one total confusion matrix produced from pooling all testing data from all populations, *output*_confusion.txt.

An array containing the accuracy scores per population (rows), per tree sequence tested (coloumns) is saved in *output*.kfold.log. 

## Painting

To use the trained classifier on real data first requires tree sequences to be constructed for each chromosome needing to be painted. This can be done in a way that is best for the data at hand, following the guidelines and suggestions in the RELATE documentation. Once tree sequences have been constructed, each chromosome can be painted with the command: 

```
python paint.py 
	-ts <chromosome.trees> 
	-model <path_to_output_files/model.h5> 
	-nn <number_of_nodes> 
	-poplab <poplabels_file> 
	-sample <sample_order> 
	-out <output_prefix> 
	-map <genetic_map> 
```

The <poplabels_file> is the same as that used in the RELATE inference and must have the same population labels as in the <population_ages_file> used for classifier training. The <sample_order> is a file with the names of the samples in the correct order. This can be produced from the VCF file using bcftools or from the .sample file used in RELATE inference. It is needed to determine the correct sample order and should be one sample name per line. The map file is the genetic map file for the corresponding chromosome.  

The output of this is 1. a *output*_painted.npz file. This is an array with dimensions N x S x 3, where N = number of trees and S = number of **haploid** samples. For each tree and each samples there is a path label, the softmax value for that assignment and the right genomic position of the tree intervals across the chromosome. 2. 

The *output*_intervals.npz file is an N x S x 4 array. For each sample and each tree, the interval in bp the tree spans in columns 1 and 2, the path label in column 3, the genetic distance of the right hand edge the tree reaches along the chromosome in column 3 and the softmax value in column 4. This file is used as input to the admixture analysis.

The painting.R script can be used to visualize individual painted haploid chromosomes by passing the *output*_intervals.npz file loaded as an array in R and the haploid sample number to plot.

## Admixture time and fraction analysis

The admixture time analysis estimates both the admixture time and fractions of the admixed populations *Neo, Yam, Bronze_Age* and *present_day*. Run the command. 

```
python admixture.py 
	-poplab <poplabels_file> 
	-out <output_prefix> 
	-samples <sample_order> 
	-ints <list_of_interval_filenames> | -intsf <file_of interval_filenames>
```

The <poplabels_file> and <sample_order> are the same as for painting above. There is a choice between options -ints where a list of intervals filenames are passed or -intsf where a text file containing the interval filenames, one on each line, is passed. 

The output of this command is two files. A *output*.times.npz file and a *output*.fractions.npz file. Both are arrays of S x 2 x 2 dimensions where S is the number of diploid samples. For each diploid sample there are two estimates of the admixture time/fraction, one from each of the two admixing ancestries and an associated standard error.    
 
