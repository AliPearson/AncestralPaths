# AncestralPaths

AncestralPaths is a local ancestry painting technique based on genealogies produced by [RELATE](https://myersgroup.github.io/relate/). The method aims to assign 'path' ancestries to segments of chromosomes the decribe the path back through a population history that segment has taken. Steps of the method as as follows:

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

#Population ages file:
To specify the number and ages of each sample to simulate and train over a sample ages file is required. This should be a two column file, tab separated, with the population label and age in units of generations, one line per sample. The population labels should be one of *present_day, Bronze_Age, Neo, Yam, WHG, EHG, CHG, Ana, BAA*, caps-sensitive. The are populations from the model of European population structure. This ![figure](Model_schematic.pdf) describes the relationships of these ancient groups and present day Europeans. Ancient genomes from all of these groups are publically available.

## Usage

The run the simulation, through training and testing of a classifier over a speicified number of samples:

```
./pipeline.sh <output_prefix> <population_ages_file> <path_to_relate>
```
where:

- <output_prefix> specifies the prefix for all the output files and output directory
- <population_ages_file> is the file specifying the population and ages of each sample described above
- <path_to_relate> is the relative path the the RELATE module

### Ouput 

Two diretories, *simulated_files* and *relate_files* contain the output tree sequences, VCFs and associated files for simulating and constructing RELATE genealogies. 

The model_*output*.h5 is the trained classifier that is used to paint real chromosomes in the next step. 

The results from testing the classifier in the form of confusion matrices are saved. One for each population as *output*_*population*_confusion.txt, and one total confusion matrix when pooling all testing data from all populations, *output*_confusion.txt.

An array containing the accuracy scores per population (rows), per tree sequence tested (coloumns) is saved in *output*.kfold.log. 

#Painting



 
