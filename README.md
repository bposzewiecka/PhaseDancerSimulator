# Tree-Seg-Dup

Simulations 

##  Dependencies

To run Tree-Seg-Dup you should have above software installed:

* python3
* [pbsim2](https://github.com/yukiteruono/pbsim2)
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* bedtools

## Configuration

### Step 1: Installation of the required software

The following dependencies should be installed: pbsim2, snakemake, bedtools.

### Step 2: Cloning the Tree-Seg-Dup repository

To clone Tree-Seg-Dup repository, the following line should be executed:

```
git clone https://github.com/bposzewiecka/tree-seg-dup.git
```

### Step 4: Creating the configuration file

Configuration file **config.yaml** list parameters of each simulation(s) and coordinates of regions to simulate from.

Each entry of *simulations* dictionary lists parameters of simulation(s).

| Property | Description | Value |
|---|---|---|
| topology | | See: topologies |  
| simulations-number | | integer |
| reference | | string - key from *references* dictionary |
| muation-rates | | list of probabilities |
| chemistry | | |
| coverages |     | list of integers |
| accuracies | | list fo numbers from range [70-100] |
| length-mean | | integer |
| length-sd | | integer |
| type | | leaves, all  | 

#### Topologies

* flat number_of_leaves, 
* bifurcating list_of_numbers, 
* cascading number_of_leaves, 
* random number_of_leaves 

| Property | Description | Value |
|---|---|---|
| region | | string |  
| reference | | string |

The following *yaml* file can be created to simulate duplication of two different topology types. First set of parameters *myflat* simulate 
one topology tree with root and 10 leaves. Second set of parameters *myrandom* simulate 4 random topology trees with 10 leaves.

In both simulations the same input region is used. For each simulation *reference* parameter is equal to *"chr1_30"*. 
In the dictianary *references* the name (*'chr1_30'*) is used to specify the region from the *hg38* genome build used for the simulations.

```yaml
simulations:
    myflat:
        topology: 'flat 10'
        simulations-number: 1
        reference: 'chr1_30'
        mutation-rates: [0.001, 0.005]
        chemistry: 'P6C4'
        coverages: [40, 60, 80]
        accuracies: [85, 90, 95]
        length-mean: 18000
        length-sd: 3000
        type: 'leaves'
    myrandom:
        topology: 'random 10'
        simulations-number: 4
        reference: 'chr1_30'
        mutation-rates: [0.001, 0.005]
        chemistry: 'P6C4'
        coverages: [40, 60, 80]
        accuracies: [85, 90, 95]
        length-mean: 18000
        length-sd: 3000
        type: 'all'
references:
    chr1_30:
       region: 'chr1:30000000-30200000'
       reference: 'hg38'
```

### Step 4: Copy the reference genomes used for simulation(s) to appropiate ditectory (or create symbolic link to it)

Reference genomes used for sinmulations should be placed in the *data/refs* directory. Fasta file shoud have *.fa* extension.

The following command can be exectued to download hg38 build of the human genome.

```
mkdir -p data/refs
cd data/refs
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

### Step 5: Starting the Snakemake workflow

To start the Snakemake workflow, the following line of code should be executed with the *number_of_threads* replaced by the maximum number of threads that workflow can use.

```
snakemake --cores number_of_threads 
```

## References

* Our article
* [Resolving multicopy duplications *de novo* using polyploid phasing](https://pubmed.ncbi.nlm.nih.gov/28808695/)
