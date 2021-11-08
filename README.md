# Tree-Seg-Dup

##  Dependencies

To run Tree-Seg-Dup you should have above software installed:

* python3
* [pbsim2](https://github.com/yukiteruono/pbsim2)
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* bedtools

## Configuration

### Step 1: Installation of the required software

The following dependencies should be installed: pbsim2, snakemake, bedtools to run Tree-Seg-Dup simulator.

### Step 2: Cloning the Tree-Seg-Dup repository

To clone Tree-Seg-Dup repository, the following line should be executed:

```
git clone https://github.com/bposzewiecka/tree-seg-dup.git
```

### Step 4: Creating the configuration file

Configuration file **config.yaml** must include two dictianaries:

* *simulations*, that contains parameters of each simulation(s) by the name of each simulation.
* *regions*, that contains coordinates of regions to simulate from by the name of each region. 

Each entry from the *simulations* dictionary must have following properities:

| Property | Description | Values |
|---|---|---|
| topology | Type of the tree topology | See: topologies |  
| simulations-number | Number of the simulations | integer |
| region | Name of the region from the *regions* dictianary | string - key from *regions* dictionary |
| mutation-rates | List of probabilities  | list of probabilities |
| chemistry | HMM model of quality code for chemistry | P4C2, P5C, P6C4, R103, R94, R95 |
| coverages | List of the simulated coverages | list of integers |
| accuracies | List of the simulated reads accuracies | list of numbers from range 70-100 |
| length-mean | Mean of the simulated reads length  | integer |
| length-sd | Standard deviation of simulated reads length | integer |
| type | Simulate all region from tree topoloy or only leaves? | all, leaves |


### Topologies

* flat number_of_leaves,
* bifurcating list_of_numbers,
* cascading number_of_leaves,
* random number_of_leaves

Each entry from the *regions* dictionary must have following properties:

| Property | Description | Values |
|---|---|---|
| coordinates | Region coordinates in format chromosome:from-to | string |  
| reference | Reference name, from which region will be extracted. If reference name is *ref* the file *ref.fa* should be in the data/refs directory.  | string |


### Example configuration file

The following *yaml* file can be created to simulate duplication of two different topology types. First set of parameters *myflat* simulate
one topology tree with root and 10 leaves. Second set of parameters *myrandom* simulate 4 random topology trees with 10 leaves.

In both simulations the same input region is used. For each simulation *reference* parameter is equal to *"chr1_30"*.
In the dictionary *references* the name (*'chr1_30'*) is used to specify the region from the *hg38* genome build used for the simulations.

```yaml
simulations:
    myflat:
        topology: 'flat 10'
        simulations-number: 1
        region: 'chr1_30'
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
        region: 'chr1_30'
        mutation-rates: [0.001, 0.005]
        chemistry: 'P6C4'
        coverages: [40, 60, 80]
        accuracies: [85, 90, 95]
        length-mean: 18000
        length-sd: 3000
        type: 'all'
regions:
    chr1_30:
       coordinates: 'chr1:30000000-30200000'
       reference: 'hg38'
```

### Step 4: Copying the reference genomes used for simulation(s) to the appropriate directory

Reference genomes used for simulations should be placed in the *data/refs* directory.
Symbolic links to reference genomes can be also used.
Fasta file or symbolic link to it should have a *.fa* extension.

The following command can be executed to download hg38 human genome build.

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

## Output


## References

* Our article
* [Resolving multicopy duplications *de novo* using polyploid phasing](https://pubmed.ncbi.nlm.nih.gov/28808695/)
