# Tree-Seg-Dup

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

### Step 3: Creating the configuration file

```yaml
simulations:
    my_flat:
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
    my_random:
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

### Step 4: Starting the snakemake workflow

To start the main algorithm, the following line of code should be executed with the *number_of_threads* replaced by the maximum number of threads that workflow can use.

```
snakemake --cores number_of_threads 
```
