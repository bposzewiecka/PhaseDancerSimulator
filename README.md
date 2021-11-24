# Tree-Seg-Dup - Simulator of collapsed segmental duplications

To account for a complex history of the formation of segmental duplications **Tree-Seg-Dup** simulates such regions based on the tree topology. Fragment of the reference genome is assigned to the root of the tree. Child sequences are created by copying a parent node sequence and mutating each base at a fixed per base rate. Then, long-reads (PacBio or Oxford Nanopore) are simulated using [PBSIM2](https://github.com/yukiteruono/pbsim2).

**Tree-Seg-Dup** outputs:

* image presenting tree topology.
* XML file storing tree topology ([gexf](https://gephi.org/gexf/format/)).
* information about the mutations for every node of the tree (in *vcf* format).
* reads simulated for leaves or all nodes (depending on the *type* setting).
* reads mapped on the root sequence grouped by node (RG tag)
 
##  Dependencies

To run Tree-Seg-Dup you should have above software installed:

* python3
* [snakemake](https://snakemake.readthedocs.io/en/stable/)


## Configuration

### Step 1: Installation of the required software

The following dependencies should be installed: pbsim2, snakemake, bedtools to run Tree-Seg-Dup simulator.

### Step 2: Cloning the Tree-Seg-Dup repository and installing dependencies

To clone Tree-Seg-Dup repository and installing dependencies, the following lines should be executed:

```
git clone https://github.com/bposzewiecka/tree-seg-dup.git
./install.sh
```

### Step 3: Creating the configuration file

Configuration file **config.yaml** must include two dictionaries:

* *simulations*, that contains parameters of each simulation(s) by the name of each simulation.
* *regions*, that contains coordinates of regions to simulate from by the name of each region.

Each entry from the *simulations* dictionary must have the following properties:

| Property | Description | Values |
|---|---|---|
| topology | Type of the tree topology | See: topologies |  
| simulations-number | Number of the simulations | integer |
| region | Name of the region from the *regions* dictionary | string - key from *regions* dictionary |
| mutation-rates | List of probabilities  | list of probabilities |
| chemistry | HMM model of quality code for chemistry. See: chemistries |  P4C2, P5C3, P6C4, R103, R94, R95 |
| coverages | List of the simulated coverages | list of integers |
| accuracies | List of the simulated reads accuracies | list of numbers from range 70-100 |
| length-mean | Mean of the simulated reads length  | integer |
| length-sd | Standard deviation of simulated reads length | integer |
| type | Simulate  from all tree nodes or only leaves? | all, leaves |

### Chemistries

Reads are simulated using FIC-HMM model for  chemistries of PacBio and Nanopore. 

| Name |  Technology  | 
|---|---|
| P6C4 | PacBio | 
| P5C3| PacBio | 
| P4C2 | PacBio | 
| R103 | Nanopore | 
| R95 | Nanopore | 
| R94 | Nanopore | 

### Topologies
 
Tree topologies can be simulated using 4 topology types: flat, cascading, bifurcating and random.
 
| |  Flat  |  Cascading |
|---|---|---|
| Topology | ![Flat topology](/images/flat.png?raw=true "Flat topology") | ![Cascading topology](/images/cascading.png?raw=true "Cascading topology") |
| Property format | **flat** number_of_leaves | **cascading** number_of_leaves |
| Example value | flat 6 | cascading 6 |

| |  Bifurcating  |  Random |
|---|---|---|
| Topology | ![Bifurcating topology](/images/bifurcating.png?raw=true "Bifurcating topology") | ![Random topology](/images/random.png?raw=true "Random topology") |
| Property format |  **bifurcating** nodes_multiplier_level_1 nodes_multiplier_level_2 nodes_multiplier_level_3 ... | **random** number_of_leaves |
| Example value | bifurcating 2 3 2 | random 6 |

Each entry from the *regions* dictionary must have the following properties:

| Property | Description | Values |
|---|---|---|
| reference | Reference name, from which region will be extracted. If reference name is *ref* the file *ref.fa* should be in the data/refs directory.  | string |
| chrom | Region chromosome | string |  
| start | Region start coordinate | number |
| end | Region end coordinate | number |

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
        chemistries: ['P6C4', 'R103']
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
        chemistries: ['P6C4', 'R103']
        coverages: [40, 60, 80]
        accuracies: [85, 90, 95]
        length-mean: 18000
        length-sd: 3000
        type: 'all'
regions:
    chr1_30:
       chrom: 'chr1'
       start: 30000000
       end: 30200000
       reference: 'hg38'
```

### Step 4: Placing the reference genomes (or symbolic link to them) used for simulations in the appropriate directory

Reference genomes used for simulations should be placed in the *data/refs* directory.
Symbolic links to reference genomes can be also used.
Fasta file or symbolic link to it should have a *.fa* extension.

The following command can be executed to download *hg38* human genome build.

```
mkdir -p data/refs
cd data/refs
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

### Step 5: Starting the Snakemake workflow

To start the Snakemake workflow, the following line of code should be executed with the *number_of_threads* replaced by the maximum number of threads that can be used for simulations.

```
snakemake --cores number_of_threads
```

## Note

The **Tree-Seg-Dup** simulator is based on the idea of simulating segmental duplications derived from tree topologies proposed in  [^2].
Simulator was used to validate and benchmark PhaseDancer [^1] assembler.

Embdeded software: 
* [PBSIM2](https://github.com/yukiteruono/pbsim2) - [GNU General Public License v2.0](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
* [Samtools](http://www.htslib.org/) - [Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0) Licence](https://creativecommons.org/licenses/by-nc-nd/4.0/)
* [Bedtools](https://bedtools.readthedocs.io/en/latest/) - [GNU Public License (Version 2)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

[^1]: [PhaseDancer](https://github.com/bposzewiecka/phaseDancer)
[^2]: [Resolving multicopy duplications *de novo* using polyploid phasing](https://pubmed.ncbi.nlm.nih.gov/28808695/)
