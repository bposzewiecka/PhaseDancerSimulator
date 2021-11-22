from src.trees import generate_tree, get_topology_and_sizes, save_tree
from src.tree_images import plot_tree
from src.mutation_simulator import mutate_sequences
from src.utils import  get_fasta_record 

configfile: "config.yaml"

seed = 785

PBSIM_MODELS_PATH = '/home/basia/bin/pbsim2/data/'

def get_output_files_contig_simulation(name, pattern):
    parameters = config['simulations'][name]

    return [ pattern.format(**parameters, name=name, sim_number=i)  for i in range(parameters['simulations-number'])]

def get_output_files_reads_simulation(name, pattern):
    parameters = config['simulations'][name]
    coverages = parameters['coverages']
    accuracies = parameters['accuracies']
    
    return [ pattern.format(**parameters, name=name, sim_number=i, accuracy=accuracy, coverage=coverage)  for i in range(parameters['simulations-number']) for coverage in coverages for accuracy in accuracies ]
       
rule main:
    input:
        [ get_output_files_contig_simulation(name, 'data/simulations/{name}/sim-{sim_number}/tree.{region}.{name}.sim-{sim_number}.png') for name in config['simulations'].keys() ],
        [ get_output_files_reads_simulation(name, 'data/simulations/{name}/sim-{sim_number}/{region}.{name}.sim-{sim_number}.{type}.{chemistry}.{accuracy}.{coverage}x/{region}.{name}.sim-{sim_number}.{type}.{chemistry}.{accuracy}.{coverage}x.fastq') for name in config['simulations'].keys() ]

rule simulate_trees_and_mutate:
    input:
        region = 'data/input/{region}.fasta'
    output:
        tree_png = 'data/simulations/{name}/sim-{sim_number}/tree.{region}.{name}.sim-{sim_number}.png',
        tree_xml = 'data/simulations/{name}/sim-{sim_number}/tree.{region}.{name}.sim-{sim_number}.xml',
        fasta_all = 'data/simulations/{name}/sim-{sim_number}/mutated-sequences.{region}.{name}.sim-{sim_number}.all.fasta',
        fasta_leaves = 'data/simulations/{name}/sim-{sim_number}/mutated-sequences.{region}.{name}.sim-{sim_number}.leaves.fasta'
    params:
        vcf_pattern_fn = lambda wildcards: 'data/simulations/{name}/sim-{sim_number}/SEQ-{{seq_number}}.{region}.{name}.sim-{sim_number}.vcf'.format(**wildcards)
    run:	
        parameters = config['simulations'][wildcards.name]
        topology_string = parameters['topology']
        probabilities = parameters['mutation-rates']     

        tree = generate_tree(**get_topology_and_sizes(topology_string))
        region = str(get_fasta_record(input.region).seq).upper()
        
        mutate_sequences(region, tree, probabilities, output.fasta_all, output.fasta_leaves, params.vcf_pattern_fn)

        save_tree(tree, output.tree_xml)
        plot_tree(tree, output.tree_png)        

rule simulate_reads:
    input:
        ref = 'data/simulations/{name}/sim-{sim_number}/mutated-sequences.{region}.{name}.sim-{sim_number}.{type}.fasta',
        model = PBSIM_MODELS_PATH + '{chemistry}.model'
    output:
        'data/simulations/{name}/sim-{sim_number}/{region}.{name}.sim-{sim_number}.{type}.{chemistry}.{accuracy}.{coverage}x/{region}.{name}.sim-{sim_number}.{type}.{chemistry}.{accuracy}.{coverage}x.fastq' 
    params:
        sim_dir = 'data/simulations/{name}/sim-{sim_number}/{region}.{name}.sim-{sim_number}.{type}.{chemistry}.{accuracy}.{coverage}x/simulated_reads',
        length_mean = lambda wildcards:  str(config['simulations'][wildcards.name]['length-mean']),
        length_sd = lambda wildcards:  str(config['simulations'][wildcards.name]['length-sd'])
    log:
        pbsim = 'data/simulations/{name}/sim-{sim_number}/{region}.{name}.sim-{sim_number}.{type}.{chemistry}.{accuracy}.{coverage}x/logs/pbsim/pbsim.log'	
    shell:
        " mkdir -p  {params.sim_dir}; "
        " cd {params.sim_dir}; "
        " pbsim --hmm_model {input.model} "
        " --length-mean {params.length_mean} " 
        " --length-sd {params.length_sd} "
        " --prefix sim --id-prefix sim "
        " --accuracy-mean 0.{wildcards.accuracy} "
        " --depth {wildcards.coverage} --seed " + str(seed) + " ../../../../../../{input.ref} 2> ../../../../../../{log.pbsim}; "
        " cd ../../../../../.. ; "
        " cat {params.sim_dir}/*.fastq > {output} "

rule get_bed_reference:
    output:
        bed = 'data/input/{region}/{region}.bed'
    params:
        chromosome = lambda wildcards: config['regions'][wildcards.region]['chrom'],
        start = lambda wildcards: config['regions'][wildcards.region]['start'],
        end = lambda wildcards: config['regions'][wildcards.region]['end']
    shell:
        "echo '{params.chromosome}\t{params.start}\t{params.end}' > {output.bed}"

rule get_fasta_from_region:
    input:
        region = lambda wildcards: 'data/refs/{reference}.fa'.format(reference = config['regions'][wildcards.region]['reference']),
        bed = 'data/input/{region}/{region}-{name}.bed'
    output:
        fasta = 'data/input/{region}/{region}-{name}.fasta'
    shell:
        'bedtools getfasta -fi {input.region} -bed {input.bed} > {output.fasta}'
