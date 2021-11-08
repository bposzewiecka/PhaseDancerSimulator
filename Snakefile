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
        get_output_files_contig_simulation('flat', 'data/simulations/{name}/sim-{sim_number}/tree.{reference}.{name}.sim-{sim_number}.png'),
        get_output_files_reads_simulation('flat', 'data/simulations/{name}/sim-{sim_number}/{reference}.{name}.sim-{sim_number}.{type}.{chemistry}.{accuracy}.{coverage}x/regions.{reference}.{name}.sim-{sim_number}.{type}.{chemistry}.{accuracy}.{coverage}x.fastq')

rule simulate_trees_and_mutate:
    input:
        reference = 'data/input/{reference}.fasta'
    output:
        tree_png = 'data/simulations/{name}/sim-{sim_number}/tree.{reference}.{name}.sim-{sim_number}.png',
        tree_xml = 'data/simulations/{name}/sim-{sim_number}/tree.{reference}.{name}.sim-{sim_number}.xml',
        fasta_all = 'data/simulations/{name}/sim-{sim_number}/regions.{reference}.{name}.sim-{sim_number}.all.fasta',
        fasta_leaves = 'data/simulations/{name}/sim-{sim_number}/regions.{reference}.{name}.sim-{sim_number}.leaves.fasta'
    params:
        vcf_pattern_fn = lambda wildcards: 'data/simulations/{name}/sim-{sim_number}/SEQ-{{seq_number}}.{reference}.{name}.sim-{sim_number}.vcf'.format(**wildcards)
    run:	
        parameters = config['simulations'][wildcards.name]
        topology_string = parameters['topology']
        probabilities = parameters['mutation-rates']     

        tree = generate_tree(**get_topology_and_sizes(topology_string))
        reference = str(get_fasta_record(input.reference).seq).upper()
        
        mutate_sequences(reference, tree, probabilities, output.fasta_all, output.fasta_leaves, params.vcf_pattern_fn)

        save_tree(tree, output.tree_xml)
        plot_tree(tree, output.tree_png)        

rule simulate_reads:
    input:
        ref = 'data/simulations/{name}/sim-{sim_number}/regions.{reference}.{name}.sim-{sim_number}.{type}.fasta',
        model = PBSIM_MODELS_PATH + '{chemistry}.model'
    output:
        'data/simulations/{name}/sim-{sim_number}/{reference}.{name}.sim-{sim_number}.{type}.{chemistry}.{accuracy}.{coverage}x/regions.{reference}.{name}.sim-{sim_number}.{type}.{chemistry}.{accuracy}.{coverage}x.fastq' 
    params:
        sim_dir = 'data/simulations/{name}/sim-{sim_number}/{reference}.{name}.sim-{sim_number}.{type}.{chemistry}.{accuracy}.{coverage}x/simulated_reads',
        length_mean = lambda wildcards:  str(config['simulations'][wildcards.name]['length-mean']),
        length_sd = lambda wildcards:  str(config['simulations'][wildcards.name]['length-sd'])
    log:
        pbsim = 'data/simulations/{name}/sim-{sim_number}/{reference}.{name}.sim-{sim_number}.{type}.{chemistry}.{accuracy}.{coverage}x/logs/pbsim/pbsim.log'	
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

rule get_bed:
    output:
        bed = 'data/input/{reference}.bed'
    params:
        region = lambda wildcards: config['references'][wildcards.reference]['region'].replace(':', '\t').replace('-', '\t')
    shell:
        "echo '{params.region}' > {output.bed}"

rule get_fasta_from_reference:
    input:
        reference = lambda wildcards: 'data/refs/{reference}.fa'.format(reference = config['references'][wildcards.reference]['reference']),
        bed = 'data/input/{reference}.bed'
    output:
        fasta = 'data/input/{reference}.fasta'
    shell:
        "bedtools getfasta -fi {input.reference} -bed {input.bed} > {output.fasta}"       
