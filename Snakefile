from src.trees import generate_tree, get_topology_and_sizes, save_tree
from src.tree_images import plot_tree
from src.mutation_simulator import mutate_sequences
from src.utils import  get_fasta_record 

configfile: "config.yaml"

seed = 785

PBSIM_MODELS_PATH = '/home/basia/bin/pbsim2/data/'

def get_output_files(name, pattern):
    parameters = config['simulations'][name]

    return [ pattern.format(**parameters, name=name, sim_number=i)  for i in range(parameters['simulations-number'])]
       
rule main:
   input:
       get_output_files('flat', 'data/simulations/{name}/sim-{sim_number}/tree.{reference}.{name}.sim-{sim_number}.png')

rule simulate_trees_and_mutate:
    input:
        reference = 'data/input/{reference}.fasta'
    output:
        tree_png = 'data/simulations/{name}/sim-{sim_number}/tree.{reference}.{name}.sim-{sim_number}.png',
        tree_xml = 'data/simulations/{name}/sim-{sim_number}/tree.{reference}.{name}.sim-{sim_number}.xml',
        fasta = 'data/simulations/{name}/sim-{sim_number}/regions.{reference}.{name}.sim-{sim_number}.fasta'
    params:
        vcf_pattern_fn = lambda wildcards: 'data/simulations/{name}/sim-{sim_number}/SEQ-{{seq_number}}.{reference}.{name}.sim-{sim_number}.vcf'.format(**wildcards)
    run:	
        parameters = config['simulations'][wildcards.name]
        topology_string = parameters['topology']
        probabilities = parameters['mutation-rates']     

        tree = generate_tree(**get_topology_and_sizes(topology_string))
        reference = get_fasta_record(input.reference)
        
        mutate_sequences(reference, tree, probabilities, output.fasta, params.vcf_pattern_fn)

        save_tree(tree, output.tree_xml)
        plot_tree(tree, output.tree_png)        

rule simulate_reads:
    input:
        ref = 'data/simulations/{name}/sim-{sim_number}/regions.{reference}.{name}.sim-{sim_number}.fasta',
        model = PBSIM_MODELS_PATH + '{chemistry}.model'
    output:
        'data/simulations/{name}/sim-{sim_number}/{reference}.{name}.{chemistry}.{accuracy}.{coverage}x/regions.{reference}.{name}.sim-{sim_number}.{chemistry}.{accuracy}.{coverage}x.fastq' 
    params:
        sim_dir = 'data/simulations/{name}/sim-{sim_number}/{reference}.{name}.{chemistry}.{accuracy}.{coverage}x/simulated_reads'
    log:
        pbsim = 'data/simulations/{name}/sim-{sim_number}/{reference}.{name}.{chemistry}.{accuracy}.{coverage}x/logs/pbsim/pbsim.log'
    params:
       length_mean = lambda wildcards:  str(config['simulations'][wildcards.name]['length-mean']),
       length_sd = lambda wildcards:  str(config['simulations'][wildcards.name]['length-sd']),
    shell:
        " mkdir -p  {params.sim_dir}; "
        " cd {params.sim_dir}; "
        " pbsim --hmm_model {input.model} "
        " --length-mean {params.length_mean} " +
        " --length-sd {params.length_sd} "
        " --prefix sim --id-prefix sim "
        " --accuracy-mean 0.{wildcards.accuracy} "
        " --depth {wildcards.coverage} --seed " + str(seed) + " ../../../{input.ref} 2> ../../../{log.pbsim}; "
        " cd ../../.. ; "
        " cat {params.sim_dir}/*.fastq > {output} "
