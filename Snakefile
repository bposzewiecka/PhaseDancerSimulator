from src.trees import generate_tree, get_topology_and_sizes, save_tree
from src.tree_images import plot_tree
from src.mutation_simulator import mutate_sequences
from src.utils import  get_fasta_record, get_fasta_record_by_name, save_yaml
from src.validation import get_alignment_results, get_variants_validation_table, hamming_distance

import pysam
from Bio import SeqIO

configfile: "config.yaml"

seed = 785

PBSIM_MODELS_PATH = 'pbsim2/data/'
PBSIM_BIN_PATH = 'pbsim2/bin/pbsim'
SAMTOOLS_BIN_PATH = 'samtools-1.14/samtools'
BEDTOOLS_BIN_PATH = 'bedtools2/bin/bedtools'
MINIMAP_BIN_PATH = 'minimap2/minimap2'

OUTPUT_DIR = ''

PACBIO = 'pacbio'
NANOPORE = 'nanopore'


def get_technology(chemistry):
    if chemistry in ['P4C2', 'P5C3', 'P6C4']:
        return PACBIO
    if chemistry in ['R103', 'R94', 'R95']:
        return NANOPORE
    raise

def get_output_files_contig_simulation(name, pattern):
    parameters = config['simulations'][name]
    regions =  parameters['regions']

    return [ pattern.format(**parameters, name=name, sim_number=i, region=region)  for i in range(parameters['simulations-number']) for region in regions]

def get_output_files_reads_simulation(name, pattern):
    parameters = config['simulations'][name]
    coverages = parameters['coverages']
    accuracies = parameters['accuracies']
    chemistries = parameters['chemistries']
    regions = parameters['regions']
   
    return [ pattern.format(**parameters, name=name, sim_number=i, accuracy=accuracy, coverage=coverage, chemistry=chemistry, region=region)  for i in range(parameters['simulations-number']) for coverage in coverages for accuracy in accuracies for chemistry in chemistries for region in regions]

ruleorder: map_reads > index_reads
           
rule main:
    input:
        [ get_output_files_contig_simulation(name, OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/tree-{region}-{name}-sim{sim_number}.png') for name in config['simulations'].keys()],
        [ get_output_files_contig_simulation(name, OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/starts/contigroot-start-mutatedseq-{region}-{name}-sim{sim_number}.fasta') for name in config['simulations'].keys() if 'start-contig' in config['simulations'][name]],
        [ get_output_files_reads_simulation(name, OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x.bam.bai') for name in config['simulations'].keys() ],
        expand(OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{ttype}-{chemistry}-{accuracy}-{coverage}x/contig{seq_number}.{assembler_name}.{stats}.yaml', name='vsmallflat', sim_number=0, seq_number=0, region ='region0new', ttype= 'leaves', accuracy= 85, chemistry='P6C4', coverage=40, assembler_name='phaseDancer', stats = ['stats', 'variants_stats' ])

rule simulate_trees_and_mutate:
    input:
        region = OUTPUT_DIR + 'data/input/{region}/{region}.fasta'
    output:
        tree_png = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/tree-{region}-{name}-sim{sim_number}.png',
        tree_xml = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/tree-{region}-{name}-sim{sim_number}.xml',
        fasta_all = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/mutatedseq-{region}-{name}-sim{sim_number}-all.fasta',
        fasta_leaves = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/mutatedseq-{region}-{name}-sim{sim_number}-leaves.fasta'
    params:
        vcf_pattern_fn = lambda wildcards: 'data/simulations/{name}/sim{sim_number}/contig{{seq_number}}-{region}-{name}-sim{sim_number}.vcf'.format(**wildcards)
    run:	
        parameters = config['simulations'][wildcards.name]
        topology_string = parameters['topology']
        probabilities = parameters['mutation-rates']     
        random_prefix_size = parameters.get('random-prefix-size', 0)
	random_suffix_size = parameters.get('random-suffix-size', 0)
        
        tree = generate_tree(**get_topology_and_sizes(topology_string))
        region = str(get_fasta_record(input.region).seq).upper()
        
        mutate_sequences(region, tree, probabilities, output.fasta_all, output.fasta_leaves, params.vcf_pattern_fn, wildcards.region, random_prefix_size, random_suffix_size)

        save_tree(tree, output.tree_xml)
        plot_tree(tree, output.tree_png)        

rule simulate_reads:
    input:
        ref = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/mutatedseq-{region}-{name}-sim{sim_number}-{type}.fasta',
        model = PBSIM_MODELS_PATH + '{chemistry}.model'
    output:
        OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x.fastq' 
    params:
        sim_dir = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/simulated_reads',
	difference_ratio = lambda wildcards: '6:50:54' if get_technology(wildcards.chemistry) == PACBIO else '23:31:46',
        length_mean = lambda wildcards:  str(config['simulations'][wildcards.name]['length-mean']),
        length_sd = lambda wildcards:  str(config['simulations'][wildcards.name]['length-sd'])
    log:
        pbsim = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/logs/pbsim/pbsim.log'	
    shell:
        " mkdir -p  {params.sim_dir}; "
        " export SIM_PWD=`pwd`;"
        " cd {params.sim_dir}; "
        " $SIM_PWD/" + PBSIM_BIN_PATH  + " --hmm_model $SIM_PWD/{input.model} "
        " --length-mean {params.length_mean} " 
        " --length-sd {params.length_sd} "
        " --prefix sim --id-prefix sim "
        " --accuracy-mean 0.{wildcards.accuracy} "
	" --difference-ratio {params.difference_ratio}  "
        " --depth {wildcards.coverage} --seed " + str(seed) + " $SIM_PWD/{input.ref} 2> $SIM_PWD/{log.pbsim}; "
        " cd $SIM_PWD ; "
        " cat {params.sim_dir}/*.fastq > {output} "

rule get_starts_sequences:
    input:
        ref = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/mutatedseq-{region}-{name}-sim{sim_number}-all.fasta'
    output:
        starts = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/starts/contigroot-start-mutatedseq-{region}-{name}-sim{sim_number}.fasta'
    params:
        fn_pattern = lambda wildcards: OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/starts/{{contig}}-start-mutatedseq-{region}-{name}-sim{sim_number}.fasta'.format(**wildcards)
    run:
        for record in SeqIO.parse(input.ref, 'fasta'):
    
           fn = params.fn_pattern.format(contig=record.name)
            
           start_contig = config['simulations'][wildcards.name]['start-contig']

           start = start_contig['start-coordinate']
           length = start_contig['contig-size']
    
           with open(fn, 'w') as f:
               f.write(f'>{record.name} start={start} length={length}\n')
               f.write(str(record.seq[start:start + length]))
               f.write('\n')
             
rule get_bed_reference:
    output:
        bed = OUTPUT_DIR + 'data/input/{region}/{region}.bed'
    params:
        chromosome = lambda wildcards: config['regions'][wildcards.region]['chrom'],
        start = lambda wildcards: config['regions'][wildcards.region]['start'],
        end = lambda wildcards: config['regions'][wildcards.region]['end']
    shell:
        "echo '{params.chromosome}\t{params.start}\t{params.end}' > {output.bed}"

rule get_fasta_from_region:
    input:
        region = lambda wildcards: OUTPUT_DIR + 'data/refs/{reference}.fa'.format(reference = config['regions'][wildcards.region]['reference']),
        bed = OUTPUT_DIR + 'data/input/{region}/{region}.bed'
    output:
        fasta = OUTPUT_DIR + 'data/input/{region}/{region}.fasta'
    shell:
        './' + BEDTOOLS_BIN_PATH + ' getfasta -fi {input.region} -bed {input.bed} > {output.fasta}; sed -i "1s/.*/>{wildcards.region}/" {output.fasta} '

rule map_reads:
    input:
        ref = OUTPUT_DIR + 'data/input/{region}/{region}.fasta',
	reads = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x.fastq'
    output:
        bam = temp(OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x.notgrouped.bam'),
        bai = temp(OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x.notgrouped.bam.bai')
    params:
         minimap2_preset = lambda wildcards: 'map-pb' if get_technology(wildcards.chemistry) == PACBIO else 'map-ont'
    shell:
        './' + MINIMAP_BIN_PATH + ' -ax {params.minimap2_preset} {input.ref} {input.reads} | samtools sort - > {output.bam}; '
	'./' + SAMTOOLS_BIN_PATH + ' index {output.bam} '

rule add_group_type:
    input:
        bam = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x.notgrouped.bam',
        bai = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x.notgrouped.bam.bai'
    output:
        bam = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x.bam'
    run:        
        with pysam.AlignmentFile(input.bam, 'rb') as bam:  

            with pysam.AlignmentFile(output.bam, 'wb',  template=bam) as grouped_bam:

                for read in bam.fetch():
                    simulation_number = int(read.query_name[3:].split('_')[0])
                    read.set_tag('RG', f"sim{simulation_number:03d}")
                    grouped_bam.write(read)
          
rule index_reads:
    input:
        '{name}.bam'	 
    output:
       	'{name}.bam.bai'    
    shell:	   
        './'+ SAMTOOLS_BIN_PATH + ' index {input}'

rule identity_statistics_contig:
    input:
        reference = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/mutatedseq-{region}-{name}-sim{sim_number}-all.fasta',
        assembly = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/contig{seq_number}.{assembler_name}.fasta'
    output:
        stats = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/contig{seq_number}.{assembler_name}.stats.yaml'
    run:
         reference = get_fasta_record_by_name(input.reference, 'contig' + wildcards.seq_number) 
         assembly = get_fasta_record(input.assembly)

         results = get_alignment_results(reference, assembly)
 
         save_yaml(output.stats, { f"contig{wildcards.seq_number}.{wildcards.assembler_name}" : results } )


rule variants_statistics_contig:
    input:
        reference = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/mutatedseq-{region}-{name}-sim{sim_number}-all.fasta',
        assembly = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/contig{seq_number}.{assembler_name}.fasta',
        vcf = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/contig{seq_number}-{region}-{name}-sim{sim_number}.vcf'
    output:
        stats = OUTPUT_DIR + 'data/simulations/{name}/sim{sim_number}/{region}-{name}-sim{sim_number}-{type}-{chemistry}-{accuracy}-{coverage}x/contig{seq_number}.{assembler_name}.variants_stats.yaml'
    run:
        contig = 'contig' + wildcards.seq_number	    

        reference = str(get_fasta_record_by_name(input.reference, 'contig' + wildcards.seq_number).seq)
        assembly = str(get_fasta_record(input.assembly).seq)

        validation_table =  get_variants_validation_table(reference, assembly, contig, input.vcf)
        save_yaml(output.stats, { f"contig{wildcards.seq_number}.{wildcards.assembler_name}": { "validation-table": validation_table, "hamming-distance": hamming_distance(validation_table) }})
