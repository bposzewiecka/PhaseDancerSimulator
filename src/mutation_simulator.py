from src.trees import ROOT_NAME
import random
import vcfpy
from src.utils import save_fasta_records

def get_name(node):

    return 'contig' + str(node)

def mutate_sequences(reference, tree, probabilities, fasta_all_fn, fasta_leaves_fn, vcf_pattern_fn, chrom_name, random_prefix_size = 0, random_suffix_size = 0):

    mutated_sequences = {}

    def simulate_mutations(probability, cummulated_mutations):
    
        nucleotides = {
            'A' :('T', 'G', 'C'),
            'T' :('A', 'G', 'C'),
            'G' :('A', 'T', 'C'),
            'C' :('A', 'T', 'G'),
        }

        mutations = cummulated_mutations.copy()
           
        for i, base in enumerate(reference):

            if random.uniform(0, 1) < probability:

                base = random.choice(nucleotides[base])
                mutations[i + random_prefix_size] = base

        return  mutations 
            
    def mutate_sequences_rec(node, parent_node, cummulated_mutations):
              
       probability = 0 if node == ROOT_NAME else  random.choice(probabilities)   
       tree.nodes[node]['mutation_rate'] = probability
       mutations = simulate_mutations(probability, cummulated_mutations)

       save_as_vcf(reference, mutations, node, vcf_pattern_fn, chrom_name, random_prefix_size)
       mutated_sequences[node] = get_mutated_sequence(reference, mutations, random_prefix_size, random_suffix_size)

       for child in tree[node]:      
           if child != parent_node:
               mutate_sequences_rec(child, node, mutations)
    
    mutate_sequences_rec(ROOT_NAME, None, {})
    save_mutated_sequences(mutated_sequences, tree, fasta_all_fn, fasta_leaves_fn)
                

def get_vcf_header(reference, sample_name, chrom_name):

    header = vcfpy.Header()
    header.add_line(vcfpy.HeaderLine('fileformat', 'VCFv4.2'))
    header.add_contig_line(vcfpy.OrderedDict([('ID', chrom_name), ('length', len(reference))]))
    header.add_format_line(vcfpy.OrderedDict([('ID', 'GT'), ('Number', 1), ('Type', 'String'), ('Description', 'Genotype')]))
    header.samples = vcfpy.SamplesInfos(sample_names=[sample_name])
    
    return header

def get_vcf_record(reference, mutations, node,  i, chrom_name, random_prefix_size):
    
    substitution = vcfpy.Substitution('SNV', mutations[i]),
    
    d = { 
        'CHROM':  chrom_name,
        'POS': i + 1,
        'ID': ['variant' +  str(i + 1)], 
        'REF': reference[i-random_prefix_size],
        'ALT': substitution,
        'QUAL': '.',
        'FILTER': ['PASS'], 
        'INFO': {},
        'FORMAT': ['GT'],
        'calls': [vcfpy.Call(get_name(node),{'GT': '1/1'})]
    }
    
    return  vcfpy.Record(**d)

def save_as_vcf(reference, mutations, node, seq_pattern_fn, chrom_name, random_prefix_size):

    vcf_fn = seq_pattern_fn.format(seq_number = node) 
    
    header = get_vcf_header(reference,  get_name(node), chrom_name)
    
    with vcfpy.Writer.from_path(vcf_fn, header) as f_vcf:
        for i in  sorted(mutations):
            record = get_vcf_record(reference, mutations,  node,  i, chrom_name, random_prefix_size)
            f_vcf.write_record(record)

def generate_random_sequence(size):

    nucleotides =  ['A','T', 'G', 'C']

    return ''.join([ random.choice(nucleotides)  for i in range(size)])


def get_mutated_sequence(reference, mutations, random_prefix_size, random_suffix_size):

    mutated_sequence = list(reference)

    for i, base in mutations.items():
         mutated_sequence[i - random_prefix_size] = base

    random_prefix = generate_random_sequence(random_prefix_size)
    random_suffix = generate_random_sequence(random_suffix_size)
    mutated_sequence = ''.join(mutated_sequence)

    return random_prefix + mutated_sequence + random_suffix 


def save_mutated_sequences(mutated_sequences, tree, fasta_all_fn, fasta_leaves_fn):

    all_sequences = {}
    leaves_sequences = {}

    for node, sequence in  mutated_sequences.items():
        seq_name = get_name(node)
        all_sequences[seq_name] = sequence
     
        if node != ROOT_NAME and len(tree[node])==1:  
            leaves_sequences[seq_name] = sequence

    save_fasta_records(all_sequences, fasta_all_fn)
    save_fasta_records(leaves_sequences, fasta_leaves_fn)  
