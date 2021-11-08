from src.trees import ROOT_NAME
import random
import vcfpy
from src.utils import save_fasta_records

CHROM_ID  = 1

def get_sample_name(node):

    return 'sample' + str(node)

def mutate_sequences(reference, tree, probabilities, fasta_all_fn, fasta_leaves_fn, vcf_pattern_fn):

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
                mutations[i] = base

        return  mutations 
            
    def mutate_sequences_rec(node, parent_node, cummulated_mutations):
              
       probability = 0 if node == ROOT_NAME else  random.choice(probabilities)   
       tree.nodes[node]['mutation_rate'] = probability
       mutations = simulate_mutations(probability, cummulated_mutations)

       save_as_vcf(reference, mutations, node, vcf_pattern_fn)
       mutated_sequences[node] = get_mutated_sequence(reference, mutations)

       for child in tree[node]:      
           if child != parent_node:
               mutate_sequences_rec(child, node, mutations)
    
    mutate_sequences_rec(ROOT_NAME, None, {})
    save_mutated_sequences(mutated_sequences, tree, fasta_all_fn, fasta_leaves_fn)
                

def get_vcf_header(reference, sample_name):

    header = vcfpy.Header()
    header.add_line(vcfpy.HeaderLine('fileformat', 'VCFv4.2'))
    header.add_contig_line(vcfpy.OrderedDict([('ID', CHROM_ID), ('length', len(reference))]))
    header.add_format_line(vcfpy.OrderedDict([('ID', 'GT'), ('Number', 1), ('Type', 'String'), ('Description', 'Genotype')]))
    header.samples = vcfpy.SamplesInfos(sample_names=[sample_name])
    
    return header

def get_vcf_record(reference, mutations, node,  i):
    
    substitution = vcfpy.Substitution('SNV', mutations[i]),
    
    d = { 
        'CHROM':  CHROM_ID,
        'POS': i + 1,
        'ID': ['variant' +  str(i + 1)], 
        'REF': reference[i],
        'ALT': substitution,
        'QUAL': '.',
        'FILTER': ['PASS'], 
        'INFO': {},
        'FORMAT': ['GT'],
        'calls': [vcfpy.Call(get_sample_name(node),{'GT': '1/1'})]
    }
    
    return  vcfpy.Record(**d)

def save_as_vcf(reference, mutations, node, seq_pattern_fn):

    vcf_fn = seq_pattern_fn.format(seq_number = node) 
    
    header = get_vcf_header(reference,  get_sample_name(node))
    
    with vcfpy.Writer.from_path(vcf_fn, header) as f_vcf:
        for i in  sorted(mutations):
            record = get_vcf_record(reference, mutations,  node,  i)
            f_vcf.write_record(record)


def get_mutated_sequence(reference, mutations):

    mutated_sequence = list(reference)

    for i, base in mutations.items():
         mutated_sequence[i] = base

    return  ''.join(mutated_sequence)


def save_mutated_sequences(mutated_sequences, tree, fasta_all_fn, fasta_leaves_fn):

    all_sequences = {}
    leaves_sequences = {}

    for node, sequence in  mutated_sequences.items():
        seq_name = get_sample_name(node)
        all_sequences[seq_name] = sequence
     
        if node != ROOT_NAME and len(tree[node])==1:  
            leaves_sequences[seq_name] = sequence

    save_fasta_records(all_sequences, fasta_all_fn)
    save_fasta_records(leaves_sequences, fasta_leaves_fn)  
