import math
from Bio import SeqIO 
import edlib
import re
import vcfpy

def get_identity_by_cigar(cigar):

    alignment_size = 0
    errors = 0

    for bases_sign, sign in re.findall('(\d+(=|D|I|X))', cigar):
        bases = int(bases_sign[:-1])

        alignment_size += bases

        if sign != '=':
           errors += bases

    return {
        'alignment_size': alignment_size,
        'edit_distance': errors,
        'identity': 1 - errors / alignment_size,
        'phred_quality': math.log10(errors/alignment_size) * -10
    }

def get_alignment_results(reference, assembly):

    reference_seq = str(reference.seq)
    assembly_seq = str(assembly.seq)

    result = edlib.align(assembly_seq[100:-100], reference_seq, mode="HW", task="path")

    return get_identity_by_cigar(result['cigar'])

def get_variants_validation_table(reference, assembly, contig, vcf_fn):

    result = edlib.align(assembly, reference, mode="HW", task="path")
    start_loci = result['locations'][0][0]
    alignment = edlib.getNiceAlignment(result,  assembly, reference)
    query_aligned = alignment['query_aligned']
    target_aligned = alignment['target_aligned']

    reader = vcfpy.Reader.from_path(vcf_fn)

    variants = { record.POS - 1: (no, record.ALT[0].value) for no, record in enumerate(reader) }  

    i = 0

    validation_table =  [None] * len(variants)
    
    for base_query, base_target in zip(query_aligned, target_aligned):
        
        if base_target == '-':
            continue
        
        if start_loci + i in variants:
            no, base_alt = variants[start_loci + i]
            correct = int(base_alt == base_query)
            validation_table[no] = correct

        i+=1
        
    return validation_table

def hamming_distance(validation_table):
    
    return sum([1 for correct in validation_table if correct == 0 or correct is None])/len(validation_table)
