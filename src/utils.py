from Bio import SeqIO
import yaml

def get_fasta_record(fasta_fn):

    with open(fasta_fn , 'r') as f:

        for record in SeqIO.parse(f, 'fasta'):
            return record

def get_fasta_record_by_name(fasta_fn, name):

    with open(fasta_fn , 'r') as f:

        for record in SeqIO.parse(f, 'fasta'):
            if name==record.id: 
                return record


def save_fasta_records(sequences_by_read_name, fasta_fn):

    with open(fasta_fn, 'w') as f:
        for read_name, sequence in sequences_by_read_name.items():
            f.write('>')
            f.write(read_name)
            f.write('\n')
            f.write(sequence)
            f.write('\n')

def save_yaml(yaml_fn, data):

    if not yaml_fn.endswith('.yaml'):
        yaml_fn += '.yaml'

    with open(yaml_fn, 'w') as f:
        yaml.dump(data, f, default_flow_style = None)

def get_number_of_reads(fasta_fn):

    with open(fasta) in f:
        n = 0

        for line in f:
            if line.startswith(">"):
                n += 1
       
        return n
