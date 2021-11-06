from Bio import SeqIO

def get_fasta_record(fasta_fn):

    with open(fasta_fn , 'r') as f:

        for record in SeqIO.parse(f, 'fasta'):
            return record

def save_fasta_records(sequences_by_read_name, fasta_fn):

    with open(fasta_fn, 'w') as f:
        for read_name, sequence in sequences_by_read_name.items():
            f.write('>')
            f.write(read_name)
            f.write('\n')
            f.write(sequence)
            f.write('\n')
