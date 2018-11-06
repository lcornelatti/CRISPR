# import sys


def get_string(filename):
    fasta = {}
    active_sequence_name = ''
    seq_name = 'ref|NW_018734404.1| Aedes aegypti strain LVP_AGWG chromosome 1 genomic scaffold, AaegL5.0 Primary ' + \
               'Assembly, whole genome shotgun sequence'
    with open(filename) as file_one:

        for line in file_one:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line[1:]
                if active_sequence_name not in fasta:
                    fasta[active_sequence_name] = []
                continue
            sequence = line
            fasta[active_sequence_name].append(sequence)

    print(fasta[seq_name][0])

    return fasta[seq_name]


genome = get_string(u"C:\\Users\\leoco\\Downloads\\aae_ref_AaegL5.0_chr1.fa")
print(len(genome))
