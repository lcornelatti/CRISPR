# import sys


def read_fasta(filename):
    fasta = {}
    active_sequence_name = ''
    # seq_name = '22 dna:chromosome chromosome:GRCh38:22:1:50818468:1 REF'
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

    for key in fasta.keys():
        return fasta[key]


def remove_ns(sequence_list):
    start = 0
    end = 0
    for i in range(0, len(sequence_list)):
        if not (sequence_list[i].startswith("NNN") and sequence_list[i].endswith("NNN")):
            start = i
            break

    for i in range(start, len(sequence_list)):
        if sequence_list[i].startswith("NNN") and sequence_list[i].endswith("NNN"):
            end = i
            break

    assert end != 0 and start < end

    return sequence_list[start:end]


def remove_ns2(sequence_list):
    return [x for x in sequence_list if list(set(x)) != ['N']]


def list_to_string(list):
    nstring = ''
    for item in list:
        nstring += item

    return nstring


def file_to_string(filepath):
    fasta_sequence = read_fasta(filepath)
    fasta_sequence2 = remove_ns2(fasta_sequence)
    string_genome = ''.join(fasta_sequence2)
    print(len(string_genome))
    return string_genome


if __name__ == "__main__":
    fp = r"C:\Users\leoco\Desktop\Project\Homo_sapiens.GRCh38.dna.chromosome.22.fa"
    file_to_string(fp)
