from pprint import pprint
import fasta_reader
import time


M = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]
max_mismatches = 5

# Validate IUPAC code to see if 2 strings are the same
# Strcheck MUST be composed of ACGT only, otherwise this will not work
def validate_iupac(strref, strcheck):

    assert len(strref) == len(strcheck)

    for cref, ccheck in zip(strref, strcheck):
        if ccheck == 'A':
            if cref not in ('AWMRDHVN'): return False
        elif ccheck == 'C':
            if cref not in ('CSMYBHVN'): return False
        elif ccheck == 'G':
            if cref not in ('GSKRBDVN'): return False
        elif ccheck == 'T':
            if cref not in ('TWKYBDHN'): return False
        else: return False
    return True

# Gets complementary DNA strand
def get_complement(sequence):
    lower = sequence.lower()
    comp = lower.replace('t', 'A').replace('a', 'T').replace('g', 'C').replace('c', 'G')
    return comp


# Finds potential guides within a sequence (must contain an appropriate PAM)
def find_guides(sequence, guide_length, pam_ref):
    guides = []

    for i in range(0, (len(sequence) - (len(pam_ref) + guide_length) + 1)):

        i_guide = sequence[i:i + guide_length]
        i_pam = sequence[i+guide_length:i+guide_length+len(pam_ref)]

        if validate_iupac(pam_ref, i_pam):
            guides.append({'sequence':i_guide, 'pam':i_pam})

    return guides


# Calculates position of differences between two strings
def calculate_differences(str1, str2):
    assert len(str1) == len(str2)
    diffs = []
    for i in range(0, len(str1)):
        if str1[i] != str2[i]:
            diffs.append(i)
    return diffs


# Get potential offtargets
def get_potential_offtargets(guide, genome, max_mismatches, pam_ref):

    potential_offtargets = []
    guide_length = len(guide['sequence'])
    pam_length = len(pam_ref)

    for i in range(0, (len(genome) - pam_length - guide_length + 1)):

        this_sequence = genome[i:i + guide_length]
        this_pam = genome[(i+guide_length):(i+guide_length+pam_length)]

        if validate_iupac(pam_ref, this_pam):
            diffs = calculate_differences(this_sequence, guide['sequence'])

            if len(diffs) <= max_mismatches:
                potential_offtargets.append({'sequence': this_sequence,
                                             'differences': diffs,
                                             'start': i,
                                             'pam': this_pam
                                             })

    return potential_offtargets


# Return the sum of all pairs in an array
def sum_pairs(arr):
    n = len(arr)
    sum = 0
    for i in range(n - 1, -1, -1):
        sum += i * arr[i] - (n - 1 - i) * arr[i]
    return sum


# Return single hit score of an offtarget hit
def single_hit_score(offtarget):

    num_mis = len(offtarget)

    # Calculate weight of each mismatch
    t1 = 1
    for mismatch in offtarget:
        t1 *= (1-M[mismatch])

    # Calculate penalty for mismatch distance
    dsum = sum_pairs(offtarget)
    pairs = num_mis*(num_mis-1)*0.5
    d = dsum/pairs if pairs != 0 else 0
    t2 = 1/(((19-d)/19)*4 + 1)

    # Calculate penalty for high mismatch
    t3 = 1/(num_mis**2)
    return 100*t1*t2*t3


def main():
    start = time.time()
    # Get genome from file and calculate complement and reverse
    genome_fp = r"C:\Users\leoco\Desktop\Project\Homo_sapiens.GRCh38.dna.chromosome.22.fa"
    genome = fasta_reader.file_to_string(genome_fp)
    # genome = 'gctcattacgacccgagaccgacgcagacgtggtaataccaacgaactcggatgttctaaaggtgtttgcttcctacgcgtcggtgcaatcctggtagtg'.upper()
    genome_comp = get_complement(genome)
    genome_reverse = genome[::-1]
    genome_comp_reverse = genome_comp[::-1]

    # Get all possible guides
    sequence = 'atccacctaatccagaatccaacagggctt'.upper()
    pam_ref = 'NGG'
    guide_length = 20
    complement = get_complement(sequence)
    guides = find_guides(sequence, guide_length, pam_ref)
    comp_guides = find_guides(complement, guide_length, pam_ref)
    rev_guides = find_guides(sequence[::-1], guide_length, pam_ref)
    crev_guides = find_guides(complement[::-1], guide_length, pam_ref)
    guides = guides + comp_guides + rev_guides + crev_guides
    # pprint(guides)

    # For each guide calculate score
    for guide in guides:

        offtargets = get_potential_offtargets(guide, genome, max_mismatches, pam_ref)
        offtargets += get_potential_offtargets(guide, genome_reverse, max_mismatches, pam_ref)
        offtargets += get_potential_offtargets(guide, genome_comp, max_mismatches, pam_ref)
        offtargets += get_potential_offtargets(guide, genome_comp_reverse, max_mismatches, pam_ref)

        scores_sum = 0
        for item in offtargets:
            shs = single_hit_score(item['differences'])
            item['shs'] = shs
            scores_sum += shs
        guide['score'] = 100 * (100/(100 + scores_sum))
        guide['potential_offtargets'] = offtargets

    pprint(guides)
    end = time.time()
    print(end-start)
if __name__ == "__main__":
    main()

