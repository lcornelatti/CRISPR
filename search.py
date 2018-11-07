M = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]

def get_complement(sequence):
    lower = sequence.lower()
    comp = lower.replace('t', 'A').replace('a', 'T').replace('g', 'C').replace('c', 'G')
    return comp


def find_guides(sequence):
    guides = []

    for i in range(21, len(sequence)-1):
        j = len(sequence) - i
        if sequence[i] == 'G' and sequence[i+1] == 'G':
            guides.append([sequence[(i-21):(i-1)],sequence[(i-1):(i+2)]])
        if sequence[j] == 'G' and sequence[j-1] == 'G':
            guides.append([sequence[(j+2):(j+22)][::-1],sequence[(j-1):(j+2)][::-1]])

    return guides


def calculate_differences(str1, str2):
    assert len(str1) == len(str2)
    diffs = []
    for i in range(0, len(str1)):
        if str1[i] != str2[i]:
            diffs.append(i)
    return diffs


# Get potential offtargets. Only works if PAM is exactly equal
def get_potential_offtargets(guide, genome):
    potential_offtargets = []
    for i in range(0, len(genome)-22):
        this_guide = genome[(i+20):(i+23)]
        this_sequence = genome[i:i+20]
        if genome[(i+20):(i+23)] == guide[1]:
            diffs = calculate_differences(this_sequence, guide[0])
            if len(diffs) < 5:
                potential_offtargets.append([this_sequence, diffs, i])
    return potential_offtargets


def sumPairs(arr):
    n = len(arr)
    sum = 0
    for i in range(n - 1, -1, -1):
        sum += i * arr[i] - (n - 1 - i) * arr[i]
    return sum


def single_hit_score(offtarget):
    print('gf is cute')
    num_mis = len(offtarget)

    # Calculate weight of each mismatch
    t1 = 1
    for mismatch in offtarget:
        t1 *= (1-M[mismatch])

    # Calculate penalty for mismatch distance
    dsum = sumPairs(offtarget)
    pairs = num_mis*(num_mis-1)*0.5
    d = dsum/pairs if pairs != 0 else 0
    t2 = 1/(((19-d)/19)*4 + 1)

    # Calculate penalty for high mismatch
    t3 = 1/(num_mis**2)
    return t1*t2*t3

def main():

    genome = ('attcagaTCCTAATCCAGAATCCAACAGGGaattatataggaaaaACCACTTCCAGAATCCAACAGGGatatttcatctagactggtccgtccccagcacccaactgactccccgggatccgaggctccctgagttg' +
             'caatgccgtggtaccgaaatcaACCGGATCCAGAATCCAACAGGGgtagagtgtccgccgtgttgtccgcgactagacgt').upper()
    genome_comp = get_complement(genome)
    genome_reverse = genome[::-1]
    genome_comp_reverse = genome_comp[::-1]
    sequence = 'atccacctaatccagaatccaacagggctt'.upper()
    complement = get_complement(sequence)
    guides = find_guides(sequence)
    comp_guides = find_guides(complement)
    guides = guides + comp_guides
    scores = []

    for guide in guides:
        offtargets = get_potential_offtargets(guide, genome)
        offtargets += get_potential_offtargets(guide, genome_reverse)
        offtargets += get_potential_offtargets(guide, genome_comp)
        offtargets += get_potential_offtargets(guide, genome_comp_reverse)
        print('xd')
        print(offtargets)
        for item in offtargets:
            print(guide, item)
            shs = single_hit_score(item[1])
            print(shs)
if __name__ == "__main__":
    main()

