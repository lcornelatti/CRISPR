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


print("Hello world")

def main():

    genome = 'attcagaaattatataggaaaaatatttcatctagactggtccgtccccagcacccaactgactccccgggatccgaggctccctgagttgcaatgccgtggtaccgaaatcagtagagtgtccgccgtgttgtccgcgactagacgt'.upper()
    sequence = 'atccacctaatccagaatccaacagggctt'.upper()
    complement = get_complement(sequence)
    guides = find_guides(sequence)
    comp_guides = find_guides(complement)
    guides = guides + comp_guides

    scores = []
    for guide in guides:
        score = score_guide(guide)
        scores.append(score)

if __name__ == "__main__":
    main()

