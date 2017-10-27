import random
from numpy.random import choice
import names
import argparse
import os

# 0 denotes 5'-UTR, 1 denotes CDS region,
# 2 denotes intronic region, 3 denotes 3'-UTR region.
# and 4 denotes the end of the gene.
UTR_5 = 0
CDS = 1
INTRON = 2
UTR_3 = 3
END = 4
GENE_STATES = [UTR_5, CDS, INTRON, UTR_3, END]

# Gene Transition Probability
GENE_TP = [[0.97, 0.03, 0, 0, 0],
            [0, 0.99, 0.009, 0.001, 0],
            [0, 0.005, 0.995, 0, 0],
            [0, 0, 0, 0.96, 0.04],
            [0, 0, 0, 0, 1]]

def random_base():
    return 'ATGC'[random.randint(0, 3)]

def state_transition(state, transitionProbability):
    nState = len(transitionProbability[0])
    return choice(range(nState), 1, p=transitionProbability[state])[0]

def simulate_gene(start):
    name = names.get_last_name() + '_%d' % start
    state = UTR_5
    currPos = start

    mRNAStart, exonStart, cdsStart = start, start, None
    seq, features = [], []
    while state != END:
        seq.append(random_base())
        newState = state_transition(state, GENE_TP)
        if state == UTR_5 and newState == CDS:
            cdsStart = currPos + 1

        elif state == CDS and newState == INTRON:
            features.append((name, 'CDS', cdsStart, currPos))
            features.append((name, 'exon', exonStart, currPos))

        elif state == CDS and newState == UTR_3:
            features.append((name, 'CDS', cdsStart, currPos))

        elif state == INTRON and newState == CDS:
            cdsStart = currPos + 1
            exonStart = currPos + 1

        elif state == UTR_3 and newState == END:
            features.append((name, 'exon', exonStart, currPos))
            features.append((name, 'mRNA', mRNAStart, currPos))
        state = newState
        currPos += 1
    
    return ''.join(seq), features

def simulate_genome(maxLength, gene_frequency):
    currLength = 1
    features = []
    seq = []
    while currLength < maxLength:
        if random.random() < gene_frequency:
            s, feature = simulate_gene(currLength)
            seq.append(s)
            features.extend(feature)
            currLength += len(s)
        else:
            seq.append(random_base())
            currLength += 1

    return ''.join(seq), features

def parse_argument():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name', required=True, help='Genome name.')
    parser.add_argument('-f', '--freq', type=float, default=1e-2, help='Gene frequency.')
    parser.add_argument('-l', '--length', type=int, default=1e5, help='(Approximately) Max genome length.')
    parser.add_argument('-o', '--output', required=True, help='Output directory.')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_argument()
    seq, features = simulate_genome(maxLength=args.length, gene_frequency=args.freq)

    with open(os.path.join(args.output, '%s.fasta' % args.name), 'w') as outFile:
        print('>%s Total length %d Gene count %d' % (args.name, len(seq), sum(feature[2] == 'mRNA' for feature in features)), file=outFile)
        for i in range(0, len(seq), 70):
            print(seq[i:i+70], file=outFile)

    with open(os.path.join(args.output, '%s.gff' % args.name), 'w') as outFile:
        print('##gff-version 3.2.1', file=outFile)
        for feature in features:
            print('%s\t.\t%s\t%d\t%d\t.\t+\t0\tID=%s' % (args.name, feature[1], feature[2], feature[3], feature[0]), file=outFile)

