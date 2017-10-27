import argparse
import random

def mutate(seq, rate=1e-6):
    newSeq = []
    for base in seq:
        if random.random() < rate:
            newSeq.append('ATGC'[random.randint(0, 3)])
        else:
            newSeq.append(base)

    return ''.join(newSeq)

def get_mRNA_sequences_from_gtf(referenceSeqPath, gtfPath):
    with open(referenceSeqPath) as inFile:
        inFile.readline()
        referenceSeq = ''.join([line.strip() for line in inFile.readlines()])

    seqs = []
    with open(gtfPath) as inFile:
        seq = []
        for line in inFile.readlines():
            if line.startswith('##'):
                continue

            tokens = line.strip().split('\t')
            if tokens[2] == 'CDS':
                seq.append(referenceSeq[int(tokens[3]):int(tokens[4])+1])

            if tokens[2] == 'mRNA':
                seqs.append(''.join(seq))

    return seqs

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reference', required=True, help='Reference genome sequence file.')
    parser.add_argument('-g', '--gff', required=True, help='Annotation file.')
    parser.add_argument('-o', '--output', required=True, help='Output fastq file')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    p = 0.98
    readLength = 75

    reads = []
    sequences = get_mRNA_sequences_from_gtf(args.reference, args.gff)
    for sequence in sequences:
        dice = random.random()
        while dice < p:
            start = random.randint(0, len(sequence) - readLength)
            reads.append(mutate(sequence[start:start+readLength]))
            dice = random.random()

    with open(args.output, 'w') as outFile:
        for i, read in enumerate(reads, 1):
            print('@ read id=%d length=%d' % (i, len(read)), file=outFile)
            print(read, file=outFile)
            print('+', file=outFile)
            print('I' * len(read), file=outFile)
