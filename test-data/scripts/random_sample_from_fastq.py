import os
import argparse
import random
import gzip
import time

random.seed(12345)

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--fastq1', required=True, help='Fastq file 1')
    parser.add_argument('-b', '--fastq2', help='Fastq file 2')
    parser.add_argument('-r', '--ratio', required=True, type=float, help='Sampling ratio.')
    parser.add_argument('-p', '--prefix', required=True, help='Output prefix.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase verbosity')
    
    return parser.parse_args()

def random_sample_paired_end(fastq1, fastq2, ratio, prefix):
    outFileName1 = os.path.join(os.path.dirname(fastq1), prefix + '_1.fastq')
    outFileName2 = os.path.join(os.path.dirname(fastq2), prefix + '_2.fastq')
    inFile1 = open(fastq1)
    inFile2 = open(fastq2)
    outFile1 = open(outFileName1, 'w')
    outFile2 = open(outFileName2, 'w')

    i = 0
    chunkSize = 1000
    chunk = 0
    l1, l2 = [], []
    eof = False
    while not eof:
        if random.random() < ratio:
            # Sample this read.
            chunk += 1

            for _ in range(4):
                i += 1
                if i % 10000000 == 0:
                    print('[%s] Processed %dth line.' % (time.ctime(), i))
                m1 = inFile1.readline()
                m2 = inFile2.readline()
                l1.append(m1)
                l2.append(m2)
                if len(m1.strip()) == 0 or len(m2.strip()) == 0:
                    eof = True
                
            if chunk == chunkSize:
                print('\n'.join(l1), file=outFile1)
                print('\n'.join(l2), file=outFile2)
                chunk = 0
                l1, l2 = [], []
        else:
            # Ignore this read.
            for _ in range(4):
                i += 1
                if i % 10000000 == 0:
                    print('[%s] Processed %dth line.' % (time.ctime(), i))
                m1 = inFile1.readline()
                m2 = inFile2.readline()
                if len(m1.strip()) == 0 or len(m2.strip()) == 0:
                    eof = True

    inFile1.close()
    inFile2.close()
    outFile1.close()
    outFile2.close()

def random_sample_single_end(fastq1, N):
    pass
    
if __name__ == '__main__':
    args = parse_arguments()
    random_sample_paired_end(args.fastq1, args.fastq2, ratio=args.ratio, prefix=args.prefix) 
