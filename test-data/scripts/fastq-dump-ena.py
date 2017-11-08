#!/usr/bin/python

import argparse
import subprocess
import os

# specify arguments
parser = argparse.ArgumentParser(description="Dump fastq file from ENA.")
parser.add_argument("ACCESSION", help="Run accession.")
parser.add_argument("-p", "--paired", action="store_true", help="Specify this option if the reads are paired.")
parser.add_argument("-o", "--output", default=".", help="Output directory")

args = parser.parse_args()

# make directory if it does not exist
if not os.path.exists(args.output):
    os.makedirs(args.output)

accession = args.ACCESSION
dir1 = accession[:6]
dir2 = ["", "/%03d" % int(accession[9:])][len(accession) > 9]

if not args.paired:
    path = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s%s/%s/%s.fastq.gz" % (dir1, dir2, accession)
    print("[%s] Calling wget %s" % (__file__, path))
    subprocess.call(["wget", "-P", args.output, path])

else:
    path1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s%s/%s/%s_1.fastq.gz" % (dir1, dir2, accession, accession)
    path2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s%s/%s/%s_2.fastq.gz" % (dir1, dir2, accession, accession)
    print("[%s] Calling wget %s" % (__file__, path1))
    subprocess.call(["wget", "-P", args.output, path1])
    print("[%s] Calling wget %s" % (__file__, path2))
    subprocess.call(["wget", "-P", args.output, path2])