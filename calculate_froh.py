#!/usr/bin/env python3


import sys
import argparse
import gzip
from functools import partial
import pysam


# # This file was produced by: bcftools roh(1.9+htslib-1.9)
# The command line was: bcftools roh -O rz --AF-file VCF_sorted.ragoo.beagle.frq.AF_file.txt.gz VCF_sorted.ragoo.vcf.gz
#
# RG    [2]Sample       [3]Chromosome   [4]Start        [5]End  [6]Length (bp)  [7]Number of markers    [8]Quality (average fwd-bwd phred score)
# RG      Q41_adapt_filt_sbfi_filt_6-22748_R1adapt_filt_sbfi_filt_6-22748_R2rmdupread1.sorted.bam 1_RaGOO 2025286 3730652 1705367 60      14.1
# RG      Q41_adapt_filt_sbfi_filt_6-22748_R1adapt_filt_sbfi_filt_6-22748_R2rmdupread1.sorted.bam 1_RaGOO 6454687 7897686 1443000 58      48.6
# RG      Q41_adapt_filt_sbfi_filt_6-22748_R1adapt_filt_sbfi_filt_6-22748_R2rmdupread1.sorted.bam 1_RaGOO 44337934        48625379        4287446 196     59.3
# RG      Q41_adapt_filt_sbfi_filt_6-22748_R1adapt_filt_sbfi_filt_6-22748_R2rmdupread1.sorted.bam 1_RaGOO 62138218        71826524        9688307 396     63.2
# RG      Q41_adapt_filt_sbfi_filt_6-22748_R1adapt_filt_sbfi_filt_6-22748_R2rmdupread1.sorted.bam 1_RaGOO 111727448       114092861       2365414 66      32.3
# RG      Q41_adapt_filt_sbfi_filt_6-22807_R1adapt_filt_sbfi_filt_6-22807_R2rmdupread1.sorted.bam 1_RaGOO 14394366        15629061        1234696 34      28.5


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference", type=pysam.FastaFile, required=True)
    parser.add_argument("-m", "--min-score", default=30, type=float)
    parser.add_argument("roh", type=partial(gzip.open, mode="rt"))
    parser.add_argument("out", nargs="?", default=sys.stdout, type=argparse.FileType("wt"))
    args = parser.parse_args()

    individuals = dict()
    
    for line in args.roh:
        if line[0] == "#":
            continue
        group, sample, chrom, start, end, length, markers, qual = line.rstrip().split()
        if chrom.startswith("Z") or chrom.startswith("Chr0"):
            continue
        length, qual = int(length), float(qual)
        if qual < args.min_score:
            continue
        individuals[sample] = individuals.get(sample, 0) + length

    tot = sum([length for (chrom, length) in zip(args.reference.references, args.reference.lengths) if
               not (chrom.startswith("Z") or chrom.startswith("Chr0"))])
        
    # tot = sum(args.reference.lengths)
    for sam in individuals:
        print(sam, individuals[sam], individuals[sam] / tot, sep="\t", file=args.out)
    args.roh.close()
    args.out.close()

    return


main()
