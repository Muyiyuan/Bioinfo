#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse


def start_equal_end(bed_dict, chrom, start, end, equal):
    if equal == 'True':
        bed_dict[chrom].append('{}\t{}'.format(start, end))
    else:
        if start != end:
            bed_dict[chrom].append('{}\t{}'.format(start, end))

def pos_to_bed(pos_file, bed_name, equal):
    pos_dict = {}
    bed_dict = {}
    with open(pos_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            chrom = line.split('\t')[0]
            pos = int(line.split('\t')[1])
            if chrom not in pos_dict:
                pos_dict[chrom] = []
            pos_dict[chrom].append(pos)
            if chrom not in bed_dict:
                bed_dict[chrom] = []
    for chrom in sorted(pos_dict.keys()):
        pos_dict[chrom] = sorted(pos_dict[chrom])
        start = 0
        end = 0
        for i, pos in enumerate(pos_dict[chrom]):
            if i == 0:
                start = pos
                end = pos
            if i > 0 and i < len(pos_dict[chrom]) - 1:
                if pos - end == 1:
                    end = pos
                else:
                    start_equal_end(bed_dict, chrom, start, end, equal)
                    start = pos
                    end = pos
            if i == len(pos_dict[chrom]) - 1:
                if pos - end == 1:
                    end = pos
                    start_equal_end(bed_dict, chrom, start, end, equal)
                else:
                    start_equal_end(bed_dict, chrom, start, end, equal)
                    start = pos
                    end = pos
                    start_equal_end(bed_dict, chrom, start, end, equal)
    with open('{}.bed'.format(bed_name), 'w') as f:    
        for chrom in sorted(bed_dict.keys()):
            for region in bed_dict[chrom]:
                f.write('{}\t{}\n'.format(chrom, region))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='generate bed file according pos file')
    parser.add_argument('-p', '--pos', help='pos file', required=True)
    parser.add_argument('-n', '--name', help='bed name', required=True)
    parser.add_argument('-e', '--equal', help='start = end', choices=['True', 'False'], default='True')
    args = parser.parse_args()

    pos_file = args.pos
    bed_name = args.bed
    equal = args.equal

    pos_to_bed(pos_file, bed_name, equal)