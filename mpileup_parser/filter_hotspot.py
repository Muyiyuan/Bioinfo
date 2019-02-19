#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import re
import argparse


def get_bed_dict(bed_file):
    bed_dict = {}
    with open(bed_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if re.match(r'^chr.*', line):
                line = line.strip()
                chrom = line.split('\t')[0]
                start = int(line.split('\t')[1])
                end = int(line.split('\t')[2]) + 1
                if chrom not in bed_dict:
                    bed_dict[chrom] = []
                for i in range(start, end):
                    bed_dict[chrom].append(i)
    return bed_dict

def get_hotspot(hotspot_file, bed_dict, out_prefix):
    r = open('{}.hotspot.txt'.format(out_prefix), 'w')
    r.write('Chr\tPos\tRef\tAlt\tGene\tAA\n')
    with open(hotspot_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if re.match(r'^chr.*', line):
                chrom = line.split('\t')[0]
                pos = int(line.split('\t')[1])
                if chrom in bed_dict and pos in bed_dict[chrom]:
                    r.write(line)
    r.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter hotspot based on bed file')
    parser.add_argument('-s', '--spot', help='hotspot file')
    parser.add_argument('-b', '--bed', help='bed file')
    parser.add_argument('-o', '--out', help='out prefix')
    args = parser.parse_args()

    hotspot_file = args.spot
    bed_file = args.bed
    out_prefix = args.out

    bed_dict = get_bed_dict(bed_file)
    get_hotspot(hotspot_file, bed_dict, out_prefix)