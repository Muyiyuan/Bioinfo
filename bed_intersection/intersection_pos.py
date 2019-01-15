#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import argparse
from multiprocessing import Process


chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
            'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 
            'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

def split_bed(bed, out_path, chr_list):
    bed_name = os.path.splitext(bed)[0]
    for chrom in chr_list:
        os.system("awk '{if ($1 == \"%s\") print $0}' %s > %s/%s_%s.bed" % (chrom, bed, out_path, bed_name, chrom))

def read_bed(bed):
    bed_dict = {}
    with open(bed, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            chrom = line.split('\t')[0]
            start = int(line.split('\t')[1])
            end = int(line.split('\t')[2]) + 1
            if chrom not in bed_dict:
                bed_dict[chrom] = []
            for pos in range(start, end):
                bed_dict[chrom].append(pos)
    return bed_dict

def make_result_dir(out_path):
    if not os.path.exists(out_path):
        os.mkdir(out_path)

def intersection(bed1, bed2, out_path, chrom):
    bed1_name = os.path.splitext(bed1)[0]
    bed2_name = os.path.splitext(bed2)[0]
    r = open('{}/intersection_{}.pos'.format(out_path, chrom), 'w')
    bed1_dict = read_bed('{}/{}_{}.bed'.format(out_path, bed1_name, chrom))
    bed2_dict = read_bed('{}/{}_{}.bed'.format(out_path, bed2_name, chrom))
    if chrom in bed1_dict and chrom in bed2_dict:
        intersection_pos = list(set(bed1_dict[chrom]).intersection(set(bed2_dict[chrom])))
        for pos in intersection_pos:
            r.write('{}\t{}\n'.format(chrom, pos))
    r.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='bed region intersection')
    parser.add_argument('-b1', '--bed1', help='bed1 path', required=True)
    parser.add_argument('-b2', '--bed2', help='bed2 path', required=True)
    parser.add_argument('-o', '--out', help='out path', required=True)
    args = parser.parse_args()

    bed1 = args.bed1
    bed2 = args.bed2
    out_path = args.out
    make_result_dir(out_path)
    split_bed(bed1, out_path, chr_list)
    split_bed(bed2, out_path, chr_list)
    procs = []
    for chrom in chr_list:
        proc = Process(target=intersection, args=(bed1, bed2, out_path, chrom))
        procs.append(proc)
        proc.start()
    for proc in procs:
        proc.join()
    os.system('cat {}/intersection_chr*.pos > {}/intersection.pos'.format(out_path, out_path))
