#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import re
import argparse
from multiprocessing import Process


def read_config_file(config_file):
    config_dict = {}
    with open(config_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            config_dict[line.split(' = ')[0]] = line.split(' = ')[1] 
    return config_dict

def get_sample_list(sample_file):
    sample_list = []
    with open(sample_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            sample = line.strip()
            sample_list.append(sample)
    sample_list = sorted(sample_list)
    return sample_list

def make_result_dir(out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
def mpileup_parser(sample, out_dir):
    f = open('{}/{}.mpileup.txt'.format(out_dir, sample), 'r')
    r = open('{}/{}.parser.txt'.format(out_dir, sample), 'w')
    e = open('{}/{}.error.txt'.format(out_dir, sample), 'w')
    r.write('Chr\tPos\tRef\tAlt\tDP\tVD\tAF\n')
    e.write('Chr\tPos\tRef\tAlt\tDP\tVD\tAF\n')
    lines = f.readlines()
    for line in lines:
        chrom = line.strip().split('\t')[0]
        pos = line.strip().split('\t')[1]
        ref = line.strip().split('\t')[2].upper()
        depth = int(line.strip().split('\t')[3])
        base = line.strip().split('\t')[4].upper()
        if ref in ['A', 'T', 'C', 'G'] and depth > 0:            
            indel_temp = re.findall(r'[-\+]\d+\w+', base)
            indel_list = []
            for indel in indel_temp:
                regexp = re.match(r'([-\+])(\d+)(\w+)', indel)
                if int(regexp.group(2)) == len(regexp.group(3)):
                    indel_list.append(indel)
                else:
                    real_indel = regexp.group(1) + regexp.group(2) + regexp.group(3)[:int(regexp.group(2))]
                    indel_list.append(real_indel)
            indel_dict = {}
            total_indel = 0
            for indel in list(set(indel_list)):
                indel_dict[indel] = base.count(indel)
                total_indel += base.count(indel)

            match = base.count('.') + base.count(',') - total_indel - base.count('^.') - base.count('^,')
            base_a = base.count('A') - base.count('^A')
            base_t = base.count('T') - base.count('^T')
            base_c = base.count('C') - base.count('^C')
            base_g = base.count('G') - base.count('^G')
            base_n = base.count('N') - base.count('^N')
            base_star = base.count('*') - base.count('^*')
            for key in indel_dict.keys():
                base_a -= key.count('A') * indel_dict[key]
                base_t -= key.count('T') * indel_dict[key]
                base_c -= key.count('C') * indel_dict[key]
                base_g -= key.count('G') * indel_dict[key]
                base_n -= key.count('N') * indel_dict[key]
                base_star -= key.count('*') * indel_dict[key]

            read = ''
            count_depth = base_a + base_t + base_c + base_g + base_n + base_star + match
            for key, value in indel_dict.items():
                count_depth += value
                if 'N' not in key:
                    read += '{}:{};'.format(key, value)
            if ref == 'A':
                read = read + 'A:{};T:{};C:{};G:{};DEL:{}'.format(match, base_t, base_c, base_g, base_star)
            elif ref == 'T':
                read = read + 'A:{};T:{};C:{};G:{};DEL:{}'.format(base_a, match, base_c, base_g, base_star)
            elif ref == 'C':            
                read = read + 'A:{};T:{};C:{};G:{};DEL:{}'.format(base_a, base_t, match, base_g, base_star)
            elif ref == 'G':
                read = read + 'A:{};T:{};C:{};G:{};DEL:{}'.format(base_a, base_t, base_c, match, base_star)

            for key in read.split(';'):
                var = key.split(':')[0]
                if var != ref:
                    af = int(key.split(':')[1]) / float(depth)
                    vd = int(key.split(':')[1])
                    if re.match(r'\+\d+(\w+)', var):
                        new_alt = ref + re.match(r'[-\+]\d+(\w+)', var).group(1)
                        r.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, ref, new_alt, depth, vd, af))
                        if af < 0.15 and af > 0:
                            e.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, ref, new_alt, depth, vd, af))
                    elif re.match(r'-\d+(\w+)', var):
                        new_ref = ref + re.match(r'[-\+]\d+(\w+)', var).group(1)
                        r.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, new_ref, ref, depth, vd, af))
                        if af < 0.15 and af > 0:
                            e.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, new_ref, ref, depth, vd, af))
                    else:
                        r.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, ref, var, depth, vd, af))
                        if af < 0.15 and af > 0 and var != 'DEL':
                            e.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, ref, var, depth, vd, af))
            if count_depth != int(depth):
                print('{}:{}:{}:{}:{}:{}'.format(sample, chrom, pos, depth, count_depth, read))
    f.close()
    r.close()
    e.close()

def get_error_file(root_dir, sub_dir, bam_suffix, sample, ref_genome, samtools, bed_file, out_dir):
    os.chdir('{}/{}/{}'.format(root_dir, sample, sub_dir))
    bam_file = '{}.{}'.format(sample, bam_suffix)
    os.system('{} mpileup -d 1000000 -A -B -q 1 -Q 1 -l {} -f {} {} > {}/{}.mpileup.txt'.format(samtools, bed_file, ref_genome, bam_file, out_dir, sample))
    mpileup_parser(sample, out_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='mpileup parser')
    parser.add_argument('-c', '--config', help='config file')
    args = parser.parse_args()

    config_file = args.config
    config_dict = read_config_file(config_file)
    root_dir = config_dict['root_dir']
    sub_dir = config_dict['sub_dir']
    bam_suffix = config_dict['bam_suffix']
    sample_file = config_dict['sample_file']
    ref_genome = config_dict['ref_genome']
    samtools = config_dict['samtools']
    bed_file = config_dict['bed_file']
    out_dir = config_dict['out_dir']
    sample_list = get_sample_list(sample_file)
    make_result_dir(out_dir)
    procs = []
    for sample in sample_list:
        proc = Process(target=get_error_file, args=(root_dir, sub_dir, bam_suffix, sample, ref_genome, samtools, bed_file, out_dir))
        procs.append(proc)
        proc.start()
    for proc in procs:
        proc.join()