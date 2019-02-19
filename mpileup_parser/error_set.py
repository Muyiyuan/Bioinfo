#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import numpy as np
from scipy.stats import weibull_min


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

def get_error_set(root_dir, error_suffix, sample_list, out_prefix):
    var_list = []
    dp_0_dict = {}  # include sample with vd >= 0
    dp_1_dict = {}  # include sample with vd >= 1
    vd_dict = {}
    sample_0_dict = {}  # include sample with vd >= 0
    sample_1_dict = {}  # include sample with vd >= 1
    af_list_dict = {}
    for sample in sample_list:
        with open('{}/{}.{}'.format(root_dir, sample, error_suffix)) as f:
            for line in f:
                if re.match(r'^chr.*', line):
                    line = line.strip()
                    chrom, pos, ref, alt, dp, vd, af = line.split('\t')
                    var = '{}_{}_{}_{}'.format(chrom, pos, ref, alt)
                    var_list.append(var)
    var_list = list(set(var_list))
    for var in var_list:
        af_list_dict[var] = []
        dp_0_dict[var] = 0
        dp_1_dict[var] = 0
        vd_dict[var] = 0
        sample_0_dict[var] = 0
        sample_1_dict[var] = 0       
    for sample in sample_list:
        with open('{}/{}.{}'.format(root_dir, sample, error_suffix)) as f:
            for line in f:
                if re.match(r'^chr.*', line):
                    line = line.strip()
                    chrom, pos, ref, alt, dp, vd, af = line.split('\t')
                    var = '{}_{}_{}_{}'.format(chrom, pos, ref, alt)
                    dp_0_dict[var] += int(dp)
                    vd_dict[var] +=  int(vd)
                    sample_0_dict[var] += 1
                    af_list_dict[var].append(af)
                    if float(af) > 0:
                        dp_1_dict[var] += int(dp)
                        sample_1_dict[var] += 1
    r = open('{}/{}.error.txt'.format(root_dir, out_prefix), 'w')
    r.write('Chr\tPos\tRef\tAlt\tSample(T|0|1)\tDP(0|1)\tVD\tAF(0|1)\tAF_List\tDistribution\tMean/Shape\tSD/Scale\n')
    sample = len(sample_list)
    for var in sorted(var_list):
        chrom, pos, ref, alt = var.split('_')
        err_0 = sample_0_dict[var]
        err_1 = sample_1_dict[var]
        dp_0 = dp_0_dict[var]
        dp_1 = dp_1_dict[var]
        vd = vd_dict[var]
        af_0 = vd / float(dp_0)
        af_1 = 0.0
        if dp_1 != 0:
            af_1 = vd / float(dp_1)
        af_list = ','.join(af_list_dict[var])
        float_af_list = []
        for af in af_list_dict[var]:
            float_af_list.append(float(af))
        af_max = max(float_af_list)
        gas_list = []
        wei_list = []
        for i in float_af_list:
            if i != af_max:
                gas_list.append(i)
                if i != 0.0:
                    wei_list.append(i)
        if len(float_af_list) == 1 or len(float_af_list) == float_af_list.count(0.0):
                r.write('{}\t{}\t{}\t{}\t{}|{}|{}\t{}|{}\t{}\t{}|{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, ref, alt, sample, err_0, err_1, dp_0, dp_1, vd, af_0, af_1, af_list, 'NA', 'NA', 'NA'))
        else:
            if len(wei_list) < 5:
                dis = "Gaussian"
                mean = np.mean(gas_list)
                sd = np.std(gas_list)
                r.write('{}\t{}\t{}\t{}\t{}|{}|{}\t{}|{}\t{}\t{}|{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, ref, alt, sample, err_0, err_1, dp_0, dp_1, vd, af_0, af_1, af_list, dis, mean, sd))
            else:
                dis = 'Weibull'
                shape, loc, scale = weibull_min.fit(wei_list, floc=0)
                r.write('{}\t{}\t{}\t{}\t{}|{}|{}\t{}|{}\t{}\t{}|{}\t{}\t{}\t{}\t{}\n'.format(chrom, pos, ref, alt, sample, err_0, err_1, dp_0, dp_1, vd, af_0, af_1, af_list, dis, shape, scale))
    r.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Error set')
    parser.add_argument('-c', '--config', help='config file')
    args = parser.parse_args()

    config_file = args.config
    config_dict = read_config_file(config_file)
    root_dir = config_dict['root_dir']
    error_suffix = config_dict['error_suffix']
    sample_file = config_dict['sample_file']
    out_prefix = config_dict['out_prefix']
    sample_list = get_sample_list(sample_file)
    get_error_set(root_dir, error_suffix, sample_list, out_prefix)