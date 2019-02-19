#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import operator
import argparse
import collections


def read_config_file(config_file):
    config_dict = {}
    with open(config_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            config_dict[line.split(' = ')[0]] = line.split(' = ')[1]
    return config_dict
    
def get_hotspot_list(hotspot_file):
    hotspot_dict = {}
    with open(hotspot_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if re.match(r'^chr.*', line):
                line = line.strip()
                chrom, pos, ref, alt, gene, aa = line.split('\t')
                hotspot_dict['{}_{}_{}_{}'.format(chrom, pos, ref, alt)] = [gene, aa]
    return hotspot_dict

def error_statistics(error_file, hotspot_dict):
    bin_dict, err_dict, hot_dict = {}, {}, {}
    type_list = ['T>G', 'A>C', 'T>C','A>G', 'T>A', 'A>T', 'G>T', 'C>A', 'G>C', 'C>G', 'G>A', 'C>T']
    type_dict = {}.fromkeys(type_list, 0)
    bin_chr, bin_start, bin_end = '', 0, 0
    with open(error_file, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            line = line.strip()
            chrom, pos, ref, alt, sample, dp, vd, af, af_list = line.split('\t')
            if i == 1:
                bin_chr, bin_start, bin_end = chrom, int(pos), int(pos) + 19
                bin_region = '{}_{}_{}'.format(bin_chr, bin_start, bin_end)
                if bin_region not in bin_dict:
                    bin_dict[bin_region] = 0
                bin_dict[bin_region] += int(vd)
            if i > 1:
                if chrom == bin_chr and int(pos) >= bin_start and int(pos) <= bin_end:
                    bin_dict[bin_region] += int(vd)
                else:
                    bin_chr, bin_start, bin_end = chrom, int(pos), int(pos) + 19
                    bin_region = '{}_{}_{}'.format(bin_chr, bin_start, bin_end)
                    if bin_region not in bin_dict:
                        bin_dict[bin_region] = 0
                    bin_dict[bin_region] += int(vd)
            if i > 0:
                if len(ref) == 1 and len(alt) == 1:
                    type_dict['{}>{}'.format(ref, alt)] += int(vd)
                if int(sample.split('|')[2]) not in err_dict:
                    err_dict[int(sample.split('|')[2])] = 0
                err_dict[int(sample.split('|')[2])] += 1
                var = '{}_{}_{}_{}'.format(chrom, pos, ref, alt)
                if var in hotspot_dict:
                    gene, aa = hotspot_dict[var][0], hotspot_dict[var][1]
                    new_af_list = []
                    for str_af in af_list.split(','):
                        new_af_list.append(float(str_af))
                    hot_dict['{}_{}_{}'.format(var, gene, aa)] = new_af_list
    return bin_dict, err_dict, hot_dict, type_dict

def generate_html_report(raw_bin_dict, raw_err_dict, raw_hot_dict, raw_type_dict, dedup_bin_dict, dedup_err_dict, dedup_hot_dict, dedup_type_dict, out_dir, echarts_js, datatool_js, baseline_html):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if not os.path.exists('{}/baseline.html'.format(out_dir)):
        os.system('cp {} {}'.format(baseline_html, out_dir))
    if not os.path.exists('{}/echats.js'.format(out_dir)):
        os.system('cp {} {}'.format(echarts_js, out_dir))
    if not os.path.exists('{}/dataTool.js'.format(out_dir)):
        os.system('cp {} {}'.format(datatool_js, out_dir))
    r = open('{}/data.js'.format(out_dir), 'w')

    bin_err_xaxis = sorted(list(set(list(raw_bin_dict.keys()) + list(dedup_bin_dict.keys()))))
    raw_bin_yaxis = []
    dedup_bin_yaxis = []
    for bin_region in bin_err_xaxis:
        if bin_region in raw_bin_dict:
            raw_bin_yaxis.append(raw_bin_dict[bin_region])
        else:
            raw_bin_yaxis.append(0)
        if bin_region in dedup_bin_dict:
            dedup_bin_yaxis.append(dedup_bin_dict[bin_region])
        else:
            dedup_bin_yaxis.append(0)
    
    err_dis_xaxis = sorted(list(set(list(raw_err_dict.keys()) + list(dedup_err_dict.keys()))))
    raw_err_yaxis = []
    dedup_err_yaxis = []
    for error in err_dis_xaxis:
        if error in raw_err_dict:
            raw_err_yaxis.append(raw_err_dict[error])
        else:
            raw_err_yaxis.append(0)
        if error in dedup_err_dict:
            dedup_err_yaxis.append(dedup_err_dict[error])
        else:
            dedup_err_yaxis.append(0)

    hot_af_xaxis = sorted(list(set(list(raw_hot_dict.keys()) + list(dedup_hot_dict.keys()))))
    raw_hot_yaxis = []
    dedup_hot_yaxis = []
    for var in hot_af_xaxis:
        if var in raw_hot_dict:
            raw_hot_yaxis.append(raw_hot_dict[var])
        else:
            raw_hot_yaxis.append([])
        if var in dedup_hot_dict:
            dedup_hot_yaxis.append(dedup_hot_dict[var])
        else:
            dedup_hot_yaxis.append([])

    type_list = ['T>G', 'A>C', 'T>C','A>G', 'T>A', 'A>T', 'G>T', 'C>A', 'G>C', 'C>G', 'G>A', 'C>T']
    raw_type_data = []
    dedup_type_data = []
    for i, err_type in enumerate(type_list):
        raw_type_data.append([0, i, raw_type_dict[err_type]])
        dedup_type_data.append([0, i, dedup_type_dict[err_type]])
    raw_max = max(list(raw_type_dict.values()))
    raw_min = min(list(raw_type_dict.values()))
    dedup_max = max(list(dedup_type_dict.values()))
    dedup_min = min(list(dedup_type_dict.values()))
    
    r.write('var bin_err_xaxis = {}\n'.format(bin_err_xaxis))
    r.write('var raw_bin_yaxis = {}\n'.format(raw_bin_yaxis))
    r.write('var dedup_bin_yaxis = {}\n'.format(dedup_bin_yaxis))
    r.write('var err_dis_xaxis = {}\n'.format(err_dis_xaxis))
    r.write('var raw_err_yaxis = {}\n'.format(raw_err_yaxis))
    r.write('var dedup_err_yaxis = {}\n'.format(dedup_err_yaxis))
    r.write('var hot_af_xaxis = {}\n'.format(hot_af_xaxis))
    r.write('var raw_hot_yaxis = {}\n'.format(raw_hot_yaxis))
    r.write('var dedup_hot_yaxis = {}\n'.format(dedup_hot_yaxis))
    r.write('var type_list = {}\n'.format(type_list))
    r.write('var raw_type_data = {}\n'.format(raw_type_data))
    r.write('var raw_min = {}\n'.format(raw_min))
    r.write('var raw_max = {}\n'.format(raw_max))
    r.write('var dedup_type_data = {}\n'.format(dedup_type_data))
    r.write('var dedup_min = {}\n'.format(dedup_min))
    r.write('var dedup_max = {}\n'.format(dedup_max))
    
    r.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='error set statistics')
    parser.add_argument('-c', '--config', help='config file')
    args = parser.parse_args()

    config_file = args.config
    config_dict = read_config_file(config_file)
    hotspot_file = config_dict['hotspot_file']
    raw_error_file = config_dict['raw_error_file']
    dedup_error_file = config_dict['dedup_error_file']
    out_dir = config_dict['out_dir']
    baseline_html = config_dict['baseline_html']
    echarts_js = config_dict['echarts_js']
    datatool_js = config_dict['datatool_js']

    hotspot_dict = get_hotspot_list(hotspot_file)
    raw_bin_dict, raw_err_dict, raw_hot_dict, raw_type_dict = error_statistics(raw_error_file, hotspot_dict)
    dedup_bin_dict, dedup_err_dict, dedup_hot_dict, dedup_type_dict = error_statistics(dedup_error_file, hotspot_dict)
    generate_html_report(raw_bin_dict, raw_err_dict, raw_hot_dict, raw_type_dict, dedup_bin_dict, dedup_err_dict, dedup_hot_dict, dedup_type_dict, out_dir, echarts_js, datatool_js, baseline_html)