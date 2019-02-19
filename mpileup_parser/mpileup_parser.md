#### Description

<p>Samtools mpileup and parser</p>

#### Usage

```shell
usage: mpileup_parser.py [-h] [-c CONFIG]

mpileup parser

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
                        config file
```

#### Note
1. Config file
```shell
root_dir = /datapool/Analysis/ResDev/RD03_Pan70/NA12878_Error_Set
sub_dir = 04.DEDUP
bam_suffix = dedup.bam
sample_file = 
ref_genome = /datapool/RefData/Genome/hg19/ucsc.hg19.fa
samtools = /datapool/Apps/Production/miniconda2/bin/samtools
bed_file = /datapool/RefData/Research/panel/RD03_Pan70/Sim_V7_3_primary_targets_expand20.bed
out_dir = /datapool/Analysis/ResDev/RD03_Pan70/NA12878_Error_Set/mpileup_parser
```
2. Sample file
```shell
RD2019012101CF
RD2019012102CF
...
RD2019012119CF
RD2019012120CF
```
