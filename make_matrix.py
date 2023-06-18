#!/usr/bin/env python

import sys
import re
import gzip
from collections import defaultdict

vcfs = sys.argv[1:]

#print(vcfs)

def parse_variant(line, giab=False):
    content = line.strip().split('\t')
    variant_id = f'{content[0]}:{content[1]}:{content[3]}>{content[4]}'
    if not giab:
        return(variant_id, '1', 'NA')
    else:
        if 'CALL=TP' in line or 'BASE=FN' in line:
            truth = 'YES'
        elif 'CALL=FP' in line:
            truth = 'NO'
        else:
            return None
        if 'CALL=TP' in line or 'CALL=FP' in line:
            return (variant_id, '1', truth)
        else:
            return (variant_id, '0', truth)


fields = ['REGIONS', 'TRUTH', 'CLAIR3_STANDART', 'DV_STANDART', 'FB_STANDART', 'HC_1DCNN', 'HC_2DCNN',
        'HC_HARDFILTER', 'OCTOPUS_FOREST', 'OCTOPUS_STANDART', 'STRELKA_STANDART']

variant_data = defaultdict(dict)
for vcf in vcfs:
#    print(vcf)
    giab = 'HG00' in vcf
    regions_set = re.findall('nonmaster_[^\.]+', vcf)[0] if 'nonmaster' in vcf else 'master'
    if giab:
        vcf_prefix = re.findall('HG00\d_[EXGEN]+OME_[A-Z0-9]+_[A-Z0-9]+_[A-Z0-9]+', vcf)[0]
    else:
        vcf_prefix = re.findall('[NARUSZ1245780]+.*_[A-Z0-9]+_[A-Z0-9]+', vcf)[0]
#    else:
#        vcf_prefix = re.findall('wes_\d+\.sample_Almaz\d+.*_[A-Z0-9]+_[A-Z0-9]+', vcf)[0]
    caller_filter = '_'.join(vcf_prefix.split('_')[-2:])
    sample_name = '_'.join(vcf_prefix.split('_')[:-2]).replace('_BWA', '')
    with gzip.open(vcf, 'rt') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue
            if not giab and 'PASS' not in line:
                continue
            parsed = parse_variant(line, giab=giab)
            if parsed is None:
                continue
            if parsed[0] not in variant_data[sample_name]:
                variant_data[sample_name][parsed[0]] = {k: '0' for k in fields}
                variant_data[sample_name][parsed[0]]['TRUTH'] = parsed[2]
                variant_data[sample_name][parsed[0]]['REGIONS'] = regions_set
            variant_data[sample_name][parsed[0]][caller_filter] = parsed[1]


fields_cat = '\t'.join(fields)
print(f'sample\tgiab\tvar_id\t{fields_cat}')

for sn, sdata in variant_data.items():
    for var_id, var_data in sdata.items():
        fmt_data = '\t'.join([var_data[k] for k in fields])
        giab_flag = 'HG00' in sn
        out_line = f'{sn}\t{giab_flag}\t{var_id}\t{fmt_data}'
        print(out_line)



