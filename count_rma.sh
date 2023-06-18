#!/bin/bash

grep -P '^#|AF=(1.00000e+00$|[5-9]\.\d+e-01)' gnomad_all.vcf > gnomad_hg38_rma.vcf
grep -P '^#|AF=(1.00000e+00$|9\.9\d+e-01$)' gnomad_hg38_rma.vcf  > gnomad_hg38_rma_1pct.vcf
bedtools intersect -header -wa -a gnomad_hg38_rma.vcf -b gencode_v43.cds.bed | uniq - > gnomad_hg38_rma_coding.vcf

# Count pathogenics
bedtools intersect -wa -a gnomad_hg38_rma_1pct.vcf -b <( zcat ./clinvars/clinvar_20230326.vcf.gz | perl -pe 's|^([\dXYM])|chr\1|' |  grep -P '^#|CLNSIG=Pathogen'  ) | wc -

# Coding pathogenic vars
bedtools intersect -header -wa -a gnomad_hg38_rma_1pct.vcf -b <( zcat ./clinvars/clinvar_20230326.vcf.gz | perl -pe 's|^([\dXYM])|chr\1|' |  grep -P '^#|CLNSIG=Pathogen'  ) | bedtools intersect -wa -a - -b gencode_v43.cds.bed

# Coding P and LP vars
bedtools intersect -header -wa -a gnomad_hg38_rma_1pct.vcf -b <( zcat ./clinvars/clinvar_20230326.vcf.gz | perl -pe 's|^([\dXYM])|chr\1|' |  grep -P '^#|CLNSIG=(Pathogen|Likely_path)'  ) | bedtools intersect -wa -a - -b gencode_v43.cds.bedl


# T2T
grep -P '^#|AF=(1.00000e+00$|[5-9]\.\d+e-01)' gnomad_all_t2t.vcf > gnomad_t2t_rmaKept.vcf &
grep -P '^#|MismatchedRefAllele\tAF=([1-4]\.\d+e-01|\d\.\d+e-0[2-9])' gnomad_all_t2tReject.vcf > gnomad_t2t_rmaNovel.vcf
cat <( grep -P '^#' gnomad_t2t_rmaKept.vcf ) <( cat <( grep -vP '^#' gnomad_t2t_rmaKept.vcf ) <( grep -vP '^#' gnomad_t2t_rmaNovel.vcf ) | sort -k1,1 -k2,2n ) > gnomad_t2t_rma.vcf

bedtools intersect -header -wa -a gnomad_t2t_rma.vcf -b ../assembly_mf/gencode_v43_t2t.cds.bed | uniq - > gnomad_t2t_rma.coding.vcf
bedtools intersect -header -wa -a gnomad_t2t_rmaKept.vcf -b ../assembly_mf/gencode_v43_t2t.cds.bed | uniq - > gnomad_t2t_rmaKept.coding.vcf
bedtools intersect -header -wa -a gnomad_t2t_rmaNovel.vcf -b ../assembly_mf/gencode_v43_t2t.cds.bed | uniq - > gnomad_t2t_rmaNovel.coding.vcf
