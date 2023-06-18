#!/bin/bash

for i in *.vcf.gz
do 
	smp=$( echo $i | grep -oP 'HG00\d' )
	/media/ironwolf/callers_proj/rtg-tools-3.12/rtg vcfeval -b /media/ironwolf/callers_proj/concordance_analysis/${smp}*gz --bed-regions /media/ironwolf/callers_proj/concordance_analysis/GRCh37_WES_CDS.bed -c $i -e /media/ironwolf/callers_proj/concordance_analysis/giab/nonmaster_XY.bed -o ${i%%.vcf.gz}.nonmaster_XY.rtg/ -m combine -t /media/ironwolf/callers_proj/GRCh37.primay_assembly.genome.rtg/
done
