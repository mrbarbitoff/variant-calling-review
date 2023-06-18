#!/bin/bash

# docker run -v ./:/data -it broadinstitute/gatk gatk LiftoverVcf -I /data/gnomad_hg38_rma.vcf -O /data/gnomad_t2t_rma.vcf --REJECT /data/gnomad_t2t_rmaReject.vcf -C /data/grch38-chm13v2.chain -R /data/chm13v2.0.fa

docker run -v ./:/data broadinstitute/gatk gatk LiftoverVcf -I /data/gnomad_all.vcf -O /data/gnomad_all_t2t.vcf --REJECT /data/gnomad_all_t2tReject.vcf -C /data/grch38-chm13v2.chain -R /data/chm13v2.0.fa
