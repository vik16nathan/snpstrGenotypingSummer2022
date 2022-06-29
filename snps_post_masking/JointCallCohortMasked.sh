#!/bin/bash -i

# conda activate gatk

Ref=$1

 gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R $Ref \
   -V gendb://P-pyrhulla-masked \
   -O P_pyrhulla_masked.vcf.gz \

# conda deactivate gatk
