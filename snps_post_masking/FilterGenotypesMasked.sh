#!/bin/bash -i

#Question - shoud I mask the reference fasta?
   gatk VariantFiltration \
   -R Pyrhulla_pyrhulla_PrimerRef.fasta \
   -V P_pyrhulla_masked_post.vcf.gz \
   -O P_pyrhulla_masked_filteredG.vcf.gz \
   --genotype-filter-expression "GQ < 20" \
   --genotype-filter-name "my_filters"
