#!/bin/bash -i
gatk VariantFiltration \
   -R PrimerRef.fasta \
   -V MergedVariantsBeforeFiltering_post.vcf.gz \
   -O MergedVariants_filteredG.vcf.gz \
   --genotype-filter-expression "GQ < 20" \
   --genotype-filter-name "my_filters"