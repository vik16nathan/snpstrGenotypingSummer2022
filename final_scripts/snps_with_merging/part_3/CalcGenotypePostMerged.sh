#!/bin/bash -i

 gatk --java-options "-Xmx4g" CalculateGenotypePosteriors \
   -V MergedVariantsBeforeFiltering.vcf  \
   -O MergedVariantsBeforeFiltering_post.vcf.gz