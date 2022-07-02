#!/bin/bash -i

 gatk --java-options "-Xmx4g" CalculateGenotypePosteriors \
   -V P_pyrhulla_masked.vcf \
   -O P_pyrhulla_masked_post.vcf.gz


