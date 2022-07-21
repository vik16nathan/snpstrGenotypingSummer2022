#!/bin/bash -i

#IMPORTANT! Determine (by hand) the line number of the first non-header line in the output .vcf table.
#This will be used as an input for making the final output table
tail -n+74 P_pyrhulla_annotated.vcf > P_pyrhulla_annotated_table_only_no_masking.vcf