#!/bin/bash -i

gatk VariantsToTable \
-V MergedVariantsBeforeFiltering.vcf \
-F CHROM -F POS -F EVENTLENGTH -F TYPE -F REF -F ALT -F MULTI-ALLELIC -F VAR -F NSAMPLES -F NCALLED -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum \
-O filtering.table \
--show-filtered

mv filtering.table filter_test.txt

awk ' $11 >= 2 {print $0 } ' filter_test.txt > filter_test2.txt
awk ' $12 <= 60 {print $0 } ' filter_test2.txt > filter_test3.txt
awk ' $13 <= 3 {print $0 } ' filter_test3.txt > filter_test4.txt
awk ' $14 >= 40 {print $0 } ' filter_test4.txt > filter_test5.txt
awk ' $15 >= -12.5 {print $0 } ' filter_test5.txt > filter_test6.txt
awk ' BEGIN{print"CHROM     POS     EVENTLENGTH     TYPE    REF     ALT     MULTI-ALLELIC   VAR     NSAMPLES        NCALLED   QD      FS      SOR     MQ      MQRankSum       ReadPosRankSum"} $16 >= -8.0 {print $0 } ' filter_test6.txt > final_filtering.table