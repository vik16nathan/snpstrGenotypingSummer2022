# SNPSTR Genotyping

The scripts in this repository represent the work I've done in the past two and a half months to streamline the process of jointly genotyping SNPs and STRs contained with roughly 30 forward and reverse primers, where each set of forward and reverse primers represents a 300 b.p. region where a STR and possibly multiple SNPs are found.

The ultimate goal of this pipeline is to automate the detection of SNPSTRs for any number of samples of a particular animal species in Annika Mozer's FOGS (Forensic Genetics for Species Protection) project, which will hopefully save Annika and other researchers invaluable time. Although the current pipeline isn't perfect, it has proven to be ~80% accurate at STR genotyping for *Pyrrhula pyrrhula* and quite good at assigning SNP alternative alleles to the STR primers, samples, and alleles to which they belong, particularly for primers with lots of functioning reads. 

There are three main parts of the pipeline:
1. Processing reads prior to STRait Razor
2. From the processed and separated reads in (1), making, improving, and using .config files containing the proper flanks to genotype STRs based on the most frequent copy numbers out of all processed reads for a primer-sample pair
3. Merge .fastq files separated by primer, sample, and the alleles identified in (2), find SNPs that correspond to each primer, sample, and allele, and then consolidate and filter results into the desired SNPSTR table format

All scripts in the final_scripts directory should be **copied into one directory** (with NO subdirectories); this working directory should already contain the following inputs (also not within subdirectories):
* A list of sample IDs (should be in the form of speciesPrefix_sample_names.lst or any .lst file; one sample name per line and no trailing whitespace lines)
* A set of unmerged, untrimmed Illumina MySeq reads corresponding to the PCR products in the format speciesPrefix-${sample_number}_S1_L001_R1_001.fastq and speciesPrefix-${sample_number}_S1_L001_R2_001.fastq (make sure to also remember what speciesPrefix is!!)
* A reference file containing a sequence for each primer (in .fasta file format)
* The .conda environment fullSNPSTREnv (to aid in downloading packages/programs needed to run everything)
    * Note: .conda was very problematic with GATK - a manual installation and alias may be needed (see instructions in Appendix below)

# To run the full SNPSTR Pipeline:
1. ./processReadsBeforeSTRaitRazor.sh $speciesPrefix $referenceGenome $sampleListFile (look within the process_reads_before_strait_razor subdirectory within GitLab)
2. ./makeConfigAndRunAllSTRaitRazor.sh $referenceGenome $sampleListFile (within config_and_strait_razor directory)
    * Note: there are some leftover R scripts and Jupyter Notebooks in the main directory of the repository that allow an organized comparison of STR genotypes and excel table genotypes for samples with manually-computed genotypes. The order of these scripts is:
        * a. Rscript getExcelTableGenotypes.R
        * b. Rscript compareSTRaitWithExcel.r
        * c. STRait_Razor_Allele_Mismatch_Evaluation.ipynb (need to copy files and modify path names)
3. Within the snps_with_merging directory, run four separate scripts, since there's some things you need to do by hand prior to running the next step (look within these files for more information):
    * a. ./findSNPsWithMergingPart1.sh $sampleListFile (within part_1)
    * b. ./findSNPsWithMergingPart2.sh $sampleListFile $genomics-db_workspace_path
    * c.  ./findSNPsWithMergingPart3.sh
    * d. ./organizeSNPSTRResults.sh

The final output tables should be stored in the finalExcelOutputs subdirectory within the main working directory. Each table contains SNPs and STRs for one primer and all samples that were able to be genotyped (in .tsv format). These tables closely resemble the format that needs to be uploaded to the FOGS database.

# BIG ISSUES
* Determining STR and SNP zygosity using a certain threshold for the proportion of reads
* There's barely any indel handling within the .config file generation process
* Not all reads are the same length after merging, meaning that GATK finds many STR pieces along with SNPs

# SMALLER ISSUES

