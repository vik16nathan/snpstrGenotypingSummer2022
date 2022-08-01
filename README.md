# SNPSTR Genotyping

The scripts in this repository represent the work I've done in the past two months to streamline the process of jointly genotyping SNPs and STRs contained with roughly 30 forward and reverse primers, where each set of forward and reverse primers represents a 300 b.p. region where a STR and possibly multiple SNPs are found.

The ultimate goal of this pipeline is to automate the detection of SNPSTRs for any number of samples of a particular animal species in Annika Mozer's FOGS (Forensic Genetics for Species Protection) project, which will hopefully save Annika and other researchers invaluable time. Although the current pipeline isn't perfect, it has proven to be ~80% accurate at STR genotyping for *Pyrrhula pyrrhula* and quite good at assigning SNP alternative alleles to the STR primers, samples, and alleles to which they belong, particularly for primers with lots of functioning reads. 

There are three main parts of the pipeline:
1. Processing reads prior to STRait Razor
2. From the processed and separated reads in (1), making, improving, and using .config files containing the proper flanks to genotype STRs based on the most frequent copy numbers out of all processed reads for a primer-sample pair
3. Merge .fastq files separated by primer, sample, and the alleles identified in (2), find SNPs that correspond to each primer, sample, and allele, and then consolidate and filter results into the desired SNPSTR table format

All scripts in the final_scripts directory should be copied into one directory (with NO subdirectories); this directory should already contain the following information (also not within subdirectories):
* A list of sample IDs (should be in the form of Species-name_sample_names.lst or any .lst file; one sample name per line and no trailing whitespace lines)
* A set of unmerged, untrimmed Illumina MySeq reads corresponding to the PCR products in the format Species-name-${sample_number}_S1_L001_R1_001.fastq and Species-name-${sample_number}_S1_L001_R2_001.fastq
# Steps to STR genotyping pipeline
1. Initial data processing/cleaning (done by Annika)
2. Fasta/Fastq Input File Generation, separated by primer and sample
3. Config file creation (see link: https://drive.google.com/file/d/1xQ2A38eSKoq75_ttMmsK4Yl4AuyFVFVQ/view?usp=sharing)
    * .config files created using multiple flanking region options often have too many flanks/incorrect flanks that contain the repeat motif
    * Created four files within fix_config_files to help with this, with one primer's multiple-flank config file at a time
    * Full description/order of analyses: https://docs.google.com/document/d/1SrBHJA-7HQTlyU8EoovvQSa5ho_orZcSchW9jNAJpww/edit?usp=sharing

4. Run STRait razor/filter output
    * runAllSTRaitRazor.sh
    * processMultipleFlanksSTRaitRazorResults.r
5. Gather results into one table of STRait Razor Genotypes (homozygous or heterozygous, with full genotypes of alleles)
    * getSTRaitRazorGenotypes.r
6. Compare to Excel table
    * getExcelTableGenotypes.r
    * compareSTRaitWithExcel.r

7. Convert genotyping errors to a .txt file for easy viewing/visualize copy number errors/post-processing
    * STRait_Razor_Allele_Mismatch_Evaluation.ipynb

8. Compare old results using one-line .config files (PERF_tsv_to_STRaitRazor_config.r) with new results using complex config file creation - see which primers are problematic and try to fix issues
    * Compare_Old_and_New_STRait_Razor_Results.ipynb

##BIG ISSUES
* Determining STR and SNP zygosity using a certain threshold for the proportion of reads
* There's barely any indel handling within the .config file generation process
* Not all reads are the same length after merging, meaning that GATK finds many STR pieces along with SNPs
