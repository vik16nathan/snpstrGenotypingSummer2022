# SNPSTR Genotyping

The scripts in this repository represent the work I've done in the past two and a half months to streamline the process of jointly genotyping SNPs and STRs contained with roughly 30 forward and reverse primers, where each set of forward and reverse primers represents a 300 b.p. region where a STR and possibly multiple SNPs are found.

The ultimate goal of this pipeline is to automate the detection of SNPSTRs for any number of samples of a particular animal species in Annika Mozer's FOGS (Forensic Genetics for Species Protection) project, which will hopefully save Annika and other researchers invaluable time. Although the current pipeline isn't perfect, it has proven to be ~80% accurate at STR genotyping for *Pyrrhula pyrrhula* and quite good at assigning SNP alternative alleles to the STR primers, samples, and alleles to which they belong, particularly for primers with lots of functioning reads. 

There are three main parts of the pipeline:
1. Processing reads prior to STRait Razor (trimming, QC, mapping, separating by primer and sample)
2. From the processed and separated reads in (1), making, improving, and using .config files containing the proper flanks to genotype STRs based on the most frequent copy numbers out of all processed reads for a primer-sample pair
3. Merge .fastq files separated by primer, sample, and the alleles identified in (2), find SNPs that correspond to each primer, sample, and allele, and then consolidate and filter results into the desired SNPSTR table format

All scripts in the final_scripts directory should be **copied into one directory** (with NO subdirectories); this working directory should already contain the following inputs (also not within subdirectories):
* A list of sample IDs (should be in the form of speciesPrefix_sample_names.lst or any .lst file; one sample name per line and no trailing whitespace lines)
* A set of unmerged, untrimmed Illumina MySeq reads corresponding to the PCR products in the format speciesPrefix-${sample_number}_S1_L001_R1_001.fastq and speciesPrefix-${sample_number}_S1_L001_R2_001.fastq (make sure to also remember what speciesPrefix is!!)
* A reference file containing a sequence for each primer (in .fasta file format)
* The .conda environment fullSNPSTREnv (to aid in downloading packages/programs needed to run everything)
    * Note: .conda was very problematic with GATK - a manual installation and alias may be needed (see instructions in Appendix below)


# To run the full SNPSTR Pipeline:
Before starting, activate the conda environment to make sure all software tools required for the analysis are installed. Note that the script may nonetheless say that some tools are missing - if an entire tool is missing, try conda install nameOfTool, or if an R package is missing, see the last bullet point under "Issues."

Type in conda env create -f fullSNPSTREnv.yml, and then conda activate fullSNPSTREnv.

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

# Issues

Let . be the working directory.

* PERF on the reference genome sometimes produces false flanks
    * Problematic output files to look at: ./config/PrimerX_flank10.config
    * Files to be edited: PrimerRef_perf.tsv (maybe trying to fix this by hand based on the 30 lines in the reference genome would eliminate some errors that sliding window could cause. Sliding window should be fixing these issues, but it's nice to have an extra built-in check beforehand)
    * Need to create PrimerRef_fixed.bed anyways for makeFinalOutputTable.r (might as well fix PrimerRef_perf.tsv and then use this to make PrimerRef_fixed.bed)

* There's no indel handling within sliding window, and no guarantee that the first read that contains the flank of interest doesn't have missing data (see logic in findReferencesForIncorrectFlanks.r - the moment we find a read with an incorrect flank, we use it as our reference for sliding window without considering any other options; we could be choosing a bad read)
    * Problematic output files: ./config/PrimerX_multiple_flanks.config
    * Scripts to be edited: findReferencesForIncorrectFlanks.r (to find multiple possible reads that contain a flank); findTrueFlanks.r (with indel handling - see old/incorrect work for some functions that could maybe help with this)
    
* Determining STR and SNP zygosity using a certain threshold for the proportion of reads
    * Ideally - determine this proportion (in place of 0.15) based on what works best for *many* different samples (copy the script you wrote from waay back in here)
    * Problematic output files: strait_razor_genotypes_0.15_multiple_flanks.tsv
    * Files to be edited: 
        * STR zygosity: last line of makeConfigAndRunSTRaitRazor.sh
            * Rscript getSTRaitRazorGenotypes.r 0.15 - replace 0.15 with percentage of all STRait razor genotypes required for an allele to be identified as "significant" (we identified 0.15 as the ideal proportion based on one dataset, but this may not be consistent with other datasets)
            * Keep running excel comparison with STRait razor genotypes generated with various thresholds (see note under part (2) of running the main pipeline)


* Not all reads are the same length after merging, meaning that GATK finds many STR pieces along with SNPs
    * Problematic output files: ./FastqInputs/${sample}_Primer${primer}_allele${allele}_merged.fastq
    * Ways to fix this: Filter out all lines in the above file(s) that are below a certain length, or find the most common read length in these files and extract only reads with this read length 

* Conda and R - Might need to manually install/update R within the conda environment, and install R packages within R console within the .conda environment
    * 1. conda activate fullSNPSTREnv
    * 2. Type "R" into the terminal
    * 3. Within the R console: type install.packages("package_name") for any packages that claim to not be installed


# Old/Incorrect Work
* In fix_config_files: findTrueIncorrectFlanksWithReferences.r	- this was my first attempt at **indel handling**, and there are functions that involve a five base-pair window to determine whether or not our sliding window is the 4 b.p. motif plus or minus a base pair from a 1 b.p. indel. The issue is that four base pairs is so small that there are many windows that are, by chance, one base pair away from being the correct window. Also, findTrueFlanks.r is where sliding window is *actually* used to fix the flanks, so this is what needed to be fixed.
* masking directory: A complete, utter, and irredeemable disaster. Didn't merge reads beforehand, and didn't realize that the 250 b.p. cutoff for Illumina MiSeq reads meant that the forward and reverse reads have different STR starting/ending positions (basically this is all wrong and also way too time-consuming since I used .fasta inputs and searched for headers rather than using .fastq inputs in the FastqInputs directory.)
* snps_post_masking - didn't work because masking didn't work. However, a lot of the logic was reused in the final_scripts/snps_with_merging subdirectories

# Installing GATK by hand
1. https://github.com/broadinstitute/gatk/releases - download .zip file, unzip it, copy unzipped gatk-4.2.6.1 into home directory (~) 
    * Note: check for updates!
2. Go to home ~ directory and edit .bashrc. Add the following line: alias gatk='~/gatk-4.2.6.1/gatk'
3. Close and reopen your terminal.
4. Calling gatk will now invoke the version on your machine

