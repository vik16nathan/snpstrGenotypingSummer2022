# SNPSTR Genotyping

Goal: to integrate and improve upon existing bioinformatics software to accurately genotype and detect STRs and SNPs within NGS data from PCR-amplified DNA, with the ultimate goal of using SNPSTRs to detect illegal animal trade.


# Steps to STR genotyping pipeline
1. Initial data processing/cleaning (see Analysis.sh)
2. Fasta/Fastq Input File Generation, separated by primer and sample
3. Config file creation (see link: https://drive.google.com/file/d/1xQ2A38eSKoq75_ttMmsK4Yl4AuyFVFVQ/view?usp=sharing)
    * .config files created using multiple flanking region options often have too many flanks/incorrect flanks that contain the repeat motif
    * Created four files within fix_config_files to help with this
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

