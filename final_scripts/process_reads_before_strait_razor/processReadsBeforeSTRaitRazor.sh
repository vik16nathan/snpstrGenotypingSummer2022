#!/bin/bash -i

#Trim, QC, and map reads

referenceGenome=$1
sampleListFile=$2
#conda activate qc
./trimming.sh
./fastqc.sh
./multiqc.sh
#conda deactivate qc

#conda activate mapping
./mapping.sh $referenceGenome
./mapping_post.sh
#conda deactivate mapping

#conda activate straitRazorEnv
./makeAllSTRaitRazorDirs.sh
./separateBamsByPrimerAndSample.sh $sampleListFile
./filterPrimerSamplePairsWithTooFewReads.sh $sampleListFile $threshold
./createFastqAllPrimersAllSamples.sh $sampleListFile

#Install STRaitRazor if it hasn't already been installed
#right now I'm editing .bashrc within this script... that's dangerous
./installSTRaitRazor.sh

#Important - make sure to add the following line to your ~/.bashrc file:
#alias str8rzr='~/STRaitRazor/str8rzr'

#conda deactivate straitRazorEnv

