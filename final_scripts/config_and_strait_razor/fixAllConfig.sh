#!/bin/bash

sampleListFile=$1
for primer in {1..30}
do

    Rscript findIncorrectFlanks.r $primer
    Rscript findReferencesForIncorrectFlanks.r $primer $sampleListFile
    Rscript findTrueFlanks.r $primer
    Rscript fixConfigFile.r $primer

done
