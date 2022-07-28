#!/bin/bash

for primer in {1..30}
do

    Rscript findIncorrectFlanks.r $primer
    Rscript findReferencesForIncorrectFlanks.r $primer
    Rscript findTrueFlanks.r $primer
    Rscript fixConfigFile.r $primer

done
