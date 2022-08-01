#!/bin/bash -i
./VariantsToTableMerged.sh
./filteringMerged.sh
./CalcGenotypePostMerged.sh
./FilterGenotypesMerged.sh
./VariantAnnotatorMerged.sh
