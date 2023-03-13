#!/usr/bin/env zsh


mapdir=results/hisat2

########################################## SECTION IV ###############################################
# remove all .sam files
echo "trashing .sam files..."
rm $mapdir/*.sam