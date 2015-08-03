#!/bin/bash
# this script takes removes high coverage regions and regions containing existing transposons from the regions to be available for RSVSim
# USE: list_for_RSVSim.py

chromsomes=/lscr2/andersenlab/kml436/git_repos2/Transposons/files/chromosomes.bed
existing_transposons=/lscr2/andersenlab/kml436/ID_transposable_element.bed
high_coverage_regions=/lscr2/andersenlab/kml436/sv_files/N2_high_coverage_regions


bedtools subtract -a $chromsomes -b $existing_transposons > intermediate.bed
bedtools subtract -a intermediate.bed -b $high_coverage_regions > RSVSim.bed

rm intermediate.bed