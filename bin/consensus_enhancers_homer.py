#!/usr/bin/python3
#
# consensus_enhancers_homer.py
#
# ####################################################################
#
#      CopyLeft (C) 2022 - Computational Genomics Lab @ UB
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# ####################################################################
#

import sys
import csv

MIN_CHIP_EVIDENCES = 1
MIN_ATAC_EVIDENCES = 1
MIN_TOTAL_EVIDENCES = 2
ATAC_EVIDENCE_FILES = set([
    "peaks/atac.anterior-12h-1.narrowPeak.bed",
    "peaks/atac.anterior-12h-2.narrowPeak.bed",
    "peaks/atac.anterior-12h-notum.narrowPeak.bed",
    "peaks/atac.posterior-0h.narrowPeak.bed",
    "peaks/atac.posterior-12h-1.narrowPeak.bed",
    "peaks/atac.posterior-12h-2.narrowPeak.bed",
    "peaks/atac.posterior-12h-wnt1.narrowPeak.bed",
    "peaks/atac.posterior-48h.narrowPeak.bed"
])

CHIP_EVIDENCE_FILES = set([
    "peaks/chip.anterior-1.narrowPeak.bed",
    "peaks/chip.anterior-2.narrowPeak.bed",
    "peaks/chip.posterior-1.narrowPeak.bed",
    "peaks/chip.posterior-2.narrowPeak.bed"
])


def count_evidences(row, evidence_files):
    evidences = 0
    files_for_peak = row["Parent files"]
    files_for_peak = set(files_for_peak.split("|"))

    evidences = len(evidence_files.intersection(files_for_peak))
    return evidences

def classify_peak(atac_evidences, chip_evidences):
    category = "Unclassified"
    total_evidences = chip_evidences + atac_evidences
    if chip_evidences < MIN_CHIP_EVIDENCES:
        # Not an enhancer
        if atac_evidences >= MIN_ATAC_EVIDENCES:
            category = "ATAC_ONLY"
    else:
        # Possible enhancer
        if atac_evidences >= MIN_ATAC_EVIDENCES:
            category = "ATAC+CHIP"
        else:
            category = "CHIP_ONLY"
    
    if total_evidences < MIN_TOTAL_EVIDENCES:
        category = "UNKNOWN"
    return category

def get_intersection(row):
    return row["Parent files"]


def main():
    peak_file = sys.argv[1]

    category_count = {
        "UNKNOWN": 1,
        "ATAC_ONLY": 1,
        "ATAC+CHIP": 1,
        "CHIP_ONLY": 1
    }
    with open(peak_file) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            atac_evidences = count_evidences(row, ATAC_EVIDENCE_FILES)
            chip_evidences = count_evidences(row, CHIP_EVIDENCE_FILES)
            
            peak_category = classify_peak(atac_evidences, chip_evidences)
            peak_intersection = get_intersection(row)
            peak_id = peak_category + "_" + str(category_count[peak_category])
            category_count[peak_category] += 1
            print("{}\t{}\t{}\t{}\t{}\t+\t{}\n".format(row["chr"], row["start"], row["end"], peak_id, peak_category, peak_intersection))
    

main()
