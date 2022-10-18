#
# compute_motif_location.py
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
from collections import defaultdict

def compute_center(start, end):
    return int(int(start) + ((int(end) - int(start)) / 2))

def read_peaks_file(peaks_file):
    peak_coordinates = {}
    peak_counter = {}
    with open(peaks_file, "r") as peaks_fh:
        for line in peaks_fh:
            line = line.strip()
            if not line:
                continue
            chromosome, start, end, name, *_ = line.split("\t")
            dup_suffix = ""
            if name not in peak_counter:
                peak_counter[name] = 1
            else:
                dup_suffix = "-Dup" +  str(peak_counter[name])
                peak_counter[name] += 1
            name = name + dup_suffix
            peak_coordinates[name] = {}
            peak_coordinates[name]["chromosome"] = chromosome
            peak_coordinates[name]["start"] = start
            peak_coordinates[name]["end"] = end
            peak_coordinates[name]["center"] = compute_center(start, end)
    return peak_coordinates

def compute_motif_coordinates(peak_info, offset, strand, seq_length):
    motif_start = int(peak_info["center"]) + int(offset)
    if strand == "-":
        motif_start -= int(seq_length) - 1
    return motif_start, int(motif_start) + int(seq_length)

def read_and_print_locations(prefix, motif_file, peak_coordinates):
    motif_counter = defaultdict(int)
    with open(motif_file, "r") as motif_fh:
        next(motif_fh) # Skip header
        for line in motif_fh:
            line = line.strip()
            if not line:
                continue
            name, offset, seq, motif, strand, score = line.split("\t")
            fullname = motif
            motif, *_ = motif.split("/") 

            motif_counter[motif] += 1
            seq_length = len(seq)
            motif_start, motif_end = compute_motif_coordinates(peak_coordinates[name], offset, strand, seq_length)
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tID={};ELEMENT={};SEQ={};FULLNAME={}".format(
                peak_coordinates[name]["chromosome"], 
                "HOMER", 
                "TF_binding_site", 
                motif_start, 
                motif_end, 
                score, 
                strand, 
                ".", 
                "tf_" + motif +  "_" + prefix + "_" + str(motif_counter[motif]),
                name,
                seq,
                fullname
                )
            )
            

def main():
    if len(sys.argv) != 4:
        sys.stderr.write("Must provide 3 arguments: experiment prefix, peaks file, motif location file.")
        sys.exit(1)
    
    prefix = sys.argv[1]
    peaks_file = sys.argv[2]
    motif_file = sys.argv[3]

    peak_coordinates = read_peaks_file(peaks_file)
    read_and_print_locations(prefix, motif_file, peak_coordinates)


if __name__ == "__main__":
    main()
