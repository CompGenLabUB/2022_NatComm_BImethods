#!/bin/bash
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

#SET THE PATHS TO TOOLS WITHIN THE PIPELINE...
PATH=$PATH:$HOME/bin
export PATH=$ATACBIN/python/bin/python:$PATH
export PYTHONPATH=$ATACBIN/python/bin/python
export PATH=$HOME/.local/bin:$PATH
export PATH=$ATACBIN/python/bin:$PATH
export PATH=$ATACBIN/bowtie-1.0.0:$PATH
export PATH=$ATACBIN/samtools-0.1.18:$PATH
export PATH=$ATACBIN/bedtools2/bin:$PATH
export PATH=$ATACBIN/fseq/bin:$PATH
export PATH=$ATACBIN:$PATH
export PATH=$ATACBIN/ghostscript-9.20/bin:$PATH
export PATH=$ATACBIN/weblogo:$PATH
export PATH=$ATACBIN/homer/bin:$PATH
export PATH=$ATACBIN/ucsc_tools:$PATH
export PATH=$ATACBIN/symlinks:$PATH

#SET UP DATA LABELS...
RUN_NAME=""
RUN_DIRECTORY="" 
R1=""
R2=""
BOWTIE_REFERENCE=""
FASTA_REFERENCE=""
CHROM_SIZES=""
MACS2_GENOME_SIZE=""

#SET UP OPTIONS...
echo >&2
echo "Running ATAC-seq Pipeline version 1.0.1 started `date`" >&2
echo "Settings:" >&2
echo  >&2
while getopts ":n:o:1:2:b:g:s:chm:u:" opt;
 do {
  case $opt in
    n)
        echo "-n $OPTARG" >&2
        RUN_NAME=$OPTARG
        ;;
    o)
        echo "-o $OPTARG" >&2
        RUN_DIRECTORY=$OPTARG
        ;;
    1)
        echo "-1 $OPTARG" >&2
        R1=$OPTARG
        ;;
    2)
        echo "-2 $OPTARG" >&2
        R2=$OPTARG
        ;;
    b)
        echo "-b $OPTARG" >&2
        BOWTIE_REFERENCE=$OPTARG
        ;;
    g)
        echo "-g $OPTARG" >&2
        FASTA_REFERENCE=$OPTARG
        ;;
    s)
        echo "-s $OPTARG" >&2
        CHROM_SIZES=$OPTARG
        ;;
    m)
        echo "-m $OPTARG" >&2
        MACS2_GENOME_SIZE=$OPTARG
        ;;

    h)
        echo "ATAC Pipeline Usage:" >&2
        echo >&2
        echo "pipeline.sh -n run_name -o output_dir -1 r1.fastq -2 r2.fastq -b bowtie_index -g reference.fa -s chrom.sizes -m mm" >&2
        echo >&2
        echo "-h help (this print out)" >&2
        echo "-n run name e.g. treatmentXYZ" >&2
        echo "-o output directory e.g. /path/myoutputdirectory" >&2
        echo "-1 fastq file containing read1 e.g. /path/R1.fastq" >&2
        echo "-2 fastq file containing read2 e.g. /path/R2.fastq" >&2
        echo "-b bowtie index e.g. /path/index (index must be the prefix only)" >&2
        echo "-g reference genome in fasta format e.g. /path/reference.fa (Note: thish should be the reference used to build the bowtie-index)" >&2
        echo "-s chromosome sizes e.g. /path/chrom.sizes" >&2
        echo "-m macs2 prebuild genome size. Use either hs, mm, ce or dm where hs=2.7e9, mm=1.87e9, ce=9e7 and dm=1.2e8" >&2
        exit 0
        ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :)
        echo "Option -$OPTARG requires an argument." >&2
        exit 1
        ;;
    esac
 }; done
echo >&2

#START PROCESSING DATA...

#SET UP RUNNING DIRECTORY...
mkdir ${RUN_DIRECTORY} #make the run directory.

cd ${RUN_DIRECTORY} #move into the run directory.
export RESULTS_DIR=${RUN_DIRECTORY}/results

#STEP 3 run alignment using bowtie
echo "STEP 3 - running bowtie" >&2
mkdir ${RESULTS_DIR}
bowtie --chunkmbs 1024 -y -v 2 -p 7 \
       --best --strata -m 3 -k 1 --sam-nohead \
       --sam ${BOWTIE_REFERENCE} -1 ${R1} -2 ${R2} \
     > ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam
echo "done" >&2
echo >&2
wait

#STEP 4a convert bowtie sam to bam
echo "STEP 4a - converting sam to bam" >&2
samtools "view" -b -T ${FASTA_REFERENCE} \
         -o ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam.bam \
         -S ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam
echo "done" >&2
echo >&2
wait

#STEP 4b sort the bam
echo "STEP 4b - sorting the bam" >&2
samtools "sort" ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam.bam \
                ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam.bam.sorted
echo "done" >&2
echo >&2
wait

#STEP 4c count aligned reads in sorted bam
# (> "${RESULTS_DIR}/samtools_flagstat.txt")
echo "STEP 4c - counting reads" >&2
samtools "flagstat" ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam.bam.sorted.bam
echo "done" >&2
echo >&2
wait


#STEP 6 filter to mapped reads only
echo "STEP 6 - filtering to mapped reads only" >&2
samtools "view" -b -F 4 \
      -o ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam.bam.sorted.filtered.bam \
         ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam.bam.sorted.bam
echo "done" >&2
echo >&2
wait


#STEP 7 filter nucleosome free region (NF)
echo "STEP 7 - filtering to NF region" >&2
bamtools filter -in ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam.bam.sorted.filtered.bam \
               -out ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam.bam.sorted.filtered_NF.bam \
               -insertSize "<=100"
echo "done" >&2
echo >&2
wait


#STEP 8 Convert processed bam to bed so formated for fseq peak caller
echo "STEP 8 - converting from bam to bed for use with fseq" >&2
bedtools bamtobed -i ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam.bam.sorted.filtered_NF.bam \
                   > ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam.bam.sorted.filtered_NF.bam.bed
echo "done" >&2
echo >&2
wait

#Step 8a shift coordinates within bed as per ENCODE project (and Marta) for tn5
echo "Step 8a shift coordinates within bed as per ENCODE project (tn5 shift)" >&2
awk '{ if ($6=="+"){  $2=$2+4 }
       else if ($6=="-") { $3=$3-5 };
       print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; }
    ' ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam.bam.sorted.filtered_NF.bam.bed \
    > ${RESULTS_DIR}/r1_r2_${RUN_NAME}.sam.bam.sorted.filtered_NF.bam.bed.tn5mofidied.bed
echo "done" >&2
echo >&2
wait


#Extra step using MACS2
echo "STEP EXTRA - peak calling using macs2" >&2
docker run -v ${RESULTS_DIR}:/data \
    fooliu/macs2 callpeak                           \
              -t /data/r1_r2_${RUN_NAME}.sam.bam.sorted.filtered_NF.bam.bed.tn5mofidied.bed \
              -n /data/${RUN_NAME}peaks_macs2 -f BED -g 773939492   \
              --outdir /data/                       \
           2> ${RESULTS_DIR}/${RUN_NAME}_peaks_macs2.docker_macs2.log
echo "done" >&2
echo >&2

#PIPELINE END
echo "Pipeline complete - `date`" >&2
exit

# cd $ATACBIN
# ./ENCODE_ATAC_PIPELINE_svr4_tn5_shift.sh -n [name of job] -o [path:output directory] \
#                                          -1 [path:fastq R1] -2 [path:fastq R2] \
#                                          -b [path:bowtie1 genome index] -g [path:fasta genome] -s [path:chrom sizes file] \
#                                          -m mm -u http://leightonucsc.cmp.uea.ac.uk 2>&1 | tee [path to log file]
