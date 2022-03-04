#!/bin/bash
#
# TRIMMO_TruSeq.sh <set> <filename> <idir> <odir>
#

PX=$1;
FL=$2;
ID=$3;
OD=$4;

cat <<EOF 1>&2;
###
### Running trimmomatic on $PX
###
#
#   SET: $PX
#  FILE: $FL
#  Idir: $ID
#  Odir: $OD
#
EOF

export TMC=/usr/local/install/Trimmomatic-0.32/;
export TRMPAR="LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:30 TOPHRED33";
export TRMPECLP="ILLUMINACLIP:$TMC/adapters/TruSeq2-PE.fa:2:30:10";
export TRMSGCLP="ILLUMINACLIP:$TMC/adapters/TruSeq2-SE.fa:2:30:10";

java -jar $TMC/trimmomatic-0.32.jar PE \
      ${ID}/${FL}_1.fastq.gz           \
      ${ID}/${FL}_2.fastq.gz           \
      ${OD}/${PX}.r1.trimmope.fastq.gz \
      ${OD}/${PX}.r1.trimmosg.fastq.gz \
      ${OD}/${PX}.r2.trimmope.fastq.gz \
      ${OD}/${PX}.r2.trimmosg.fastq.gz \
      $TRMPECLP $TRMPAR                \
   2> ${OD}/${PX}.pe.trimmo.log 1>&2 ;

exit 0;
