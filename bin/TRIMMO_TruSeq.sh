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
# USAGE:
#
#   TRIMMO_TruSeq.sh <set> <filename> <idir> <odir>
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
