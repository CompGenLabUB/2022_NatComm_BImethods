# 2022_NatComm_BImethods

This page contains the Bioinformatics supplementary material for the
manuscript __Pascual-Carreras et al, Nature Communications, 2022__ (in press);
current publicly available preprint version was posted at
[bioRxiv](https://www.biorxiv.org/content/10.1101/2020.12.08.416008v1).

We really apreciate if you cite this paper when using any of the referenced materials.


## Abstract

For successful regeneration, the identity of the missing tissue must
be specified according to the pre-existing tissue. Planarians are
ideal for the study of the mechanisms underlying this process because
the same field of cells can regrow a head or a tail according to the
missing body part. After amputation, the differential activation of
the WNT/beta-catenin signal specifies anterior versus posterior
identity. Initially, both ''wnt1'' and ''notum'' (Wnt inhibitor) are
expressed in all wounds, but 48 hours later they are restricted to
posterior or anterior facing wounds, respectively, by an unknown
mechanism. Here we show that 12 hours after amputation the chromatin
accessibility of cells in the wound region changes according to the
polarity of the pre-existing tissue in a WNT/beta-catenin-dependent
manner. Genomic analyses suggest that homeobox transcription factors
and chromatin-remodeling proteins are direct WNT/betacatenin targets,
which trigger the expression of posterior effectors. Finally, we
identified FoxG as a ''wnt1'' up-stream regulator, probably via binding
to its first intron enhancer.


## Supplementary Material

Bioinformatics protocols for that manuscript are described, along with
relevant figures and tables, in the main PDF file:

<p align="center">
 <a href="https://compgen.bio.ub.edu/datasets/2022_NatComm/2022_NatComm_BImethods.pdf"
    alt="Computational Supplementary Methods for Eudald et al, 2022.">
    Computational Supplementary Methods</a>
    (184.5MB PDF)
</p>

It contains links to all the scripts and files described there that
are also distributed through this repository. The only exception due
to its file size is the custom <a
href="https://bioconductor.org/packages/release/bioc/html/BSgenome.html"
alt="BSgenome package at Bioconductor site"><tt>BSgenome
Bioconductor</tt></a> package for the <tt>dd_Smes_g4</tt> genome
annotation dataset, which is available as a direct link from our web
server to the gzipped tarball:

<p align="center">
 <a href="https://compgen.bio.ub.edu/datasets/2022_NatComm/BSgenome.Smed.PlanMine.ddSmesg4_4.0.tar.gz"
    alt="ddSmesg4_4.0 curated BSgenome package">
    <tt>BSgenome.Smed.PlanMine.ddSmesg4_4.0.tar.gz</tt></a>
    (169M gzipped tarball)
</p>

You can easily install it after downloading on your local <tt>R</tt>
distribution; just run the commands below on the <tt>shell</tt>
command-line:

```{.sh}
R CMD build BSgenome.Smed.PlanMine.ddSmesg4;
R CMD check BSgenome.Smed.PlanMine.ddSmesg4_4.0.tar.gz;
R CMD INSTALL BSgenome.Smed.PlanMine.ddSmesg4_4.0.tar.gz;
```

We wish you find those materials useful for your planarian research.

    
## CopyLeft

<table width="100%"
       cellspacing="0" cellpadding="0"
       style="border: none;">
<tr>
<td width=150px">
    <img alt="GNU-GPLv3 license logo"
         src="./bin/gplv3.png" width="150px"></td>
<td><p>
Scripts developed at the lab are distributed under
<a href="./bin/gpl-3.0.txt"
   alt="GNU-GPLv3 software license">
   GNU-GPLv3 license</a>.
They are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU General Public License for more details.
</p></td>
</tr>
</table>

<hr />

<p align="center">
2022 - <a href="https://compgen.bio.ub.edu"
          alt="Computational Genomics Lab home page">
          Computational Genomics Lab</a><br />
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"
   alt="This data is licensed under a Creative Commons Attribution 4.0 International License.">
   <img alt="Creative Commons License"
        style="border-width:0" height="30px"
        src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a>
</p>
