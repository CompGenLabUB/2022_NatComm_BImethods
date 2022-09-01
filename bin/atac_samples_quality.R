#
#  atac_samples_quality.R
#
#      running all quality analyses provided by ATACseqQC
#      using planarian genome annotation.
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

library(ATACseqQC);
source(system.file("extdata", "IGVSnapshot.R", package = "ATACseqQC"));
library(Rsamtools);
library(GenomicFeatures);
library(BSgenome);
library(BSgenome.Smed.PlanMine.ddSmesg4);
library(ChIPpeakAnno);

genome <- Smed;
seqlev <- names(genome);

possibleTag <- list(
    "integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                "TC", "UQ"), 
  "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                "U2"));

# most confident annotation set
txdb.hgh.filename <- "./refseqs/annotation/smes_v2_hconf_SMESG.sqlite";
txdb.hgh <- loadDb(txdb.hgh.filename);
# larger annotation set
txdb.med.filename <- "./refseqs/annotation/smes_v2_repeatfilt_SMESG.sqlite";
txdb.med <- loadDb(txdb.med.filename);

runatacseqQC <- function(set, bamfile, bamfile.dir, image.dir, stats.dir) {

    bamfile.set   <- set
  # bamfile.dir   <- dirname(bamfile);
    bamFILE       <- paste0(bamfile.dir,"/",bamfile);
    
    bamfile.clean <- paste0(bamfile.dir,"/",
                            sub(".bam", ".clean.bam",
                                basename(bamfile)));

    image.prefix <- paste0(image.dir,'/',set);
    
    #
    cat("# ", set, " : bamQC stats\n", sep="");
    
    bamfile.stats <- bamQC(bamFILE,
                           outPath = bamfile.clean);
    write.table(bamfile.stats$idxstats,
                paste0(stats.dir, '/', set, '_bammappingstats.tbl'),
                sep="\t");
    #     bamQC(
    #        bamfile,
    #        index = bamfile,
    #        mitochondria = "chrM",
    #        outPath = sub(".bam", ".clean.bam", basename(bamfile)),
    #        doubleCheckDup = FALSE
    #      )
    # ---------------------------> must find which sequence is mitochrondria
    #  in our case, sequences of genome set have no significant blastn match

    #
    cat("# ", set, " : library complexity estimation\n", sep="");

    png(paste0(image.prefix,"_librarycomplexity.png"),
        height=6,width=8,units="in",res=600);
    estimateLibComplexity(readsDupFreq(bamFILE));
    dev.off();

    #
    cat("# ", set, " : fragment size distribution\n", sep="");
    
    png(paste0(image.prefix,"_fragmentsizesdist.png"),
        height=6,width=8,units="in",res=600);
    fragSize <- fragSizeDist(bamFILE, bamfile.set);
    dev.off();

    #
    cat("# ", set, " : bamfile tag scanning\n", sep="");
    bamTop100 <- scanBam(BamFile(bamFILE, yieldSize = 100),
                         param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag;
    tags <- names(bamTop100)[lengths(bamTop100)>0];
    tags;

    #
    for (TR.set in c("hgh","med")) {
        
        #
        cat("# ", set, " [ ",TR.set," ] : create dir\n", sep="");
        
        outPath <- paste0(bamfile.dir,"/",
                          sub(".bam", ".splited", basename(bamfile)),
                          "_", TR.set);
        dir.create(outPath);
        
        bamfiles <- file.path(outPath,
                              c("NucleosomeFree.bam",
                                "mononucleosome.bam",
                                "dinucleosome.bam",
                                "trinucleosome.bam"));

        if (TR.set == "hgh") {
            # most confident annotation set
            txdb <- txdb.hgh; 
        } else {
            # larger annotation set
            txdb <- txdb.med;
        };
        
        image.set <- paste0('_',TR.set);
        
        #
        cat("# ", set, " [ ",TR.set," ] : shifting bam coords\n", sep="");
        seqinformation <- seqinfo(txdb);
        which <- as(seqinformation, "GRanges");
        gal <- readBamFile(bamfile.clean,
                           tag=tags, which=which,
                           asMates=TRUE, bigFile=TRUE);
        bamfile.shifted <- file.path(outPath,
                                     sub(".bam", ".shifted.bam",
                                         basename(bamfile.clean)));
        gal1 <- shiftGAlignmentsList(gal, outbam=bamfile.shifted);


        #
        cat("# ", set, " [ ",TR.set," ] : getting transcripts and TSSs\n", sep="");
        txs <- transcripts(txdb);
        TSS <- promoters(txs, upstream=0, downstream=1);
        TSS <- unique(TSS);

        #
        cat("# ", set, " [ ",TR.set," ] : PT score\n", sep="");
        pt <- PTscore(gal1, txs);
        png(paste0(image.prefix,"_PTscore_plot",image.set,".png"),
            height=6,width=8,units="in",res=600);
        plot(pt$log2meanCoverage, pt$PT_score, 
             xlab="log2 mean coverage",
             ylab="Promoter vs Transcript");
        dev.off();

        #
        cat("# ", set, " [ ",TR.set," ] : NFR score\n", sep="");
        nfr <- NFRscore(gal1, txs);
        png(paste0(image.prefix,"_NFRscore_plot",image.set,".png"),
            height=6,width=8,units="in",res=600);
        plot(nfr$log2meanCoverage, nfr$NFR_score, 
             xlab="log2 mean coverage",
             ylab="Nucleosome Free Regions score",
             main="NFRscore for 200bp flanking TSSs",
             xlim=c(-10, 0), ylim=c(-5, 5));
        dev.off();

        #
        cat("# ", set, " [ ",TR.set," ] : TSSE score\n", sep="");
        tsse <- TSSEscore(gal1, txs);
        cat("# ", set, " [ ",TR.set," ] : TSSE score = ",
            tsse$TSSEscore,"\n", sep="");
        png(paste0(image.prefix,"_TSSEscore_plot",image.set,".png"),
            height=6,width=8,units="in",res=600);
        plot(100*(-9:10-.5), tsse$values, type="b", 
             xlab="distance to TSS",
             ylab="aggregate TSS score");
        dev.off();

        #
        cat("# ", set, " [ ",TR.set," ] : shifting, splitting and saving\n", sep="");
        objs <- splitBam(bamfile.clean, tags=tags, outPath=outPath,
                         txs=txs, genome=genome, seqlev=seqlev);

        #
        cat("# ", set, " [ ",TR.set," ] : cum.pct of tag allocation in nucleosome-free and mononucleosome\n", sep="");
        png(paste0(image.prefix,"_cumulativecoverage_nf-mn-dn",image.set,".png"),
            height=6,width=10,units="in",res=600);
        cumulativePercentage(bamfiles[1:3], which); # as(seqinformation["chr1"], "GRanges"))
        dev.off();

        #
        cat("# ", set, " [ ",TR.set," ] : calculating signal enrichment\n", sep="");
        ## estimate the library size for normalization
        (librarySize <- estLibSize(bamfiles));
        ## calculate the signals around TSSs.
        NTILE <- 101;
        dws <- ups <- 1010;
        sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", 
                                             "mononucleosome",
                                             "dinucleosome",
                                             "trinucleosome")], 
                                  TSS=TSS,
                                  librarySize=librarySize,
                                  seqlev=seqlev,
                                  TSS.filter=0.5,
                                  n.tile = NTILE,
                                  upstream = ups,
                                  downstream = dws);
        ## log2 transformed signals
        sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1));

        #
        cat("# ", set, " [ ",TR.set," ] : TSS heatmap\n", sep="");
        png(paste0(image.prefix,"_TSS_heatmap",image.set,".png"),
            height=10,width=5,units="in",res=600);
        featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                              zeroAt=.5, n.tile=NTILE);
        dev.off();
        
        #
        cat("# ", set, " [ ",TR.set," ] : TSS normalized signals\n", sep="");
        ## get signals normalized for nucleosome-free and nucleosome-bound regions.
        png(paste0(image.prefix,"_TSS_featurealndist",image.set,".png"),
            height=6,width=8,units="in",res=600);
        out <- featureAlignedDistribution(sigs, 
                                          reCenterPeaks(TSS, width=ups+dws),
                                          zeroAt=.5, n.tile=NTILE, type="l", 
                                          ylab="Averaged coverage");
        dev.off();
        ## rescale the nucleosome-free and nucleosome signals to 0~1
        range01 <- function(x){(x-min(x))/(max(x)-min(x))};
        out <- apply(out, 2, range01);
        png(paste0(image.prefix,"_TSS_signalfractiondist",image.set,".png"),
            height=6,width=8,units="in",res=600);
        matplot(out, type="l", xaxt="n", 
                xlab="Position (bp)", 
                ylab="Fraction of signal");
        axis(1, at=seq(0, 100, by=10)+1, 
             labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2);
        abline(v=seq(0, 100, by=10)+1, lty=2, col="gray");
        dev.off();

        cat("# ", set, " [ ",TR.set," ] : DONE...\n", sep="");
        
   }; # for TR.set

} # runatacseqQC
