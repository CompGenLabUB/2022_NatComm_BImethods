#
# prm_vs_enh_count_plots.R
#
#   Functions to plot different filterings
#   made on promoter/enhancer foxG ratios.
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

prm_vs_enh_plots <- function(TT,LBL,IGDR,mxx=0,mxy=0) {
    
  TT$PRM.foxRatio <- ifelse(TT$PRM==0, 0, TT$PRM.FOXG/TT$PRM);
  TT$ENH.foxRatio <- ifelse(TT$ENH==0, 0, TT$ENH.FOXG/TT$ENH);

  t.mxx <- max(mxx,TT$PRM.foxRatio);
  t.mxy <- max(mxy,TT$ENH.foxRatio);
  
  ggplot(TT,
         aes(x=PRM.foxRatio,
             y=ENH.foxRatio,
             color=as.factor(ENH.FOXG),
             size=as.factor(PRM.FOXG))) + 
    geom_point() +
    geom_point(data=TT[ TT$SMESG == "SMESG000024191.1",  ],
               aes(x=PRM.foxRatio,
                   y=ENH.foxRatio),
               color="red", fill="red", shape=1, size=2, stroke = 2) +
    theme_bw() + theme(legend.title.align=0.5) +
    scale_x_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05)),limits=c(0, t.mxx)) +
    scale_y_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05)),limits=c(0, t.mxy)) +
    scale_color_viridis_d() +
    labs(x="Ratio #FoxG vs #all motifs in promoter x gene",
         y="Ratio #FoxG vs #all motifs in enhancer x gene",
         color="# Enhancers\nx FOXG",
         size="# Promoters\nx FOXG") +
    geom_hline(aes(yintercept=0.015),color="red",linetype="dashed") + 
    geom_vline(aes(xintercept=0.015),color="red",linetype="dashed");

  ggsave(file=paste0(IGDR,"_",LBL,".png"), width=7, height=6, dpi=300);

} # prm_vs_enh_plots

prm_vs_enh_count_plots <- function(TT,LBL,IGDR) {
    
  TT$PRM.foxRatio <- ifelse(TT$PRM==0, 0, TT$PRM.FOXG/TT$PRM);
  TT$ENH.foxRatio <- ifelse(TT$ENH==0, 0, TT$ENH.FOXG/TT$ENH);

  # saving the whole data frame first
  write.table(TT,file=paste0(IGDR,"_",LBL,".wholeset.tbl"),
              sep="\t",col.names=TRUE,quote=FALSE);

  ggplot(TT,
         aes(x=PRM.foxRatio,
             y=ENH.foxRatio,
             color=as.factor(ENH.FOXG),
             size=as.factor(PRM.FOXG))) + 
    geom_point() +
    geom_point(data=TT[ TT$SMESG == "SMESG000024191.1",  ],
               aes(x=PRM.foxRatio,
                   y=ENH.foxRatio),
               color="red", fill="red", shape=1, size=2, stroke = 2) +
    theme_bw() + theme(legend.title.align=0.5) +
    facet_wrap( ~ as.factor(DGESTR)) +
    scale_x_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05))) +
    scale_y_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05))) +
    scale_color_viridis_d() +
    labs(x="Ratio #FoxG vs #all motifs in promoter x gene",
         y="Ratio #FoxG vs #all motifs in enhancer x gene",
         color="# Enhancers\nx FOXG",
         size="# Promoters\nx FOXG") +
    geom_hline(aes(yintercept=0.015),color="red",linetype="dashed") + 
    geom_vline(aes(xintercept=0.015),color="red",linetype="dashed");
  ggsave(file=paste0(IGDR,"_",LBL,".png"),
         width=12, height=6, dpi=300);
  write.table(TT,file=paste0(IGDR,"_",LBL,".tbl"),
              sep="\t",col.names=TRUE,quote=FALSE);

  tt <- TT[ TT$PRM.foxRatio > 0.015 & TT$ENH.foxRatio > 0.015, ]
  ggplot(tt,
         aes(x=PRM.foxRatio,
             y=ENH.foxRatio,
             color=as.factor(ENH.FOXG),
             size=as.factor(PRM.FOXG))) + 
    geom_point() +
    geom_point(data=tt[ tt$SMESG == "SMESG000024191.1",  ],
               aes(x=PRM.foxRatio,
                   y=ENH.foxRatio),
               color="red", fill="red", shape=1, size=2, stroke = 2) +
    theme_bw() + theme(legend.title.align=0.5) +
    facet_wrap( ~ as.factor(DGESTR)) +
    scale_x_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05))) +
    scale_y_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05))) +
    scale_color_viridis_d() +
    labs(x="Ratio #FoxG vs #all motifs in promoter x gene",
         y="Ratio #FoxG vs #all motifs in enhancer x gene",
         color="# Enhancers\nx FOXG",
         size="# Promoters\nx FOXG") +
    geom_hline(aes(yintercept=0.015),color="red",linetype="dashed") + 
    geom_vline(aes(xintercept=0.015),color="red",linetype="dashed");
  ggsave(file=paste0(IGDR,"_",LBL,".Pfxr-and-Efxr_gt0015.png"),
         width=12, height=6, dpi=300);
  write.table(tt,file=paste0(IGDR,"_",LBL,".Pfxr-and-Efxr_gt0015.tbl"),
              sep="\t",col.names=TRUE,quote=FALSE);
       
  tt <- TT[ TT$PRM.foxRatio > 0.015 | TT$ENH.foxRatio > 0.015, ]
  ggplot(tt,
         aes(x=PRM.foxRatio,
             y=ENH.foxRatio,
             color=as.factor(ENH.FOXG),
             size=as.factor(PRM.FOXG))) + 
    geom_point() +
    geom_point(data=tt[ tt$SMESG == "SMESG000024191.1",  ],
               aes(x=PRM.foxRatio,
                   y=ENH.foxRatio),
               color="red", fill="red", shape=1, size=2, stroke = 2) +
    theme_bw() + theme(legend.title.align=0.5) +
    facet_wrap( ~ as.factor(DGESTR)) +
    scale_x_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05))) +
    scale_y_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05))) +
    scale_color_viridis_d() +
    labs(x="Ratio #FoxG vs #all motifs in promoter x gene",
         y="Ratio #FoxG vs #all motifs in enhancer x gene",
         color="# Enhancers\nx FOXG",
         size="# Promoters\nx FOXG") +
    geom_hline(aes(yintercept=0.015),color="red",linetype="dashed") + 
    geom_vline(aes(xintercept=0.015),color="red",linetype="dashed");
  ggsave(file=paste0(IGDR,"_",LBL,".Pfxr-or-Efxr_gt0015.png"),
         width=12, height=6, dpi=300);
  write.table(tt,file=paste0(IGDR,"_",LBL,".Pfxr-or-Efxr_gt0015.tbl"),
              sep="\t",col.names=TRUE,quote=FALSE);
       
  tt <- TT[ xor( TT$PRM.foxRatio > 0.015 , TT$ENH.foxRatio > 0.015 ), ]
  ggplot(tt,
         aes(x=PRM.foxRatio,
             y=ENH.foxRatio,
             color=as.factor(ENH.FOXG),
             size=as.factor(PRM.FOXG))) + 
    geom_point() +
    geom_point(data=tt[ tt$SMESG == "SMESG000024191.1",  ],
               aes(x=PRM.foxRatio,
                   y=ENH.foxRatio),
               color="red", fill="red", shape=1, size=2, stroke = 2) +
    theme_bw() + theme(legend.title.align=0.5) +
    facet_wrap( ~ as.factor(DGESTR)) +
    scale_x_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05))) +
    scale_y_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05))) +
    scale_color_viridis_d() +
    labs(x="Ratio #FoxG vs #all motifs in promoter x gene",
         y="Ratio #FoxG vs #all motifs in enhancer x gene",
         color="# Enhancers\nx FOXG",
         size="# Promoters\nx FOXG") +
    geom_hline(aes(yintercept=0.015),color="red",linetype="dashed") + 
    geom_vline(aes(xintercept=0.015),color="red",linetype="dashed");
  ggsave(file=paste0(IGDR,"_",LBL,".Pfxr-xor-Efxr_gt0015.png"),
         width=12, height=6, dpi=300);
  write.table(tt,file=paste0(IGDR,"_",LBL,".Pfxr-xor-Efxr_gt0015.tbl"),
              sep="\t",col.names=TRUE,quote=FALSE);
       
  l <- function(x) 0.015 - x;
  f <- function(x) sin(pi/2) * 0.000075/x;
  # ggplot() +  xlim(-5, 5) + ylim(-5, 5) + geom_function(fun=f)
  t.mxx <- max(tt$PRM.foxRatio);
  t.mxy <- max(tt$ENH.foxRatio);
  tt <- TT[l(TT$PRM.foxRatio) < TT$ENH.foxRatio, ];
  ggplot(tt,
         aes(x=PRM.foxRatio,
             y=ENH.foxRatio,
             color=as.factor(ENH.FOXG),
             size=as.factor(PRM.FOXG))) + 
    geom_point() +
    geom_point(data=tt[ tt$SMESG == "SMESG000024191.1",  ],
               aes(x=PRM.foxRatio,
                   y=ENH.foxRatio),
               color="red", fill="red", shape=1, size=2, stroke = 2) +
    theme_bw() + theme(legend.title.align=0.5) +
    facet_wrap( ~ as.factor(DGESTR)) +
    scale_x_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05)),limits=c(0, t.mxx)) +
    scale_y_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05)),limits=c(0, t.mxy)) +
    scale_color_viridis_d() +
    labs(x="Ratio #FoxG vs #all motifs in promoter x gene",
         y="Ratio #FoxG vs #all motifs in enhancer x gene",
         color="# Enhancers\nx FOXG",
         size="# Promoters\nx FOXG",
         title="l(PRM FoxG ratio) < ENH FoxG ratio") +
    geom_hline(aes(yintercept=0.015),color="red",linetype="dashed") + 
    geom_vline(aes(xintercept=0.015),color="red",linetype="dashed") +
    # geom_abline(aes(intercept=sqrt(0.015), slope = -1), color="red",linetype="dotted") +
    stat_function(fun = l, color="red", linetype="dotted",
                  n=250, size=0.5, na.rm=TRUE, inherit.aes=FALSE) +
    stat_function(fun = f, color="blue", linetype="dotted",
                  n=250, size=0.5, na.rm=TRUE, inherit.aes=FALSE);
  ggsave(file=paste0(IGDR,"_",LBL,".lPfxr_lt_Efxr.png"),
       width=12, height=6, dpi=300);
  write.table(tt,file=paste0(IGDR,"_",LBL,".lPfxr_lt_Efxr.tbl"),
              sep="\t",col.names=TRUE,quote=FALSE);

  tt <- TT[f(TT$PRM.foxRatio) < TT$ENH.foxRatio, ];
  ggplot(tt,
         aes(x=PRM.foxRatio,
             y=ENH.foxRatio,
             color=as.factor(ENH.FOXG),
             size=as.factor(PRM.FOXG))) + 
    geom_point() +
    geom_point(data=tt[ tt$SMESG == "SMESG000024191.1",  ],
               aes(x=PRM.foxRatio,
                   y=ENH.foxRatio),
               color="red", fill="red", shape=1, size=2, stroke = 2) +
    theme_bw() + theme(legend.title.align=0.5) +
    facet_wrap( ~ as.factor(DGESTR)) +
    scale_x_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05)),limits=c(0, t.mxx)) +
    scale_y_sqrt(breaks=c(0,0.01,0.025,seq(0.05,1,0.05)),limits=c(0, t.mxy)) +
    scale_color_viridis_d() +
    labs(x="Ratio #FoxG vs #all motifs in promoter x gene",
         y="Ratio #FoxG vs #all motifs in enhancer x gene",
         color="# Enhancers\nx FOXG",
         size="# Promoters\nx FOXG",
         title="f(PRM FoxG ratio) < ENH FoxG ratio") +
    geom_hline(aes(yintercept=0.015),color="red",linetype="dashed") + 
    geom_vline(aes(xintercept=0.015),color="red",linetype="dashed") +
    stat_function(fun = l, color="red", linetype="dotted",
                  n=250, size=0.5, na.rm=TRUE, inherit.aes=FALSE) +
    stat_function(fun = f, color="blue", linetype="dotted",
                  n=250, size=0.5, na.rm=TRUE, inherit.aes=FALSE);
  ggsave(file=paste0(IGDR,"_",LBL,".fPfxr_lt_Efxr.png"),
         width=12, height=6, dpi=300);
  write.table(tt,file=paste0(IGDR,"_",LBL,".fPfxr_lt_Efxr.tbl"),
              sep="\t",col.names=TRUE,quote=FALSE);
} # prm_vs_enh_count_plots
