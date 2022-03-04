## IMPROVED CORRELATION Plots

corrplotfunct <- function(CE, cepx, thy.cols, VX.lim=NULL, VY.lim=NULL,
                          VX.logflg=NULL, VY.logflg=NULL, XYratio=TRUE, binflg=FALSE,
                          cep.title=NULL, cep.subtitle=NULL,
                          thy.digits=6, thy.cex.cor=0.75) {

  CEP <- list();
  KK <- 1;
  #
  cepy <- cepx;
  #
  max.n <- length(cepx);
  max.m <- length(cepy);

  for (n in 1:max.n) {
    
    N <- cepx[[n]];
    # CEP[[N]] <- list();
    if (missing(cep.title) || is.null(cep.title)) {
      J.t <- cep.Title <- N;
    } else {
      if (is.list(cep.title)) {
        J.t <- cep.title[[N]];
      } else {
        if (is.vector(cep.title) && length(cep.title)>1) {
          J.t <- cep.title[[n]];
        } else {
          J.t <- cep.title;
        };
      };
    };
    #    
    if (missing(cep.subtitle) || is.null(cep.subtitle)) {
      J.s <- cep.SubTitle <- "";
    } else {
      if (is.list(cep.subtitle)) {
        J.s <- cep.subtitle[[N]];
      } else {
        if (is.vector(cep.subtitle) && length(cep.subtitle)>1) {
          J.s <- cep.subtitle[[n]];
        } else {
          J.s <- cep.SubTitle;
        };
      };
    };

    for (m in 1:max.m) {

      M <- cepy[[m]];
      #
      if (n == m) { # plot varname

          vx <- CE[ , N];
          ##
          if (missing(VX.lim) || is.null(VX.lim)) {
            vx.brk <- pretty(c(min(vx, na.rm=TRUE),
                               max(vx, na.rm=TRUE)),
                             n=6);
            vx.lim <- c(vx.brk[1], vx.brk[length(vx.brk)]);
          } else {
            if (is.list(VX.lim)) {
              vx.brk <- pretty(c(VX.lim[[N]][1],
                                 VX.lim[[N]][2]),
                               n=6);
              vx.lim <- VX.lim;
            } else {
              vx.brk <- pretty(c(VX.lim[1],
                                 VX.lim[2]),
                               n=6);
              vx.lim <- VX.lim;
            };
          };
          ##
          if (missing(thy.cols) || is.null(thy.cols)) {
             vx.col <- "blue";
          } else {
             vx.col <- thy.cols[[n]];
          }
          #
          if (missing(VX.logflg) || is.null(VX.logflg)) {
            vx.logflg <- FALSE;
          } else {
            if (is.list(VX.logflg)) {
              vx.logflg <- VX.logflg[[N]];
            } else {
              vx.logflg <- VX.logflg;
            }
          }
          #
          cep <- ggplot(data=data.frame(YVIOL = vx,
                                        XVIOL = as.factor(rep(1,length(vx)))),
                        aes(x=XVIOL, y=YVIOL)
                        ) +
                   geom_violin(fill=vx.col) + coord_flip() +
                   theme(panel.background = element_rect(fill = "#00000000", colour = "black"),
                         panel.grid.major = element_line(colour = "grey", linetype = "dotted", size=0.5),
                         panel.grid.minor = element_line(colour = "grey", linetype = "dotted", size=0.5),
                         plot.margin      = unit(x = c(20,  2,     ifelse(m < max.m, 2, 20), ifelse(n > 1, 2, 20)), units = "mm"),
                                                     # top, right, bottom,                   and left
                         plot.title       = element_text(size=35, vjust=2, face="bold"), #element_blank(),
                         axis.ticks.y     = element_blank(),
                         panel.grid.major.y = element_blank(),
                         panel.grid.minor.y = element_blank(),
                         axis.title.y = element_blank(),
                         axis.title.x = element_blank(),
                         axis.text.y  = element_blank(),
                         axis.text.x  = element_text(vjust=1, hjust=0, angle=-35, color="black", face="bold")
                         # axis.text.x  = element_text(color="black")
                         ) +
                   scale_x_discrete(limits=as.factor(c(0,1,2))) + 
                   scale_y_continuous(limits=vx.lim,
                                      breaks=vx.brk,
                                      labels=if(vx.logflg) { sprintf("10e+%02d", vx.brk) } else { vx.brk }) +
                   ggtitle(J.t) +
                   annotate("text", label = J.s,
                            x = as.factor(2), y = vx.lim[2],
                            vjust = 0, hjust = 1, size = 5, colour = "black")

                   # ggtitle(N)

      } else {

        if (n > m) { # write correlations on upper diagonal

          x <- CE[ , N];
          y <- CE[ , M];
          #
          # corerlation parameters
          # k.pe   <- cor(x=x, y=y, method="pearson");  # parametric correlation
          k.pe.tst <- cor.test(x=x, y=y, method="pearson");
          # k.sp   <- cor(x=x, y=y, method="spearman"); # non-parametric correlation
          k.sp.tst <- cor.test(x=x, y=y, method="spearman", alternative="two.sided");
          #
          k.pe.txt  <- format(k.pe.tst$estimate, digits=thy.digits);
          k.pe.ptxt <- format(k.pe.tst$p.value,  digits=thy.digits);
          k.sp.txt  <- format(k.sp.tst$estimate, digits=thy.digits);
          k.sp.ptxt <- format(k.sp.tst$p.value,  digits=thy.digits);
          #
          k.pe.Signif <- symnum(k.pe.tst$p.value, corr = FALSE, na = FALSE,
                                cutpoints = c(    0, 0.001,  0.01,  0.05,  0.1,   1),
                                  symbols = c("***",  "**",   "*",   ".",  " "));
          k.sp.Signif <- symnum(k.sp.tst$p.value, corr = FALSE, na = FALSE,
                                cutpoints = c(    0, 0.001,  0.01,  0.05,  0.1,   1),
                                  symbols = c("***",  "**",   "*",   ".",  " "));
          #
          # linear model parameters
          k.lm      <- lm(y ~ x);
          # k.lm.coefs     <- coef(k.lm)
          # k.lm.intercept <- k.lm.coefs[1]; # k.lm$coefficients[1];
          # k.lm.slope     <- k.lm.coefs[2]; # k.lm$coefficients[2];
          # k.lm.df        <- k.lm$df.residual;
          # k.lm.se        <- sqrt( sum( residuals(k.lm)^2 ) / k.lm.df );
          # k.lm.confint   <- confint(k.lm, level=0.99)
          k.lm.data  <- summary(k.lm);
          k.lm.coefs <- k.lm.data$coefficients;
          k.lm.intercept.Signif <- symnum(k.lm.coefs[7], corr = FALSE, na = FALSE,
                                          cutpoints = c(    0, 0.001,  0.01,  0.05,  0.1,   1),
                                            symbols = c("***",  "**",   "*",   ".",  " "));
          k.lm.slope.Signif     <- symnum(k.lm.coefs[8], corr = FALSE, na = FALSE,
                                          cutpoints = c(    0, 0.001,  0.01,  0.05,  0.1,   1),
                                            symbols = c("***",  "**",   "*",   ".",  " "));
          #
          cep <- ggplot(data=data.frame()) +
                       geom_point() +
                       theme(panel.background = element_rect(fill = "#00000000", colour = "black"),
                             plot.margin  = unit(x = c(2, 2, 2, 2), units = "mm"), # top, right, bottom, and left
                             plot.title   = element_blank(),
                             axis.ticks   = element_blank(),
                             panel.grid   = element_blank(),
                             axis.title.x = element_blank(),
                             axis.title.y = element_blank(),
                             axis.text.x  = element_blank(),
                             axis.text.y  = element_blank()
                             ) +
                       xlim(0,2) + ylim(0, 1.2) +
                       # scale_x_discrete(name=NULL, limits=c(0, 1), breaks=NULL) +
                       # scale_y_discrete(name=NULL, limits=c(0, 1), breaks=NULL) +
                       annotate("text", label = paste(M, "x", N, sep=" "), x = 0.05, y = 1.1,
                                hjust = 0, size = 6, colour = "black", face="bold") +
                       annotate("text", label = "Pearson's product-moment correlation",                       x = 0.05, y = 0.925,
                                hjust = 0, size = 4, colour = "black", face="bold") +
                       annotate("text", label = paste("corr = ", k.pe.txt, sep=""),                           x = 0.15, y = 0.850,
                                hjust = 0, size = 3.5, colour = "black") +
                       annotate("text", label = paste("p-val = ", k.pe.ptxt, " ", k.pe.Signif, sep=""), x = 0.15, y = 0.775,
                                hjust = 0, size = 3.5, colour = "black") +
                       annotate("text", label = "Spearman's rank correlation rho (two sided)",                x = 0.05, y = 0.675,
                                hjust = 0, size = 4, colour = "black", face="bold") +
                       annotate("text", label = paste("Rho = ", k.sp.txt, sep=""),                            x = 0.15, y = 0.600,
                                hjust = 0, size = 3.5, colour = "black") +
                       annotate("text", label = paste("p-val = ", k.sp.ptxt, " ", k.sp.Signif, sep=""), x = 0.15, y = 0.525,
                                hjust = 0, size = 3.5, colour = "black") +
                       annotate("text", label = "Linear model estimates (lm)",                                x = 0.05, y = 0.425,
                                hjust = 0, size = 4, colour = "black", face="bold") +
                       annotate("text", label = sprintf("%12s  %12s  %12s  %12s  %12s",
                                                        "Coefficients", "Estimate", "Std.Error", "t value", "Pr(>|t|)"),
                                x = 0.15, y = 0.350, hjust = 0, size = 2.5, colour = "black", face="bold") +
                       annotate("text", label = sprintf("%12s  %12.4g  %12.4g  %12.4g  %12.4g %s",
                                                        "Intercept", k.lm.coefs[1], k.lm.coefs[3],
                                                                     k.lm.coefs[5], k.lm.coefs[7], k.lm.intercept.Signif),
                                x = 0.15, y = 0.300, hjust = 0, size = 2.5, colour = "black", face="bold") +
                       annotate("text", label = sprintf("%12s  %12.4g  %12.4g  %12.4g  %12.4g %s",
                                                        "Slope", k.lm.coefs[2], k.lm.coefs[4],
                                                                 k.lm.coefs[6], k.lm.coefs[8], k.lm.slope.Signif),
                                x = 0.15, y = 0.250, hjust = 0, size = 2.5, colour = "black", face="bold") +
                       annotate("text", label = sprintf("Residual standard error: %f on %d degrees of freedom",
                                                        k.lm.data$sigma, k.lm.data$df[2]),
                                x = 0.15, y = 0.200, hjust = 0, size = 2.5, colour = "black", face="bold") +
                       annotate("text", label = sprintf("Multiple R-squared: %7.5f  /  Adjusted R-squared: %7.5f",
                                                        k.lm.data$r.squared, k.lm.data$adj.r.squared),
                                x = 0.15, y = 0.150, hjust = 0, size = 2.5, colour = "black", face="bold") +
                       annotate("text", label = sprintf("F-statistic: %.2f on %d and %d DF,  p-value: %12.4g",
                                                        k.lm.data$fstatistic[1], k.lm.data$fstatistic[2],
                                                        k.lm.data$fstatistic[3], k.lm.coefs[8]),
                                x = 0.15, y = 0.100, hjust = 0, size = 2.5, colour = "black", face="bold")
 
        } else {  # n > m  # plot correlations on lower diagonal

          # vx <- CE[ , N];
          # vy <- CE[ , M];
          vvxy <- CE[ , c(N, M)];
          ##
          if (missing(VX.lim) || is.null(VX.lim)) {
            vx.brk <- pretty(c(min(vx, na.rm=TRUE),
                               max(vx, na.rm=TRUE)),
                             n=6);
            vx.lim <- c(vx.brk[1], vx.brk[length(vx.brk)]);
          } else {
            if (is.list(VX.lim)) {
              vx.brk <- pretty(c(VX.lim[[N]][1],
                                 VX.lim[[N]][2]),
                               n=6);
              vx.lim <- VX.lim[[N]];
            } else {
              vx.brk <- pretty(c(VX.lim[1],
                                 VX.lim[2]),
                               n=6);
              vx.lim <- VX.lim;
            };
          };
          ##
          if (missing(VX.logflg) || is.null(VX.logflg)) {
            vx.logflg <- FALSE;
          } else {
            if (is.list(VX.logflg)) {
              vx.logflg <- VX.logflg[[N]];
            } else {
              vx.logflg <- VX.logflg;
            }
          };
          #
          if (missing(VY.logflg) || is.null(VY.logflg)) {
            vy.logflg <- vx.logflg;
          } else {
            if (is.list(VY.logflg)) {
              vy.logflg <- VY.logflg[[N]];
            } else {
              vy.logflg <- VY.logflg;
            }
          };
          ##
          if (missing(VY.lim) || is.null(VY.lim)) {
            if (XYratio) {
              vy.brk <- vx.brk;
              vy.lim <- vx.lim;
              VY.logflg <- VX.logflg;
              vy.logflg <- vx.logflg;
            } else {
              vy.brk <- pretty(c(min(vy, na.rm=TRUE),
                                 max(vy, na.rm=TRUE)),
                               n=6);
              vy.lim <- c(vy.brk[1], vy.brk[length(vy.brk)]);
            };
          } else {
            if (is.list(VY.lim)) {
              vy.brk <- pretty(c(VY.lim[[M]][1],
                                 VY.lim[[M]][2]),
                               n=6);
              vy.lim <- VY.lim[[M]];
            } else {
              vy.brk <- pretty(c(VY.lim[1],
                                 VY.lim[2]),
                               n=6);
              vy.lim <- VY.lim;
            };
          };
          ##
          if (missing(thy.cols) || is.null(thy.cols)) {
             vx.col <- "blue";
             vy.col <- "blue";
          } else {
             vx.col <- thy.cols[[n]];
             vy.col <- thy.cols[[m]];
          }
          #
          # cep <- ggplot(data = data.frame(VXVAR = vx,
          #                                 VYVAR = vy),
          cepi <- ggplot(data = data.frame(VXVAR = vvxy[ , N], # vx,
                                           VYVAR = vvxy[ , M]),# vy),
                         aes(x = VXVAR,
                             y = VYVAR)
                         )
          if (binflg) {
            cpeo <- cepi +
                      stat_density2d(geom = "tile", contour = F, aes(fill = ..density..)) +
                      scale_fill_gradient2(low="darkblue",
                                           high="darkred",
                                           mid="white",
                                           midpoint=0.5)
                       # + stat_binhex(bins = if(vx.logflg) { 0.1 } else { 5 }, na.rm = TRUE)
          } else {
            cpeo <- cepi +
                      scale_fill_manual(values=vx.col,
                                        guide="none") +
                      scale_colour_manual(values=vy.col,
                                          guide="none")
                   
          }
          cep <- cpeo +
                   theme(panel.background = element_rect(fill = "#00000000"),
                         plot.margin  = unit(x = c(2, 2, ifelse(m < max.m, 2, 10), ifelse(n > 1, 2, 10)), units = "mm"), # top, right, bottom, and left
                         plot.title   = element_blank(),
                         panel.grid.major = element_line(colour = "grey", linetype = "dotted", size=0.5),
                         panel.grid.minor = element_line(colour = "grey", linetype = "dotted", size=0.5),
                         axis.title.x = element_blank(), # element_text(size=30, vjust=-1.2, face="bold"),
                         axis.title.y = element_blank(), # element_text(size=30, vjust=-0.1, face="bold", angle = 90),
                         strip.text.x = if (m < max.m) { element_blank() } else { element_text(size=20, face="bold") },
                         strip.text.y = if (n > 1    ) { element_blank() } else { element_text(size=20, face="bold", angle = -90) },
                         axis.text.x  = if (m < max.m) { element_blank() } else { element_text(vjust=1, hjust=0, angle=-35, color="black", face="bold") },
                         axis.text.y  = if (n > 1    ) { element_blank() } else { element_text(color="black") }
                         ) +
                   geom_point(colour=vx.col, size=2, alpha=0.9, na.rm=TRUE) +
                   scale_x_continuous(limits=vx.lim,
                                      breaks=vx.brk,
                                      labels=if(vx.logflg) { sprintf("10e+%02d", vx.brk) } else { vx.brk }) +
                   scale_y_continuous(limits=vy.lim,
                                      breaks=vy.brk,
                                      labels=if(vy.logflg) { sprintf("10e+%02d", vy.brk) } else { vy.brk }) +
                   geom_smooth(method="lm",
                               alpha=1,
                               size=0.5,
                               se=FALSE,
                               color="black",
                               linetype="dotted",
                               fullrange=TRUE);

        }

      }

      CEP[[KK]] <- cep;
      KK <- as.character(as.numeric(KK) + 1);

    }; # for m

  }; # for n

  CEP

} # corrplotfunct

## CORRELATION PLOT With color factor


corrplotfunctCOL <- function(CE, cepx, thy.cols, thy.factor, factor.cols=NULL, VX.lim=NULL, VY.lim=NULL,
                             VX.logflg=NULL, VY.logflg=NULL, XYratio=TRUE, binflg=FALSE,
                             cep.title=NULL, cep.subtitle=NULL,
                             thy.digits=6, thy.cex.cor=0.75, thy.alpha=0.9) {
  
  CEP <- list();
  KK <- 1;
  #
  cepy <- cepx;
  #
  max.n <- length(cepx);
  max.m <- length(cepy);

  if (missing(factor.cols) || is.null(factor.cols)) {
    cols.factor <- "blue";
  } else {
    cols.factor <- c(factor.cols);
  };
  
  for (n in 1:max.n) {
    
    N <- cepx[[n]];
    # CEP[[N]] <- list();
    if (missing(cep.title) || is.null(cep.title)) {
      J.t <- cep.Title <- N;
    } else {
      if (is.list(cep.title)) {
        J.t <- cep.title[[N]];
      } else {
        if (is.vector(cep.title) && length(cep.title)>1) {
          J.t <- cep.title[[n]];
        } else {
          J.t <- cep.title;
        };
      };
    };
    #    
    if (missing(cep.subtitle) || is.null(cep.subtitle)) {
      J.s <- cep.SubTitle <- "";
    } else {
      if (is.list(cep.subtitle)) {
        J.s <- cep.subtitle[[N]];
      } else {
        if (is.vector(cep.subtitle) && length(cep.subtitle)>1) {
          J.s <- cep.subtitle[[n]];
        } else {
          J.s <- cep.SubTitle;
        };
      };
    };

    for (m in 1:max.m) {

      M <- cepy[[m]];
      #
      if (n == m) { # plot varname

          vx <- CE[ , N];
          ##
          if (missing(VX.lim) || is.null(VX.lim)) {
            vx.brk <- pretty(c(min(vx, na.rm=TRUE),
                               max(vx, na.rm=TRUE)),
                             n=6);
            vx.lim <- c(vx.brk[1], vx.brk[length(vx.brk)]);
          } else {
            if (is.list(VX.lim)) {
              vx.brk <- pretty(c(VX.lim[[N]][1],
                                 VX.lim[[N]][2]),
                               n=6);
              vx.lim <- VX.lim;
            } else {
              vx.brk <- pretty(c(VX.lim[1],
                                 VX.lim[2]),
                               n=6);
              vx.lim <- VX.lim;
            };
          };
          #
          if (missing(VX.logflg) || is.null(VX.logflg)) {
            vx.logflg <- FALSE;
          } else {
            if (is.list(VX.logflg)) {
              vx.logflg <- VX.logflg[[N]];
            } else {
              vx.logflg <- VX.logflg;
            }
          }
          ##
          if (missing(thy.cols) || is.null(thy.cols)) {
             vx.col <- "blue";
          } else {
            if (is.list(thy.cols)) {
               vx.col <- thy.cols[[N]];
             } else {
               if (is.vector(thy.cols) && length(thy.cols)>1) {
                 vx.col <- thy.cols[[n]];
               } else {
                 vx.col <- thy.cols;
               };
             };
          }
          #
          cep <- ggplot(data=data.frame(YVIOL = vx,
                                        XVIOL = as.factor(rep(1,length(vx)))),
                        aes(x=XVIOL, y=YVIOL)
                        ) +
                   geom_violin(fill=vx.col) + coord_flip() +
                   theme(panel.background = element_rect(fill = "#00000000", colour = "black"),
                         panel.grid.major = element_line(colour = "grey", linetype = "dotted", size=0.5),
                         panel.grid.minor = element_line(colour = "grey", linetype = "dotted", size=0.5),
                         plot.margin      = unit(x = c(20,  2,     ifelse(m < max.m, 2, 20), ifelse(n > 1, 2, 20)), units = "mm"),
                                                     # top, right, bottom,                   and left
                         plot.title       = element_text(size=15, vjust=2, face="bold"), #element_blank(),
                         axis.ticks.y     = element_blank(),
                         panel.grid.major.y = element_blank(),
                         panel.grid.minor.y = element_blank(),
                         axis.title.y = element_blank(),
                         axis.title.x = element_blank(),
                         axis.text.y  = element_blank(),
                         axis.text.x  = element_text(vjust=1, hjust=0, angle=-35, color="black", face="bold")
                         # axis.text.x  = element_text(color="black")
                         ) +
                   scale_x_discrete(limits=as.factor(c(0,1,2))) + 
                   scale_y_continuous(limits=vx.lim,
                                      breaks=vx.brk,
                                      labels=if(vx.logflg) { sprintf("10e+%02d", vx.brk) } else { vx.brk }) +
                   ggtitle(J.t) +
                   annotate("text", label = J.s,
                            x = as.factor(2), y = vx.lim[2],
                            vjust = 0, hjust = 1, size = 5, colour = "black")

                   # ggtitle(N)

      } else {

        if (n > m) { # write correlations on upper diagonal

          x <- CE[ , N];
          y <- CE[ , M];
          #
          # corerlation parameters
          # k.pe   <- cor(x=x, y=y, method="pearson");  # parametric correlation
          k.pe.tst <- cor.test(x=x, y=y, method="pearson");
          # k.sp   <- cor(x=x, y=y, method="spearman"); # non-parametric correlation
          k.sp.tst <- cor.test(x=x, y=y, method="spearman", alternative="two.sided");
          #
          k.pe.txt  <- format(k.pe.tst$estimate, digits=thy.digits);
          k.pe.ptxt <- format(k.pe.tst$p.value,  digits=thy.digits);
          k.sp.txt  <- format(k.sp.tst$estimate, digits=thy.digits);
          k.sp.ptxt <- format(k.sp.tst$p.value,  digits=thy.digits);
          #
          k.pe.Signif <- symnum(k.pe.tst$p.value, corr = FALSE, na = FALSE,
                                cutpoints = c(    0, 0.001,  0.01,  0.05,  0.1,   1),
                                  symbols = c("***",  "**",   "*",   ".",  " "));
          k.sp.Signif <- symnum(k.sp.tst$p.value, corr = FALSE, na = FALSE,
                                cutpoints = c(    0, 0.001,  0.01,  0.05,  0.1,   1),
                                  symbols = c("***",  "**",   "*",   ".",  " "));
          #
          # linear model parameters
          k.lm      <- lm(y ~ x);
          # k.lm.coefs     <- coef(k.lm)
          # k.lm.intercept <- k.lm.coefs[1]; # k.lm$coefficients[1];
          # k.lm.slope     <- k.lm.coefs[2]; # k.lm$coefficients[2];
          # k.lm.df        <- k.lm$df.residual;
          # k.lm.se        <- sqrt( sum( residuals(k.lm)^2 ) / k.lm.df );
          # k.lm.confint   <- confint(k.lm, level=0.99)
          k.lm.data  <- summary(k.lm);
          k.lm.coefs <- k.lm.data$coefficients;
          k.lm.intercept.Signif <- symnum(k.lm.coefs[7], corr = FALSE, na = FALSE,
                                          cutpoints = c(    0, 0.001,  0.01,  0.05,  0.1,   1),
                                            symbols = c("***",  "**",   "*",   ".",  " "));
          k.lm.slope.Signif     <- symnum(k.lm.coefs[8], corr = FALSE, na = FALSE,
                                          cutpoints = c(    0, 0.001,  0.01,  0.05,  0.1,   1),
                                            symbols = c("***",  "**",   "*",   ".",  " "));
          #
          cep <- ggplot(data=data.frame()) +
                       geom_point() +
                       theme(panel.background = element_rect(fill = "#00000000", colour = "black"),
                             plot.margin  = unit(x = c(2, 2, 2, 2), units = "mm"), # top, right, bottom, and left
                             plot.title   = element_blank(),
                             axis.ticks   = element_blank(),
                             panel.grid   = element_blank(),
                             axis.title.x = element_blank(),
                             axis.title.y = element_blank(),
                             axis.text.x  = element_blank(),
                             axis.text.y  = element_blank()
                             ) +
                       xlim(0,2) + ylim(0, 1.2) +
                       # scale_x_discrete(name=NULL, limits=c(0, 1), breaks=NULL) +
                       # scale_y_discrete(name=NULL, limits=c(0, 1), breaks=NULL) +
                       annotate("text", label = paste(M, "x\n", N, sep=" "), x = 0.05, y = 1.1,
                                hjust = 0, size = 5, colour = "black", face="bold") +
                       annotate("text", label = "Pearson's product-moment correlation",                       x = 0.05, y = 0.925,
                                hjust = 0, size = 4, colour = "black", face="bold") +
                       annotate("text", label = paste("corr = ", k.pe.txt, sep=""),                           x = 0.15, y = 0.850,
                                hjust = 0, size = 3.5, colour = "black") +
                       annotate("text", label = paste("p-val = ", k.pe.ptxt, " ", k.pe.Signif, sep=""), x = 0.15, y = 0.775,
                                hjust = 0, size = 3.5, colour = "black") +
                       annotate("text", label = "Spearman's rank correlation rho (two sided)",                x = 0.05, y = 0.675,
                                hjust = 0, size = 4, colour = "black", face="bold") +
                       annotate("text", label = paste("Rho = ", k.sp.txt, sep=""),                            x = 0.15, y = 0.600,
                                hjust = 0, size = 3.5, colour = "black") +
                       annotate("text", label = paste("p-val = ", k.sp.ptxt, " ", k.sp.Signif, sep=""), x = 0.15, y = 0.525,
                                hjust = 0, size = 3.5, colour = "black") +
                       annotate("text", label = "Linear model estimates (lm)",                                x = 0.05, y = 0.425,
                                hjust = 0, size = 4, colour = "black", face="bold") +
                       annotate("text", label = sprintf("%12s  %12s  %12s  %12s  %12s",
                                                        "Coefficients", "Estimate", "Std.Error", "t value", "Pr(>|t|)"),
                                x = 0.15, y = 0.350, hjust = 0, size = 2.5, colour = "black", face="bold") +
                       annotate("text", label = sprintf("%12s  %12.4g  %12.4g  %12.4g  %12.4g %s",
                                                        "Intercept", k.lm.coefs[1], k.lm.coefs[3],
                                                                     k.lm.coefs[5], k.lm.coefs[7], k.lm.intercept.Signif),
                                x = 0.15, y = 0.300, hjust = 0, size = 2.5, colour = "black", face="bold") +
                       annotate("text", label = sprintf("%12s  %12.4g  %12.4g  %12.4g  %12.4g %s",
                                                        "Slope", k.lm.coefs[2], k.lm.coefs[4],
                                                                 k.lm.coefs[6], k.lm.coefs[8], k.lm.slope.Signif),
                                x = 0.15, y = 0.250, hjust = 0, size = 2.5, colour = "black", face="bold") +
                       annotate("text", label = sprintf("Residual standard error: %f on %d degrees of freedom",
                                                        k.lm.data$sigma, k.lm.data$df[2]),
                                x = 0.15, y = 0.200, hjust = 0, size = 2.5, colour = "black", face="bold") +
                       annotate("text", label = sprintf("Multiple R-squared: %7.5f  /  Adjusted R-squared: %7.5f",
                                                        k.lm.data$r.squared, k.lm.data$adj.r.squared),
                                x = 0.15, y = 0.150, hjust = 0, size = 2.5, colour = "black", face="bold") +
                       annotate("text", label = sprintf("F-statistic: %.2f on %d and %d DF,  p-value: %12.4g",
                                                        k.lm.data$fstatistic[1], k.lm.data$fstatistic[2],
                                                        k.lm.data$fstatistic[3], k.lm.coefs[8]),
                                x = 0.15, y = 0.100, hjust = 0, size = 2.5, colour = "black", face="bold")
 
        } else {  # n > m  # plot correlations on lower diagonal

          # vx <- CE[ , N];
          # vy <- CE[ , M];
          # vg <- CE[ , thy.factor];
          vvxyg <- CE[ , c(N, M, thy.factor)];
          ##
          if (missing(VX.logflg) || is.null(VX.logflg)) {
            vx.logflg <- FALSE;
          } else {
            if (is.list(VX.logflg)) {
              vx.logflg <- VX.logflg[[N]];
            } else {
              vx.logflg <- VX.logflg;
            }
          };
          #
          if (missing(VY.logflg) || is.null(VY.logflg)) {
            vy.logflg <- vx.logflg;
          } else {
            if (is.list(VY.logflg)) {
              vy.logflg <- VY.logflg[[N]];
            } else {
              vy.logflg <- VY.logflg;
            }
          };
          ##
          if (missing(VX.lim) || is.null(VX.lim)) {
            vx.brk <- pretty(c(min(vx, na.rm=TRUE),
                               max(vx, na.rm=TRUE)),
                             n=6);
            vx.lim <- c(vx.brk[1], vx.brk[length(vx.brk)]);
          } else {
            if (is.list(VX.lim)) {
              vx.brk <- pretty(c(VX.lim[[N]][1],
                                 VX.lim[[N]][2]),
                               n=6);
              vx.lim <- VX.lim[[N]];
            } else {
              vx.brk <- pretty(c(VX.lim[1],
                                 VX.lim[2]),
                               n=6);
              vx.lim <- VX.lim;
            };
          };
          #
          if (missing(VY.lim) || is.null(VY.lim)) {
            if (XYratio) {
              vy.brk <- vx.brk;
              vy.lim <- vx.lim;
              VY.logflg <- VX.logflg;
              vy.logflg <- vx.logflg;
            } else {
              vy.brk <- pretty(c(min(vy, na.rm=TRUE),
                                 max(vy, na.rm=TRUE)),
                               n=6);
              vy.lim <- c(vy.brk[1], vy.brk[length(vy.brk)]);
            };
          } else {
            if (is.list(VY.lim)) {
              vy.brk <- pretty(c(VY.lim[[M]][1],
                                 VY.lim[[M]][2]),
                               n=6);
              vy.lim <- VY.lim[[M]];
            } else {
              vy.brk <- pretty(c(VY.lim[1],
                                 VY.lim[2]),
                               n=6);
              vy.lim <- VY.lim;
            };
          };
          ## cols.factor 
          if (missing(thy.cols) || is.null(thy.cols)) {
             vx.col <- "blue";
             vy.col <- "blue";
          } else {
            if (is.list(thy.cols)) {
              vx.col <- thy.cols[[N]];
              vy.col <- thy.cols[[M]];
            } else {
              if (is.vector(thy.cols) && length(thy.cols)>1) {
                vx.col <- thy.cols[[n]];
                vy.col <- thy.cols[[m]];
              } else {
                vx.col <- vy.col <- thy.cols;
              };
            };
          }
          #
          VVXYG <- data.frame(VXVAR = vvxyg[, N],
                              VYVAR = vvxyg[, M],
                              GPVAR = vvxyg[, thy.factor]);
          cepi <- ggplot() +
                   geom_point(data=VVXYG,
                              aes(x = VXVAR,
                                  y = VYVAR,
                                  group=as.factor(GPVAR),
                                  colour=as.factor(GPVAR),
                                  fill=as.factor(GPVAR)),
                              size=2,
                              alpha=thy.alpha,
                              # colour=cols.factor,
                              na.rm=TRUE)
          if (binflg) {
            cpeo <- cepi +
                      stat_density2d(geom = "tile", contour = F, aes(fill = ..density..)) +
                      scale_fill_gradient2(low="darkblue",
                                           high="darkred",
                                           mid="white",
                                           midpoint=0.5)
                       # + stat_binhex(bins = if(vx.logflg) { 0.1 } else { 5 }, na.rm = TRUE)
          } else {
            cpeo <- cepi +
                      scale_colour_manual(values=cols.factor,
                                          guide="none") +
                      scale_fill_manual(values=cols.factor,
                                        guide="none")
          }
          cep <- cpeo +
                   theme(panel.background = element_rect(fill = "#00000000"),
                         plot.margin  = unit(x = c(2, 2, ifelse(m < max.m, 2, 10), ifelse(n > 1, 2, 10)), units = "mm"), # top, right, bottom, and left
                         plot.title   = element_blank(),
                         panel.grid.major = element_line(colour = "grey", linetype = "dotted", size=0.5),
                         panel.grid.minor = element_line(colour = "grey", linetype = "dotted", size=0.5),
                         axis.title.x = element_blank(), # element_text(size=30, vjust=-1.2, face="bold"),
                         axis.title.y = element_blank(), # element_text(size=30, vjust=-0.1, face="bold", angle = 90),
                         strip.text.x = if (m < max.m) { element_blank() } else { element_text(size=20, face="bold") },
                         strip.text.y = if (n > 1    ) { element_blank() } else { element_text(size=20, face="bold", angle = -90) },
                         axis.text.x  = if (m < max.m) { element_blank() } else { element_text(vjust=1, hjust=0, angle=-35, color="black", face="bold") },
                         axis.text.y  = if (n > 1    ) { element_blank() } else { element_text(color="black") }
                         ) +
                   scale_x_continuous(limits=vx.lim,
                                      breaks=vx.brk,
                                      labels=if(vx.logflg) { sprintf("10e+%02d", vx.brk) } else { vx.brk }) +
                   scale_y_continuous(limits=vy.lim,
                                      breaks=vy.brk,
                                      labels=if(vy.logflg) { sprintf("10e+%02d", vy.brk) } else { vy.brk }) +
                   geom_smooth(data=VVXYG,
                               aes(x = VXVAR,
                                   y = VYVAR,
                                   group=as.factor(GPVAR),
                                   colour=as.factor(GPVAR),
                                   fill=as.factor(GPVAR)),
                               method="lm",
                               alpha=.2,
                               size=0.5,
                               linetype="dotted",
                               fullrange=TRUE) +
                   geom_smooth(data=VVXYG,
                               aes(x = VXVAR,
                                   y = VYVAR),
                               method="lm",
                               alpha=1,
                               size=0.5,
                               se=FALSE,
                               color="black",
                               linetype="dashed",
                               fullrange=TRUE);
          
        }

      }

      CEP[[KK]] <- cep;
      KK <- as.character(as.numeric(KK) + 1);

    }; # for m

  }; # for n

  CEP

} # corrplotfunctCOL


        
# Multiple plot function (Modified)
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplotpairs <- function(..., plotlist=NULL, file, cellsize=NULL, cols=1, title=NULL, layout=NULL) {
  require(grid);

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist);

  numPlots = length(plots);
  
  if (missing(cellsize) || is.null(cellsize)) {
    scols <- 10; # size of cols
  } else {
    scols <- cellsize;
  }
  fct <- 1.2;
  thyrows <- ceiling(numPlots/cols);

  # def.par <- par();
  # par(oma=c(0,0,0,0));
  
  if (is.null(title)) {
      # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout.mtx <- matrix(seq(1, cols * thyrows),
                             ncol = cols, nrow = thyrows);
        layout.lst <- list(widths  = unit(x = c(scols * fct, rep(scols, cols - 1)), units = "cm"),
                           heights = unit(x = c(rep(scols, thyrows - 1), scols * fct), units = "cm"));
      };
  } else {
     titleplot <- list(
                    ggplot(data=data.frame()) +
                         geom_point() +
                         theme(panel.background = element_rect(fill = "#00000000", colour = "white"),
                               plot.margin  = unit(x = c(2, 2, 2, 2), units = "mm"), # top, right, bottom, and left
                               axis.ticks   = element_blank(),
                               panel.grid   = element_blank(),
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               axis.text.x  = element_blank(),
                               axis.text.y  = element_blank()
                               ) +
                         xlim(0,1) + ylim(0, 1) +
                         annotate("text", label = title, x = 0.5, y = 0,
                                  hjust = 0.5, vjust = 0, size = 8, colour = "black", face="bold")
                    );
                       
     plots <- c(titleplot, plots);
     numPlots <- numPlots + 1;
     if (is.null(layout)) {
         layout.mtx <- rbind(rep(1, cols),
                             matrix(seq(1, cols * thyrows) + 1,
                                    ncol = cols, nrow = thyrows));
         layout.lst <- list(widths  = unit(x = c(scols * fct, rep(scols, cols - 1)), units = "cm"),
                            heights = unit(x = c(scols / 2, rep(scols, thyrows - 1), scols * fct), units = "cm"));
     };
  };

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout.mtx), ncol(layout.mtx),
                                               widths  = layout.lst$widths,
                                               heights = layout.lst$heights)))

    # grid.show.layout(grid.layout(nrow(layout.mtx), ncol(layout.mtx), widths  = layout.lst$widths, heights = layout.lst$heights))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      # if (is.null(title)) {
        matchidx <- as.data.frame(which(layout.mtx == i, arr.ind = TRUE));
      # } else {
      #   matchidx <- as.data.frame(which(layout.mtx == i + 1, arr.ind = TRUE));
      # };

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }

  # par(def.par);

} # multiplotpairs
