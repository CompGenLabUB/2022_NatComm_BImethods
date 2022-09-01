#
# #######################################################################
#
#   Several R functions provided by the topGO package
#   have been adapted/extended on this file by Josep F Abril @ CGL,
#   to provide specific visualizations for the analyses on this project.
#
# #######################################################################
#

###
# https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html
library(DOSE);
library(data.table);
library(topGO);
## 
## 
## ------------------------------------------------------------------------------------------------ ##
## IMPORTANT NOTE: "_fix" functions have been adapted/extended from the original topGO R functions. ##
## ------------------------------------------------------------------------------------------------ ##
## 

THYSG.NODES <- 5;

# readMappings from topGO lib
# GODIR must be predefined
go_mappings.GN <- readMappings(paste0(GODIR,"/smesg-go.onlysmesg.tsv"));
go_mappings.TR <- readMappings(paste0(GODIR,"/smest-go.onlysmest.tsv"));
all_gene_names.GN <- names(go_mappings.GN);
all_gene_names.TR <- names(go_mappings.TR);

library(GOplot);
library(ggpubr);
library(grid);
library(gridExtra);
library(cowplot);   # get_legend() & plot_grid() functions
library(patchwork);
library(gtable);
library(stringr);

#a1c57f99 - GREEN
#b1caf699 - BLUE  #6495ed
#f69c9c99 - RED   #ff3030
      # annotation for pheatmap
      # annotation_col <- data.frame(Condition = samples.info$condition,
      #                                   Time = samples.info$time,
      #                              Replicate = samples.info$replicate), annotation_col = annotation_col

# Draw adjacent table for GOBubble and GOCircle
draw_table <- function(data, col){
  id <- term <- NULL;
  colnames(data) <- tolower(colnames(data));
  if (missing(col)){
      tt1 <- ttheme_default()
  } else {
      text.col <- c(rep(col[1], sum(data$category == 'BP')),
                    rep(col[2], sum(data$category == 'CC')),
                    rep(col[3], sum(data$category == 'MF')));
      tt1 <- ttheme_minimal(core = list(bg_params = list(fill = text.col, col=NA, alpha= 1/3)), 
                            colhead = list(fg_params = list(col = "black")));
  };
  table <- tableGrob(subset(data, select = c(id, term)),
                     cols = c('ID', 'Description'), rows = NULL, theme = tt1);
  return(table);
} # draw_table was taken from GOplot/R/Helper.R

# theme_blank was also taken from GOplot/R/Helper.R
theme_blank <- theme(axis.line = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     panel.background = element_blank(),
                     panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     plot.background = element_blank());

GOBar_fix <- function (data, display, order.by.zscore = T, title, zsc.col,
                       zscolabs = c()) 
{
    id <- adj_pval <- zscore <- NULL;
    if (missing(display)) 
        display <- "single";
    if (missing(title)) 
        title <- "";
    if (missing(zsc.col)) 
        zsc.col <- c("firebrick1", "white", "dodgerblue1");
    if (length(zscolabs) == 0) {
        zscolabs <- c("decreasing", "", "increasing");
    };

    colnames(data) <- tolower(colnames(data)); ### TOLOWER!!!
    data$adj_pval <- -log(data$adj_pval, 10);
    data <- data %>%
        filter(!is.na(adj_pval)) %>%
        filter(!is.na(zscore)) %>%
        filter(!is.na(logfc));
    
    sub <- data[!duplicated(data$term), ];
    mxsubsco <- max(abs(sub$zscore));
    if (order.by.zscore == T) {
        sub <- sub[order(sub$zscore, decreasing = T), ];
        leg <- theme(legend.position = "bottom");
        g <- ggplot(sub, aes(x = factor(id, levels = stats::reorder(id, adj_pval)),
                             y = adj_pval,
                             fill = zscore)) +
            geom_bar(stat = "identity", colour = "black") +
            scale_fill_gradient2("z-score", space = "Lab",
                                 low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1],
                                 guide = guide_colourbar(title.position = "top", title.hjust = 0.5),
                                 limits = c(-mxsubsco,mxsubsco),
                                 breaks = c(-mxsubsco, 0, mxsubsco),
                                 labels = zscolabs) + 
            labs(title = title,
                 x = "",
                 y = expression(-log[10]("Adjusted p-value"))) + 
            leg;

    }
    else {
        sub <- sub[order(sub$adj_pval, decreasing = T), ];
        leg <- theme(legend.justification = c(1, 1),
                     legend.position = c(0.98, 0.995),
                     legend.background = element_rect(fill = "transparent"), 
                     legend.box = "vertical",
                     legend.direction = "horizontal");
        g <- ggplot(sub, aes(x = factor(id, levels = reorder(id, adj_pval)),
                             y = zscore, fill = adj_pval)) +
            geom_bar(stat = "identity", colour = "black") +
            scale_fill_gradient2("Significance", space = "Lab",
                                 low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1],
                                 guide = guide_colourbar(title.position = "top", title.hjust = 0.5),
                                 breaks = c(min(sub$adj_pval), 
                                            max(sub$adj_pval)),
                                 labels = c("low", "high")) + 
            labs(title = title, x = "", y = "z-score") +
            leg;
    };
    if (display == "single") {
        g + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                  axis.line = element_line(colour = "grey80"), 
                  axis.ticks = element_line(colour = "grey80"),
                  axis.title = element_text(size = 14, face = "bold"),
                  axis.text = element_text(size = 14), 
                  panel.background = element_blank(),
                  panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  plot.background = element_blank());
    }
    else {
        g + facet_grid(. ~ category, space = "free_x", scales = "free_x") + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                  axis.line = element_line(colour = "grey80"), 
                  axis.ticks = element_line(colour = "grey80"), 
                  axis.title = element_text(size = 14, face = "bold"), 
                  axis.text = element_text(size = 14),
                  panel.background = element_blank(), 
                  panel.border = element_blank(),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.background = element_blank());
    };
} # GOBar_fix

GOHeat_fix <- function (data, nlfc, fill.col, thygo = NULL) 
{
  x <- y <- z <- NULL;
  if (missing(nlfc)) 
    nlfc <- 0
  else nlfc <- nlfc;
  if (missing(fill.col)) 
    fill.col <- c("firebrick", "white", "dodgerblue")
  else fill.col <- fill.col;
  if (is.null(thygo)) {
    thygo <- "";
  };
  distance <- dist(data);
  cluster <- hclust(distance);
  M <- dim(data)[2];
  nterm <- M - nlfc;
  if (nlfc == 0) {
    s <- rowSums(data[, 1:nterm]);
    tmp <- NULL;
    for (r in 1:nrow(data)) {
      tmp <- c(tmp, as.numeric(gsub(1, s[r], data[r, 1:nterm])));
    };
  }
  else {
    tmp <- NULL;
    for (r in 1:nrow(data)) {
      tmp <- c(tmp, as.numeric(gsub(1, data[r, (nterm + 1)], data[r, 1:nterm])));
    };
  };
  df <- data.frame(x = factor(rep(cluster$order, each = nterm)),
                   y = str_trunc(rep(colnames(data[,1:nterm]), length(rownames(data))),
                                 60, side="right", ellipsis = "…"),
                   z = tmp, lab = rep(rownames(data), each = nterm));
  df_o <- df[order(df$x), ];
  g <- ggplot() +
    geom_tile(data = df_o, aes(x = x, y = y, fill = z)) +
    scale_x_discrete(breaks = 1:length(unique(df_o$x)), labels = unique(df_o$lab)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_text(size = 10),
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    coord_fixed(clip = 'off') +
    annotate("text", label=thygo, x=0, y=0, hjust=2, vjust=2, size = rel(12));
  if (nlfc == 0) {
    g + scale_fill_gradient2("Count", space = "Lab", low = fill.col[2], 
                             mid = fill.col[3], high = fill.col[1]);
  }
  else {
    g + scale_fill_gradient2("logFC", space = "Lab", low = fill.col[3], 
                             mid = fill.col[2], high = fill.col[1]);
  };
} # GOHeat_fix

GOBubble_fix <- function(data, display, title, colour,
                         labels, ID = T, table.legend = T, 
                         table.col = T, bg.col = F,
                         pval.label = NULL) 
{
    zscore <- adj_pval <- category <- count <- id <- term <- NULL;
    if (missing(display)) 
        display <- "single";
    if (missing(title)) 
        title <- "";
    if (missing(colour)) 
        cols <- c("chartreuse4", "brown2", "cornflowerblue")
    else cols <- colour;
    if (missing(labels)) 
        labels <- 5;
    if (bg.col == T & display == "single") 
        cat("Parameter bg.col will be ignored. To use the parameter change display to 'multiple'");
    if (is.null(pval.label))
        pval.label <- parse(text= expression(-log[10]("adj. p-value")));
    
    colnames(data) <- tolower(colnames(data));
    if (!"count" %in% colnames(data)) {
        rang <- c(5, 5);
        data$count <- rep(1, dim(data)[1]);
    }
    else {
        rang <- c(1, 30);
    };
    data$adj_pval <- -log(data$adj_pval, 10);

    sub <- data[!duplicated(data$term), ] %>%
                   filter(!is.na(adj_pval)) %>%
                   filter(!is.na(zscore));
    if (!is.character(labels)) {
        sub2 <- subset(sub, subset = sub$adj_pval >= labels);
        lab.threshold <- labels;
    } else {
        sub2 <- subset(sub, sub$id %in% labels | sub$term %in% labels);
        lab.threshold <- min(sub2$adj_pval);
    };

    g <- ggplot(sub, aes(zscore, adj_pval, fill = category, size = count)) + 
            labs(title = title, x = "z-score", y = pval.label) +
            geom_point(shape = 21, col = "black", alpha = 1/2) + 
            geom_hline(yintercept = lab.threshold, col = "orange") + # yintercept = 1.3
            scale_size(range = rang, guide = "none");
    if (display == "single") {
        g <- g + scale_fill_manual("Category", values = cols, 
                                   labels = c("Biological Process", "Cellular Component", "Molecular Function")) +
            theme(legend.position = "bottom") + 
            annotate("text", x = min(sub$zscore), y = lab.threshold * 0.99, # x = min(sub$zscore) + 0.2, y = 1.4, 
                     label = paste0("Threshold adj.pval=",sprintf("%.2f",10^-lab.threshold)),
                     colour = "orange", size = 4, hjust=0, vjust=1)
        if (ID) 
            g <- g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, label = id), size = 5)
        else
            g <- g + geom_text(data = sub2, aes(x = zscore, y = adj_pval, label = term), size = 4);
        if (table.legend) {
            if (table.col) 
                table <- draw_table(sub2, col = cols)
            else table <- draw_table(sub2);
            g <- g + theme(axis.text = element_text(size = 14), 
                           axis.line = element_line(colour = "grey80"), 
                           axis.ticks = element_line(colour = "grey80"), 
                           axis.title = element_text(size = 14, face = "bold"), 
                           panel.background = element_blank(),
                           panel.grid.minor = element_blank(), 
                           panel.grid.major = element_line(colour = "grey80"), 
                           plot.background = element_blank());
            graphics::par(mar = c(0.1, 0.1, 0.1, 0.1));
            grid.arrange(g, table, ncol = 2);
        }
        else {
            g + theme(axis.text = element_text(size = 14),
                      axis.line = element_line(colour = "grey80"), 
                      axis.ticks = element_line(colour = "grey80"), 
                      axis.title = element_text(size = 14, face = "bold"), 
                      panel.background = element_blank(),
                      panel.grid.minor = element_blank(), 
                      panel.grid.major = element_line(colour = "grey80"), 
                      plot.background = element_blank());
        };
    }
    else {
        if (bg.col) {
            dummy_col <- data.frame(category = c("BP", "CC", "MF"),
                                    adj_pval = sub$adj_pval[1:3], zscore = sub$zscore[1:3], 
                                    size = 1:3, count = 1:3);
            g <- g + geom_rect(data = dummy_col, aes(fill = category), 
                               xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
                               alpha = 0.1) +
                     facet_grid(. ~ category, space = "free_x", scales = "free_x") +
                     scale_fill_manual(values = cols, guide = "none");
        }
        else {
            g <- g + facet_grid(. ~ category, space = "free_x", scales = "free_x") +
                     scale_fill_manual(values = cols, guide = "none");
        }
        if (ID) {
            g + geom_text(data = sub2,
                          aes(x = zscore, y = adj_pval, label = id),
                          size = 5) +
                theme(axis.title = element_text(size = 14, face = "bold"),
                      axis.text = element_text(size = 14), 
                      axis.line = element_line(colour = "grey80"), 
                      axis.ticks = element_line(colour = "grey80"), 
                      panel.border = element_rect(fill = "transparent", colour = "grey80"),
                      panel.background = element_blank(), 
                      panel.grid = element_blank(),
                      plot.background = element_blank());
        }
        else {
            g + geom_text(data = sub2,
                          aes(x = zscore, y = adj_pval, label = term),
                          size = 5) +
                theme(axis.title = element_text(size = 14, face = "bold"),
                      axis.text = element_text(size = 14), 
                      axis.line = element_line(colour = "grey80"), 
                      axis.ticks = element_line(colour = "grey80"), 
                      panel.border = element_rect(fill = "transparent", colour = "grey80"),
                      panel.background = element_blank(), 
                      panel.grid = element_blank(),
                      plot.background = element_blank());
        };
    };
} # GOBubble_fix

GOCircle_fix <- function (data, title, nsub, rad1, rad2, table.legend = T, zsc.col, 
                          lfc.col, label.size, label.fontface, thygo = NULL, thygocol = NULL,
                          updnlabel = NULL, updnlabels = c("upregulated", "downregulated"),
                          zscolabs = c()) 
{
    xmax <- y1 <- zscore <- y2 <- ID <- logx <- logy2 <- logy <- logFC <- NULL;

    if (missing(title)) 
        title <- "";
    if (missing(nsub)) 
        if (dim(data)[1] > 10) 
            nsub <- 10
        else nsub <- dim(data)[1];
    if (missing(rad1)) 
        rad1 <- 2;
    if (missing(rad2)) 
        rad2 <- 3;
    if (missing(zsc.col)) 
        zsc.col <- c("red", "white", "blue");
    if (missing(lfc.col))
        lfc.col <- c("cornflowerblue", "firebrick1")
    else lfc.col <- rev(lfc.col);
    if (missing(label.size)) 
        label.size = 6;
    if (missing(label.fontface)) 
        label.fontface = "bold";
    
    if (is.null(thygo)) {
        thygo <- "";
    };
    if (is.null(thygocol)) {
        thygocol <- "#000000";
    };
    if (is.null(updnlabel)) {
        updnlabel <- parse(text= expression(log[2]("Fold Change")));
    };
    if (length(updnlabels) == 0) {
        updnlabels <- c("upregulated", "downregulated");
                   #  c(parse(text= expression(log[2](FC)>=0)),
                   #    parse(text= expression(log[2](FC)<0)));
    };
    if (length(zscolabs) == 0) {
        zscolabs <- c("decreasing", "", "increasing");
    };

    colnames(data) <- toupper(colnames(data)); ### TOUPPER!!!
    data$ADJ_PVAL <- -log(data$ADJ_PVAL, 10);
    data <- data %>%
        filter(!is.na(ADJ_PVAL)) %>%
        filter(!is.na(ZSCORE)) %>%
        filter(!is.na(LOGFC));
    
    suby <- data[!duplicated(data$TERM), ] %>%
              filter(!is.na(ADJ_PVAL)) %>%
              filter(!is.na(ZSCORE)) %>%
              filter(!is.na(LOGFC));
    if (is.numeric(nsub) == T) {
        suby <- suby[1:nsub, ];
    }
    else {
        if (strsplit(nsub[1], ":")[[1]][1] == "GO") {
            suby <- suby[suby$ID %in% nsub, ];
        }
        else {
            suby <- suby[suby$TERM %in% nsub, ];
        };
        nsub <- length(nsub);
    };
    N <- dim(suby)[1];
    r_pval <- round(range(suby$ADJ_PVAL), 0) + c(-2, 2);
    ymax <- c();
    for (i in 1:length(suby$ADJ_PVAL)) {
        val <- (suby$ADJ_PVAL[i] - r_pval[1])/(r_pval[2] - r_pval[1]);
        ymax <- c(ymax, val);
    };
    df <- data.frame(x = seq(0, 10 - (10/N), length = N),
                     xmax = rep(10/N - 0.2, N),
                     y1 = rep(rad1, N),
                     y2 = rep(rad2, N),
                     ymax = ymax, 
                     zscore = suby$ZSCORE, ID = suby$ID);
    scount <- data[!duplicated(data$TERM), which(colnames(data) == "COUNT")][1:nsub];
    idx_term <- which(!duplicated(data$TERM) == T);
    xm <- c();
    logs <- c();
    for (sc in 1:length(scount)) {
        idx <- c(idx_term[sc], idx_term[sc] + scount[sc] - 1);
        val <- stats::runif(scount[sc], df$x[sc] + 0.06, (df$x[sc] + df$xmax[sc] - 0.06));
        xm <- c(xm, val);
        r_logFC <- round(range(data$LOGFC[idx[1]:idx[2]]), 0) + c(-1, 1);
        for (lfc in idx[1]:idx[2]) {
            val <- (data$LOGFC[lfc] - r_logFC[1])/(r_logFC[2] - r_logFC[1]);
            logs <- c(logs, val);
        };
    };
    cols <- c();
    for (ys in 1:length(logs))
        cols <- c(cols, ifelse(data$LOGFC[ys] > 0, "upregulated", "downregulated"));
    
    dfp <- data.frame(logx = xm, logy = logs,
                      logFC = factor(cols), 
                      logy2 = rep(rad2, length(logs)));

    mxsubsco <- max(abs(df$zscore));

    c0 <- ggplot() +  # this is the main plot, without legends
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y1, ymax = y1 + ymax, fill = zscore),
                  colour = "black") + 
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y2, ymax = y2 + 1),
                  fill = "gray70") +
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y2 + 0.5, ymax = y2 + 0.5),
                  colour = "white", linetype= "dotted") + # colour = "yellow", linetype= "dotted") +
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y2 + 0.25, ymax = y2 + 0.25), 
                  colour = "white", linetype= "dotted") +
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y2 + 0.75, ymax = y2 + 0.75), 
                  colour = "white", linetype= "dotted") +
        geom_text(data = df, aes(x = x + (xmax/2), y = y2 + 1.15,
                                 label = ID, angle = 360 - (x = x + (xmax/2))/(10/360)),
                  size = label.size, fontface = label.fontface) + 
        coord_polar() + labs(title = title) +
        ylim(1, rad2 + 1.6) + xlim(0, 10) +
        theme_blank +
        geom_point(data = dfp, aes(x = logx, y = logy2 + logy), 
                   pch = 21, fill = "transparent", colour = "black", size = 3) +
        geom_point(data = dfp, aes(x = logx, y = logy2 + logy, colour = logFC),
                   size = 2.5)+
        scale_colour_manual(updnlabel,
                            values = lfc.col,
                            guide = guide_legend(title.position = "top",
                                                 title.hjust = 0.5, nrow=2,
                                                 title.size=rel(0.7)),
                            labels = updnlabels) +
        scale_fill_gradient2("z-score", space = "Lab",
                             low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1], 
                             guide = guide_colourbar(title.position = "top",
                                                     title.hjust = 0.5,
                                                     title.size=rel(1.5)), 
                             limits = c(-mxsubsco,mxsubsco),
                             breaks = c(-mxsubsco,0,mxsubsco), # c(min(df$zscore), 0, max(df$zscore)),
                             labels = zscolabs)  + #, legend.spacing.x = unit(1.0, 'cm')) +
        theme(plot.margin = ggplot2::margin(0,0,0,0,"pt"),
              legend.position = "none", # legend.position = "bottom", 
              legend.background = element_rect(fill = "transparent"), 
              legend.box = "horizontal",
              legend.direction = "horizontal");
    c1 <- ggplot() + # this plot to generate first legend grob
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y1, ymax = y1 + ymax, fill = zscore),
                  colour = "black") + 
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y2, ymax = y2 + 1),
                  fill = "gray70") +
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y2 + 0.5, ymax = y2 + 0.5),
                  colour = "white") +
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y2 + 0.25, ymax = y2 + 0.25), 
                  colour = "white") +
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y2 + 0.75, ymax = y2 + 0.75), 
                  colour = "white") +
        geom_text(data = df, aes(x = x + (xmax/2), y = y2 + 1.15,
                                 label = ID, angle = 360 - (x = x + (xmax/2))/(10/360)),
                  size = label.size, fontface = label.fontface) + 
        coord_polar() + labs(title = title) +
        ylim(1, rad2 + 1.6) + xlim(0, 10) +
        theme_blank +
        geom_point(data = dfp, aes(x = logx, y = logy2 + logy), 
                   pch = 21, fill = "transparent", colour = "black", size = 3) +
        geom_point(data = dfp, aes(x = logx, y = logy2 + logy, colour = logFC),
                   size = 2.5) + 
        scale_colour_manual(updnlabel,
                            values = lfc.col,
                            guide = guide_legend(title.position = "top",
                                                 title.hjust = 0.5, nrow=2,
                                                 title.size=rel(0.7)),
                            labels = updnlabels) +
        scale_fill_continuous(guide = 'none') + #, legend.spacing.x = unit(1.0, 'cm')) +
        theme(plot.margin = ggplot2::margin(0,0,0,0,"pt"),
              legend.position = "bottom", 
              legend.background = element_rect(fill = "transparent"), 
              legend.box = "horizontal",
              legend.direction = "horizontal");
    c2 <- ggplot() + # this plot to generate second legend grob
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y1, ymax = y1 + ymax, fill = zscore),
                  colour = "black") + 
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y2, ymax = y2 + 1),
                  fill = "gray70") +
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y2 + 0.5, ymax = y2 + 0.5),
                  colour = "white") +
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y2 + 0.25, ymax = y2 + 0.25), 
                  colour = "white") +
        geom_rect(data = df, aes(xmin = x, xmax = x + xmax,
                                 ymin = y2 + 0.75, ymax = y2 + 0.75), 
                  colour = "white") +
        geom_text(data = df, aes(x = x + (xmax/2), y = y2 + 1.15,
                                 label = ID, angle = 360 - (x = x + (xmax/2))/(10/360)),
                  size = label.size, fontface = label.fontface) + 
        coord_polar() + labs(title = title) +
        ylim(1, rad2 + 1.6) + xlim(0, 10) +
        theme_blank+
        geom_point(data = dfp, aes(x = logx, y = logy2 + logy), 
                   pch = 21, fill = "transparent", colour = "black", size = 3) +
        geom_point(data = dfp, aes(x = logx, y = logy2 + logy, colour = logFC),
                   size = 2.5) +
        scale_colour_discrete(guide = 'none') +
        scale_fill_gradient2("z-score", space = "Lab",
                             low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1], 
                             guide = guide_colourbar(title.position = "top",
                                                     title.hjust = 0.5,
                                                     title.size=rel(1.5)), 
                             limits = c(-mxsubsco,mxsubsco),
                             breaks = c(-mxsubsco,0,mxsubsco), # c(min(df$zscore), 0, max(df$zscore)),
                             labels = zscolabs) + #, legend.spacing.x = unit(1.0, 'cm')) +
        theme(plot.margin = ggplot2::margin(0,0,0,0,"pt"),
              legend.position = "bottom", 
              legend.background = element_rect(fill = "transparent"), 
              legend.box = "horizontal",
              legend.direction = "horizontal");

    # merging all grobs to create radial plot with properly placed legends
    leg1 <- get_legend(c1);
    leg2 <- get_legend(c2);
    blank_p <- plot_spacer() + theme_void();
    leg12 <- plot_grid(blank_p, leg1, blank_p, leg2, blank_p, nrow = 1,
                       rel_widths = c(0.5, 1, 0.5, 1, 0.5));
    # c <- plot_grid(c0, leg12, nrow = 2,
    #                align = "h",
    #                axis = "t",
    #                rel_heights = c(2.75, 0.25));
    c <- ggdraw(blank_p + xlim(-1,1) + ylim(-1,1) +
                draw_plot(c0,-1.325,-1.275,2.75,2.75) +
                draw_plot(leg1,-1,-1.45,1,1) +
                draw_plot(leg2,0.05,-1.425,1,1)) +
                geom_text(aes(label=thygo, x=0.5275, y=0.535),
                          hjust=0.5, vjust=0.5, size = rel(14), color = thygocol);
    
    if (table.legend) {
        table <- draw_table(suby);
        graphics::par(mar = c(0.1, 0.1, 0.1, 0.1));
        grid.arrange(c, table, ncol = 2);
    }
    else {
        c + theme(plot.background = element_rect(fill = "aliceblue"), 
                  panel.background = element_rect(fill = "white"));
    };
} # GOCircle_fix

darken <- function(color, factor=1.4){
    col <- col2rgb(color);
    col <- col/factor;
    col <- rgb(t(col), maxColorValue=255);
    col
}

lighten <- function(color, factor=1.4){
    col <- col2rgb(color);
    col <- col*factor;
    col <- rgb(t(col), maxColorValue=255);
    col
}

GOfoldchange <- function(circ, title, ofile="foldchange+zscore_GOplot", oextns="png",
                         col = NULL, thygo = NULL, nsub, zsc.col, lfc.col, fclim = NULL,
                         updnlabel = NULL, updnlabels = c(), zscolabs = c(), 
                         debugflg=FALSE, saveflg=FALSE) 
{        
        if (missing(title)) 
            title <- "";
        if (missing(zsc.col)) 
            zsc.col <- c("red", "white", "blue");
        if (missing(lfc.col)) {
            lfc.col <- c("cornflowerblue", "firebrick1");
        } else {
            lfc.col <- rev(lfc.col);
        };

        if (is.null(thygo)) {
            thygo <- "";
        };

        if (is.null(updnlabel)) {
           updnlabel <- parse(text= expression(log[2]("Fold Change")));
        };
        if (length(updnlabels) == 0) {
            updnlabels <- c(parse(text= expression(log[2](FC)>=0)),
                            parse(text= expression(log[2](FC)<0)));
        };
        if (length(zscolabs) == 0) {
            zscolabs <- c("decreasing", "", "increasing");
        };
        
        colnames(circ) <- toupper(colnames(circ)); ### TOUPPER!!!
        dcir <- circ[order(-circ$COUNT, circ$ID), ];
        dcir$TERM <- str_trunc(dcir$TERM, 50, side="right", ellipsis = "…") # "\u2026" "..." "…";
        suby <- dcir[!duplicated(dcir$ID), c("CATEGORY","ID","TERM","COUNT","ADJ_PVAL","ZSCORE")] %>%
                   filter(!is.na(ADJ_PVAL)) %>%
                   filter(!is.na(ZSCORE));
        suby.brks <- suby$ID;
        suby.lbls <- paste0(suby$ID," ",suby$TERM);
        suby.nrow <- nrow(suby);
        
        if (missing(nsub)) { 
            nsub <- suby.nrow;
        } else {
            if (nsub > suby.nrow) { nsub <- suby.nrow };
            suby <- suby[ c(1:nsub), ];
            suby.brks <- suby$ID;
            suby.lbls <- paste0(suby$ID," ",suby$TERM);
            suby.nrow <- nrow(suby);
        };

        if (is.null(col)) {
            tt1 <- ttheme_default(base_size = 8);
        } else {
            text.col <- rep(c(col, darken(col)), ceiling(suby.nrow)/2);
            tt1 <- ttheme_minimal(core = list(bg_params = list(fill = text.col, col=NA, alpha= 1/3)), 
                                  colhead = list(fg_params = list(col = "black")),
                                  base_size = 8 );
        };
        dcir$ID <- factor(dcir$ID, levels=suby.brks);
        dcir$logFCcol <- factor(ifelse(dcir$LOGFC < 0, "cornflowerblue", "firebrick1"),
                                levels=c("cornflowerblue", "firebrick1"));
        table <- tableGrob(subset(suby, select = c(ID, TERM)),
                           cols = NULL, # cols = c('ID', 'Description'),
                           rows = NULL,
                           theme = tt1);
        table$heights <- unit(rep(1,nrow(table)), "null");

        # margins -> unit(c(t, r, b, l), "cm")
        thyrmar <- 0; # floor(log10(suby.nrow+0.1))*2
        fmax <- max(abs(min(dcir$LOGFC)),abs(max(dcir$LOGFC)),5);
        if (is.null(fclim)) {
            fclim <- c(fmax,-fmax);
        };
        fcplot <- ggplot(dcir,
                         aes(x=factor(ID, levels=rev(suby.brks)), y=LOGFC, group=ID)) +
                    geom_point(color="black", stroke=1.5) +
                    geom_point(aes(color=logFCcol)) + theme_bw() +
                    geom_text(aes(x=factor(ID, levels=rev(suby.brks)), y=0, label=COUNT),
                              size=1.5, nudge_x=0.1,nudge_y=0.1,hjust=0,vjust=1,color="green2") +
                    scale_y_reverse(limits=fclim) +
                    ylab(updnlabel) +
                    scale_x_discrete(breaks=rev(suby.brks),
                                     labels=rev(suby.brks), position = "top") +
                    coord_flip() +
                    ggtitle(title) +
                    scale_color_manual(name="",
                                       values=c("firebrick1", "cornflowerblue"),
                                       breaks=c("firebrick1", "cornflowerblue"),
                                       labels=updnlabels,
                                       drop = FALSE) +
                    guides(color = guide_legend(override.aes = list(size = 3), nrow=2)) +
                    theme(plot.margin = unit(c(25, thyrmar, 25, 5), "pt"), # t, r, b, l
                          plot.background = element_blank(),
                          legend.position="bottom",
                          legend.direction="vertical",
                          legend.spacing.x = unit(0.1, 'cm'),
                          legend.title.align = 1,
                          legend.text=element_text(size=rel(1),vjust=0.5,hjust=0),
                          legend.title=element_blank(), # element_text(size=rel(1.5),family="bold",vjust=0.5,hjust=1.5),
                          axis.title.y=element_blank(),
                          # axis.title.x=element_blank(), # axis.text.x=element_blank(),
                          panel.grid.major.x=element_line(linetype="dashed",color="black",size=0.35),
                          panel.grid.minor.x=element_line(linetype="dotted",color="black",size=0.25));

         pval <- ifelse(suby$ADJ_PVAL>0 , -log10(suby$ADJ_PVAL), 0);
         pmax <- max(min(pval),max(pval),5);
         zmax <- max(abs(min(suby$ZSCORE)),abs(max(suby$ZSCORE)));
         suby$ADJ_PVALfix <- -log10(suby$ADJ_PVAL) + 0.025;

         zplot <- ggplot(suby,
                         aes(x=factor(ID, levels=rev(suby.brks)), y=ADJ_PVALfix, fill=ZSCORE)) +
                    geom_bar(stat="identity", colour="black") +
                    geom_hline(yintercept=c(-log10(0.05),-log10(0.0101)),
                               linetype=c('dashed','dotted'), color="red", size=0.25) + theme_bw() +
                    scale_y_continuous(limits=c(0,pmax), expand = c(0, 0)) +
                    scale_x_discrete(breaks=suby.brks, labels=suby.brks, position = "bottom") +
                    scale_fill_gradient2("z-score", 
                                         space = "Lab", low = zsc.col[3], mid = zsc.col[2], high = zsc.col[1], 
                                         guide = guide_colourbar(title.position = "left", direction = "horizontal",
                                                                 title.vjust = 1, title.hjust = 1), 
                                         limits = c(-zmax, zmax),
                                         breaks = c(-zmax, 0, zmax),
                                         labels = zscolabs) +
                    ylab(expression(-log[10]("Adjusted p-value"))) +
                    coord_flip() +
                    ggtitle("") +
                    theme(plot.margin = unit(c(25, 5, 25, thyrmar), "pt"), # t, r, b, l
                          plot.background = element_blank(),
                          legend.position="bottom",
                          legend.direction="horizontal",
                          legend.text=element_text(size=rel(.75),vjust=1.5,hjust=0.5),
                          legend.title=element_text(size=rel(1.5),family="bold",vjust=0.5,hjust=1.5),
                          axis.title.y=element_blank(), # axis.text.y=element_blank(),
                          panel.grid.major.x=element_line(linetype="dashed",color="black",size=0.35),
                          panel.grid.minor.x=element_line(linetype="dotted",color="black",size=0.25));
                          # ,axis.text.x=element_blank())
        
         # gtable_show_layout(gfcplot)
         # tutorial at https://stackoverflow.com/questions/37984000/how-to-manage-the-t-b-l-r-coordinates-of-gtable-to-plot-the-secondary-y-axi

         if (debugflg) {
                gfcplot <- ggplotGrob(fcplot);
                gzplot  <- ggplotGrob(zplot);
                gfcplot <- gtable::gtable_add_cols(gfcplot, unit(1, "null"), -1);
                gfcplot <- gtable::gtable_add_grob(gfcplot, table, t=7, l=10, r=10);
                ggarrange(gfcplot,gzplot,   ncol = 2, nrow = 1, widths=c(8 ,4), align = "h");
                if (saveflg) {
                    cat("Save combo plot to: ",paste0(ofile,".check.",oextns),"\n");
                    ggsave(paste0(ofile,".check.",oextns),
                           width=12, height=4+(0.1*suby.nrow), dpi=300, units=c("in"),
                           bg="white");
                };
         };

         gfcplot2 <- ggplotGrob(fcplot + theme(axis.text.y=element_blank()));
         gzplot2  <- ggplotGrob(zplot + theme(axis.text.y=element_blank()));
         gfcplot2 <- gtable::gtable_add_cols(gfcplot2, unit(1, "null"), -1);
         gfcplot2 <- gtable::gtable_add_grob(gfcplot2, table, t=7, l=10, r=10);
         ggarrange(gfcplot2,gzplot2, ncol = 2, nrow = 1, widths=c(8 ,4), align = "h");
         if (saveflg) {
             cat("Save combo plot to: ",paste0(ofile,".",oextns),"\n");
             ggsave(paste0(ofile,".",oextns),
                    width=12, height=4+(0.1*suby.nrow), dpi=300, units=c("in"),
                    bg="white");
         };
         # plot(gfcplot) + zplot + plot_layout(ncol=2, nrow=1, widths = c(6,4), heights = c(10,10))
         # p3 <- ggplot() + table
         # fcplot + table + zplot + plot_layout(ncol=3, nrow=1, widths = c(4,8,4), heights = c(10,10,10))
         # fcplot + zplot + plot_layout(ncol=2, nrow=1, widths = c(6,4), heights = c(10,10))
         # ggarrange(fcplot, table, ncol = 2, nrow = 1, widths = c(4,8), heights = c(10,10))
        
} # GOfoldchange

runGOsUPDN_ext <- function(udflg="N", contrast.n, contrast.lbl, set.lbl,
                           THYFUNGO, table_exp, go_mapping, annot_names) {
    
    if (identical(udflg,"U")) {
       udlbl = "UP";
       lfcCOL <- c("#ff3030");
       DE_genes <- subset(table_exp, P_adj <= MIN.PV & log_FC >= MIN.FC)$Gene_ID;
    } else {
       if (identical(udflg,"D")) {
         udlbl = "DOWN";
         lfcCOL <- c("#6495ed");
         DE_genes <- subset(table_exp, P_adj <= MIN.PV & -log_FC >= MIN.FC)$Gene_ID;
       } else  {
         udlbl = "BOTH";
         lfcCOL <- c("#ff3030", "#6495ed");
                     # RED-UP BLUE-DOWN
         DE_genes <- subset(table_exp, P_adj <= MIN.PV & abs(log_FC) >= MIN.FC)$Gene_ID;
       };
    };
    THYGOCOL <- ifelse(identical(THYFUNGO,"BP"),
                       "#a1c57f99",
                       ifelse(identical(THYFUNGO,"MF"),
                              "#b1caf699", "#f69c9c99"));

    cat("# Working on ", contrast.n," ",udlbl, ": ", contrast.lbl, " / ", set.lbl,"\n", sep="");

    gene_list <- factor(as.integer(annot_names %in% DE_genes));
    names(gene_list) <- annot_names;
    cat("# GENE LIST DONE...\n");
    
    cat("# Analysis of ",THYFUNGO,"...\n");
    top_GO_data <- new("topGOdata", ontology = THYFUNGO, allGenes = gene_list,       # all_genes,
                       nodeSize = 10, annot = annFUN.gene2GO, gene2GO = go_mapping); # affyLib = "drosophila2.db")
    cat("# ",THYFUNGO,": topGO DONE...\n");
    
    result_top_GO_elim    <- runTest(top_GO_data, algorithm = "elim",    statistic = "Fisher");
    result_top_GO_classic <- runTest(top_GO_data, algorithm = "classic", statistic = "Fisher");

    rv <- c(as.integer(geneData(result_top_GO_elim)),
            as.integer(geneData(result_top_GO_classic)));
    cat(do.call(sprintf,
                c(fmt=paste0('topGO %s %s %s signif.nodes %d\n',
                             'F.elim: Annot. %d  Signif. %d  NodeSize %d  SigTerms %d\n',
                             'F.clas: Annot. %d  Signif. %d  NodeSize %d  SigTerms %d\n'),
                udlbl, contrast.lbl, THYFUNGO, THYSG.NODES, as.list(rv))));
    cat("# ",THYFUNGO,": runTest DONE...\n");

    res_top_GO <- GenTable(top_GO_data,
                           Fisher.elim    = result_top_GO_elim,
                           Fisher.classic = result_top_GO_classic,
                           orderBy = "Fisher.elim", topNodes = 100, numChar=1000);
    res_top_GO$Category    <- rep(THYFUNGO, nrow(res_top_GO));
    res_top_GO$AnnotatedGenes   <- lapply(res_top_GO$GO.ID,
                                     function(x) as.character(unlist(genesInTerm(object = top_GO_data, whichGO = x))));
    res_top_GO$SignificantGenes <- lapply(res_top_GO$AnnotatedGenes, function(x) intersect(x,  DE_genes));
    res_top_GO$AnnotatedGenes   <- sapply(res_top_GO$AnnotatedGenes, paste, collapse=", ");
    res_top_GO$SignificantGenes <- sapply(res_top_GO$SignificantGenes, paste, collapse=", ");

    write.table(res_top_GO,
                file=paste0(statdir,"topGOdata.extgenetable.",contrast.lbl,"_",set.lbl,'_',udlbl,"_",THYFUNGO,'_sn',THYSG.NODES,".tbl"),
                sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE);
    cat("# ",THYFUNGO,": Extended GenTable DONE...\n");

    # Generate the GOplot plotting object
    cat("# GOplot object...\n");
    EC$goa  <- res_top_GO[ ,c("Category","GO.ID","Term","SignificantGenes","Fisher.elim") ];
    colnames(EC$goa)  <- c("Category","ID","Term","Genes","adj_pval");
    EC$glst  <- table_exp[ DE_genes, c("Gene_ID","log_FC","Avg_Expr","stat","P_value","P_adj","B") ];
    colnames(EC$glst) <- c("ID","logFC","AveExpr","t","P.Value","adj.P.Val","B");
    EC$genes <- table_exp[ DE_genes, c("Gene_ID","log_FC") ];
    colnames(EC$genes) <- c("ID","logFC");
    EC$process <- head(EC$goa$Term,25);
    circ <- circle_dat(EC$goa,EC$glst);
    circ$adj_pval <- as.numeric(circ$adj_pval); # fix needed 
    chord <- chord_dat(circ, EC$genes, EC$process);
    write.table(circ[ !is.na(circ$logFC), ],
                file=paste0(statdir,"topGOdata.zscores.",contrast.lbl,"_",set.lbl,'_',udlbl,"_",THYFUNGO,'_sn',THYSG.NODES,".tbl"),
                sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE);
    cat("# ",THYFUNGO,": zscores GenTable DONE...\n");
    circ$term <- str_trunc(circ$term, 50, side="right", ellipsis = "…"); # "\u2026" "..." "…"
    cmax.rows <- min(nrow(chord), 100);
    cmax.cols <- min(ncol(chord), 100);

    cat("# combo-plot...\n");
    skiptonext <- FALSE;
    tryCatch(phmb <- GOfoldchange(circ,
                                  title=paste0('FoldChange + Z-score for',contrast.lbl,' ',set.lbl,' ',udlbl,' ',THYFUNGO,' sn',THYSG.NODES),
                                      ofile=paste0(imgdir,'GOplots_comboGOplot_FC+zSco_',contrast.lbl,'_',set.lbl,'_',udlbl,"_",THYFUNGO,'_sn',THYSG.NODES,'_annot'), oextns="png",
                                  col = THYGOCOL, thygo = THYFUNGO, debugflg=TRUE, saveflg=TRUE),
             error = function(e) { skiptonext <<- TRUE });

    cat("# hist-plot...\n");
    skiptonext <- FALSE;
    tryCatch(phmb <- GOBar(circ, title=paste0('Z-score barplot for',contrast.lbl,' ',set.lbl,' ',udlbl,' ',THYFUNGO,' sn',THYSG.NODES)) +
                          geom_hline(yintercept=c(-log10(0.05),-log10(0.0101)), linetype=c('dashed','dotted'), color="red", size=0.25),
             error = function(e) { skiptonext <<- TRUE });
    if (skiptonext) {
      skiptonext <- FALSE;
    } else {
      ggplot2::ggsave(plot=phmb,
                      filename=paste0(imgdir,'GOplots_histplot_',contrast.lbl,'_',set.lbl,'_',udlbl,"_",THYFUNGO,'_sn',THYSG.NODES,'_annot.png'),
                      height=max(4, 4+(16 * cmax.rows/100)), width=16, units="in",dpi=600,
                      bg="white");
    };

    cat("# bubble-plot...\n");
    skiptonext <- FALSE;
    tryCatch(phmb <- GOBubble(circ, labels = 2, colour = c(THYGOCOL)),
             error = function(e) { skiptonext <<- TRUE });
    if (skiptonext) {
      skiptonext <- FALSE;
    } else {
      ggplot2::ggsave(plot=phmb,
                      filename=paste0(imgdir,'GOplots_bubbleplot_',contrast.lbl,'_',set.lbl,'_',udlbl,"_",THYFUNGO,'_sn',THYSG.NODES,'_annot.png'),
                      width=16,height=8,units="in",dpi=600,
                      bg="white");
    };
    
    cat("# circle-plot...\n");
    tryCatch(phmc <- GOCircle(circ, nsub = min(nrow(circ),25), label.size = 2, lfc.col = lfcCOL, thygo=THYFUNGO),
             error = function(e) { skiptonext <<- TRUE });
    if (skiptonext) {
      skiptonext <- FALSE;
    } else {
      ggplot2::ggsave(plot=phmc,
                    filename=paste0(imgdir,'GOplots_circleplot_',contrast.lbl,'_',set.lbl,'_',udlbl,"_",THYFUNGO,'_sn',THYSG.NODES,'_annot.png'),
                    width=16,height=8,units="in",dpi=600,
                    bg="white");
    };
    
    cat("# GO heatmap counts...\n");
    if (cmax.rows > 0) {
      phmA <- GOHeat_fix(chord[,-ncol(chord)], nlfc = 0, fill.col = c("red", "yellow", "blue"), thygo=THYFUNGO);
      ggplot2::ggsave(plot=phmA,
                      filename=paste0(imgdir,'GOplots_heatmapCounts_',contrast.lbl,'_',set.lbl,'_',udlbl,"_",THYFUNGO,'_sn',THYSG.NODES,'_annot.png'),
                      height=max(4, 4+(16 * cmax.rows/100)),
                      width=max(8, 8+(14 * cmax.cols/100)),units="in",dpi=600,
                      bg="white");
      cat("# GO heatmap logFC...\n");
      phmB <- GOHeat_fix(chord, nlfc = 1, fill.col = c("red", "yellow", "blue"), thygo=THYFUNGO);
      ggplot2::ggsave(plot=phmB,
                      filename=paste0(imgdir,'GOplots_heatmapLogFC_',contrast.lbl,'_',set.lbl,'_',udlbl,"_",THYFUNGO,'_sn',THYSG.NODES,'_annot.png'),
                      height=max(4, 4+(16 * cmax.rows/100)),
                      width=max(8, 8+(14 * cmax.cols/100)),units="in",dpi=600,
                      bg="white");
    };
    cat("# GO pheatmap...\n");
    thychord <- t(chord[,c(1:ncol(chord)-1)]);
    cmax.rows <- min(nrow(thychord), 100);
    cmax.cols <- min(ncol(thychord), 100);
    if (cmax.rows > 0) {
      phm <- pheatmap::pheatmap(thychord,
                                main=paste(contrast.lbl,set.lbl,udlbl,THYFUNGO,collapse=" "),
                                fontsize_row = 8, fontsize_col = 7, legend = FALSE);
      png(file=paste0(imgdir,"GOplots_heatmaps_",contrast.lbl,'_',set.lbl,'_',udlbl,"_",THYFUNGO,'_sn',THYSG.NODES,'_annot.png'),
          res=600, height=max(4, 4+(10 * cmax.rows/100)),
          width=max(6, 4+(15 * cmax.cols/100)), unit="in", pointsize=10);
      grid::grid.newpage();
      grid::grid.draw(phm$gtable);
      dev.off();
    }
    cat("# PLOTS DONE...\n");

    # # Restore output to console
    # sink() 
    # sink(type="message")

    # return the GO dataframe
    res_top_GO
} # runGOsUPDN_ext

runGOsUPDN_extS <- function(udflg="N",
                            contrast.n, contrast.lbl, set.lbl,
                            THYFUNGO, table_exp, go_mapping,
                            annot_names, outpfx) {
    
    if (identical(udflg,"U")) {
       udlbl = "UP";
       lfcCOL <- c("#ff3030");
       UDlabels <- c("upregulated");
       DE_genes <- subset(table_exp,
                          P_adj <= MIN.PV & log_FC >= MIN.FC)$Gene_ID;
    } else {
       if (identical(udflg,"D")) {
         udlbl = "DOWN";
         UDlabels <- c("downregulated");
         lfcCOL <- c("#6495ed");
         DE_genes <- subset(table_exp,
                            P_adj <= MIN.PV & -log_FC >= MIN.FC)$Gene_ID;
       } else  {
         udlbl = "BOTH";
         lfcCOL <- c("#ff3030", "#6495ed");
                          # RED-UP BLUE-DOWN
         UDlabels <- c("downregulated", "upregulated");
         DE_genes <- subset(table_exp,
                            P_adj <= MIN.PV & abs(log_FC) >= MIN.FC)$Gene_ID;
       };
    };
    
    OUTPFX <- paste0(outpfx,"_",set.lbl,".",contrast.lbl,"_",udlbl,".GOA_",THYFUNGO);
   
    THYGOCOL <- ifelse(identical(THYFUNGO,"BP"),
                       "#a1c57f99",
                       ifelse(identical(THYFUNGO,"MF"),
                              "#b1caf699",
                              "#f69c9c99")); 

    cat("# Working on ", contrast.n," ",udlbl, ": ", contrast.lbl, " / ", set.lbl,"\n", sep="");

    gene_list <- factor(as.integer(annot_names %in% DE_genes));
    names(gene_list) <- annot_names;
    cat("# GENE LIST DONE...\n");
    
    cat("# Analysis of ",THYFUNGO,"...\n");
    top_GO_data <- new("topGOdata", ontology = THYFUNGO, allGenes = gene_list,       # all_genes,
                       nodeSize = 10, annot = annFUN.gene2GO, gene2GO = go_mapping); # affyLib = "drosophila2.db")
    cat("# ",THYFUNGO,": topGO DONE...\n");
    
    result_top_GO_elim    <- runTest(top_GO_data, algorithm = "elim",    statistic = "Fisher");
    result_top_GO_classic <- runTest(top_GO_data, algorithm = "classic", statistic = "Fisher");

    rv <- c(as.integer(geneData(result_top_GO_elim)),
            as.integer(geneData(result_top_GO_classic)));
    cat(do.call(sprintf,
                c(fmt=paste0('topGO %s %s %s signif.nodes %d\n',
                             'F.elim: Annot. %d  Signif. %d  NodeSize %d  SigTerms %d\n',
                             'F.clas: Annot. %d  Signif. %d  NodeSize %d  SigTerms %d\n'),
                udlbl, contrast.lbl, THYFUNGO, THYSG.NODES, as.list(rv))));
    cat("# ",THYFUNGO,": runTest DONE...\n");

    res_top_GO <- GenTable(top_GO_data,
                           Fisher.elim    = result_top_GO_elim,
                           Fisher.classic = result_top_GO_classic,
                           orderBy = "Fisher.elim", topNodes = 100, numChar=1000);
    res_top_GO$Category    <- rep(THYFUNGO, nrow(res_top_GO));
    res_top_GO$AnnotatedGenes   <- lapply(res_top_GO$GO.ID,
                                     function(x) as.character(unlist(genesInTerm(object = top_GO_data, whichGO = x))));
    res_top_GO$SignificantGenes <- lapply(res_top_GO$AnnotatedGenes, function(x) intersect(x,  DE_genes));
    res_top_GO$AnnotatedGenes   <- sapply(res_top_GO$AnnotatedGenes, paste, collapse=", ");
    res_top_GO$SignificantGenes <- sapply(res_top_GO$SignificantGenes, paste, collapse=", ");

    write.table(res_top_GO,
                file=paste0(OUTPFX,".topGOdata.extgenetable.sn",THYSG.NODES,".tbl"),
                sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE);
    cat("# ",THYFUNGO,": Extended GenTable DONE...\n");

    # Generate the GOplot plotting object
    cat("# GOplot object...\n");
    EC$goa  <- res_top_GO[ ,c("Category","GO.ID","Term","SignificantGenes","Fisher.elim") ];
    colnames(EC$goa)  <- c("Category","ID","Term","Genes","adj_pval");
    EC$glst  <- table_exp[ DE_genes, c("Gene_ID","log_FC","Avg_Expr","stat","P_value","P_adj","B") ];
    colnames(EC$glst) <- c("ID","logFC","AveExpr","t","P.Value","adj.P.Val","B");
    EC$genes <- table_exp[ DE_genes, c("Gene_ID","log_FC") ];
    colnames(EC$genes) <- c("ID","logFC");
    EC$process <- head(EC$goa$Term,25);
    circ <- circle_dat(EC$goa,EC$glst);
    circ$adj_pval <- as.numeric(circ$adj_pval); # fix needed 
    chord <- chord_dat(circ, EC$genes, EC$process);
    write.table(circ[ !is.na(circ$logFC), ],
                file=paste0(OUTPFX,".topGOdata.zscores.sn",THYSG.NODES,".tbl"),
                sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE);
    cat("# ",THYFUNGO,": zscores GenTable DONE...\n");
    # cat("# GOplot object2...\n", colnames(circ),"\n")
    circ$term <- str_trunc(circ$term, 50, side="right", ellipsis = "…"); # "\u2026" "..." "…"
    cmax.rows <- min(nrow(chord), 100);
    cmax.cols <- min(ncol(chord), 100);

    cat("# combo-plot...\n")
    skiptonext <- FALSE;
    tryCatch(phmb <- GOfoldchange(circ,
                                  title=paste0('FoldChange + Z-score for',contrast.lbl,' ',
                                               set.lbl,' ',udlbl,' ',THYFUNGO,' sn',THYSG.NODES),
                                  ofile=paste0(OUTPFX,'.GOplots_comboGOplot_FC+zSco_sn',THYSG.NODES,'_annot'),
                                  col = THYGOCOL, thygo=THYFUNGO,
                                  debugflg=TRUE, saveflg=TRUE),
             error = function(e) { skiptonext <<- TRUE });

    cat("# hist-plot...\n")
    skiptonext <- FALSE;
    tryCatch(phmb <- GOBar(circ,
                           title=paste0('Z-score barplot for',contrast.lbl,' ',
                                        set.lbl,' ',udlbl,' ',THYFUNGO,' sn',THYSG.NODES)) +
                     geom_hline(yintercept=c(-log10(0.05),-log10(0.0101)),
                                linetype=c('dashed','dotted'), color="red", size=0.25),
             error = function(e) { skiptonext <<- TRUE });
    if (skiptonext) {
      skiptonext <- FALSE;
    } else {
      ggplot2::ggsave(plot=phmb,
                      filename=paste0(OUTPFX,'.GOplots_histplot_sn',THYSG.NODES,'_annot.png'),
                      height=max(4, 4+(16 * cmax.rows/100)), width=16, units="in",dpi=600,
                      bg="white");
    };


    cat("# bubble-plot...\n");
    skiptonext <- FALSE;
    tryCatch(phmb <- GOBubble_fix(circ, labels = 2, # if labels=number then threshold is -log10(adjpval)>=number
                                  colour = c(THYGOCOL), table.col=TRUE),
             error = function(e) { skiptonext <<- TRUE });
    if (skiptonext) {
      skiptonext <- FALSE;
    } else {
      ggplot2::ggsave(plot=phmb,
                      filename=paste0(OUTPFX,'.GOplots_bubbleplot_sn',THYSG.NODES,'_annot.png'),
                      width=16,
                      height=max(6, 4+(16 * cmax.rows/100)),units="in",dpi=600,
                      bg="white");
    };
        
    cat("# circle-plot...\n");
    tryCatch(phmc <- GOCircle_fix(circ, nsub = min(nrow(circ),25), label.size = 2,
                                  lfc.col = lfcCOL, updnlabels = UDlabels,
                                  thygo = THYFUNGO, thygocol = THYGOCOL),
             error = function(e) { skiptonext <<- TRUE });
    if (skiptonext) {
      skiptonext <- FALSE;
    } else {
      ggplot2::ggsave(plot=phmc,
                    filename=paste0(OUTPFX,'.GOplots_circleplot_sn',THYSG.NODES,'_annot.png'),
                    width=16,height=8,units="in",dpi=600,
                    bg="white");
    };

    if (cmax.rows > 0) {
      cat("# GO heatmap counts...\n");
      phmA <- GOHeat_fix(chord[,-ncol(chord)], nlfc = 0,
                         fill.col = c("red", "yellow", "blue"), thygo=THYFUNGO);
      ggplot2::ggsave(plot=phmA,
                      filename=paste0(OUTPFX,'.GOplots_heatmapCounts_sn',THYSG.NODES,'_annot.png'),
                      width =max(4, 6+(12 * cmax.rows/100)),
                      height=max(4, 2+(12 * cmax.cols/100)),units="in",dpi=600,
                      bg="white");
      cat("# GO heatmap logFC...\n");
      phmB <- GOHeat_fix(chord, nlfc = 1,
                         fill.col = c("red", "yellow", "blue"), thygo=THYFUNGO);
      ggplot2::ggsave(plot=phmB,
                      filename=paste0(OUTPFX,'.GOplots_heatmapLogFC_sn',THYSG.NODES,'_annot.png'),
                      width =max(4, 6+(12 * cmax.rows/100)),
                      height=max(4, 2+(12 * cmax.cols/100)),units="in",dpi=600,
                      bg="white");
    };

    cat("# GO pheatmap...\n");
    thychord <- t(chord[,c(1:ncol(chord)-1)]);
    cmax.rows <- min(nrow(thychord), 100);
    cmax.cols <- min(ncol(thychord), 100);
    if (cmax.rows > 0) {
      phm <- pheatmap::pheatmap(thychord,
                                main=paste(contrast.lbl,set.lbl,THYFUNGO,collapse=" "),
                                fontsize_row = 8, fontsize_col = 7, legend = FALSE);
      png(file=paste0(OUTPFX,'.GOplots_heatmaps_sn',THYSG.NODES,'_annot.png'),
          height=max(4, 2+(12 * cmax.rows/100)),
          width =max(4, 6+(12 * cmax.cols/100)), unit="in", res=600, pointsize=10);
      grid::grid.newpage();
      grid::grid.draw(phm$gtable);
      # grid.text(THYFUNGO, x=0, y=0, hjust=0, vjust=0, size = rel(12));
      dev.off();
    }

    cat("# GO plotting significant nodes...\n");
    png(file=paste0(OUTPFX,'.GOplots_significantnodes_sn',THYSG.NODES,'_annot.png'),
        width=6,height=16,unit="in",pointsize=8,res=600);
    showSigOfNodes(top_GO_data, score(result_top_GO_elim),
                   firstSigNodes = THYSG.NODES, sigForAll = TRUE,
                   useInfo = 'all', putCL = 1);
    dev.off();
    pdf(file=paste0(OUTPFX,'.GOplots_significantnodes_sn',THYSG.NODES,'_annot.pdf'),
        width=6,height=16,pointsize=8,paper="special",onefile=TRUE);
    showSigOfNodes(top_GO_data, score(result_top_GO_elim),
                   firstSigNodes = THYSG.NODES, sigForAll = TRUE,
                   useInfo = 'all', putCL = 1);
    dev.off();

    cat("# PLOTS DONE...\n");

    # # Restore output to console
    # sink() 
    # sink(type="message")

    # return the GO dataframe
    res_top_GO
} # runGOsUPDN_extS


runGOsPrEn_extS <- function(contrast.lbl, set.lbl,
                            THYFUNGO, table_exp, go_mapping,
                            chosen_IDS, annot_names, outpfx) {
    lfcCOL <- c("#ff3030", "#6495ed")
               # RED  -> log(E/P)>0 more selected motives at enhancer
                                        # BLUE -> log(E/P)<0 more selected motives at promoter
    
    THYGOCOL <- ifelse(identical(THYFUNGO,"BP"),
                       "#a1c57f99",
                       ifelse(identical(THYFUNGO,"MF"),
                              "#b1caf699",
                              "#f69c9c99"));

    rownames(table_exp) <- table_exp$Gene_ID;
    OUTPFX <- paste0(outpfx,"_",set.lbl,".",contrast.lbl,".GOA_",THYFUNGO);
        
    # # avoiding undefined due to division by zero
    # table_exp$log_FC <- ifelse(table_exp$PRM.foxRatio == 0 | table_exp$ENH.foxRatio == 0,
    #                            NA,
    #                            log2(table_exp$ENH.foxRatio / table_exp$PRM.foxRatio));
    # # yet another approach
    table_exp$PRM.foxRatio[ table_exp$PRM.foxRatio == 0 ] <- 1e-3;
    table_exp$ENH.foxRatio[ table_exp$ENH.foxRatio == 0 ] <- 1e-3;
    table_exp$log_FC <- log2(table_exp$ENH.foxRatio / table_exp$PRM.foxRatio);
    
    CH_genes <- subset(table_exp, Gene_ID %in% chosen_IDS)$Gene_ID;

    cat("# Working on : ", contrast.lbl, " / ", set.lbl,"\n", sep="");

    gene_list <- factor(as.integer(annot_names %in% CH_genes));
    names(gene_list) <- annot_names;
    cat("# GENE LIST DONE...\n");

    cat("# Analysis of ",THYFUNGO,"...\n")
    top_GO_data <- new("topGOdata", ontology = THYFUNGO, allGenes = gene_list,       # all_genes,
                       nodeSize = 10, annot = annFUN.gene2GO, gene2GO = go_mapping); # affyLib = "drosophila2.db")
    cat("# ",THYFUNGO,": topGO DONE...\n");
    
    result_top_GO_elim    <- runTest(top_GO_data, algorithm = "elim",    statistic = "Fisher");
    result_top_GO_classic <- runTest(top_GO_data, algorithm = "classic", statistic = "Fisher");
    
    rv <- c(as.integer(geneData(result_top_GO_elim)),
            as.integer(geneData(result_top_GO_classic)))
    cat(do.call(sprintf,
                c(fmt=paste0('topGO %s %s signif.nodes %d\n',
                             'F.elim: Annot. %d  Signif. %d  NodeSize %d  SigTerms %d\n',
                             'F.clas: Annot. %d  Signif. %d  NodeSize %d  SigTerms %d\n'),
                contrast.lbl, THYFUNGO, THYSG.NODES, as.list(rv))));
    cat("# ",THYFUNGO,": runTest DONE...\n");

    res_top_GO <- GenTable(top_GO_data,
                           Fisher.elim    = result_top_GO_elim,
                           Fisher.classic = result_top_GO_classic,
                           orderBy = "Fisher.elim", topNodes = 100, numChar=1000);
    res_top_GO$Category    <- rep(THYFUNGO, nrow(res_top_GO));
    res_top_GO$AnnotatedGenes   <- lapply(res_top_GO$GO.ID,
                                     function(x) as.character(unlist(genesInTerm(object = top_GO_data, whichGO = x))));
    res_top_GO$SignificantGenes <- lapply(res_top_GO$AnnotatedGenes, function(x) intersect(x,  CH_genes));
    res_top_GO$AnnotatedGenes   <- sapply(res_top_GO$AnnotatedGenes, paste, collapse=", ");
    res_top_GO$SignificantGenes <- sapply(res_top_GO$SignificantGenes, paste, collapse=", ");

    write.table(res_top_GO,
                file=paste0(OUTPFX,".topGOdata.extgenetable.sn",THYSG.NODES,".tbl"),
                sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE);
    cat("# ",THYFUNGO,": Extended GenTable DONE...\n");

    # Generate the GOplot plotting object
    cat("# GOplot object...\n");
    EC$goa  <- res_top_GO[ ,c("Category","GO.ID","Term","SignificantGenes","Fisher.elim") ];
    colnames(EC$goa)  <- c("Category","ID","Term","Genes","adj_pval");
    EC$glst  <- table_exp[ CH_genes, c("Gene_ID","log_FC","PRM.foxRatio","ENH.foxRatio") ];
    colnames(EC$glst) <- c("ID","logFC","PRM.foxRatio","ENH.foxRatio");
    EC$genes <- table_exp[ CH_genes, c("Gene_ID","log_FC") ];
    colnames(EC$genes) <- c("ID","logFC");
    EC$process <- head(EC$goa$Term,25);
    circ <- circle_dat(EC$goa,EC$glst);
    circ$adj_pval <- as.numeric(circ$adj_pval); # fix needed 
    chord <- chord_dat(circ, EC$genes, EC$process);
    write.table(circ[ !is.na(circ$logFC), ],
                file=paste0(OUTPFX,".topGOdata.zscores.sn",THYSG.NODES,".tbl"),
                sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE);
    cat("# ",THYFUNGO,": zscores GenTable DONE...\n");
    
    # cat("# GOplot object2...\n", colnames(circ),"\n")
    circ$term <- str_trunc(circ$term, 50, side="right", ellipsis = "…"); # "\u2026" "..." "…"
    cmax.rows <- min(nrow(chord), 100);
    cmax.cols <- min(ncol(chord), 100);
    
    cat("# combo-plot...\n");
    skiptonext <- FALSE;
    tryCatch(phmb <- GOfoldchange(circ,
                                  title=paste0('FoldChange + Z-score for',contrast.lbl,' ',
                                               set.lbl,' ',THYFUNGO,' sn',THYSG.NODES),
                                  ofile=paste0(OUTPFX,'.GOplots_comboGOplot_FC+zSco_sn',THYSG.NODES,'_annot'),
                                  oextns = "png",
                                  col = THYGOCOL, thygo=THYFUNGO,
                                  fclim = c(5,-5),
                                  updnlabel = parse(text= expression(log[2]("FoxG"["enhancer"]/"FoxG"["promoter"]))),
                                  updnlabels = c(parse(text= expression("FoxG"["enh."]>="FoxG"["prom."])),
                                                 parse(text= expression("FoxG"["enh."]<"FoxG"["prom."]))),
                                  zscolabs=c("promoter","","enhancer"),
                                  debugflg=TRUE, saveflg=TRUE),
             error = function(e) { skiptonext <<- TRUE });

    cat("# hist-plot...\n");
    skiptonext <- FALSE;
    tryCatch(phmb <- GOBar_fix(circ,
                               title=paste0('Z-score barplot for',contrast.lbl,' ',
                                            set.lbl,' ',THYFUNGO,' sn',THYSG.NODES),
                               zscolabs=c("promoter","","enhancer")) +
                         geom_hline(yintercept=c(-log10(0.05),-log10(0.0101)),
                                    linetype=c('dashed','dotted'), color="red", size=0.25),
             error = function(e) { skiptonext <<- TRUE });
    if (skiptonext) {
      skiptonext <- FALSE;
    } else {
      ggplot2::ggsave(plot=phmb,
                      filename=paste0(OUTPFX,'.GOplots_histplot_sn',THYSG.NODES,'_annot.png'),
                      height=max(4, 4+(16 * cmax.rows/100)), width=16, units="in",dpi=600,
                      bg="white");
    };

    cat("# bubble-plot...\n");
    skiptonext <- FALSE;
    tryCatch(phmb <- GOBubble_fix(circ, labels = 2, # if labels=number then threshold is -log10(adjpval)>=number
                                  colour = c(THYGOCOL), table.col=TRUE),
             error = function(e) { skiptonext <<- TRUE });
    if (skiptonext) {
      skiptonext <- FALSE;
    } else {
      ggplot2::ggsave(plot=phmb,
                      filename=paste0(OUTPFX,'.GOplots_bubbleplot_sn',THYSG.NODES,'_annot.png'),
                      width=16,
                      height=max(6, 4+(16 * cmax.rows/100)),units="in",dpi=600,
                      bg="white");
    };
        
    cat("# circle-plot...\n");
    tryCatch(phmc <- GOCircle_fix(circ, nsub = min(nrow(circ),25), label.size = 2,
                                  lfc.col = lfcCOL, thygo = THYFUNGO, thygocol = THYGOCOL,
                                  updnlabel = parse(text= expression(log[2]("FoxG"["enh."]/"FoxG"["prom."]))),
                                  updnlabels=c(parse(text= expression("FoxG"["enh."]<"FoxG"["prom."])),
                                               parse(text= expression("FoxG"["enh."]>="FoxG"["prom."]))),
                                  zscolabs=c("promoter","","enhancer")),
             error = function(e) { skiptonext <<- TRUE });
    if (skiptonext) {
      skiptonext <- FALSE;
    } else {
      ggplot2::ggsave(plot=phmc,
                    filename=paste0(OUTPFX,'.GOplots_circleplot_sn',THYSG.NODES,'_annot.png'),
                    width=16,height=8,units="in",dpi=600,
                    bg="white");
    };
    
    if (cmax.rows > 0) {
      cat("# GO heatmap counts...\n");
      phmA <- GOHeat_fix(chord[,-ncol(chord)], nlfc = 0,
                         fill.col = c("red", "yellow", "blue"), thygo=THYFUNGO);
      ggplot2::ggsave(plot=phmA,
                      filename=paste0(OUTPFX,'.GOplots_heatmapCounts_sn',THYSG.NODES,'_annot.png'),
                      width =max(4, 6+(12 * cmax.rows/100)),
                      height=max(4, 2+(12 * cmax.cols/100)),units="in",dpi=600,
                      bg="white");
      cat("# GO heatmap logFC...\n");
      phmB <- GOHeat_fix(chord, nlfc = 1,
                         fill.col = c("red", "yellow", "blue"), thygo=THYFUNGO);
      ggplot2::ggsave(plot=phmB,
                      filename=paste0(OUTPFX,'.GOplots_heatmapLogFC_sn',THYSG.NODES,'_annot.png'),
                      width =max(4, 6+(12 * cmax.rows/100)),
                      height=max(4, 2+(12 * cmax.cols/100)),units="in",dpi=600,
                      bg="white");
    };

    cat("# GO pheatmap...\n");
    thychord <- t(chord[,c(1:ncol(chord)-1)]);
    cmax.rows <- min(nrow(thychord), 100);
    cmax.cols <- min(ncol(thychord), 100);
    if (cmax.rows > 0) {
      phm <- pheatmap::pheatmap(thychord,
                                main=paste(contrast.lbl,set.lbl,THYFUNGO,collapse=" "),
                                fontsize_row = 8, fontsize_col = 7, legend = FALSE);
      png(file=paste0(OUTPFX,'.GOplots_heatmaps_sn',THYSG.NODES,'_annot.png'),
          height=max(4, 2+(12 * cmax.rows/100)),
          width =max(4, 6+(12 * cmax.cols/100)), unit="in", res=600, pointsize=10);
      grid::grid.newpage();
      grid::grid.draw(phm$gtable);
      # grid.text(THYFUNGO, x=0, y=0, hjust=0, vjust=0, size = rel(12));
      dev.off();
    }

    cat("# GO plotting significant nodes...\n");
    png(file=paste0(OUTPFX,'.GOplots_significantnodes_sn',THYSG.NODES,'_annot.png'),
        width=6,height=16,unit="in",pointsize=8,res=600);
    showSigOfNodes(top_GO_data, score(result_top_GO_elim),
                   firstSigNodes = THYSG.NODES, sigForAll = TRUE,
                   useInfo = 'all', putCL = 1);
    dev.off();
    pdf(file=paste0(OUTPFX,'.GOplots_significantnodes_sn',THYSG.NODES,'_annot.pdf'),
        width=6,height=16,pointsize=8,paper="special",onefile=TRUE);
    showSigOfNodes(top_GO_data, score(result_top_GO_elim),
                   firstSigNodes = THYSG.NODES, sigForAll = TRUE,
                   useInfo = 'all', putCL = 1);
    dev.off();

    cat("# PLOTS DONE...\n");

} # runGOsPrEn_extS
