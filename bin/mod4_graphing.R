#!/usr/bin/env Rscript
# Import libraries
# Create figures from data
suppressPackageStartupMessages(library("optparse"))                                
suppressPackageStartupMessages(library("stats"))

suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(plyr))

# Sourcing useful functions
source("/home/jaker/pfalci/201908_telomere_lengths/github/tld/fxns/def_plotting.R")

# Import options
option_list <- list(make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
                                help = "Print extra output [default]"),
                    make_option(c("-q", "--quietly"), action = "store_false",
                                dest = "verbose", help = "Print little output"),
                    make_option(c("-r", "--res_path"), type = "character",
                                help = "path to results from previous scripts",
                                metavar = "PATH"),
                    make_option(c("-o", "--out_path"), type = "character",
                                help = "full path to folder where results should be stored"),
                    make_option(c("-p", "--prefix"), type = "character",
                                help = "prefix for naming output files, must match from module 1")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(option_list=option_list))

# print some progress messages to stderr if \"quietly\" wasn't requested
if ( opt$verbose ) {
  write(paste("\nCreating graphs from telomere ends:",
              Sys.time(), collapse = ""), stderr())
}

# Confirm all inputs
cat("\n These are the options you submitted: \n",
    paste("\tResults path:", opt$res_path, collapse = ""), "\n",
    paste("\tOut path:", opt$out_path, collapse = ""), "\n",
    paste("\tPrefix:", opt$prefix, collapse = ""), "\n",
    paste("Are these correct? (y/n)", collapse = ""), "\n"
)

# Exit if user does not confirm variables
opt_check <- readLines(con = "stdin", n = 1)
if ( opt_check == "y") {
  cat("\nOptions have been confirmed, continuing... \n\n")
} else {
  cat("\nOptions have not been confirmed, exiting. \n\n")
  q(save = "no", status = 1, runLast = FALSE)
}



# Import csv files 
result.path <- opt$res_path

df.f <- paste(c(result.path, "/", opt$prefix, ".df.Rda"), collapse = "")
df <- readRDS(file = df.f)

result.df.f <- paste(c(result.path, "/", opt$prefix, ".result.df.Rda"),
                     collapse = "")
result.df <- readRDS(file = result.df.f)

per.tel.stats.f <- paste(c(result.path, "/", opt$prefix, ".per.tel.stats.Rda"),
                         collapse = "")
per.tel.stats <- readRDS(file = per.tel.stats.f)

result.tel.stats.f <- paste(c(result.path, "/", opt$prefix,
                              ".result.tel.stats.Rda"), collapse = "")
result.tel.stats <- readRDS(file = result.tel.stats.f)

end.bam.df.f <- paste(c(result.path, "/", opt$prefix, ".end.bam.csv"),
                      collapse = "")
end.bam.df <- read.csv(file = end.bam.df.f, stringsAsFactors = FALSE)
#end.bam.df$s.name[end.bam.df$s.name=="3d7"] <- "wt"

read_count.f <- paste(c(result.path, "/", opt$prefix, ".read_count.Rda"),
                      collapse = "")
read_count <- readRDS(file = read_count.f)

read_count.Abam.f <- paste(c(result.path, "/", opt$prefix,
                             ".read_count.Abam.Rda"),
                           collapse = "")
read_count.Abam <- readRDS(file = read_count.Abam.f)

df.f <- paste(c(result.path, "/", prefix, ".df.Rda"), collapse = "")
df <- readRDS(file = df.f)

result.df.f <- paste(c(result.path, "/", prefix, ".result.df.Rda"), 
                     collapse = "")
result.df <- readRDS(file = result.df.f)

per.tel.stats.f <- paste(c(result.path, "/", prefix, ".per.tel.stats.Rda"), 
                         collapse = "")
per.tel.stats <- readRDS(file = per.tel.stats.f)

result.tel.stats.f <- paste(c(result.path, "/", prefix, 
                              ".result.tel.stats.Rda"), collapse = "") 
result.tel.stats <- readRDS(file = result.tel.stats.f)

end.bam.df.f <- paste(c(result.path, "/", prefix, ".end.bam.csv"), 
                      collapse = "")
end.bam.df <- read.csv(file = end.bam.df.f, stringsAsFactors = FALSE)

read_count.f <- paste(c(result.path, "/", prefix, ".read_count.Rda"), 
                      collapse = "")
read_count <- readRDS(file = read_count.f)

read_count.Abam.f <- paste(c(result.path, "/", prefix, 
                             ".read_count.Abam.Rda"), 
                           collapse = "")
read_count.Abam <- readRDS(file = read_count.Abam.f)


# Create table for data on chromosomes
bam_read_count <- ddply(end.bam.df,.(s.name, r.type,
                                     seqnames, chr.end), summarize,
                        reads.per.chr=length(unique(r.name)),
                        tel.length=mean(tel.length))

wo_12_13_bam_read_count <- bam_read_count[!(bam_read_count['chr.end']=="3'" & 
                                              bam_read_count['s.name']=="irr" & 
                                              (bam_read_count['seqnames']==12 | 
                                                 bam_read_count['seqnames']==13)),]

# Perform t-test
t1 <- result.df[result.df['s.name'] == "wt" & 
                  result.df['r.type'] == "nl" & 
                  result.df['norm'] == "normalized",]$tel.length
t2 <- result.df[result.df['s.name'] == "irr" & 
                  result.df['r.type'] == "nl" & 
                  result.df['norm'] == "normalized",]$tel.length
t_test.res <- t.test(t1, t2, var.equal = FALSE, alternative = "two.sided")

# Create figures on read counts
plot_bam_read_count <- bam_read_count[bam_read_count$r.type=="nl",]
attach(plot_bam_read_count)
plot_bam_read_count <- plot_bam_read_count[order(chr.end, decreasing = TRUE), ]
detach(plot_bam_read_count)

# Percent telomeric repeats in each window histogram
df_hist <- df
per.tel.stats.h <- per.tel.stats
perHistogram <- ggplot(df_hist, aes(x = per.tel.repeat, 
                                    color = s.name, 
                                    fill = s.name)) +
  geom_histogram(aes(y=..density..), 
                 position = "dodge",
                 alpha = 0.10,
                 bins = 33) +
  geom_vline(data = per.tel.stats.h,
             aes(xintercept = mode.1,
                 color = s.name), 
             linetype = "dashed",
             alpha = 1,
             size = 2) +
  geom_vline(data = per.tel.stats.h,
             aes(xintercept = mode.2,
                 color = s.name),
             linetype = "dashed",
             alpha = 1,
             size = 2) +
  geom_density(alpha=0.4)

# Saving Histogram of percent telomere distributions
hist.title <- "Supplementary 1: Telomere Repeat Distributions"
hist.x <- "Telomere Repeat Percent"
hist.y <- "Density"
hist.color <- "Sample Type"
hist.fill <- "Sample Type"

per.tel.hist <- perHistogram + 
  hist_lb(hist.title, hist.x, hist.y, hist.color, hist.fill) +
  def_th + theme(legend.key.size = unit(5, "lines")) + cl_pal

hist_tel_rep.f <- paste(c(result.path, "/", prefix, 
                          ".telomere_rep_dist.png"), collapse = "")
def.ggsave(hist_tel_rep.f, plot = per.tel.hist)

# Temp fix for modes
result.tel.stats[result.tel.stats$norm=="normalized" & result.tel.stats$s.name=="wt",]$mode.1 <-
  result.tel.stats[result.tel.stats$norm=="normalized" & result.tel.stats$s.name=="wt",]$mode.2

# Bar graph for read count after telomere processing
plot_bam_read_count$seqnames <- as.numeric(gsub("chr_", "", 
                                                plot_bam_read_count$seqnames))
plot_bam_read_count$s.id <- plot_bam_read_count$seqnames
title <- "A. Number of Reads by End"
y.lab <- "Telomere Read Count"
x.lab <- "Chromosome"
color.lab <- "Sample Type"

# Setup annotations
ann_text_9 <- data.frame(s.id = 9, reads.per.chr = 1300, s.name = "irr", 
                       lab = "*",
                       chr.end=factor("5'", levels = c("5'", "3'")))

ann_text_12 <- data.frame(s.id = 12, reads.per.chr = 600, s.name = "irr", 
                       lab = "*",
                       chr.end=factor("3'", levels = c("5'", "3'")))

ann_text_13 <- data.frame(s.id = 13, reads.per.chr = 600, s.name = "irr", 
                       lab = "*",
                       chr.end=factor("3'", levels = c("5'", "3'")))

bar_pl <- def.bar(plot_bam_read_count,
                  title, x.lab, y.lab,
                  color.lab) +
  geom_text(data = ann_text_9, label = "*  ", size = 15) +
  geom_text(data = ann_text_12, label = "*  ", size = 15) +
  geom_text(data = ann_text_13, label = "*  ", size = 15)
  
  def_th + theme(legend.key.size = unit(5, "lines")) + cl_pal


bar_pl.f <- paste(c(result.path, "/", prefix, ".telomere_read_count.png"), collapse = "")
def.ggsave(bar_pl.f, plot = bar_pl, height = 8)

# Create Histograms
h.title <- "B. Telomere Length Distribution"
h.x <- "Telomere Length (kb)"
h.y <- "Density (f/kb)"
h.color <- "Sample Type"
h.fill <- "Sample Type"

hist <- ggplot(result.df[result.df$norm=="normalized",], 
               aes(x = tel.length, color = s.name, fill = s.name)) +
  geom_histogram(aes(y = ..density..),
                 binwidth = 0.10, 
                 position = "identity",
                 alpha = 0.10) +
  geom_vline(data = result.tel.stats[result.tel.stats$norm=="normalized",],
             aes(xintercept = result.tel.stats[result.tel.stats$norm=="normalized",]$mode.1,
                 color = s.name),
             linetype = "dashed",
             alpha = 1,
             size = 2) +
  geom_density(alpha = 0.40)

# Putting histogram together
hist_pl <- hist + hist_lb(h.title, h.x, 
                          h.y, h.color, 
                          h.fill) + 
  def_th + theme(legend.key.size = unit(3, "lines"), legend.position = "top",
                 legend.text=element_text(size = 36)) + 
  cl_pal
hist.pl.f <- paste(c(result.path, "/", prefix, ".histo_60-60_200SW.png"), 
                   collapse = "")
def.ggsave(hist.pl.f, plot = hist_pl)

# Making Figure 1
bar_pl_f <- bar_pl + theme(legend.position = "none")
hist_pl_f <- hist_pl + theme(legend.position = "none")
leg <- get_legend(hist_pl)
title <- textGrob("Figure 1", x = 0.07,
                  gp = gpar(fontsize = 80, 
                            fontfamily = "Arial", 
                            fontface = "plain"))
lay <- rbind(c(1,1), c(3,4), c(2,2))

fig1 <- grid.arrange(title, leg, bar_pl_f, hist_pl_f, nrow = 3, ncol = 2,
                     layout_matrix = lay,
                     widths = c(5, 5), heights = c(1.5, 7, 1.5))

fig1.f <- paste(c(result.path, "/", prefix, ".Figure_1.png"), collapse = "")
def.ggsave(fig1.f, plot = fig1)

# Density Plots, AllSW and ccs
scat_title <- "Normalized Density Plot--60/60"
y.lab <- "Telomere Lengths (kb)"
x.lab <- "Read Lengths (kb)"
color.lab <- "s.win"
shape.lab <- "s.win"

scat_pl <- def.scatter(result.df[result.df$norm=="normalized", ],
                       scat_title, x.lab, y.lab,
                       color.lab, shape.lab) +
  facet_wrap(s.win ~ r.type ~ s.name) 

scat_pl.f <- paste(c(result.path, "/", prefix, ".60_60_DensityPlot.png"), 
                   collapse = "")
def.ggsave(scat_pl.f, plot = scat_pl, height = 16)

# Density Plot, just s.win 200
scat_title <- "Normalized Density Plot, 200SW--60/60"
y.lab <- "Telomere Lengths (kb)"
x.lab <- "Read Lengths (kb)"
color.lab <- "s.win"
shape.lab <- "s.win"

scat_pl <- def.scatter(result.df[result.df$norm=="normalized", ],
                       scat_title, x.lab, y.lab,
                       color.lab, shape.lab) + 
  theme(panel.spacing = unit(3, "lines")) +
  facet_wrap(r.type ~ s.name, 
             ncol = 2) 

scat_pl.f <- paste(c(result.path, "/", prefix, 
                     ".60_60_DensityPlot.200sw.png"), collapse = "")
def.ggsave(scat_pl.f, plot = scat_pl, height = 12)