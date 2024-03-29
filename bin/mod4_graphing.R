#!/usr/bin/env Rscript
# Import libraries
# Create figures from data
suppressPackageStartupMessages(require("optparse"))                                
suppressPackageStartupMessages(require("stats"))

suppressPackageStartupMessages(require(RColorBrewer))
suppressPackageStartupMessages(require(gridExtra))
suppressPackageStartupMessages(require(grid))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(cowplot))
suppressPackageStartupMessages(require(plyr))

# Sourcing useful functions
source('/tld/fxns/def_plotting.R')

# Import options
option_list <- list(make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
                                help = "Print extra output [default]"),
                    make_option(c("-q", "--quietly"), action = "store_false",
                                dest = "verbose", help = "Print little output"),
                    make_option(c("-r", "--res_path"), type = "character",
                                help = "path to results from previous scripts",
                                metavar = "PATH"),
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
    paste("\tPrefix:", opt$prefix, collapse = ""), "\n",
    "\n"
)

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

read_count.f <- paste(c(result.path, "/", opt$prefix, ".read_count.Rda"),
                      collapse = "")
read_count <- readRDS(file = read_count.f)

read_count.Abam.f <- paste(c(result.path, "/", opt$prefix,
                             ".read_count.Abam.Rda"),
                           collapse = "")
read_count.Abam <- readRDS(file = read_count.Abam.f)

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

read_count.f <- paste(c(result.path, "/", opt$prefix, ".read_count.Rda"), 
                      collapse = "")
read_count <- readRDS(file = read_count.f)

read_count.Abam.f <- paste(c(result.path, "/", opt$prefix, 
                             ".read_count.Abam.Rda"), 
                           collapse = "")
read_count.Abam <- readRDS(file = read_count.Abam.f)

# Create table for data on chromosomes
bam_read_count <- ddply(end.bam.df,.(s.name, r.type,
                                     seqnames, chr.end), summarize,
                        reads.per.chr=length(unique(r.name)),
                        tel.length=mean(tel.length))

# Perform t-test
t1 <- result.df[result.df['s.name'] == unique(result.df$s.name)[1] & 
                  result.df['r.type'] == "nl" & 
                  result.df['norm'] == "normalized",]$tel.length
t2 <- result.df[result.df['s.name'] == unique(result.df$s.name)[2] & 
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

hist_tel_rep.f <- paste(c(result.path, "/", opt$prefix, 
                          ".telomere_rep_dist.png"), collapse = "")
def.ggsave(hist_tel_rep.f, plot = per.tel.hist)

# Bar graph for read count after telomere processing
plot_bam_read_count$s.id <- as.factor(plot_bam_read_count$seqnames)
title <- "A. Number of Reads by End"
y.lab <- "Telomere Read Count"
x.lab <- "Chromosome"
color.lab <- "Sample Type"

bar_pl <- def.bar(plot_bam_read_count,
                  title, x.lab, y.lab,
                  color.lab) +
  def_th + theme(legend.key.size = unit(5, "lines"))


bar_pl.f <- paste(c(result.path, "/", opt$prefix, ".telomere_read_count.png"), collapse = "")
def.ggsave(bar_pl.f, plot = bar_pl, height = 10)

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
  def_th + theme(legend.key.size = unit(4, "lines"), legend.position = "top",
                 legend.text=element_text(size = 38), 
                 legend.spacing.x = unit(1.0, 'cm')) + 
  cl_pal

hist.pl.f <- paste(c(result.path, "/", opt$prefix, ".histo.png"), 
                   collapse = "")
def.ggsave(hist.pl.f, plot = hist_pl)

# Making Figure 1
bar_pl_f <- bar_pl + theme(legend.position = "none")
hist_pl_f <- hist_pl + theme(legend.position = "none")
leg <- get_legend(hist_pl)
title <- textGrob("Telomere Reads and Lengths", 
                  #x= 0.1,
                  gp = gpar(fontsize = 80, 
                            fontfamily = "HersheySans", 
                            fontface = "bold"))
lay <- rbind(c(1,1), c(3,4), c(2,2))

fig1 <- grid.arrange(title, leg, bar_pl_f, hist_pl_f, nrow = 3, ncol = 2,
                     layout_matrix = lay,
                     widths = c(5, 5), heights = c(1.5, 7, 1.5))

fig1.f <- paste(c(result.path, "/", opt$prefix, ".combined_plot.png"), collapse = "")
def.ggsave(fig1.f, plot = fig1)

# Density Plot, just Sample Type 200
scat_title <- "Normalized Density Plot"
y.lab <- "Telomere Lengths (kb)"
x.lab <- "Read Lengths (kb)"
color.lab <- "Sample Type"
shape.lab <- "Sample Type"

scat_pl <- def.scatter(result.df[result.df$norm=="normalized", ],
                       scat_title, x.lab, y.lab,
                       color.lab, shape.lab) + 
  theme(panel.spacing = unit(3, "lines")) +
  facet_wrap(r.type ~ s.name, 
             ncol = 2) 

scat_pl.f <- paste(c(result.path, "/", opt$prefix, 
                     ".density_plot.png"), collapse = "")
def.ggsave(scat_pl.f, plot = scat_pl, height = 12)