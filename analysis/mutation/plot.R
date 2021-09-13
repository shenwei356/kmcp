#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description = "", formatter_class = "argparse.RawTextHelpFormatter")
parser$add_argument("infile", type = "character", help = "usage: ./plot.R all.se150.fasta.blastn.filter.tsv.gz.stats.gz")
parser$add_argument("outfile", type = "character", default = "figure.jpg", help = "figure.jpg")

args <- parser$parse_args()

library(tidyverse)
library(dplyr)
library(scales)
library(ggthemes)
library(ggrepel)


# df <- read_delim("all.se150.fasta.blastn.filter.tsv.gz.stats.gz", "\t")
#  df <-   read_delim("all.se150.fasta.blastn.filter.tsv.gz.stats.gz.test", "\t")

df <- read_delim(args$infile, "\t")

every_nth <- function(x,
                      nth,
                      empty = TRUE,
                      inverse = FALSE) {
  if (!inverse) {
    if (empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if (empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}


custom_breaks_x <- seq(75, 100, 1)
custom_breaks_y <- seq(10, 100, 5)

p <-
  ggplot(df,
         aes(
           x = pident,
           y = qcov * 100,
           color = ifelse(acov > 1, 100, acov * 100)
           # color = mismatches
           # color = gaps
         )) +
  
  geom_point(size = 0.8, alpha = 0.5) +
  geom_vline(xintercept = 95,
             colour = "black",
             linetype = 2) +
  geom_vline(xintercept = 98,
             colour = "black",
             linetype = 2) +
  geom_hline(yintercept = 60,
             colour = "black",
             linetype = 2) +
  geom_hline(yintercept = 50,
             colour = "black",
             linetype = 2) +
  
  xlab('alignment identity%') +
  ylab('k-mer matches%') +
  labs(color = "len(alignment)/len(query)%") +
  # labs(color = "mismatches") +
  # labs(color = "gaps") + 
  
  # scale_colour_gradient_tableau("Orange")  +
  scale_colour_gradient2_tableau("Orange-Blue Diverging")  +
  scale_x_continuous(expand = c(0, 0),
                     breaks = custom_breaks_x,
                     labels =  every_nth(custom_breaks_x, 5, inverse = TRUE)) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = custom_breaks_y,
                     labels =  every_nth(custom_breaks_y, 2, inverse = TRUE)) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "grey30", size = 0.8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # plot.margin = unit(c(0.1,0.4,0.1,0.1),"cm"),
    axis.ticks.y = element_line(size = 0.6),
    axis.ticks.x = element_line(size = 0.6),
    
    strip.background = element_rect(
      colour = "white",
      fill = "grey95",
      size = 0.2
    ),
    strip.text.y = element_text(size = 13),
    
    legend.text = element_text(size = 13),
    legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm"),
    
    text = element_text(size = 13, family = "arial"),
    plot.title = element_text(size = 15)
  )

ggsave(p,
       file = args$outfile,
       width = 7,
       height = 5)