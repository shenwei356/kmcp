#!/usr/bin/env Rscript
library(argparse)

parser <-
  ArgumentParser(description = "", formatter_class = "argparse.RawTextHelpFormatter")
parser$add_argument("infile", type = "character", help = "usage: ./plot.R all.se150.fasta.blastn.filter.tsv.gz.stats.gz")
parser$add_argument("outfile",
                    type = "character",
                    default = "figure.jpg",
                    help = "figure.jpg")

args <- parser$parse_args()

library(tidyverse)
library(dplyr)
library(scales)
library(ggthemes)
library(ggrepel)
library(cowplot)


# df <- read_delim("all.se150.fasta.blastn.filter.tsv.gz.stats.uniq.gz", "\t")

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


theme1 <- theme_bw() +
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

## ------------------------------------------------------------------

custom_breaks_x0 <- seq(75, 100, 1)
custom_breaks_y0 <- seq(10, 100, 5)

addplot <- function(p) {
  p +
    xlab('alignment identity(%)') +
    ylab('k-mer similarity(%)') +
    geom_point(size = 1.2, alpha = 0.7) +
    scale_colour_gradient2_tableau("Orange-Blue Diverging")  +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = custom_breaks_x0,
      labels =  every_nth(custom_breaks_x0, 5, inverse = TRUE)
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = custom_breaks_y0,
      labels =  every_nth(custom_breaks_y0, 2, inverse = TRUE)
    ) +
    theme1
}

p1 <-
  ggplot(df,
         aes(
           x = pident,
           y = qcov * 100,
           color = ifelse(acov > 1, 100, acov * 100)
         )) +
  labs(color = "align. region(%)")

p2 <-
  ggplot(df,
         aes(x = pident,
             y = qcov * 100,
             color = mismatches)) +
  labs(color = "mismatches")

p3 <-
  ggplot(df,
         aes(x = pident,
             y = qcov * 100,
             color = gaps)) +
  labs(color = "gaps")

## ------------------------------------------------------------------

df2 <-
  df %>% filter(pident > 85 &
                  gaps <= 5 & acov >= 0.95) %>% select(-qseq,-sseq)

p4 <-
  ggplot(df2,
         aes(x = pident,
             y = qcov * 100,
             color = mismatches)) +
  labs(color = "mismatches") +
  geom_smooth(formula = y ~ poly(x, 3),
              size = 1,
              color = "#009E73")

## ------------------------------------------------------------------

# p, fpr of single bloom filter.
# k, theshold of query coverage.
# l, number of k-mers.
FPR <- function(p, k, l) {
  exp(-l * (k - p) * (k - p) / 2 / (1 - p))
}

model <- lm(pident ~ poly(qcov, 3, raw = TRUE), data = df2)

df3 <- data.frame(qcov = seq(0.2, 1, 0.01))
df3$pident <- predict(model, df3)

tmp <- c()
fpr150 <- c()
fpr100 <- c()
for (i in seq(1, length(df3$qcov))) {
  t_qcov <- df3$qcov[i]
  t_pident <- df3$pident[i]
  tmp <- c(tmp,
           dim(df2 %>% filter(qcov >= t_qcov &
                                pident > t_pident))[1] /
             dim(df2 %>% filter(pident > t_pident))[1])
  fpr150 <- c(fpr150, FPR(0.3, t_qcov, 130))
  fpr100 <- c(fpr100, FPR(0.3, t_qcov, 80))
  
}
df3$recall <- tmp
df3$SE150 <- fpr150
df3$SE100 <- fpr100

df3$recall <- df3$recall * 100
df3$qcov <- df3$qcov * 100

s150 <- -log10(min(df3$SE150)) / 100

custom_breaks_x <- seq(20, 100, 5)
custom_breaks_y <- seq(0, 100, 5)
custom_breaks_y2 <- seq(0,-log10(min(df3$SE150)), 1)

df3 <- df3 %>% gather(-qcov, -pident, -recall, key = "rlen", value = "fpr")
df3$rlen <- recode(df3$rlen, SE150="150 bp", SE100="100 bp")

p5 <- ggplot(df3, aes(qcov, recall)) +
  xlab("k-mer similarity threshold(%)") +
  # geom_point(size = 1.2, color="grey30") +
  geom_smooth(color = "#0072B2", size = 0.9) +
  geom_line(aes(y = -log10(fpr) / s150, color = rlen), size = 0.9) +
  scale_color_manual(values = c("#E69F00", "#D55e00")) +
  labs(color = "read length") +
  theme1 +
  geom_hline(
    yintercept = -log10(0.001) / s150,
    colour = "grey30",
    linetype = 2
  ) +
  geom_label(
    x = 80,
    y = -log10(0.001) / s150,
    label = "FPR=0.001",
    size = 4
  ) +
  # geom_hline(yintercept = 2 / s150,
  #            colour = "grey20",
  #            linetype = 2) +
  geom_hline(
    yintercept = -log10(0.05) / s150,
    colour = "grey20",
    linetype = 2
  ) +
  geom_label(
    x = 80,
    y = -log10(0.05) / s150,
    label = "FPR=0.05",
    size = 4
  ) +
  geom_vline(xintercept = 55,
             colour = "black",
             linetype = 2) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = custom_breaks_x,
    labels =  every_nth(custom_breaks_x, 2, inverse = TRUE)
  ) +
  scale_y_continuous(
    name = "Recall(%)",
    sec.axis = sec_axis(
      ~ . * s150,
      name = "-log10(Query FPR)",
      breaks = custom_breaks_y2,
      labels = every_nth(custom_breaks_y2, 2)
    ),
    expand = c(0, 0),
    breaks = custom_breaks_y,
    labels =  every_nth(custom_breaks_y, 2, inverse = TRUE)
  ) +
  theme(
    axis.title.y = element_text(color = "#36648e", size = 12),
    axis.title.y.right = element_text(color = "#D55E00", size = 12)
  )

## ------------------------------------------------------------------

theme_legend <-  theme(
  legend.key.size = unit(0.45, "cm"),
  legend.spacing.x = unit(0.2, "cm"),
  legend.spacing.y = unit(0.2, "cm"),
  legend.text = element_text(size = 10, angle = 0)
)

pident55 <-
  df3 %>% filter(qcov >= 55 & qcov < 56) %>% select(pident)
p <- plot_grid(
  plot_grid(
    addplot(p1) +
      theme_legend +
      theme(legend.position = "bottom"),
    addplot(p2) +
      scale_colour_gradient2_tableau(trans = "reverse") +
      theme_legend +
      theme(legend.position = "bottom"),
    addplot(p3) +
      scale_colour_gradient2_tableau(trans = "reverse") +
      theme_legend +
      theme(legend.position = "bottom"),
    ncol = 3,
    rel_widths = c(1, 1, 1),
    labels = c("a", "b", "c")
  ),
  plot_grid(
    addplot(p4) +
      geom_vline(
        xintercept = pident55$pident[1],
        colour = "black",
        linetype = 2
      ) +
      geom_hline(
        yintercept = 55,
        colour = "black",
        linetype = 2
      ) +
      theme_legend +
      theme(legend.position = "right") +
      scale_colour_gradient2_tableau(trans = "reverse"),
    # ggdraw() + draw_label(""),
    p5,
    
    ncol = 2,
    rel_widths = c(0.46, 0.5),
    labels = c("d", "e")
  ),
  nrow = 2,
  rel_heights = c(1, 1),
  labels = c("", "")
)

ggsave(p,
       file = args$outfile,
       width = 10,
       height = 7)

summary(model)