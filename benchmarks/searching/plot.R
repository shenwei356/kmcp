#!/usr/bin/env Rscript
library(tidyverse)
library(dplyr)
library(scales)
library(ggthemes)
library(cowplot)

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
        plot.title = element_text(size = 15),
        plot.margin = unit(c(0.1,0.1,0.1,0.5), "cm"),
    )

## ------------------------------------------------------------------


df <- read.csv("bench.kmcp-cobs.csv")

df$app <- factor(df$app, levels = c('KMCP', 'COBS'))
df$group <- factor(df$group, levels = unique(df$group))

n <- length(unique(df$app))

colors <- colorblind_pal()(n+1)[2:(n+1)]

p1 <- ggplot(df,
            aes(app, time, color = app)) +
    geom_boxplot(width = 0.6) + 
    geom_jitter(width = 0.2) +
    geom_point(size = 1.2) +
    coord_flip() +
    facet_grid(. ~ group, scales = "free_x") + 
    xlab(NULL) +
    ylab('time (second)') +
    expand_limits(y = 0) + 
    scale_color_manual(values = colors) +
    theme1 +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 11, color="black"),
          axis.text.x = element_text(size = 10, color="grey20"),
          axis.title.x = element_text(size = 10, color="grey20"),
          )

## ------------------------------------------------------------------

df2 <- read.csv("bench.kmcp-mash-sourmash.csv")

df2$app <- factor(df2$app, levels = c('KMCP', 'Mash', 'Sourmash'))

n <- length(unique(df2$app))

colors <- colorblind_pal()(n+1)[2:(n+1)]

p2 <- ggplot(df2,
            aes(app, time, color = app)) +
    geom_boxplot(width = 0.6) + 
    geom_jitter(width = 0.2) +
    geom_point(size = 1.2) +
    coord_flip() +
    facet_grid(. ~ group, scales = "free_x") + 
    xlab(NULL) +
    ylab('time (second)') +
    expand_limits(y = 0) + 
    scale_color_manual(values = colors) +
    theme1 +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 11, color="black"),
          axis.text.x = element_text(size = 10, color="grey20"),
          axis.title.x = element_text(size = 10, color="grey20"),
    )

## ------------------------------------------------------------------

p <- plot_grid(
    p1,
    p2,
    nrow = 2,
    labels = c("a", "b"),
    rel_heights = c(1.8, 2)
) + theme ( # fill the gap in sub figures
    panel.background = element_rect(fill = "white", colour = NA),
) 



ggsave(p,
       file = "bench.searching.jpg",
       width = 5,
       height = 3.5)
