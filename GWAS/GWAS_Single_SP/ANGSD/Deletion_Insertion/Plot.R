##############################
# Install required libraries #
##############################
#install.packages("ggfortify")
library(ggfortify)
#install.packages("qqman")
library(qqman)
library(tidyverse)
#install.packages("ggtext")
library(ggtext)
#devtools::install_github("norment/normentR")
library(normentR)
library(grid)
library(gridExtra)
library(dplyr)
#install.packages("sjlabelled")
library("sjlabelled")
library(ggplot2)
library(cowplot)
#remotes::install_github("HanjoStudy/quotidieR")
library(quotidieR)
#install.packages("ggthemes")
library("ggthemes")

###############
# Import data #
###############
dat = commandArgs(trailingOnly=TRUE)
print(-log10(0.05/as.numeric(dat[3])))
gwas = read.table(dat[1],skip=1)
bonf_threshold = -log10(0.05/as.numeric(dat[3]))

################          
# Prepare data #
################
gwas = gwas[,c(1,2,3)]
colnames(gwas) = c("CHR","BP","P")
order_contigs = read.table(dat[2])
order_contigs = as.data.frame(order_contigs)
order_contigs = as.factor(order_contigs$V1)

gwas$CHR=factor(gwas$CHR, levels=order_contigs)

sig_data <- gwas %>% 
  subset(P < 1)
notsig_data <- gwas %>%
  subset(P >= 1) %>% 
  group_by(CHR) %>% 
  sample_frac(1)
gwas_data <- bind_rows(sig_data, notsig_data)

data_cum <- gwas_data %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(BP)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(CHR, bp_add)

gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(bp_cum = BP + bp_add)

axis_set <- gwas_data %>% 
  group_by(CHR) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(P == min(P)) %>% 
  mutate(ylim = abs(floor(log10(P))) + 2) %>% 
  pull(ylim)


########
# Plot #
########
png(file="GWAS.png",width=800,height=500,type="cairo")
ggplot(gwas_data, aes(x = bp_cum, y = -log10(P), 
                                  color =CHR)) +
  geom_point(alpha = 0.75) +
  scale_x_continuous() +
  scale_color_manual(values = rep(c("blue4", "orange3"), unique(length(axis_set$CHR)))) +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(size=10),
    plot.title = element_text(hjust = 0.5,size=10)) +
    geom_hline(yintercept=bonf_threshold, linetype="dashed", 
                color = "red", size=0.8)
dev.off()
