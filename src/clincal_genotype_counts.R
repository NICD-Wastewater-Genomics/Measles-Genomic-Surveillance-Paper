#  Load required packages 
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# this generates Figure 2a. 
setwd("Measles-Genomic-Surveillance-Paper/")

#  Load Excel data 
igm_data <- read_excel("metadata/measles_paper_clinicaldata_18062025_SG.xlsx", sheet = "IgM_Positive")
seq_data <- read_excel("metadata/measles_paper_clinicaldata_18062025_SG.xlsx", sheet = "Sequencing_Outcomes")
genotype_data <- read_excel("metadata/measles_paper_clinicaldata_18062025_SG.xlsx", sheet = "Genotypes")

#  Prepare sequencing outcomes 
df_seq <- seq_data %>%
  pivot_longer(cols = c(Successful, Failed), names_to = "Outcome", values_to = "Count")

#  Prepare genotype data (D8 and B3 only) 
genotypes <- genotype_data %>%
  select(Year, D8, B3) %>%
  pivot_longer(cols = c(D8, B3), names_to = "Genotype", values_to = "Count")

df_seq_cts <- df_seq %>%
  group_by(year = Year) %>%
  summarize(total= sum(Count))

genotype_cts <- genotypes %>%
  group_by(year = Year) %>%
  summarize(total= sum(Count))

p3 <- ggplot(genotypes, aes(x = factor(Year), y = Count, fill = Genotype)) +
  geom_bar(stat = "identity") + coord_cartesian(ylim=c(0, 1.05*max(genotype_cts$total)))+
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("D8" = "#fc8d62", "B3" = "#66c2a5")) +
  labs(y = "Genotype Count", x = "Year") +
  theme_bw(base_size = 18) +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) 
ggsave("figures/clinical_sequencing_genotypes.pdf", plot = p3, width = 5, height = 5, units = "in")
