library(tidyverse)
library(data.table)
library(MetBrewer)
library(zoo)

# read freyja aggregate results
results<-read.table("aggregated_result.tsv", fill = TRUE, sep = "\t", h=T)
results<-as.data.frame(sapply(results, function(x) str_replace_all(x, "[',()\\]\\[]", ""))) # Removed the unwanted character: [], () and commas
results<-as.data.frame(sapply(results, function(x) trimws(gsub("\\s+", " ", x)))) # Removed double spaces

# read lineages file created using nextclade
lineages <- fread("lineages.tsv",header = TRUE)
# df operations
results_comb_names <- results %>%
  separate(X, into = as.character(1:6), sep = "_") %>%
  mutate(
    isolate1 = paste(`1`, `2`, sep = "_"),
    isolate2 = paste(`4`, `5`, sep = "_"),
    exp_abun1 = as.numeric(`3`),
    exp_abun2 = as.numeric(`6`)
  )

# 2. Join expected lineages from nextclade
results_comb_names <- results_comb_names %>%
  left_join(lineages %>% select(seqName, exp_lin1 = clade), by = c("isolate1" = "seqName")) %>%
  left_join(lineages %>% select(seqName, exp_lin2 = clade), by = c("isolate2" = "seqName"))

# 3. Clean up Freyja observed lineages & abundances into long format for matching
obs_long <- results_comb_names %>%
  select(isolate1, isolate2, exp_lin1, exp_lin2, exp_abun1, exp_abun2, lineages, abundances) %>%
  separate(lineages, into = c("obs_lin_a", "obs_lin_b"), sep = " ", fill = "right") %>%
  separate(abundances, into = c("obs_abun_a", "obs_abun_b"), sep = " ", fill = "right") %>%
  mutate(
    obs_lin_a = str_remove(obs_lin_a, ".*-"),
    obs_lin_b = str_remove(obs_lin_b, ".*-"),
    obs_abun_a = replace_na(as.numeric(obs_abun_a), 0),
    obs_abun_b = replace_na(as.numeric(obs_abun_b), 0)
  )

# 4. Map observed values to expected lineage 1 & 2 explicitly
df_mapped <- obs_long %>%
  mutate(
    # Isolate 1 matching
    obs_abun1 = case_when(
      exp_lin1 == obs_lin_a ~ obs_abun_a,
      exp_lin1 == obs_lin_b ~ obs_abun_b,
      TRUE ~ 0
    ),
    # Isolate 2 matching
    obs_abun2 = case_when(
      exp_lin2 == obs_lin_a ~ obs_abun_a,
      exp_lin2 == obs_lin_b ~ obs_abun_b,
      TRUE ~ 0
    )
  )

# 5. Build final plotting data frame across full 0.0 - 1.0 range
df_all <- bind_rows(
  df_mapped %>% select(exp = exp_abun1, obs = obs_abun1, lineage = exp_lin1),
  df_mapped %>% select(exp = exp_abun2, obs = obs_abun2, lineage = exp_lin2)
)
# calculate R^2
R2 <- 1 - sum((df_all$obs - df_all$exp)^2) / sum((df_all$obs - mean(df_all$obs))^2)
# Create the ggplot with combined, detailed color and size legend
plot <- df_all %>%
  filter(exp >0) %>% ggplot() + 
  geom_point(aes(exp, obs, color = lineage),alpha = 0.4, size = 3) +  # Outline color based on depth
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  theme_classic() +
  theme(
    axis.ticks.length=unit(.17, "cm"),
    legend.text = element_text(size = 18),     
    legend.title = element_text(size = 20),    
    axis.text.x = element_text(size = 16),     
    axis.text.y = element_text(size = 16),     
    axis.title.x = element_text(size = 20),    
    axis.title.y = element_text(size = 20),    
    plot.title = element_text(size = 22)
  ) +
  scale_x_continuous(
    limits = c(0,1),
    expand = c(0.02, 0)  
  ) +
  scale_y_continuous(
    limits = c(0, 1),  
    expand = c(0.02, 0)
  ) +
  xlab("Expected abundance") +
  ylab("Observed abundance")

plot
ggsave("exp-obs-plot.png", plot = plot, device = "png", width = 10, height = 10)

