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
results_comb_names <- results %>% separate(lineages, into = c("obs_lin1","obs_lin2"), sep = " ")
results_comb_names <- results_comb_names %>% 
  separate(abundances, into = c("obs_abun1","obs_abun2"), sep = " ")
results_comb_names <- results_comb_names %>% 
  separate(X, into = as.character(1:6),sep="_")
results_comb_names$isolate1 <- paste(results_comb_names$`1`,results_comb_names$`2`,sep="_")
results_comb_names$isolate2 <- paste(results_comb_names$`4`,results_comb_names$`5`,sep="_")
results_comb_names<- results_comb_names[,-c(1,2,4,5,7,12,13)]
colnames(results_comb_names)[1:2] <- c("exp_abun1","exp_abun2")
results_comb_names <- results_comb_names %>%
  mutate(obs_abun1 = replace_na(as.numeric(obs_abun1), 0),
         obs_abun2 = replace_na(as.numeric(obs_abun2), 0),
         real_exp_abun1 = pmax(exp_abun1, exp_abun2, na.rm = TRUE),
         real_exp_abun2 = pmin(exp_abun1, exp_abun2, na.rm = TRUE))

results_comb_names <- results_comb_names %>%
  mutate(res1 = abs(as.numeric(obs_abun1) - as.numeric(real_exp_abun1)))
results_comb_names <- results_comb_names %>%
  mutate(res2 = abs(as.numeric(obs_abun2) - as.numeric(real_exp_abun2)))

results_comb_names <- lineages %>% select(seqName,clade)%>%
  inner_join(results_comb_names, by =c("seqName" ="isolate1")) 
colnames(results_comb_names)[1:2] <- c("isolate1","exp_lin1")

results_comb_names <- lineages %>% select(seqName,clade)%>%
  inner_join(results_comb_names, by =c("seqName" ="isolate2")) 
colnames(results_comb_names)[1:2] <- c("isolate2","exp_lin2")

results_comb_names <- results_comb_names %>%
  separate(obs_lin1, into= c("na","obs_lin1"), sep = "-")%>%
  separate(obs_lin2, into= c("na","obs_lin2"), sep = "-")
results_comb_names <- results_comb_names %>% select(-c(exp_abun1,exp_abun2,na))
colnames(results_comb_names)[9:10] <- c("exp_abun1","exp_abun2")
results_comb_names <- results_comb_names %>% select(isolate1,isolate2,exp_lin1,
                                                    obs_lin1,exp_abun1,obs_abun1,
                                                    res1,exp_lin2,obs_lin2,
                                                    exp_abun2,obs_abun2,res2)
results_comb_names <- results_comb_names %>%
  mutate(
    sublineage_true1 = if_else(
      (coalesce(exp_lin1, "NA") == coalesce(obs_lin1, "NA") | 
         coalesce(exp_lin1, "NA") == coalesce(obs_lin2, "NA")),
      1,
      0
    ),
    sublineage_true2 = if_else(
      (coalesce(exp_lin2, "NA") == coalesce(obs_lin1, "NA") | 
         coalesce(exp_lin2, "NA") == coalesce(obs_lin2, "NA")),
      1,
      0
    )
  )

results_comb_names <- results_comb_names %>%
  mutate(
    obs_abun1_real = if_else(sublineage_true1 == "1", as.character(obs_abun1), "0"),
    obs_abun2_real = if_else(sublineage_true2 == "1", as.character(obs_abun2), "0"),
    obs_abun1_real = case_when(sublineage_true1 != sublineage_true2 ~ as.character(pmax(obs_abun1,obs_abun2)),
                               sublineage_true2 == sublineage_true1 ~ as.character(obs_abun1)
    ),
    obs_abun2_real = case_when(sublineage_true1 != sublineage_true2 ~ as.character(pmin(obs_abun1,obs_abun2)),
                               sublineage_true2 == sublineage_true1 ~ as.character(obs_abun2)
    ),
  )

results_comb_names %>% write_csv("Downloads/tb-mixed-simulations.csv")
# combine expected observed values
res1 <- results_comb_names %>% select(exp_abun1,obs_abun1_real,exp_lin1)
res2 <- results_comb_names %>% select(exp_abun2,obs_abun2_real,exp_lin2)
colnames(res1) <- c("exp","obs","lineage")
colnames(res2)<- c("exp","obs","lineage")

df_all <- rbind(res1,res2)
df_all <- df_all %>% mutate(exp = as.numeric(exp),
                            obs=as.numeric(obs))
# calculate R^2
R2 <- 1 - sum((df_all$obs - df_all$exp)^2) / sum((df_all$obs - mean(df_all$obs))^2)
# Create the ggplot with combined, detailed color and size legend
plot <- df_all %>%
  filter(exp >0) %>% ggplot() + 
  geom_point(aes(exp, obs, color = lineage),alpha = 0.7, size = 3) +  # Outline color based on depth
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
ggsave("exp-obs-plot.png", plot = plot, device = "pdf", width = 10, height = 10)

