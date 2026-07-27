library(data.table)
library(tidyverse)
library(forcats)

# read data 
df_all <- fread("amplicon_stats.csv")
sequences <- fread("primer_sequences.csv")

df_joined <- df_all %>% inner_join(sequences, by = c("genome_id" = "sequence_id",
                                                     "amplicon_number"))
df_joined <- df_joined %>%
  mutate(
    sequence_status = case_when(
      str_detect(left_aligned_sequence, "-") | str_detect(right_aligned_sequence, "-") ~ "sequence missing",
      str_detect(left_aligned_sequence, "[NRYSWKMBDHV]") | 
        str_detect(right_aligned_sequence, "[NRYSWKMBDHV]") ~ "sequence ambiguity",
      TRUE ~ "full sequence"
    )
  )
lineages <- fread("nextclade_sequence_lineages.tsv")

dropout_summary <-
  df_joined%>%
  filter(sequence_status == "full sequence")%>% 
  separate(amplicon_number, into = c("number","redundant"))%>%
  inner_join(lineages, by = c("genome_id" = "V1"))%>% 
  group_by(number, genome_id) %>% 
  summarise( successful = any(amplicon_length != 0), .groups = "drop" )%>%
  ungroup()%>% inner_join(lineages, by = c("genome_id" = "V1")) %>%
  group_by(clade,number) %>%
  summarise( n_success = sum(successful), 
             n_failure = sum(!successful),
             total = n(), dropout_ratio = n_failure / total, .groups = "drop" )

p<-dropout_summary %>%
  filter(clade != "unassigned") %>%
  mutate(
    number = factor(
      number,
      levels = sort(unique(as.numeric(number)))
    )
  ) %>%
  complete(clade, number, fill = list(dropout_ratio = NA)) %>%
  ggplot(aes(x = clade, y = number, fill = dropout_ratio)) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "darkblue",
    na.value = "darkgrey",
    name = "Dropout rate"
  ) +
  labs(x = "Clade", y = "Amplicon") +
  theme_light()+
  theme(
    axis.text.y = element_text(size = 8)
  )
ggsave(
  "dropout_heatmap.pdf",
  plot = p,
  width = 5,
  height = 5,
  units = "in"
)

