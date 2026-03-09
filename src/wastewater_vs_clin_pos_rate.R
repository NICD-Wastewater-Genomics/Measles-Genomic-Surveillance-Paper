# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(MMWRweek)
library(tidyr)
library(slider)

# The script generates Fig 1d, and SFigs 2a-d
setwd("Measles-Genomic-Surveillance-Paper/")

file_path <- "metadata/Wastewater metadata_20250218_SG.xlsx"
sheetname = "Sample Testing Results"


all_data <- read_excel(file_path, sheet = sheetname)
all_data <- all_data %>%
  rename(
    vl = `Measles Concentration (copies/uL)`
  )
all_data$vl <- as.numeric(all_data$vl)

all_data$combined_vl <- all_data$vl
all_data$combined_result <- all_data$`Measles Result`
all_data$`Sample Collection Date` = as.Date(all_data$`Sample Collection Date`)
all_data$combined_date <- all_data$`Sample Collection Date`
              
# Read data
clin_data <- read.csv("metadata/clinical_mev_pos_rate.csv")

mmwr <- MMWRweek(all_data$combined_date)
all_data$epiweek <- mmwr$MMWRweek
all_data$epiyear <- mmwr$MMWRyear
all_data <- all_data[all_data$`Measles Result` %in% c("Negative", "Positive"), ]


data <- all_data %>%
  filter(combined_result %in% c("Positive", "Negative")) %>%   # keep all valid results
  group_by(epiyear, epiweek) %>%
  summarise(
    n_total = n(),
    n_positive = sum(`Measles Result` == "Positive"),
    frac_positive = n_positive / n_total,
    vl = if (n_total >= 1)
        mean(combined_vl, na.rm = TRUE)
      else NA_real_,
    .groups = "drop"
  )

data <- data %>%
  complete(
    epiyear,
    epiweek,
    fill = list(vl = NA_real_, n_total = NA_real_, n_positive = 0, frac_positive = NA_real_)
  ) %>%
  arrange(epiyear, epiweek)

clin_data$IgM.positive[is.na(clin_data$IgM.positive)] <- 0
mmwr <- MMWRweek(clin_data$Date)
clin_data$epiweek <- mmwr$MMWRweek
clin_data$epiyear <- mmwr$MMWRyear

clin_data0 <- clin_data %>%
  group_by(epiyear, epiweek) %>%
  summarise(
    num_total_clin = sum(IgM.total.tested),
    num_positive_clin = sum(IgM.positive),
    frac_positive_clin = num_positive_clin / num_total_clin,
    .groups = "drop"
  )
  
clin_data0 <- clin_data0 %>%
  complete(
    epiyear,
    epiweek,
    fill = list(num_total_clin =  NA_real_, num_positive_clin = 0, frac_positive_clin = NA_real_)
  ) %>%
  arrange(epiyear, epiweek)

data = data %>%
  left_join(clin_data0, by = c("epiyear" = "epiyear", "epiweek" = "epiweek"))

# add epiweek binned plotting for x axis
data <- data %>%
  mutate(
    # Get start and end dates of each epiweek
    WeekStart = MMWRweek2Date(epiyear, epiweek, 1),
    WeekEnd   = MMWRweek2Date(epiyear, epiweek, 7),
    WeekLabel = if_else(
      epiweek == min(epiweek),  # first epiweek of each year
      paste0(format(WeekEnd, "%d\n%b"), "\n", epiyear),
      paste0(format(WeekEnd, "%d\n%b"))
    ),
    
    # Make sure order is preserved for plotting
    WeekLabel = factor(WeekLabel, levels = WeekLabel)
  )

binom_ci <- function(x, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half_width <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / denom
  
  tibble(
    lower = pmax(0, center - half_width),
    upper = pmin(1, center + half_width)
  )
}

data$PlotWeek <- 52*(as.integer(data$epiyear) - 2024) + data$epiweek

radius = 2
n=0.75
data <- data %>%
  arrange(epiyear, epiweek) %>%
  mutate(
    time_index = PlotWeek,
    total_roll = slide_dbl(n_total, ~mean(.x, na.rm = TRUE), .before = radius, .after = radius),
    pos_roll   = slide_dbl(n_positive, ~mean(.x, na.rm = TRUE), .before = radius, .after = radius),
    frac_pos_roll = pos_roll / total_roll,
    total_roll_clin = slide_dbl(num_total_clin, ~mean(.x, na.rm = TRUE), .before = radius, .after = radius),
    pos_roll_clin   = slide_dbl(num_positive_clin, ~mean(.x, na.rm = TRUE), .before = radius, .after = radius),
    frac_pos_roll_clin = pos_roll_clin / total_roll_clin,
    combo_clin = pos_roll_clin^n * frac_pos_roll_clin^(1-n),
    vl_roll = slide_dbl(vl, ~mean(.x, na.rm = TRUE), .before = radius, .after = radius)
  ) 

data <- data %>%
  rowwise() %>%
  mutate(
    ci = list(binom_ci(pos_roll, total_roll))
  ) %>%
  unnest(ci) %>%
  ungroup()

library(ggplot2)
library(patchwork)
scale_factor_clin <- max(data$pos_roll_clin, na.rm = TRUE) /
                     max(data$frac_pos_roll_clin, na.rm = TRUE)

## Clinical panel
p_clin <- ggplot(data, aes(x = PlotWeek)) +
  annotate(
  "rect",
  xmin = 31,
  xmax = 52,
  ymin = -Inf,
  ymax = Inf,
  fill = "grey70",
  alpha = 0.25
) +
  # Clinical positivity
  geom_line(
    aes(y = frac_pos_roll_clin, color = "Clinical IgM Positive Rate"),
    linewidth = 1.2
  ) +
  # Total clinical positives (scaled to primary axis)
  geom_line(
    aes(y = pos_roll_clin / scale_factor_clin,
        color = "Clinical IgM positives"),
    linewidth = 1.1,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    name = "Clinical IgM Positive Rate",
    expand = c(0, 0),
    limits = c(0, 0.45),
    sec.axis = sec_axis(
      ~ . * scale_factor_clin,
      name = "Clinical IgM positives"
    )
  ) +
  scale_x_continuous(
    limits = c(1,87), # through end Aug 2025.
    breaks = data$PlotWeek[seq(1, nrow(data), by = 4)],
    labels = data$WeekLabel[seq(1, nrow(data), by = 4)],
    expand = expansion(mult = c(0.0, 0.0))
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Clinical IgM Positive Rate" = "salmon",
      "Clinical IgM positives" = "firebrick"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.line = element_line(colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_line(),

    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = alpha("white", 0.8),
      color = NA
    ),
    legend.key = element_blank()
  )

# Wastewater panel

scale_factor_ww <- 
  max(data$vl_roll, na.rm = TRUE) /
  max(data$frac_pos_roll, na.rm = TRUE)

p_ww <- ggplot(data, aes(x = PlotWeek))+
  annotate(
  "rect",
  xmin = 31,
  xmax = 52,
  ymin = -Inf,
  ymax = Inf,
  fill = "grey70",
  alpha = 0.25
) +
  ## Binomial uncertainty on WW positivity
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "cornflowerblue",
    alpha = 0.25
  ) +
  
  ## WW positivity
  geom_line(
    aes(y = frac_pos_roll,color = "WW positive rate"),
    linewidth = 1.1
  ) +
  
  ## VL (log10, scaled to primary axis)
  geom_line(
    aes(y = vl_roll / scale_factor_ww, color = "Viral load(copies/uL)"),
    linewidth = 1.1,
    linetype = "dashed"
  ) +
  
  scale_y_continuous(
    name = "WW positive rate",
    expand = c(0, 0),
    
    ## Secondary axis
    sec.axis = sec_axis(
      ~ . * scale_factor_ww,
      name = expression("Viral load(copies/uL)"),
    )
  ) +
  
  scale_x_continuous(
    limits = c(1,87),
    breaks = data$PlotWeek[seq(1, nrow(data), by = 4)],
    labels = data$WeekLabel[seq(1, nrow(data), by = 4)],
    expand = expansion(mult = c(0.0, 0.0))
  ) +
  scale_color_manual(
  name = NULL,
  values = c(
    "WW positive rate" = "cornflowerblue",
    "Viral load(copies/uL)" = "deepskyblue4"
  )
) +
  theme_bw(base_size = 12) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
        legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = alpha("white", 0.8),
      color = NA
    ),
    legend.key = element_blank()
  )

p0 <- p_clin / p_ww + plot_layout(heights = c(1, 1))

ggsave(
  "figures/clin_pos_rate_and_ww_pos_stacked.pdf",
  plot = p0,
  width = 8.5,
  height = 4.5,
  units = "in",
  bg = "transparent"
)

# Generate counts figure - SFig 2a
p0 <- ggplot(data, aes(x = PlotWeek)) +
  geom_col(aes(y = n_total, fill = "Wastewater"), width = 0.7) +
  scale_y_continuous(
    name = "Wastewater Samples (Total)",
    expand=c(0,0)
  ) +
  scale_x_continuous(
    breaks = data$PlotWeek[seq(1, nrow(data), by = 4)],
    labels = data$WeekLabel[seq(1, nrow(data), by = 4)],
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_fill_manual(
    name = NULL,  
    values = c("Wastewater" = "steelblue")
  ) +
  theme_bw(base_size = 15) +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 0),
        legend.position = "top") +
  labs(
    x = "Epidemiological Week"
  )
ggsave("figures/ww_sample_counts.pdf", plot = p0, width = 8, height = 4, units = "in")


# Generate plot of clinical counts over time, SFig 2b
p0 <- ggplot(data, aes(x = PlotWeek)) +
  annotate(
  "rect",
  xmin = 31,
  xmax = 52,
  ymin = -Inf,
  ymax = Inf,
  fill = "grey70",
  alpha = 0.25
) +
  geom_col(aes(y = num_total_clin , fill = "Clinical"), linewidth = 1.2) +
  scale_y_continuous(
    name = "Clinical tests per week",
    expand=c(0,0)
  ) +
  scale_x_continuous(
    breaks = data$PlotWeek[seq(1, nrow(data), by = 4)],
    labels = data$WeekLabel[seq(1, nrow(data), by = 4)],
    expand = expansion(mult = c(0.0, 0.0))
  ) +
  scale_fill_manual(
    name = NULL,  # Removes "Indicator"
    values = c("Clinical" = "firebrick")
  ) +
  theme_bw(base_size = 15) +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 0),
        legend.position = "top") +
  labs(
    x = "Epidemiological Week"
  )
ggsave("figures/clin_tests_over_time.pdf", plot = p0, width = 8, height = 4, units = "in")


## Now investigate covariates across modes (SFig 2c,d) 
df <- data %>%
  select(frac_positive, frac_positive_clin) %>%
  filter(!is.na(frac_positive), !is.na(frac_positive_clin))

bootstrap_pred <- function(df, x_pred, B = 1000) {
  boot_mat <- replicate(B, {
    samp <- df[sample(nrow(df), replace = TRUE), ]
    fit  <- lm(frac_positive_clin ~ frac_positive, data = samp)
    predict(fit, newdata = x_pred)
  })
  
  tibble(
    frac_positive = x_pred$frac_positive,
    fit   = apply(boot_mat, 1, median),
    lower = apply(boot_mat, 1, quantile, 0.025),
    upper = apply(boot_mat, 1, quantile, 0.975)
  )
}

x_pred <- tibble(
  frac_positive = seq(min(df$frac_positive), max(df$frac_positive), length.out = 200)
)

pred <- bootstrap_pred(df, x_pred, B = 1000)
p= ggplot() +
  geom_point(
    data = df,
    aes(x = frac_positive, y = frac_positive_clin),
    size = 2, alpha = 0.5,color = "steelblue"
  ) +
  
  geom_ribbon(
    data = pred,
    aes(x = frac_positive, ymin = lower, ymax = upper),
    fill = "steelblue", alpha = 0.25
  ) +
  
  geom_line(
    data = pred,
    aes(x = frac_positive, y = fit),
    color = "steelblue", size = 1.2
  ) +
  
  theme_bw(base_size = 15) +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) +
  labs(
    x = "Positivity (WW)",
    y = "Positivity (Clinical)"
  )
ggsave("figures/covariates_positivity_only.pdf", plot = p, width = 5, height = 5, units = "in")

df <- data %>%
  select(frac_positive, num_positive_clin) %>%
  filter(!is.na(frac_positive), !is.na(num_positive_clin))

bootstrap_pred <- function(df, x_pred, B = 1000) {
  boot_mat <- replicate(B, {
    samp <- df[sample(nrow(df), replace = TRUE), ]
    fit  <- lm(num_positive_clin ~ frac_positive, data = samp)
    predict(fit, newdata = x_pred)
  })
   
  tibble(
    frac_positive = x_pred$frac_positive,
    fit = apply(boot_mat, 1, mean),
    lower    = apply(boot_mat, 1, quantile, 0.025),
    upper    = apply(boot_mat, 1, quantile, 0.975)
  )
}

x_pred <- tibble(
  frac_positive = seq(min(df$frac_positive), max(df$frac_positive), length.out = 200)
)

pred <- bootstrap_pred(df, x_pred, B = 1000)


p= ggplot() +
  geom_point(
    data = df,
    aes(x = frac_positive, y = num_positive_clin),
    size = 2, alpha = 0.5,color = "steelblue"
  ) +
  
  geom_ribbon(
    data = pred,
    aes(x = frac_positive, ymin = lower, ymax = upper),
    fill = "steelblue", alpha = 0.25
  ) +
  
  geom_line(
    data = pred,
    aes(x = frac_positive, y = fit),
    color = "steelblue", size = 1.2
  ) +
  
  theme_bw(base_size = 15) +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) +
  labs(
    x = "Positivity (WW)",
    y = "Clinical positives"
  )
ggsave("figures/covariates_positivity_and_clinical_pos.pdf", plot = p, width = 5, height = 5, units = "in")

### calculate correlation between the two (no smoothing)
h = cor.test(data$frac_positive, data$num_positive_clin, use = "complete.obs", method = "pearson")
h
print(paste('correlation results:',h))
