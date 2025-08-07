library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(patchwork)

setwd("~/GitHub/Mutation_Selection_Balance_Autos_Allos/Matlab_scripts/equal_mut_bifn")

### ===========================================================
### READ IN THE DATA, PROCESS, AND CALCULATE SUMMARY STATISTICS
### ===========================================================

# Define color mapping
dip_color <- "#F04D13"
auto_color <- "#66BED6"
allo_color <- "#7DAB5B"

# Read in all datasets
dip_rec <- read_csv("dip_rec.csv", col_names = c("s", "q", "g0", "g1", "w"))
dip_part_rec <- read_csv("dip_part_rec.csv", col_names = c("s", "q", "g0", "g1", "w"))
dip_add <- read_csv("dip_add.csv", col_names = c("s", "q", "g0", "g1", "w"))
dip_part_dom <- read_csv("dip_part_dom.csv", col_names = c("s", "q", "g0", "g1", "w"))
dip_dom <- read_csv("dip_dom.csv", col_names = c("s", "q", "g0", "g1", "w"))

auto_rec <- read_csv("auto_rec.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_part_rec <- read_csv("auto_part_rec.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_add <- read_csv("auto_add.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_part_dom <- read_csv("auto_part_rec.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_dom <- read_csv("auto_dom.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))

allo_rec <- read_csv("allo_rec.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_part_rec <- read_csv("allo_part_rec.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_add <- read_csv("allo_add.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_part_dom <- read_csv("allo_part_rec.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_dom <- read_csv("allo_dom.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))

# Add dominance labels
dip_rec$Dominance       <- "Recessive"
dip_part_rec$Dominance  <- "Partially Recessive"
dip_add$Dominance       <- "Additive"
dip_part_dom$Dominance  <- "Partially Dominant"
dip_dom$Dominance       <- "Dominant"

auto_rec$Dominance       <- "Recessive"
auto_part_rec$Dominance  <- "Partially Recessive"
auto_add$Dominance       <- "Additive"
auto_part_dom$Dominance  <- "Partially Dominant"
auto_dom$Dominance       <- "Dominant"

allo_rec$Dominance       <- "Recessive"
allo_part_rec$Dominance  <- "Partially Recessive"
allo_add$Dominance       <- "Additive"
allo_part_dom$Dominance  <- "Partially Dominant"
allo_dom$Dominance       <- "Dominant"

# Bind each ploidy group together
dip  <- bind_rows(dip_rec, dip_part_rec, dip_add, dip_part_dom, dip_dom)
auto <- bind_rows(auto_rec, auto_part_rec, auto_add, auto_part_dom, auto_dom)
allo <- bind_rows(allo_rec, allo_part_rec, allo_add, allo_part_dom, allo_dom)

# Factor by dominance
dip$Dominance <- factor(dip$Dominance, levels = c("Recessive", "Partially Recessive", "Additive", "Partially Dominant", "Dominant"))
auto$Dominance <- factor(auto$Dominance, levels = c("Recessive", "Partially Recessive", "Additive", "Partially Dominant", "Dominant"))
allo$Dominance <- factor(allo$Dominance, levels = c("Recessive", "Partially Recessive", "Additive", "Partially Dominant", "Dominant"))

# Calculate summary statistics for each dataset, including load, PD, and genotypes
dip <- dip %>% 
  mutate(G0 = g0^2, 
         G1 = 2*g0*g1, 
         G2 = g1^2, 
         load = 1-w, 
         PD = 0, 
         Model = "Diploid")

auto <- auto %>% 
  mutate(G0 = g0^2, 
         G1 = 2*g0*g1, 
         G2 = 2*g0*g2 + g1^2, 
         G3 = 2*g1*g2, 
         G4 = g2^2,
         load = 1-w, 
         PD = g2 - q^2, 
         Model = "Autotetraploid")

allo <- allo %>% 
  mutate(G0 = g00^2, 
         G1 = 2*g00*(g01+g10), 
         G2 = g01^2 + 2*(g00*g11 + g01*g10) + g10^2,
         G3 = 2*g11*(g01 + g10),
         G4 = g11^2,
         load = 1-w, 
         PD = g11 - q^2, 
         Model = "Allotetraploid")

# Combine all datasets
full_data <- bind_rows(dip, auto, allo)

# Set ploidy as a factor for consistent color and legend order
full_data$Model <- factor(full_data$Model, levels = c("Diploid", "Autotetraploid", "Allotetraploid"))

# Color mapping
ploidy_colors <- c("Diploid" = dip_color,
                   "Autotetraploid" = auto_color,
                   "Allotetraploid" = allo_color)

### ====================================================
### FIRST PLOT - BIFURCATION DIAGRAM SUBSET BY DOMINANCE
### ====================================================

# Create data frame for panel labels
panel_labels_top <- data.frame(
  Dominance = factor(c("Recessive", "Partially Recessive", "Additive", "Partially Dominant", "Dominant"),
                     levels = c("Recessive", "Partially Recessive", "Additive", "Partially Dominant", "Dominant")),
  label = c("A", "B", "C", "D", "E"),
  x = rep(3e-4, 5),  # x position for labels
  y = rep(0.95, 5)   # y position for labels
)

panel_labels_bottom <- data.frame(
  Dominance = factor(c("Recessive", "Partially Recessive", "Additive", "Partially Dominant", "Dominant"),
                     levels = c("Recessive", "Partially Recessive", "Additive", "Partially Dominant", "Dominant")),
  label = c("F", "G", "H", "I", "J"),
  x = rep(3e-4, 5),   # x position for labels
  y = rep(4.75e-8, 5) # y position for labels
)

# Top row: Allele frequency (q) plot
p_q <- ggplot(full_data, aes(x = s, y = q, color = Model)) +
  geom_line() +
  scale_color_manual(values = ploidy_colors) +
  facet_wrap(~Dominance, nrow = 1, strip.position = "top") +
  scale_x_log10(
    limits = c(1e-9, 1e-3),
    breaks = c(1e-8, 1e-6, 1e-4)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  # Add panel labels
  geom_text(data = panel_labels_top, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 4, hjust = 0) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(y = "Allele Frequency", x = NULL, color = "Model")

# Bottom row: Mutation load plot
p_load <- ggplot(full_data, aes(x = s, y = load, color = Model)) +
  geom_line() +
  scale_color_manual(values = ploidy_colors) +
  facet_wrap(~Dominance, nrow = 1) +  
  scale_x_log10(
    limits = c(1e-9, 1e-3),
    breaks = c(1e-8, 1e-6, 1e-4),
    labels = scales::label_scientific()
  ) +
  scale_y_continuous(
    limits = c(0, 5e-8),
    breaks = c(0, 1e-8, 2e-8, 3e-8, 4e-8, 5e-8),
    labels = scales::label_scientific()
  ) +
  # Add panel labels
  geom_text(data = panel_labels_bottom, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 4, hjust = 0) +
  theme_bw(base_size = 11) +
  labs(x = "s (Selection Coefficient)", y = "Mutation Load") +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  )

# Combine plots
final_plot <- p_q / p_load +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Save
ggsave("equal_mut_bifn.pdf", final_plot, width = 6.5, height = 4, units = "in")


### ================================================
### SECOND PLOT: DIFFERENCES BETWEEN AUTOS AND ALLOS
### ================================================
# First, we need to read in the additional data - note we don't need the diploid data here
setwd("~/GitHub/Mutation_Selection_Balance_Autos_Allos/Matlab_scripts/equal_mut_bifn/mu_1e-7")

auto_rec_7 <- read_csv("auto_rec.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_part_rec_7 <- read_csv("auto_part_rec.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_add_7 <- read_csv("auto_add.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_part_dom_7 <- read_csv("auto_part_rec.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_dom_7 <- read_csv("auto_dom.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))

allo_rec_7 <- read_csv("allo_rec.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_part_rec_7 <- read_csv("allo_part_rec.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_add_7 <- read_csv("allo_add.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_part_dom_7 <- read_csv("allo_part_rec.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_dom_7 <- read_csv("allo_dom.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))

# Add dominance labels
auto_rec_7$Dominance       <- "Recessive"
auto_part_rec_7$Dominance  <- "Partially Recessive"
auto_add_7$Dominance       <- "Additive"
auto_part_dom_7$Dominance  <- "Partially Dominant"
auto_dom_7$Dominance       <- "Dominant"

allo_rec_7$Dominance       <- "Recessive"
allo_part_rec_7$Dominance  <- "Partially Recessive"
allo_add_7$Dominance       <- "Additive"
allo_part_dom_7$Dominance  <- "Partially Dominant"
allo_dom_7$Dominance       <- "Dominant"

# Bind each ploidy group together
auto_7 <- bind_rows(auto_rec_7, auto_part_rec_7, auto_add_7, auto_part_dom_7, auto_dom_7)
allo_7 <- bind_rows(allo_rec_7, allo_part_rec_7, allo_add_7, allo_part_dom_7, allo_dom_7)

# Factor by dominance
auto_7$Dominance <- factor(auto_7$Dominance, levels = c("Recessive", "Partially Recessive", "Additive", "Partially Dominant", "Dominant"))
allo_7$Dominance <- factor(allo_7$Dominance, levels = c("Recessive", "Partially Recessive", "Additive", "Partially Dominant", "Dominant"))

# Calculate summary statistics for each dataset, including load, PD, and genotypes
auto_7 <- auto_7 %>% 
  mutate(G0 = g0^2, 
         G1 = 2*g0*g1, 
         G2 = 2*g0*g2 + g1^2, 
         G3 = 2*g1*g2, 
         G4 = g2^2,
         load = 1-w, 
         PD = g2 - q^2, 
         Model = "Autotetraploid")

allo_7 <- allo_7 %>% 
  mutate(G0 = g00^2, 
         G1 = 2*g00*(g01+g10), 
         G2 = g01^2 + 2*(g00*g11 + g01*g10) + g10^2,
         G3 = 2*g11*(g01 + g10),
         G4 = g11^2,
         load = 1-w, 
         PD = g11 - q^2, 
         Model = "Allotetraploid")

# Then, do the same for the 1e-6 mutation case
setwd("~/GitHub/Mutation_Selection_Balance_Autos_Allos/Matlab_scripts/equal_mut_bifn/mu_1e-6")

auto_rec_6 <- read_csv("auto_rec.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_part_rec_6 <- read_csv("auto_part_rec.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_add_6 <- read_csv("auto_add.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_part_dom_6 <- read_csv("auto_part_rec.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_dom_6 <- read_csv("auto_dom.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))

allo_rec_6 <- read_csv("allo_rec.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_part_rec_6 <- read_csv("allo_part_rec.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_add_6 <- read_csv("allo_add.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_part_dom_6 <- read_csv("allo_part_rec.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_dom_6 <- read_csv("allo_dom.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))

# Add dominance labels
auto_rec_6$Dominance       <- "Recessive"
auto_part_rec_6$Dominance  <- "Partially Recessive"
auto_add_6$Dominance       <- "Additive"
auto_part_dom_6$Dominance  <- "Partially Dominant"
auto_dom_6$Dominance       <- "Dominant"

allo_rec_6$Dominance       <- "Recessive"
allo_part_rec_6$Dominance  <- "Partially Recessive"
allo_add_6$Dominance       <- "Additive"
allo_part_dom_6$Dominance  <- "Partially Dominant"
allo_dom_6$Dominance       <- "Dominant"

# Bind each ploidy group together
auto_6 <- bind_rows(auto_rec_6, auto_part_rec_6, auto_add_6, auto_part_dom_6, auto_dom_6)
allo_6 <- bind_rows(allo_rec_6, allo_part_rec_6, allo_add_6, allo_part_dom_6, allo_dom_6)

# Factor by dominance
auto_6$Dominance <- factor(auto_6$Dominance, levels = c("Recessive", "Partially Recessive", "Additive", "Partially Dominant", "Dominant"))
allo_6$Dominance <- factor(allo_6$Dominance, levels = c("Recessive", "Partially Recessive", "Additive", "Partially Dominant", "Dominant"))

# Calculate summary statistics for each dataset, including load, PD, and genotypes
auto_6 <- auto_6 %>% 
  mutate(G0 = g0^2, 
         G1 = 2*g0*g1, 
         G2 = 2*g0*g2 + g1^2, 
         G3 = 2*g1*g2, 
         G4 = g2^2,
         load = 1-w, 
         PD = g2 - q^2, 
         Model = "Autotetraploid")

allo_6 <- allo_6 %>% 
  mutate(G0 = g00^2, 
         G1 = 2*g00*(g01+g10), 
         G2 = g01^2 + 2*(g00*g11 + g01*g10) + g10^2,
         G3 = 2*g11*(g01 + g10),
         G4 = g11^2,
         load = 1-w, 
         PD = g11 - q^2, 
         Model = "Allotetraploid")

# change the wd again to save plots in the right place
setwd("~/GitHub/Mutation_Selection_Balance_Autos_Allos/Matlab_scripts/equal_mut_bifn")

# 1e-8 data
auto_subset_8 <- auto %>% select(s, Dominance, q, load)
allo_subset_8 <- allo %>% select(s, Dominance, q, load)

diff_data_8 <- inner_join(auto_subset_8, allo_subset_8, by = c("s", "Dominance"), suffix = c("_auto", "_allo")) %>%
  mutate(
    q_diff = q_auto - q_allo,
    load_diff = load_auto - load_allo,
    mu = "1e-8"
  )

# 1e-7 data
auto_subset_7 <- auto_7 %>% select(s, Dominance, q, load)
allo_subset_7 <- allo_7 %>% select(s, Dominance, q, load)

diff_data_7 <- inner_join(auto_subset_7, allo_subset_7, by = c("s", "Dominance"), suffix = c("_auto", "_allo")) %>%
  mutate(
    q_diff = q_auto - q_allo,
    load_diff = load_auto - load_allo,
    mu = "1e-7"
  )

# 1e-6 data
auto_subset_6 <- auto_6 %>% select(s, Dominance, q, load)
allo_subset_6 <- allo_6 %>% select(s, Dominance, q, load)

diff_data_6 <- inner_join(auto_subset_6, allo_subset_6, by = c("s", "Dominance"), suffix = c("_auto", "_allo")) %>%
  mutate(
    q_diff = q_auto - q_allo,
    load_diff = load_auto - load_allo,
    mu = "1e-6"
  )

# Combine all mutation rate data
all_diff_data <- bind_rows(diff_data_8, diff_data_7, diff_data_6)
# Factor mutation rate for consistent ordering
all_diff_data$mu <- factor(all_diff_data$mu, levels = c("1e-8", "1e-7", "1e-6"))
# Define colors for mutation rates
mu_colors <- c("1e-8" = "#000000",  
               "1e-7" = "#000000",     
               "1e-6" = "#000000")   
mu_labeller <- c(
  "1e-8" = "mu == 1e-8",
  "1e-7" = "mu == 1e-7", 
  "1e-6" = "mu == 1e-6"
)

# FIRST FIGURE: Q differences with mutation rate as rows
p_q_diff <- ggplot(all_diff_data, aes(x = s, y = q_diff, color = mu)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  facet_grid(mu ~ Dominance, 
             scales = "free_y",  # Added free y-scales
             labeller = labeller(mu = as_labeller(mu_labeller, default = label_parsed))) +
  scale_color_manual(values = mu_colors, name = expression("Mutation Rate (" * mu * ")")) +
  scale_x_log10(
    limits = c(1e-9, 1e-3),
    breaks = c(1e-8, 1e-6, 1e-4),
    labels = scales::label_scientific()
  ) +
  scale_y_continuous(
    transform = "pseudo_log",
    breaks = scales::pretty_breaks(n = 4),
    labels = scales::label_scientific()
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  ) +
  labs(
    x = expression(s ~ "(Selection Coefficient)"),
    y = expression(Delta ~ "q (Auto - Allo)"),
    color = expression("Mutation Rate (" * mu * ")")
  )

# SECOND FIGURE: Load differences with mutation rate as rows  
p_load_diff <- ggplot(all_diff_data, aes(x = s, y = load_diff, color = mu)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  facet_grid(mu ~ Dominance,
             scales = "free_y",  # Added free y-scales
             labeller = labeller(mu = as_labeller(mu_labeller, default = label_parsed))) +
  scale_color_manual(values = mu_colors, name = expression("Mutation Rate (" * mu * ")")) +
  scale_x_log10(
    limits = c(1e-9, 1e-3),
    breaks = c(1e-8, 1e-6, 1e-4),
    labels = scales::label_scientific()
  ) +
  scale_y_continuous(
    transform = "pseudo_log",
    breaks = scales::pretty_breaks(n = 4),
    labels = scales::label_scientific()
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)
  ) +
  labs(
    x = expression(s ~ "(Selection Coefficient)"),
    y = expression(Delta ~ "Load (Auto - Allo)"),
    color = expression("Mutation Rate (" * mu * ")")
  )

# Save the plots
ggsave("q_differences_by_mutation_rate.pdf", p_q_diff, width = 6.5, height = 5, units = "in")
ggsave("load_differences_by_mutation_rate.pdf", p_load_diff, width = 6.5, height = 5, units = "in")


### =====================================================================
### THIRD PLOT: PD AND GENOTYPE DISTRIBUTIONS TO EXPLAIN THE BUMP IN LOAD
### =====================================================================

# plot of PD and genotype distributions to explain the bump in load
dominant_data <- full_data %>%
  filter(Dominance == "Dominant", Model %in% c("Autotetraploid", "Allotetraploid")) %>%
  select(s, Model, PD, G0, G1, G2, G3, G4) %>%
  pivot_longer(cols = c(PD, G0, G1, G2, G3, G4), 
               names_to = "Variable", 
               values_to = "Value") %>%
  mutate(Variable = factor(Variable, levels = c("PD", "G0", "G1", "G2", "G3", "G4")))

# to change PD to have \Delta
variable_labeller <- c("PD" = "Delta~~(PD)", "G0" = "G0", "G1" = "G1", "G2" = "G2", "G3" = "G3", "G4" = "G4")

p_load_bump <- ggplot(dominant_data, aes(x = s, y = Value, color = Model)) +
  geom_line(linewidth=.75) +
  facet_wrap(~Variable, scales = "free_y", ncol = 1, strip.position = "left",
             labeller = as_labeller(variable_labeller, default = label_parsed)) +
  scale_color_manual(values = c("Autotetraploid" = auto_color, "Allotetraploid" = allo_color)) +
  scale_x_log10(
    limits = c(1e-9, 1e-3),
    breaks = c(1e-8, 1e-6, 1e-4),
    labels = scales::label_scientific()
  ) +
  scale_y_continuous(labels = scales::label_scientific()) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "bottom"
  ) +
  labs(x = "s (Selection Coefficient)", y = NULL, color = "Model")

ggsave("load_bump.pdf", p_load_bump, width = 4.5, height = 7, units = "in")

