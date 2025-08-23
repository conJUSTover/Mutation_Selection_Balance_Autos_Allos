library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(patchwork)

setwd("~/GitHub/Mutation_Selection_Balance_Autos_Allos/Matlab_scripts/deleterious_bifn")


### Recreates Figure 3

# Define color mapping
dip_color <- "#F04D13"
auto_color <- "#66BED6"
allo_color <- "#7DAB5B"

# Read in datasets - Recessive
dip_rec <- read_csv("dip_rec.csv", col_names = c("s", "q", "g0", "g1", "w"))
auto_rec <- read_csv("auto_rec.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
allo_rec <- read_csv("allo_rec.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))

# Read in datasets - Additive
dip_add <- read_csv("dip_add.csv", col_names = c("s", "q", "g0", "g1", "w"))
auto_add <- read_csv("auto_add.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
allo_add <- read_csv("allo_add.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))

# Read in datasets - Dominant (three stability categories)
dip_dom_selected <- read_csv("dip_dom_selected.csv", col_names = c("s", "q", "g0", "g1", "w"))
dip_dom_neutral <- read_csv("dip_dom_neutral.csv", col_names = c("s", "q", "g0", "g1", "w"))
dip_dom_unstable <- read_csv("dip_dom_unstable.csv", col_names = c("s", "q", "g0", "g1", "w"))

auto_dom_selected <- read_csv("auto_dom_selected.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_dom_neutral <- read_csv("auto_dom_neutral.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))
auto_dom_unstable <- read_csv("auto_dom_unstable.csv", col_names = c("s", "q", "g0", "g1", "g2", "w"))

allo_dom_selected <- read_csv("allo_dom_selected.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_dom_neutral <- read_csv("allo_dom_neutral.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))
allo_dom_unstable <- read_csv("allo_dom_unstable.csv", col_names = c("s", "q", "g00", "g01", "g10", "g11", "w"))

# Add dominance and stability labels
dip_rec$Dominance <- "Recessive"
dip_rec$Stability <- "stable"
dip_add$Dominance <- "Additive"
dip_add$Stability <- "stable"
dip_dom_selected$Dominance <- "Dominant"
dip_dom_selected$Stability <- "selected"
dip_dom_neutral$Dominance <- "Dominant"
dip_dom_neutral$Stability <- "neutral"
dip_dom_unstable$Dominance <- "Dominant"
dip_dom_unstable$Stability <- "unstable"

auto_rec$Dominance <- "Recessive"
auto_rec$Stability <- "stable"
auto_add$Dominance <- "Additive"
auto_add$Stability <- "stable"
auto_dom_selected$Dominance <- "Dominant"
auto_dom_selected$Stability <- "selected"
auto_dom_neutral$Dominance <- "Dominant"
auto_dom_neutral$Stability <- "neutral"
auto_dom_unstable$Dominance <- "Dominant"
auto_dom_unstable$Stability <- "unstable"

allo_rec$Dominance <- "Recessive"
allo_rec$Stability <- "stable"
allo_add$Dominance <- "Additive"
allo_add$Stability <- "stable"
allo_dom_selected$Dominance <- "Dominant"
allo_dom_selected$Stability <- "selected"
allo_dom_neutral$Dominance <- "Dominant"
allo_dom_neutral$Stability <- "neutral"
allo_dom_unstable$Dominance <- "Dominant"
allo_dom_unstable$Stability <- "unstable"

# Bind each ploidy group together
dip <- bind_rows(dip_rec, dip_add, dip_dom_selected, dip_dom_neutral, dip_dom_unstable)
auto <- bind_rows(auto_rec, auto_add, auto_dom_selected, auto_dom_neutral, auto_dom_unstable)
allo <- bind_rows(allo_rec, allo_add, allo_dom_selected, allo_dom_neutral, allo_dom_unstable)

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
         Model = "Autotetraploid (beneath allotetraploid)")

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

# Set factor levels
full_data$Dominance <- factor(full_data$Dominance, levels = c("Recessive", "Additive", "Dominant"))
full_data$Model <- factor(full_data$Model, levels = c("Diploid", "Autotetraploid (beneath allotetraploid)", "Allotetraploid"))
full_data$Stability <- factor(full_data$Stability, levels = c("stable", "selected", "neutral", "unstable"))

# Add line type based on stability
full_data <- full_data %>%
  mutate(LineType = ifelse(Stability == "unstable", "dashed", "solid"))

# Color mapping
ploidy_colors <- c("Diploid" = dip_color,
                   "Autotetraploid (beneath allotetraploid)" = auto_color,
                   "Allotetraploid" = allo_color)

# Create data frame for panel labels
panel_labels_top <- data.frame(
  Dominance = factor(c("Recessive", "Additive", "Dominant"),
                     levels = c("Recessive", "Additive", "Dominant")),
  label = c("A", "B", "C"),
  x = rep(3e-4, 3),
  y = rep(0.95, 3)
)

panel_labels_bottom <- data.frame(
  Dominance = factor(c("Recessive", "Additive", "Dominant"),
                     levels = c("Recessive", "Additive", "Dominant")),
  label = c("D", "E", "F"),
  x = rep(3e-4, 3),
  y = rep(1e-5, 3)
)

# Top row: Allele frequency (q) plot
p_q <- ggplot(full_data, aes(x = s, y = q, color = Model, linetype = LineType, 
                             group = interaction(Model, Stability))) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = ploidy_colors) +
  scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed"), guide = "none") +
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

# Bottom row: Mutation load plot (exclude unstable for dominant case)
load_data <- full_data %>%
  filter(!(Dominance == "Dominant" & Stability == "unstable"))

p_load <- ggplot(load_data, aes(x = s, y = load, color = Model, linetype = LineType,
                                group = interaction(Model, Stability))) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = ploidy_colors) +
  scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed"), guide = "none") +
  facet_wrap(~Dominance, nrow = 1) +  
  scale_x_log10(
    limits = c(1e-9, 1e-3),
    breaks = c(1e-8, 1e-6, 1e-4),
    labels = scales::label_scientific()
  ) +
  scale_y_log10(
    limits = c(1e-10, 2.5e-5),
    breaks = c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5),
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
    strip.text.x = element_blank()
  )

# Combine plots
final_plot <- p_q / p_load +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Save
ggsave("deleterious_bifn.pdf", final_plot, width = 6.5, height = 4, units = "in")
