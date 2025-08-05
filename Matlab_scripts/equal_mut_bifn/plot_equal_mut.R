library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(patchwork)

setwd("~/GitHub/Mutation_Selection_Balance_Autos_Allos/Matlab_scripts/equal_mut_bifn")

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

# Combine all datasets (assuming you’ve already added Dominance and Ploidy columns earlier)
full_data <- bind_rows(dip, auto, allo)

# Set ploidy as a factor for consistent color and legend order
full_data$Model <- factor(full_data$Model, levels = c("Diploid", "Autotetraploid", "Allotetraploid"))

# Color mapping
ploidy_colors <- c("Diploid" = dip_color,
                   "Autotetraploid" = auto_color,
                   "Allotetraploid" = allo_color)

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
