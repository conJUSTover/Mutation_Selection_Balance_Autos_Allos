library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(patchwork)

# Define color mapping
dip_color <- "#F04D13"
auto_color <- "#66BED6"
allo_color <- "#7DAB5B"

# Read in all datasets
dip_rec <- read_csv("dip_rec.csv", col_names = FALSE)
dip_add <- read_csv("dip_add.csv", col_names = FALSE)
dip_dom <- read_csv("dip_dom.csv", col_names = FALSE)
auto_rec <- read_csv("auto_rec.csv", col_names = FALSE)
auto_add <- read_csv("auto_add.csv", col_names = FALSE)
auto_dom <- read_csv("auto_dom.csv", col_names = FALSE)

# Add identifiers and stack all data
add_labels <- function(df, type, mode, group) {
  df <- df %>%
    mutate(Model = type,
           Mode = mode,
           Group = group)
  return(df)
}

combined <- bind_rows(
  add_labels(dip_rec,  "Diploid",    "Recessive", "AlleleFreq"),
  add_labels(auto_rec, "Auto- and Allotetraploid", "Recessive", "AlleleFreq"),
  add_labels(dip_add,  "Diploid",    "Additive",  "AlleleFreq"),
  add_labels(auto_add, "Auto- and Allotetraploid", "Additive",  "AlleleFreq"),
  add_labels(dip_dom,  "Diploid",    "Dominant",  "AlleleFreq"),
  add_labels(auto_dom, "Auto- and Allotetraploid", "Dominant",  "AlleleFreq"),
  
  add_labels(dip_rec,  "Diploid",    "Recessive", "MutationLoad"),
  add_labels(auto_rec, "Auto- and Allotetraploid", "Recessive", "MutationLoad"),
  add_labels(dip_add,  "Diploid",    "Additive",  "MutationLoad"),
  add_labels(auto_add, "Auto- and Allotetraploid", "Additive",  "MutationLoad"),
  add_labels(dip_dom,  "Diploid",    "Dominant",  "MutationLoad"),
  add_labels(auto_dom, "Auto- and Allotetraploid", "Dominant",  "MutationLoad")
) %>%
  mutate(
    x = X1,
    y = case_when(
      Group == "AlleleFreq" ~ X2,
      Group == "MutationLoad" & Model == "Diploid" ~ 1 - X5,
      Group == "MutationLoad" & Model == "Auto- and Allotetraploid" ~ 1 - X6
    )
  )

# Factor for better facet labeling
combined$Mode <- factor(combined$Mode, levels = c("Recessive", "Additive", "Dominant"))
combined$Group <- factor(combined$Group, levels = c("AlleleFreq", "MutationLoad"))

allele_freq_plot <- combined %>%
  filter(Group == "AlleleFreq") %>%
  ggplot(aes(x = x, y = y, color = Model)) +
  geom_line(size = .75) +
  facet_grid(. ~ Mode, labeller = label_value) +
  scale_x_log10(
    limits = c(1e-9, 1e-3),
    breaks = c(1e-8, 1e-6, 1e-4)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  scale_color_manual(values = c("Diploid" = dip_color, "Auto- and Allotetraploid" = auto_color)) +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = "Allele Frequency (q)", color = "Ploidy") +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


mutation_load_plot <- combined %>%
  filter(Group == "MutationLoad") %>%
  ggplot(aes(x = x, y = y, color = Model)) +
  geom_line(size = .75) +
  facet_grid(. ~ Mode, labeller = label_value) +
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
  scale_color_manual(values = c("Diploid" = dip_color, "Auto- and Allotetraploid" = auto_color)) +
  theme_bw(base_size = 11) +
  labs(x = "s (Selection Coefficient)", y = "Mutation Load", color = "Ploidy") +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# Create a data frame with tags and coordinates
panel_labels <- data.frame(
  Mode = factor(c("Recessive", "Additive", "Dominant"), levels = c("Recessive", "Additive", "Dominant")),
  label_top = c("a", "b", "c"),
  label_bottom = c("d", "e", "f")
)

# Add top row labels to allele_freq_plot
allele_freq_plot <- allele_freq_plot +
  geom_text(data = panel_labels, aes(label = label_top, x = 5e-4, y = 0.95),
            inherit.aes = FALSE, size = 4, hjust = 0)

# Add bottom row labels to mutation_load_plot
mutation_load_plot <- mutation_load_plot +
  geom_text(data = panel_labels, aes(label = label_bottom, x = 5e-4, y = 4.75e-8),
            inherit.aes = FALSE, size = 4, hjust = 0)


final_plot <- (allele_freq_plot / mutation_load_plot) +
  plot_layout(heights = c(1, 1))

ggsave("equal_mut_bifn.pdf", final_plot, width = 6.5, height = 4.25, units = "in")

