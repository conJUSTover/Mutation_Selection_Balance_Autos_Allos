library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(patchwork)

# Colors
dip_color <- "#F04D13"
auto_color <- "#66BED6"

# Read in datasets
dip_rec <- read_csv("dip_rec.csv", col_names = FALSE)
auto_rec <- read_csv("auto_rec.csv", col_names = FALSE)

dip_add <- read_csv("dip_add.csv", col_names = FALSE)
auto_add <- read_csv("auto_add.csv", col_names = FALSE)

dip_dom_selected <- read_csv("dip_dom_selected.csv", col_names = FALSE)
dip_dom_neutral  <- read_csv("dip_dom_neutral.csv", col_names = FALSE)
dip_dom_unstable <- read_csv("dip_dom_unstable.csv", col_names = FALSE)

auto_dom_selected <- read_csv("auto_dom_selected.csv", col_names = FALSE)
auto_dom_neutral  <- read_csv("auto_dom_neutral.csv", col_names = FALSE)
auto_dom_unstable <- read_csv("auto_dom_unstable.csv", col_names = FALSE)

# ----- Construct data -----

# Helper function
build_df <- function(df, model, mode, group, label = "neutral", y_col = "X2") {
  df %>%
    mutate(
      x = .data[["X1"]],
      y = if (group == "MutationLoad") {
        1 - .data[[y_col]]
      } else {
        .data[[y_col]]
      },
      Model = model,
      Mode = mode,
      Group = group,
      Label = label
    ) %>%
    select(x, y, Model, Mode, Group, Label)
}

# Combine all data
df <- bind_rows(
  # Allele frequencies (column 2)
  build_df(dip_rec, "Diploid", "Recessive", "AlleleFreq", y_col = "X2"),
  build_df(auto_rec, "Auto- and Allotetraploids", "Recessive", "AlleleFreq", y_col = "X2"),
  build_df(dip_add, "Diploid", "Additive", "AlleleFreq", y_col = "X2"),
  build_df(auto_add, "Auto- and Allotetraploids", "Additive", "AlleleFreq", y_col = "X2"),
  build_df(dip_dom_neutral, "Diploid", "Dominant", "AlleleFreq", "neutral", y_col = "X2"),
  build_df(dip_dom_selected, "Diploid", "Dominant", "AlleleFreq", "selected", y_col = "X2"),
  build_df(dip_dom_unstable, "Diploid", "Dominant", "AlleleFreq", "unstable", y_col = "X2"),
  build_df(auto_dom_neutral, "Auto- and Allotetraploids", "Dominant", "AlleleFreq", "neutral", y_col = "X2"),
  build_df(auto_dom_selected, "Auto- and Allotetraploids", "Dominant", "AlleleFreq", "selected", y_col = "X2"),
  build_df(auto_dom_unstable, "Auto- and Allotetraploids", "Dominant", "AlleleFreq", "unstable", y_col = "X2"),
  
  # Mutation load (column 5 = diploid; column 6 = autotetraploid)
  build_df(dip_rec, "Diploid", "Recessive", "MutationLoad", y_col = "X5"),
  build_df(auto_rec, "Auto- and Allotetraploids", "Recessive", "MutationLoad", y_col = "X6"),
  build_df(dip_add, "Diploid", "Additive", "MutationLoad", y_col = "X5"),
  build_df(auto_add, "Auto- and Allotetraploids", "Additive", "MutationLoad", y_col = "X6"),
  build_df(dip_dom_neutral, "Diploid", "Dominant", "MutationLoad", "neutral", y_col = "X5"),
  build_df(dip_dom_selected, "Diploid", "Dominant", "MutationLoad", "selected", y_col = "X5"),
  build_df(auto_dom_neutral, "Auto- and Allotetraploids", "Dominant", "MutationLoad", "neutral", y_col = "X6"),
  build_df(auto_dom_selected, "Auto- and Allotetraploids", "Dominant", "MutationLoad", "selected", y_col = "X6")
)


# Set factor levels
df$Mode <- factor(df$Mode, levels = c("Recessive", "Additive", "Dominant"))
df$Group <- factor(df$Group, levels = c("AlleleFreq", "MutationLoad"))
df$Label <- factor(df$Label, levels = c("neutral", "selected", "unstable"))

# Set line type: dashed for unstable
df <- df %>%
  mutate(LineType = ifelse(Label == "unstable", "dashed", "solid"))

# Plot settings
color_map <- c("Diploid" = dip_color, "Auto- and Allotetraploids" = auto_color)
line_map <- c("solid" = "solid", "dashed" = "dashed")

# --- Top Row: Allele Frequencies ---
allele_freq_plot <- df %>%
  filter(Group == "AlleleFreq") %>%
  ggplot(aes(x = x, y = y, color = Model, linetype = LineType, group = interaction(Model, Label))) +
  geom_line(linewidth = .75) +
  facet_grid(. ~ Mode) +
  scale_x_log10(
    limits = c(1e-9, 1e-3),
    breaks = c(1e-8, 1e-6, 1e-4)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  scale_color_manual(values = color_map) +
  scale_linetype_manual(values = line_map, guide = "none") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank()
  ) +
  labs(y = "Allele Frequency (q)", title = NULL)

# --- Bottom Row: Mutation Load ---
mutation_load_plot <- df %>%
  filter(Group == "MutationLoad") %>%
  ggplot(aes(x = x, y = y, color = Model, linetype = LineType, group = interaction(Model, Label))) +
  geom_line(linewidth = .75) +
  facet_grid(. ~ Mode) +
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
  scale_color_manual(values = color_map) +
  scale_linetype_manual(values = line_map, guide = "none") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.placement = "outside",
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(x = "s (Selection Coefficient)", y = "Mutation Load", linetype = "Regime", color = "Ploidy")

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
  geom_text(data = panel_labels, aes(label = label_bottom, x = 5e-4, y = 1e-5),
            inherit.aes = FALSE, size = 4, hjust = 0)

# --- Combine and save ---
final_plot <- (allele_freq_plot / mutation_load_plot) +
  plot_layout(heights = c(1, 1))

ggsave("deleterious_bifn.pdf", final_plot, width = 6.5, height = 4.25, units = "in")
