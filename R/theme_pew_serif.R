# Pew · serif headline — Simcox Lab house ggplot2 theme (locked 2026-08-12).
# See /root/style/REPORT_STANDARDS.md Part E. Companion to pew_serif.mplstyle.
#
# Fonts: register the house fonts once per session (they live in ~/.fonts):
#   library(systemfonts)            # or library(showtext); showtext_auto()
#   # "Libre Franklin" (body/labels) and "Source Serif 4" (headline) resolve by family name.
# Then: p + theme_pew_serif()  and  + scale_colour_pew() / scale_fill_pew_div() etc.

suppressWarnings(suppressMessages(library(ggplot2)))

pew_ink   <- "#1a1a1a"; pew_dek <- "#6b6b6b"; pew_muted <- "#8a8a8a"
pew_rule  <- "#d6d6d6"; pew_up  <- "#e08214"; pew_down  <- "#2b7f9e"

pew_cat <- c("#e2711d","#2f7f9e","#f2b705","#46708a","#6a8d3f",
             "#a6572f","#7b5aa6","#b0555e","#7d8388")           # distinct series
pew_seq <- c("#fdecd6","#f9cd97","#f0a94f","#e2711d","#b5551a","#7f3c12")   # ordered / intensity
pew_div <- c("#b5551a","#e08a3c","#f2c79a","#f2efe9","#9fc0cc","#4f97ab","#1f6b80") # signed values

theme_pew_serif <- function(base_size = 12,
                            base_family = "Libre Franklin",
                            title_family = "Source Serif 4") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # in-plot title: neutral + descriptive, serif to match the report headline (Part D)
      plot.title      = element_text(family = title_family, face = "bold",
                                     size = rel(1.15), colour = pew_ink, hjust = 0,
                                     margin = margin(b = 6)),
      plot.subtitle   = element_text(colour = pew_dek, size = rel(0.9), hjust = 0,
                                     margin = margin(b = 8)),
      plot.caption    = element_text(colour = pew_muted, size = rel(0.72), hjust = 0),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      axis.title      = element_text(colour = pew_muted, size = rel(0.85)),
      axis.text       = element_text(colour = pew_muted),
      axis.ticks      = element_blank(),
      panel.grid      = element_blank(),      # no gridlines by default (Part E)
      panel.grid.major.y = element_blank(),
      legend.title    = element_text(colour = pew_ink, size = rel(0.85)),
      legend.text     = element_text(colour = pew_ink, size = rel(0.8)),
      legend.position = "right",
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

# scales
scale_colour_pew <- function(...) scale_colour_manual(values = pew_cat, ...)
scale_color_pew  <- scale_colour_pew
scale_fill_pew   <- function(...) scale_fill_manual(values = pew_cat, ...)
scale_fill_pew_seq <- function(...) scale_fill_gradientn(colours = pew_seq, ...)
scale_fill_pew_div <- function(midpoint = 0, ...)
  scale_fill_gradientn(colours = pew_div, rescaler = ~ scales::rescale_mid(.x, mid = midpoint), ...)

# diverging bars (up/down) helper: map a signed value to the house up/down colours
pew_updown <- function(x) ifelse(x >= 0, pew_up, pew_down)
