#!/usr/bin/env Rscript

library(argparse)
library(ggplot2)
library(scales)

# ------------------------------------------------------------
# Arguments
# ------------------------------------------------------------

parser <- ArgumentParser(
    description = "Plot a k-mer spectrum in a GenomeScope-like style"
)

parser$add_argument(
    "input",
    help = "Input k-mer histogram file: Depth NumberKmers"
)

parser$add_argument(
    "output",
    help = "Output file (e.g. spectrum.pdf or spectrum.png)"
)

parser$add_argument(
    "title",
    help = "Plot title"
)

args <- parser$parse_args()

# Settings for the heuristic detector (not extra command-line arguments).
max_plot_depth <- 250
smooth_k <- 7
persistence <- 4
flatten_fraction <- 0.20


# ------------------------------------------------------------
# Read histogram
# ------------------------------------------------------------

data <- read.table(
    args$input,
    header = FALSE,
    col.names = c("Depth", "NumberKmers")
)

# Remove invalid entries
data <- data[
    is.finite(data$Depth) &
    is.finite(data$NumberKmers) &
    data$Depth >= 1 &
    data$NumberKmers > 0,
]
data <- data[order(data$Depth), ]
if (nrow(data) < 2L || anyDuplicated(data$Depth)) {
  stop("The histogram needs at least two positive-count rows with unique depths")
}
data$WeightedKmers <- data$Depth * data$NumberKmers
if (any(!is.finite(data$WeightedKmers))) {
  stop("Depth * NumberKmers produced non-finite values")
}

# Detect the first sustained rise (valley) or, if absent, the first
# sustained flattening of the steep low-depth decline (shoulder).
# Both detectors use the multiplicity-weighted spectrum. A shoulder
# remains a plotting heuristic, not a validated biological cutoff.
first_run <- function(flag, length_needed) {
  flag[is.na(flag)] <- FALSE
  runs <- rle(flag)
  ends <- cumsum(runs$lengths)
  hits <- which(runs$values & runs$lengths >= length_needed)
  if (!length(hits)) return(NA_integer_)
  as.integer(ends[hits[1]] - runs$lengths[hits[1]] + 1L)
}

detect_error_cutoff <- function(hist, max_depth, k = 7L,
                                persistence = 4L, flatten_fraction = 0.20) {
  d <- hist[hist$Depth >= 2 & hist$Depth <= max_depth, ]
  d$SmoothLog <- log10(d$WeightedKmers)
  if (nrow(d) < 2L * persistence + 3L) {
    return(list(depth = NA_real_, method = "not detected", smoothed = d))
  }

  # The counts span several orders of magnitude, so smooth their logarithms.
  k <- min(k, if (nrow(d) %% 2L) nrow(d) else nrow(d) - 1L)
  d$SmoothLog <- runmed(log10(d$WeightedKmers), k = k, endrule = "median")

  # Centered slopes over several depth bins reduce sensitivity to noise.
  span <- min(3L, floor((nrow(d) - 1L) / 2L))
  slope <- rep(NA_real_, nrow(d))
  middle <- seq.int(1L + span, nrow(d) - span)
  slope[middle] <- (
    d$SmoothLog[middle + span] - d$SmoothLog[middle - span]
  ) / (d$Depth[middle + span] - d$Depth[middle - span])
  # Weighting can make the first few error bins rise. Begin with the
  # first persistent decline, rather than assuming depth 2 is its start.
  decline_start <- first_run(slope < 0, persistence)
  if (is.na(decline_start)) {
    return(list(depth = NA_real_, method = "not detected", smoothed = d))
  }
  decline_slope <- slope[seq.int(decline_start, nrow(d))]
  initial_slope <- median(head(decline_slope[is.finite(decline_slope)], 5L))
  if (!is.finite(initial_slope) || initial_slope >= 0) {
    return(list(depth = NA_real_, method = "not detected", smoothed = d))
  }

  # The first sustained rise must also gain at least 15% in weighted count.
  rising <- slope > 0
  rising[is.na(rising)] <- FALSE
  rising[seq_len(min(nrow(d), decline_start + persistence - 1L))] <- FALSE
  rising_runs <- rle(rising)
  rising_ends <- cumsum(rising_runs$lengths)
  rising_starts <- rising_ends - rising_runs$lengths + 1L
  valid_rises <- which(rising_runs$values &
                       rising_runs$lengths >= persistence)
  for (run in valid_rises) {
    valley_start <- rising_starts[run]
    local <- seq.int(max(1L, valley_start - span),
                     min(nrow(d), valley_start + span))
    valley_idx <- local[which.min(d$SmoothLog[local])]
    end_idx <- min(nrow(d), valley_start + persistence + span)
    rise <- max(d$SmoothLog[valley_idx:end_idx]) - d$SmoothLog[valley_idx]
    if (rise >= log10(1.15)) {
      return(list(depth = d$Depth[valley_idx], method = "valley",
                  smoothed = d))
    }
  }

  flat <- slope > initial_slope * flatten_fraction
  flat[seq_len(min(nrow(d), decline_start + persistence - 1L))] <- FALSE
  shoulder_idx <- first_run(flat, persistence)
  if (!is.na(shoulder_idx)) {
    return(list(depth = d$Depth[shoulder_idx], method = "shoulder",
                smoothed = d))
  }
  list(depth = NA_real_, method = "not detected", smoothed = d)
}

xmax <- min(max_plot_depth, max(data$Depth))
cutoff <- detect_error_cutoff(data, xmax, smooth_k,
                              persistence, flatten_fraction)

# Search for the first sustained rise followed by a sustained decline
# in the weighted spectrum, after the estimated error cutoff.
reference <- cutoff$smoothed
if (!is.na(cutoff$depth)) {
  reference <- reference[reference$Depth > cutoff$depth, ]
}
reference <- reference[reference$Depth >= 5, ]
if (!nrow(reference)) {
  reference <- data[data$Depth <= xmax, ]
  reference$SmoothLog <- log10(reference$WeightedKmers)
}
peak_idx <- NA_integer_
if (nrow(reference) >= 2L * persistence + 3L) {
  span <- 3L
  slope <- rep(NA_real_, nrow(reference))
  middle <- seq.int(1L + span, nrow(reference) - span)
  slope[middle] <- (
    reference$SmoothLog[middle + span] -
      reference$SmoothLog[middle - span]
  ) / (reference$Depth[middle + span] - reference$Depth[middle - span])

  run_starts <- function(flag, n) {
    flag[is.na(flag)] <- FALSE
    runs <- rle(flag)
    ends <- cumsum(runs$lengths)
    ends[which(runs$values & runs$lengths >= n)] -
      runs$lengths[which(runs$values & runs$lengths >= n)] + 1L
  }
  rises <- run_starts(slope > 0, persistence)
  declines <- run_starts(slope < 0, persistence)
  for (rise in rises) {
    for (decline in declines[declines >= rise + persistence]) {
      window <- seq.int(rise, decline)
      candidate <- window[which.max(reference$SmoothLog[window])]
      prominence <- reference$SmoothLog[candidate] -
        min(reference$SmoothLog[seq_len(candidate)])
      if (prominence >= log10(1.15)) {
        peak_idx <- candidate
        break
      }
    }
    if (!is.na(peak_idx)) break
  }
}
if (is.na(peak_idx)) {
  peak_idx <- which.max(reference$SmoothLog)
  reference_method <- "weighted visible maximum / no clear peak"
  # Without a peak, keep the original-count plot readable even if
  # the weighted maximum lies far to the right.
  original_y_reference <- max(reference$NumberKmers)
} else {
  reference_method <- "first weighted post-cutoff peak"
  original_y_reference <- reference$NumberKmers[peak_idx]
}
reference_depth <- reference$Depth[peak_idx]
weighted_peak_y <- 10^reference$SmoothLog[peak_idx]
ymin <- 10^(floor(log10(original_y_reference)) - 1)
ymax <- original_y_reference * 1.5
if (is.na(cutoff$depth)) {
  message("Estimated error k-mer cutoff: not detected")
} else {
  cat("Estimated error k-mer cutoff: ", cutoff$depth,
      " (", cutoff$method, ")\n", sep = "")
}
message("Weighted peak reference: depth = ", reference_depth,
        ", weighted count = ", signif(weighted_peak_y, 4),
        " (", reference_method, ")")

# Exclude counts below the visible log baseline instead of drawing
# inverted vertical lines from the baseline to those values.
plot_data <- data[data$NumberKmers >= ymin, ]
# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p <- ggplot(plot_data, aes(x = Depth, y = NumberKmers)) +

    # GenomeScope-like vertical spectrum
    geom_linerange(
        aes(ymin = ymin, ymax = NumberKmers),
        linewidth = 1,
        colour = "#56B4E9"
    ) +
  
    scale_x_continuous(
        name = "k-mer multiplicity",
        expand = expansion(mult = c(0, 0.01)),
        breaks = breaks_pretty(n = 10)
    ) +

    scale_y_log10(
        name = "Number of distinct k-mers",
        labels = scales::label_scientific(),
        expand = expansion(mult = c(0, 0.05))
    ) +

    coord_cartesian(
        xlim = c(1, xmax),
        ylim = c(ymin, ymax)
    ) +

    labs(
        title = args$title
    ) +

    theme_bw(base_size = 14) +

    theme(
        # GenomeScope-like background
        panel.background = element_rect(
            fill = "grey90",
            colour = NA
        ),

        # No grid
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        # Black border around plotting area
        panel.border = element_rect(
            colour = "black",
            fill = NA,
            linewidth = 0.8
        ),

        # Axis appearance
        axis.text = element_text(
            colour = "black"
        ),

        axis.title = element_text(
            size = 14
        ),

        # Title
        plot.title = element_text(
            hjust = 0.5,
            face = "bold",
            size = 16
        ),

        plot.margin = margin(10, 15, 10, 10)
    )

if (!is.na(cutoff$depth)) {
  p <- p +
    geom_vline(xintercept = cutoff$depth,
               linetype = "dotted", linewidth = 0.8) +
    annotate("text", x = cutoff$depth, y = ymax / 1.3,
             label = paste0("Estimated cutoff = ", cutoff$depth,
                            " (", cutoff$method, ")"),
             angle = 90, vjust = -0.5, hjust = 1, size = 4)
}
# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

ggsave(
    filename = args$output,
    plot = p,
    width = 8,
    height = 6,
    units = "in",
    dpi = 300
)

# A second, multiplicity-weighted spectrum makes higher-copy classes
# easier to see. Both the cutoff and peak above were detected on this
# weighted scale; the original plot retains distinct k-mer counts.
weighted_data <- data[data$Depth <= xmax, ]

# Do not truncate the top of this diagnostic view at the first peak.
weighted_ymax <- max(weighted_data$WeightedKmers) * 1.1
weighted_ymin <- 10^(floor(log10(max(weighted_data$WeightedKmers))) - 2)
weighted_data <- weighted_data[
    weighted_data$WeightedKmers >= weighted_ymin,
]

p_weighted <- ggplot(
    weighted_data, aes(x = Depth, y = WeightedKmers)
) +
    geom_linerange(
        aes(ymin = weighted_ymin, ymax = WeightedKmers),
        linewidth = 1, colour = "#56B4E9"
    ) +
    scale_x_continuous(
        name = "k-mer multiplicity",
        expand = expansion(mult = c(0, 0.01)),
        breaks = breaks_pretty(n = 10)
    ) +
    scale_y_log10(
        name = "k-mer instances (depth × distinct k-mers)",
        labels = scales::label_scientific(),
        expand = expansion(mult = c(0, 0.05))
    ) +
    coord_cartesian(
        xlim = c(1, xmax),
        ylim = c(weighted_ymin, weighted_ymax)
    ) +
    labs(title = paste0(args$title, " — weighted spectrum")) +
    p$theme

if (!is.na(cutoff$depth)) {
    p_weighted <- p_weighted +
        geom_vline(
            xintercept = cutoff$depth,
            linetype = "dotted", linewidth = 0.8
        ) +
        annotate(
            "text", x = cutoff$depth, y = weighted_ymax / 1.3,
            label = paste0("Estimated cutoff = ", cutoff$depth,
                           " (", cutoff$method, ")"),
            angle = 90, vjust = -0.5, hjust = 1, size = 4
        )
}

if (reference_method == "first weighted post-cutoff peak") {
    p_weighted <- p_weighted +
        geom_vline(
            xintercept = reference_depth,
            linetype = "dashed", linewidth = 0.8,
            colour = "#D55E00"
        ) +
        annotate(
            "text", x = reference_depth, y = weighted_ymax / 3,
            label = paste0("Weighted peak = ", reference_depth),
            angle = 90, vjust = -0.5, hjust = 1, size = 4,
            colour = "#D55E00"
        )
}

extension <- tools::file_ext(args$output)
if (!nzchar(extension)) {
    stop("The output filename needs an extension (e.g. .pdf or .png)")
}
weighted_output <- file.path(
    dirname(args$output),
    paste0(tools::file_path_sans_ext(basename(args$output)),
           "_weighted.", extension)
)
ggsave(
    filename = weighted_output,
    plot = p_weighted,
    width = 8,
    height = 6,
    units = "in",
    dpi = 300
)
message("Weighted spectrum saved to ", weighted_output)
