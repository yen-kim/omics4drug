#' Count Phosphosites Based on Probability Threshold
#'
#' This function processes a data frame to count unique phosphosites based on
#' a specified probability threshold. It groups the data by a "Group"
#' extracted from the `Sample` column, calculates the average phosphosite count
#' for each group, and can optionally generate a bar plot to visualize the results.
#'
#' The generated plot shows the mean count for each "Group" as a bar and
#' includes jittered points representing the individual sample counts.
#'
#' @param data A data frame containing the data. It must have columns named
#'   `Sample`, `phosphosites`, and `Probability`. The `Sample` column is
#'   expected to be in the format 'sampleID_groupName' (e.g., 'p1_control',
#'   'p3_dose-1').
#' @param prob_threshold A numeric value between 0 and 1. This is the
#'   probability cutoff; the function will only consider phosphosites with
#'   a `Probability` value greater than or equal to this threshold.
#' @param plot A character string. Specifies whether to generate a plot.
#'   Accepted values are `"bar_plot"` to generate the plot or `"no"` to
#'   return only the summary table. Defaults to `"bar_plot"`.
#'
#' @returns If `plot` is set to `"no"`, the function returns a `tibble`
#'   summarizing the mean `Number_of_phosphosites` for each `Group`.
#'   If `plot` is `"bar_plot"`, the function returns a `list` with two
#'   elements:
#'   \itemize{
#'     \item `table`: The summary `tibble` with mean counts per group.
#'     \item `plot`: A `ggplot` object of the bar plot, ready for printing.
#'   }
#'   The function will stop and return an error if the required columns are
#'   missing or the `prob_threshold` is not a valid number.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a data frame for demonstration
#' data <- data.frame(
#'   Sample = c("p1_control", "p2_control", "p3_dose-1", "p4_dose-1", "p5_dose-2", "p6_dose-2"),
#'   phosphosites = c("amino_acid_S1", "amino_acid_T2", "amino_acid_Y4",
#'                    "amino_acid_S1", "amino_acid_T2", "amino_acid_S1"),
#'   Probability = c(0.95, 0.82, 0.65, 0.91, 0.73, 0.99)
#' )
#'
#' # Get the summary table and the bar plot for a probability threshold of 0.75
#' results <- get_count_phosphosite(data, prob_threshold = 0.75,
#' plot = "bar_plot")
#' print(results$table)
#' print(results$plot)
#'
#' # Get only the summary table for a probability threshold of 0.9
#' summary_table <- get_count_phosphosite(data, prob_threshold = 0.9,
#' plot = "no")
#' print(summary_table)
#' }
get_count_phosphosite <- function(data,
                                  prob_threshold,
                                  plot = c("bar_plot", "no")) {
        plot <- match.arg(plot)

        required_cols <- c("Sample", "phosphosites", "Probability")
        if (!all(required_cols %in% colnames(data))) {
                stop(paste("Input data must contain the following columns:",
                           paste(required_cols, collapse = ", ")))
        }

        if (!is.numeric(prob_threshold) || prob_threshold < 0 || prob_threshold > 1) {
                stop("The 'prob_threshold' argument must be a numeric value between 0 and 1.")
        }

        # Filter, extract Group, and summarize
        counted_sites <- data %>%
                dplyr::filter(Probability >= prob_threshold) %>%
                dplyr::group_by(Sample) %>%
                dplyr::summarise(
                        Number_of_phosphosites = dplyr::n_distinct(phosphosites),
                        .groups = "drop") %>%
                dplyr::mutate(Group = sub("^[^_]+_", "", Sample))

        # Add the threshold column to the result
        counted_sites$prob_threshold <- prob_threshold

        if (plot == "no") {
                return(counted_sites)
        }

        # Generate bar plot if requested
        if (plot == "bar_plot") {
                summary_data <- counted_sites %>%
                        dplyr::group_by(Group) %>%
                        dplyr::summarise(Number_of_phosphosites = round(mean(Number_of_phosphosites, 0)))

                my_color <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#E76F51")

                p <- ggplot2::ggplot(data = counted_sites, ggplot2::aes(x = Group,
                                                                        y = Number_of_phosphosites,
                                                                        color = Group,
                                                                        fill = Group,
                                                                        label = Number_of_phosphosites)) +
                        ggplot2::geom_bar(data = summary_data,
                                          stat = 'identity',
                                          alpha = 0.4) +
                        ggplot2::geom_text(data = summary_data,
                                           vjust = 1.5,
                                           fontface = 'bold',
                                           size = 7) +
                        ggplot2::geom_jitter(size = 3) +
                        ggprism::theme_prism(base_size = 20,
                                             base_family = "sans",
                                             base_fontface = "bold",
                                             axis_text_angle = 45)  +
                        ggplot2::scale_color_manual(values = my_color) +
                        ggplot2::scale_fill_manual(values = my_color) +
                        ggplot2::theme(legend.position = 'none') +
                        ggplot2::labs(
                                title = paste("Phosphosite Counts (Probability >= ", prob_threshold, ")"),
                                x = "Group",
                                y = "Number of Phosphosites"
                        )

                return(list(table = counted_sites, plot = p))
        }
}
