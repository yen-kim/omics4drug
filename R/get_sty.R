#' Calculate and Visualize Phosphorylation Site (S/T/Y) Distribution
#'
#' This function calculates the absolute counts and percentages of phosphorylation
#' sites (Serine (S), Threonine (T), Tyrosine (Y)) for each sample, based on a
#' user-defined phosphorylation site probability cutoff. It can then optionally
#' generate bar plots to visualize these counts or percentages across different
#' sample groups.
#'
#' @param data A data frame containing phosphoproteomics data. This data frame
#'   is expected to have columns for sample names, phosphorylation site
#'   probabilities, amino acid type of the phosphorylation site, and a grouping
#'   variable.
#' @param plot A character string specifying the type of plot to generate.
#'   \itemize{
#'     \item "no": No plot is generated; only the summary table is returned.
#'     \item "count": Generates bar plots showing the absolute counts of S, T, and Y
#'       phosphosites for each group, faceted by amino acid type.
#'     \item "percentage": Generates bar plots showing the percentage distribution
#'       of S, T, and Y phosphosites within each group, faceted by amino acid type.
#'   }
#' @param sample_name [`data-masking`][rlang::args_data_masking]> The unquoted name
#'   of the column in `data` that contains unique sample identifiers.
#'   e.g., `SampleID`.
#' @param group <[`data-masking`][rlang::args_data_masking]> The unquoted name of
#'   the column in `data` that defines the experimental groups (e.g., "control", "treat").
#' @param ptm_site_prob_col <[`data-masking`][rlang::args_data_masking]> The unquoted
#'   name of the column in `data` that contains phosphorylation site probabilities
#'   (e.g., "Localization_Probabilities").
#' @param ptm_site_prob_val Numeric, a cutoff value for `ptm_site_prob_col`. Only
#'   sites with a probability greater than or equal to this value will be included
#'   in the analysis. Default is `0.75`.
#' @param ptm_site_amino_acid <[`data-masking`][rlang::args_data_masking]> The unquoted
#'   name of the column in `data` that contains the amino acid type of the
#'   phosphorylation site (expected to be "S", "T", or "Y").
#'   e.g., `Amino_Acid`.
#' @param selected_group A character vector specifying which groups from the `group`
#'   column should be included in the plots. If omitted or `NULL`, all groups will be plotted.
#'   This parameter is ignored if `plot = "no"`.
#'
#' @returns #' If `plot = "no"`, returns a `data.frame` (tibble) summarizing the counts and
#' percentages of S, T, Y phosphosites per sample.
#' If `plot = "count"` or `plot = "percentage"`, returns a `list` containing:
#' \itemize{
#'   \item `summary_table`: The `data.frame` of S/T/Y counts and percentages.
#'   \item `plot`: A `ggplot` object of the generated bar plot.
#' }
#' @export
#'
#'
#'
#' @examples
#' # Create data for demonstration
#' set.seed(123)
#' phospho_data <- data.frame(
#'   SampleID = rep(c("S1_ctrl", "S2_ctrl", "S3_treat", "S4_treat"), each = 25),
#'   Group = rep(c("Control", "Control", "Treatment", "Treatment"), each = 25),
#'   Protein = paste0("Prot", 1:100),
#'   Localization_Probabilities = runif(100, 0.5, 0.99),
#'   Amino_Acid = sample(c("S", "T", "Y"), 100, replace = TRUE, prob = c(0.7, 0.2, 0.1))
#' )
#'
#' # Case 1: Calculate S/T/Y distribution and return data only
#' sty_counts_df <- get_sty(
#'   data = phospho_data,
#'   plot = "no",
#'   sample_name = SampleID,
#'   group = Group,
#'   ptm_site_prob_col = Localization_Probabilities,
#'   ptm_site_prob_val = 0.75,
#'   ptm_site_amino_acid = Amino_Acid
#' )
#' print(sty_counts_df)
#'
#' # Case 2: Calculate S/T/Y distribution and generate count plots
#' # Define selected groups for plotting
#' selected_groups_for_plot <- c("Control", "Treatment")
#'
#' sty_counts_plot <- get_sty(
#'   data = phospho_data,
#'   plot = "count",
#'   sample_name = SampleID,
#'   group = Group,
#'   ptm_site_prob_col = Localization_Probabilities,
#'   ptm_site_prob_val = 0.75,
#'   ptm_site_amino_acid = Amino_Acid,
#'   selected_group = selected_groups_for_plot
#' )
#' print(sty_counts_plot$summary_table)
#' sty_counts_plot$plot # Display the count plot
#'
#' # Case 3: Calculate S/T/Y distribution and generate percentage plots
#' sty_percentage_plot <- get_sty(
#'   data = phospho_data,
#'   plot = "percentage",
#'   sample_name = SampleID,
#'   group = Group,
#'   ptm_site_prob_col = Localization_Probabilities,
#'   ptm_site_prob_val = 0.75,
#'   ptm_site_amino_acid = Amino_Acid,
#'   selected_group = selected_groups_for_plot
#' )
#' print(sty_percentage_plot$summary_table)
#' sty_percentage_plot$plot # Display the percentage plot
#'

get_sty <- function(data,
                         plot = c("no", "count", "percentage"),
                         sample_name,
                         group,
                         ptm_site_prob_col,
                         ptm_site_prob_val = 0.75,
                         ptm_site_amino_acid,
                         selected_group = NULL){

        plot <- match.arg(plot)

        # Initial calculation of S/T/Y counts and percentages per sample
        df <- data %>%
                dplyr::group_by({{sample_name}}, {{group}}) %>% # Group by sample_name and group
                dplyr::filter({{ptm_site_prob_col}} >= ptm_site_prob_val) %>%
                dplyr::count({{ptm_site_amino_acid}}, name = "count") %>% # Count S,T,Y
                tidyr::pivot_wider(names_from = {{ptm_site_amino_acid}}, values_from = count, values_fill = list(count = 0)) %>% # Fill NA with 0
                dplyr::ungroup() %>%
                # Ensure S, T, Y columns exist, create if not (e.g., if no Y sites in any sample)
                dplyr::mutate(
                        S = dplyr::coalesce(S, 0),
                        T = dplyr::coalesce(T, 0),
                        Y = dplyr::coalesce(Y, 0),
                        STY = S + T + Y,
                        pct_S = S / STY * 100,
                        pct_T = T / STY * 100,
                        pct_Y = Y / STY * 100
                ) %>%
                # Handle cases where STY might be 0 to prevent NaN for percentages
                dplyr::mutate(
                        pct_S = ifelse(STY == 0, 0, pct_S),
                        pct_T = ifelse(STY == 0, 0, pct_T),
                        pct_Y = ifelse(STY == 0, 0, pct_Y)
                )


        # Return summary table if no plot requested
        if(plot == "no") {
                return(df)
        }

        # Common data preparation for plotting
        plot_data_long <- df %>%
                dplyr::select({{sample_name}}, S, T, Y, pct_S, pct_T, pct_Y, {{group}}) %>%
                tidyr::pivot_longer(
                        cols = c(S, T, Y, pct_S, pct_T, pct_Y),
                        names_to = "amino_acid_type",
                        values_to = "value"
                ) %>%
                dplyr::mutate(
                        plot_group = paste0(rlang::as_label(rlang::enquo(group)), "-", amino_acid_type), # Use the actual group name
                        amino_acid_facet = dplyr::case_when(
                                amino_acid_type %in% c("S", "pct_S") ~ "Serine (S)",
                                amino_acid_type %in% c("T", "pct_T") ~ "Threonine (T)",
                                amino_acid_type %in% c("Y", "pct_Y") ~ "Tyrosine (Y)",
                                TRUE ~ NA_character_
                        ),
                        amino_acid_facet = factor(amino_acid_facet, levels = c("Serine (S)", "Threonine (T)", "Tyrosine (Y)"))
                )

        # Filter by selected_group if provided
        if (!is.null(selected_group)) {
                plot_data_long <- plot_data_long %>%
                        dplyr::filter({{group}} %in% selected_group)
        }

        # Summarize data for plotting (mean values per group and amino acid type)
        summary_plot_data <- plot_data_long %>%
                dplyr::group_by({{group}}, amino_acid_type, amino_acid_facet) %>%
                dplyr::summarise(mean_value = round(mean(.data$value, na.rm = TRUE), digits = 2), .groups = "drop") %>%
                dplyr::ungroup()


        # Generate plot based on 'plot' argument
        plot_obj <- NULL
        if (plot == "count") {
                plot_obj <- plot_data_long %>%
                        dplyr::filter(amino_acid_type %in% c("S", "T", "Y")) %>%
                        ggplot2::ggplot(ggplot2::aes(x = {{group}}, y = .data$value)) +
                        ggplot2::geom_col(
                                data = summary_plot_data %>% dplyr::filter(amino_acid_type %in% c("S", "T", "Y")),
                                ggplot2::aes(y = .data$mean_value, fill = {{group}}),
                                position = ggplot2::position_stack(vjust = 0.5), # Adjust position to center labels
                                color = "black", linewidth = 0.2
                        ) +
                        ggplot2::geom_point(ggplot2::aes(color = {{group}}),
                                            position = ggplot2::position_jitter(width = 0.1), alpha = 0.6, size = 2) +
                        ggplot2::geom_text(
                                data = summary_plot_data %>% dplyr::filter(amino_acid_type %in% c("S", "T", "Y")),
                                ggplot2::aes(y = .data$mean_value, label = .data$mean_value), # Label with mean values
                                vjust = -0.5, size = 3, color = "black"
                        ) +
                        ggplot2::labs(
                                x = "Group",
                                y = "Number of Phosphosites",
                                title = "Absolute Counts of S/T/Y Phosphosites"
                        ) +
                        ggplot2::facet_wrap(~ amino_acid_facet, scales = 'free_y') +
                        ggprism::theme_prism() +
                        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                                       legend.position = "top",
                                       axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

        } else if (plot == "percentage") {
                plot_obj <- plot_data_long %>%
                        dplyr::filter(amino_acid_type %in% c("pct_S", "pct_T", "pct_Y")) %>%
                        ggplot2::ggplot(ggplot2::aes(x = {{group}}, y = .data$value)) +
                        ggplot2::geom_col(
                                data = summary_plot_data %>% dplyr::filter(amino_acid_type %in% c("pct_S", "pct_T", "pct_Y")),
                                ggplot2::aes(y = .data$mean_value, fill = {{group}}),
                                position = ggplot2::position_stack(vjust = 0.5), # Adjust position to center labels
                                color = "black", linewidth = 0.2
                        ) +
                        ggplot2::geom_point(ggplot2::aes(color = {{group}}),
                                            position = ggplot2::position_jitter(width = 0.1), alpha = 0.6, size = 2) +
                        ggplot2::geom_text(
                                data = summary_plot_data %>% dplyr::filter(amino_acid_type %in% c("pct_S", "pct_T", "pct_Y")),
                                ggplot2::aes(y = .data$mean_value, label = paste0(.data$mean_value, "%")), # Label with mean percentages
                                vjust = -0.5, size = 3, color = "black"
                        ) +
                        ggplot2::labs(
                                x = "Group",
                                y = "Percentage of Phosphosites (%)",
                                title = "Percentage Distribution of S/T/Y Phosphosites"
                        ) +
                        ggplot2::facet_wrap(~ amino_acid_facet, scales = 'free_y') +
                        ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%")) + # Format y-axis as percentage
                        ggprism::theme_prism() +
                        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                                       legend.position = "top",
                                       axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
        }

        return(list(summary_table = df, plot = plot_obj))
}
