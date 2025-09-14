#' Perform Kinase Substrate Enrichment Analysis (KSEA) and generate visualizations.
#'
#' This function leverages the `KSEAapp` package to perform KSEA, either to generate
#' a kinase-substrate interaction table or a table of kinase activity scores.
#' Optionally, it can produce various plots to visualize the results, including
#' chord diagrams, sankey diagrams, bar plots, or scatter plots.
#'
#' @param data A data frame containing phosphoproteomics data. This typically includes
#'   columns for protein identifiers, phosphorylation sites, and quantitative values
#'   (e.g., log2 fold changes or intensities). The exact format should conform to
#'   the requirements of `KSEAapp::KSEA.KS_table` or `KSEAapp::KSEA.Scores`.
#' @param table_type A character string specifying the type of table to generate.
#'   Must be one of "kinase_substrate" (to get a table of kinase-substrate
#'   interactions) or "kinase_score" (to get a table of inferred kinase
#'   activity scores).
#' @param plot A character string specifying the type of plot to generate.
#'   Options include:
#'   \itemize{
#'     \item "no": No plot is generated, only the table is returned.
#'     \item "chord_diagram": Generates a chord diagram for kinase-substrate interactions
#'       (only applicable when `table_type = "kinase_substrate"`).
#'     \item "sankey_diagram": Generates a Sankey diagram for kinase-substrate interactions
#'       (only applicable when `table_type = "kinase_substrate"`).
#'     \item "bar_plot": Generates a bar plot of kinase Z-scores (only applicable
#'       when `table_type = "kinase_score"`).
#'     \item "scatter_plot": Generates a scatter plot (volcano plot-like) of
#'       kinase Z-scores vs. -log10(p-value) (only applicable when `table_type = "kinase_score"`).
#'   }
#' @param NetworKIN A logical value indicating whether to use NetworKIN scores
#'   for kinase-substrate prediction. Default is `TRUE`.
#' @param NetworKIN.cutoff A numeric value specifying the NetworKIN cutoff score.
#'   Only applicable if `NetworKIN = TRUE`. Default is `1`.
#'
#' @returns
#' A list containing:
#' \itemize{
#'   \item `table`: The generated data frame (either kinase-substrate table or
#'     kinase score table).
#'   \item `plot`: A `ggplot` object or `circlize` plot if a plot type other
#'     than "no" is specified. If `plot = "no"`, only the `table` is returned.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'data_example' is your phosphoproteomics data and 'KSData' is loaded
#' # from KSEAapp.
#' # Example 1: Get kinase-substrate table without a plot
#' ks_table <- get_KSEA(data = data_example, table_type = "kinase_substrate",
#' plot = "no")
#' print(ks_table)
#'
#' # Example 2: Get kinase scores and a bar plot
#' kinase_results <- get_KSEA(data = data_example, table_type = "kinase_score",
#' plot = "bar_plot")
#' print(kinase_results$table)
#' print(kinase_results$plot)
#'
#' # Example 3: Get kinase-substrate table and a chord diagram
#' ks_chord_results <- get_KSEA(data = data_example, table_type = "kinase_substrate",
#' plot = "chord_diagram")
#' # The chord diagram is plotted directly.
#'
#' # Example 4: Get kinase scores and a scatter plot
#' kinase_scatter_results <- get_KSEA(data = data_example, table_type = "kinase_score",
#' plot = "scatter_plot")
#' print(kinase_scatter_results$table)
#' print(kinase_scatter_results$plot)
#' }

get_KSEA <- function(data,
                      table_type = c("kinase_substrate", "kinase_score"),
                      plot = c("no", "chord_diagram", "sankey_diagram", "bar_plot", "scatter_plot"),
                      NetworKIN = TRUE,
                      NetworKIN.cutoff = 1) {

        table_type <- match.arg(table_type)
        plot <- match.arg(plot)

        # Kinase-Substrate Table and Plots
        if (table_type == "kinase_substrate") {
                kinase_substrate_table <- KSEAapp::KSEA.KS_table(
                        KSData,
                        data,
                        NetworKIN = NetworKIN,
                        NetworKIN.cutoff = NetworKIN.cutoff
                )

                if (plot == "no") {
                        return(kinase_substrate_table)
                }

                df <- kinase_substrate_table %>%
                        dplyr::mutate(Kinases = Kinase.Gene, Substrates = paste0(Substrate.Gene, "_", Substrate.Mod))

                if (plot == "chord_diagram") {
                        # Adjacency matrix for chord diagram
                        adjacencydf <- with(df, table(Kinases, Substrates))
                        set.seed(12345) # for reproducibility of layout
                        circlize::chordDiagram(adjacencydf, annotationTrack = "grid", preAllocateTracks = 1, transparency = 0.5)
                        circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
                                circlize::circos.text(circlize::CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                                            facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
                        }, bg.border = NA) # Set bg.border to NA to remove default borders around sectors
                        return(list(table = kinase_substrate_table, plot = NULL)) # Chord diagram is drawn directly
                }

                if (plot == "sankey_diagram") {
                        p <- ggplot2::ggplot(df, ggplot2::aes(axis1 = Kinases, axis2 = Substrates)) +
                                ggalluvial::geom_alluvium(ggplot2::aes(fill = Kinases), width = 1/3) +
                                ggalluvial::geom_stratum(width = 1/3, fill = "grey", color = "black") +
                                ggplot2::geom_text(stat = "stratum", ggplot2::aes(label = ggplot2::after_stat(stratum)), size = 3) +
                                ggplot2::scale_x_discrete(limits = c("Kinase", "Phosphosite"), expand = c(.05, .05)) +
                                ggplot2::theme_minimal() +
                                ggplot2::ggtitle("Kinase \u2192 Phosphosite Sankey Diagram") +
                                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
                        return(list(table = kinase_substrate_table, plot = p))
                }
                # If plot type is not 'no' and not recognized for kinase_substrate table,
                # or if the necessary packages are not loaded, it might fall through.
                # Adding a specific return for unhandled cases if needed, but 'match.arg'
                # already limits the choices.
                stop("Invalid plot type for 'kinase_substrate' table type, or missing packages for selected plot.")
        }

        # Kinase Score Table and Plots
        if (table_type == "kinase_score") {
                kinase_score_table <- KSEAapp::KSEA.Scores(
                        KSData,
                        data,
                        NetworKIN = NetworKIN,
                        NetworKIN.cutoff = NetworKIN.cutoff
                )

                if (plot == "no") {
                        return(kinase_score_table)
                }
                if (plot == "bar_plot") {
                        df <- kinase_score_table %>%
                                dplyr::select(Kinase.Gene, m, z.score, p.value) %>%
                                dplyr::mutate(point_color = dplyr::case_when(
                                        p.value < 0.05 & z.score > 0 ~ "lightpink",
                                        p.value < 0.05 & z.score < 0 ~ "lightblue",
                                        TRUE ~ "black"
                                ))

                        p <- ggplot2::ggplot(df, ggplot2::aes(x = z.score, y = stats::reorder(Kinase.Gene, z.score), fill = point_color)) +
                                ggplot2::geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
                                ggplot2::scale_fill_identity() +
                                ggplot2::labs(
                                        title = "Kinase Z-scores", # Title can be generalized if 'across drug doses' is not always true
                                        subtitle = "Light pink: p < 0.05 & Z > 0; Light blue: p < 0.05 & Z < 0",
                                        x = "Z-score", y = "Kinases"
                                ) +
                                ggprism::theme_prism() +
                                ggplot2::theme(
                                        axis.text.x = ggplot2::element_text(size = 8),
                                        axis.text.y = ggplot2::element_text(size = 9),
                                        axis.title = ggplot2::element_text(size = 12, face = "bold"),
                                        plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
                                        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, face = "italic"),
                                        legend.position = "none",
                                        panel.grid.major.x = ggplot2::element_blank(),
                                        panel.grid.minor.x = ggplot2::element_blank(),
                                        panel.background = ggplot2::element_rect(fill = "white", color = NA),
                                        plot.background = ggplot2::element_rect(fill = "white", color = NA)
                                )
                        return(list(table = kinase_score_table, plot = p))
                }

                if (plot == "scatter_plot") {
                        p <- ggplot2::ggplot(kinase_score_table, ggplot2::aes(x = z.score, y = -log10(p.value))) +
                                ggplot2::geom_point(size = 3, color = "#2a9d8f") +
                                ggplot2::geom_text(ggplot2::aes(label = Kinase.Gene), color = "#2a9d8f",
                                                   nudge_x = 0.01, nudge_y = 0.04, check_overlap = TRUE) + # Changed check_overlap to TRUE for better label placement
                                ggplot2::geom_hline(yintercept = -log10(0.05), colour = "grey", linewidth = 0.5, linetype = "dashed") + # Added dashed line
                                ggplot2::labs(
                                        title = "Kinase Activity Scatter Plot",
                                        x = "Z-score",
                                        y = "-log10(p-value)"
                                ) +
                                ggprism::theme_prism() +
                                ggplot2::theme(
                                        plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
                                        axis.title = ggplot2::element_text(size = 12, face = "bold")
                                )
                        return(list(table = kinase_score_table, plot = p))
                }
                # If plot type is not 'no' and not recognized for kinase_score table,
                # or if the necessary packages are not loaded, it might fall through.
                stop("Invalid plot type for 'kinase_score' table type, or missing packages for selected plot.")
        }
}
