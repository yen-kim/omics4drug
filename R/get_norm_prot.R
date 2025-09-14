#' Normalize Protein Data and Perform Differential Expression Analysis
#'
#' This function processes protein data by performing normalization, imputation,
#' and differential expression analysis using methods from the PhosR and limma packages.
#' It can also generate a volcano plot to visualize the results.
#'
#' @param data A data frame containing protein data. The function assumes that
#'   columns representing sample intensities have names containing suffixed with
#'   "_control" and "_treat" followed by a sample identifier (e.g., "p1_control",
#'   "p7_control", "p5_treat", "p6_treat").
#' @param alpha A numeric value between 0 and 1. This parameter is used by
#'   `PhosR::selectGrps` to filter out proteins with too many missing
#'   values. A protein is kept if it has at least `n` valid values in at least
#'   `alpha * number_of_groups` groups. Default is `0.5`.
#' @param beta A numeric value between 0 and 1. This parameter is used by
#'   `PhosR::scImpute`, which controls the number of nearest neighbors for imputation.
#'   Default is `0.7`.
#' @param plot A character string. Specifies whether to generate a volcano plot.
#'   Must be one of "no" (default) or "volcano_plot".
#' @param value A character string. Specifies which p-value to use for the
#'   volcano plot's y-axis. Must be one of "p_value" (default) or "adj_p_value".
#' @returns If `plot` is "no", the function returns a tibble with the results
#'   of the differential expression analysis. This table includes columns for
#'   log fold change (`logFC`), average expression, t-statistic, p-value, and
#'   adjusted p-value.
#'   If `plot` is "volcano_plot", the function returns a list containing
#'   two elements: `table` (the results tibble) and `plot` (a `ggplot` object
#'   of the volcano plot).
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a data frame that mimics omics data
#' # The column names must contain suffixed "_control" and "_treat"
#' data <- tibble(
#'         protein_group = LETTERS[1:10],
#'         p1_control = rnorm(10, mean = 100, sd = 10),
#'         p2_control = rnorm(10, mean = 95, sd = 15),
#'         p3_control = rnorm(10, mean = 105, sd = 12),
#'         p4_treat = c(rnorm(5, mean = 150, sd = 20), rnorm(5, mean = 50, sd = 10)),
#'         p5_treat = c(rnorm(5, mean = 145, sd = 18), rnorm(5, mean = 55, sd = 11)),
#'         p6_treat = c(rnorm(5, mean = 160, sd = 22), rnorm(5, mean = 60, sd = 13))
#' )
#'
#' # Run the function to get the normalized data table
#' # We set plot = "no" to suppress the volcano plot
#' results_table <- get_norm_prot(
#'         data = data,
#'         plot = "no"
#' )
#'# Print the head of the results table
#' print(head(results_table))
#'# Run the function to get the results and the volcano plot
#'# We set plot = "volcano_plot" and value = "adj_p_value"
#' results_with_plot <- get_norm_prot(
#'         data = data,
#'         plot = "volcano_plot",
#'         value = "adj_p_value"
#' )
#'
#' # The function returns a list, so we can access the plot and table separately
#' results_table <- results_with_plot$table
#' volcano_plot <- results_with_plot$plot
#'
#' # Print the plot
#' print(volcano_plot)
#'
#' # To save the plot to a file
#' # ggplot2::ggsave("volcano_plot.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)
#' }
get_norm_prot <- function(
                data,
                alpha = 0.5,
                beta = 0.7,
                plot = c("no", "volcano_plot"),
                value = c("p_value", "adj_p_value")) {
        plot <- match.arg(plot)
        value <- match.arg(value)
        data <- data %>% dplyr::mutate(amino_acid = NULL,
                                site = NULL,
                                modification_sites = NULL)
        # Extract sample intensity columns
        samples_id <- c("control", "treat")
        protein_intensity <- data %>% dplyr::select(dplyr::contains(samples_id))

        # Create PhosR ProteinExperiment object
        # PhosR::ProteinExperiment is used instead of PhosR::PhosphoExperiment
        ppe <- PhosR::PhosphoExperiment(assays = list(Quantification = as.matrix(protein_intensity)))

        ppe@GeneSymbol <- as.character(data[["protein_group"]])
        ppe@Residue <- as.character(data[["amino_acid"]])
        ppe@Site <- as.numeric(data[["site"]])

        # Assign metadata to the object. We use protein_group for NAMES.
        if ("protein_group" %in% colnames(data)) {
                ppe@NAMES <- as.character(data[["protein_group"]])
        } else {
                stop("Column 'protein_group' not found in input data. It is required for protein identifiers.")
        }

        # Determine sample groups for design matrix based on column names
        grps <- sub("^[^_]+_", "", colnames(ppe))

        # Normalization and Imputation ---
        # Log2 transformation and handling infinite values
        logmat <- log2(SummarizedExperiment::assay(ppe, "Quantification"))
        logmat[is.infinite(logmat)] <- NA
        SummarizedExperiment::assay(ppe, "Quantification") <- logmat

        # Filter out proteins with too many missing values
        ppe_filtered <- PhosR::selectGrps(ppe, grps, alpha, n = 1)

        # Imputation and Scaling
        set.seed(123)
        ppe_imputed_scI <- PhosR::scImpute(ppe_filtered, beta, grps)
        ppe_imputed <- PhosR::tImpute(ppe_imputed_scI, assay = "imputed")
        ppe_imputed_scaled <- PhosR::medianScaling(ppe_imputed, scale = FALSE, assay = "imputed")

        # Differential Expression Analysis with limma
        design <- stats::model.matrix(~ grps - 1)
        fit <- limma::lmFit(ppe_imputed_scaled@assays@data$scaled, design)
        contrast.matrix <- limma::makeContrasts(grpstreat-grpscontrol, levels=design)
        fit2 <- limma::contrasts.fit(fit, contrast.matrix)
        fit2 <- limma::eBayes(fit2)

        # Extract results and convert to a data frame
        norm_data <- limma::topTable(fit2, coef = "grpstreat - grpscontrol", sort.by = "none", number = Inf) %>%
                dplyr::rename(Protein = ID)

        # Generate Volcano Plot (if requested)
        if (plot == "no") {
                return(norm_data)
        }

        if (plot == "volcano_plot") {
                y_var <- if (value == "p_value") { "P.Value" } else { "adj.P.Val" }

                plot_df <- norm_data %>%
                        dplyr::mutate(point_color = dplyr::case_when(
                                .data[[y_var]] < 0.05 & logFC < -0.2 ~ "down",
                                .data[[y_var]] < 0.05 & logFC > 0.2 ~ "up",
                                TRUE ~ "NS")
                        )

                y_lab <- if (value == "p_value") { expression("-log"[10]*"P.Value") } else { expression("-log"[10]*"adj.P.Val") }

                p <- ggplot2::ggplot(data = plot_df, ggplot2::aes_string(x = "logFC", y = paste0("-log10(", y_var, ")"), col = "point_color")) +
                        ggplot2::geom_vline(xintercept = c(-0.2, 0.2), col = 'grey', linetype = 'dashed') +
                        ggplot2::geom_hline(yintercept = -log10(0.05), col = 'grey', linetype = 'dashed') +
                        ggplot2::geom_point(size = 1) +
                        ggplot2::scale_color_manual(values = c("down" = "blue", "NS" = "grey", "up" = "red"),
                                                    labels = c("Downregulated", "Not significant", "Upregulated"),
                                                    name = 'Protein Expression State') +
                        ggplot2::labs(x = 'log2 fold change', y = y_lab) +
                        ggplot2::ggtitle('Volcano plot of Differential Protein Expression') +
                        ggprism::theme_prism() +
                        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

                return(list(table = norm_data, plot = p))
        }
}
