#' Normalize Phosphoproteomics Data and Perform Differential Analysis with `limma`
#'
#' This function processes raw phosphoproteomics intensity data, performs
#' filtering, imputation, and median scaling using `PhosR` functions, and then
#' conducts differential phosphorylation analysis using the `limma` package.
#' Optionally, it can generate a volcano plot to visualize the results.
#'
#' @param data A data frame containing phosphoproteomics quantification data.
#'   It is expected to have columns:
#'   \itemize{
#'     \item `protein_group`: Character, representing protein identifiers.
#'     \item `amino_acid`: Character, representing the phosphorylated amino acid (e.g., "S", "T", "Y").
#'     \item `site`: Numeric, representing the position of the phosphorylation site.
#'     \item `modification_sites`: Character, a unique identifier for each phosphosite
#'       (e.g., "protein_S123"). This will be used as the row names for the output table.
#'     \item Columns suffixed with "_control" and "_treat" followed by a sample
#'       identifier (e.g., "p1_control", "p7_control", "p5_treat", "p6_treat") containing quantitative values.
#'   }
#' @param alpha Numeric, a parameter for `PhosR::selectGrps` function, controlling
#'   the proportion of missing values allowed. Default is `0.5`.
#' @param beta Numeric, a parameter for `PhosR::scImpute` function, controlling
#'   the maximum number of neighbors for imputation. Default is `0.7`.
#' @param plot  character string specifying whether to generate a plot.
#'   \itemize{
#'     \item "no": No plot is generated; only the differential phosphorylation table is returned.
#'     \item "volcano_plot": Generates a volcano plot to visualize differential phosphorylation.
#'   }
#' @param value A character string indicating which p-value to use for the y-axis
#'   of the volcano plot when `plot = "volcano_plot"`.
#'   \itemize{
#'     \item "p_value": Uses nominal p-values (`P.Value`) for the y-axis.
#'     \item "adj_p_value": Uses adjusted p-values (`adj.P.Val`) for the y-axis.
#'   }
#'   This parameter is ignored if `plot = "no"`.
#'
#' @returns
#' If `plot = "no"`, returns a `data.frame` (tibble) containing the results of
#' the `limma` differential analysis (e.g., `logFC`, `P.Value`, `adj.P.Val`).
#' If `plot = "volcano_plot"`, returns a `list` containing:
#' \itemize{
#'   \item `table`: The `data.frame` of differential phosphorylation results.
#'   \item `plot`: A `ggplot` object of the generated volcano plot.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a data frame that mimics omics data
#' # The column names must contain "control" and "treat"
#' data <- tibble(
#'         modification_sites = paste0("site_", 1:10),
#'         protein_group = rep(LETTERS[1:5], 2),
#'         amino_acid = c(rep("S", 5), rep("T", 5)),
#'         site = 1:10,
#'         PTM.Group = paste0("group_", 1:10),
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
#' results_table <- get_norm_phos(
#'         data = data,
#'         plot = "no"
#' )
#'# Print the head of the results table
#' print(head(results_table))
#'# Run the function to get the results and the volcano plot
#'# We set plot = "volcano_plot" and value = "adj_p_value"
#' results_with_plot <- get_norm_phos(
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
get_norm_phos <- function(data,
                          alpha = 0.5,
                          beta = 0.7,
                          plot = c("no","volcano_plot"),
                          value = c("p_value", "adj_p_value")) {
        plot <- match.arg(plot)
        value <- match.arg(value)

        # Add a PTM.Group column if it doesn't exist, as PhosR::PhosphoExperiment expects it
        if (!"PTM.Group" %in% colnames(data)) {
                data <- data %>% dplyr::mutate(PTM.Group = "empty")
        }

        # Extract sample intensity columns
        # Assumes column names contain "control" or "treat"
        samples_id <- c("control", "treat") # Expected sample group identifiers
        ptm_intensity <- data %>% dplyr::select(dplyr::contains(samples_id))

        # Create PhosR PhosphoExperiment object
        ppe <- PhosR::PhosphoExperiment(assays = list(Quantification = as.matrix(ptm_intensity)))

        # Assign metadata to ppe object
        # Ensure these columns exist in the input 'data' dataframe
        if ("protein_group" %in% colnames(data)) {
                ppe@GeneSymbol <- as.character(data[["protein_group"]])
        } else {
                warning("Column 'protein_group' not found in input data. GeneSymbol will be left empty in PhosphoExperiment object.")
        }
        if ("amino_acid" %in% colnames(data)) {
                ppe@Residue <- as.character(data[["amino_acid"]])
        } else {
                warning("Column 'amino_acid' not found in input data. Residue will be left empty in PhosphoExperiment object.")
        }
        if ("site" %in% colnames(data)) {
                ppe@Site <- as.numeric(data[["site"]])
        } else {
                warning("Column 'site' not found in input data. Site will be left empty in PhosphoExperiment object.")
        }
        if ("PTM.Group" %in% colnames(data)) {
                ppe@Sequence <- as.character(data[["PTM.Group"]])
        } else {
                warning("Column 'PTM.Group' not found in input data. Sequence will be left empty in PhosphoExperiment object.")
        }
        if ("modification_sites" %in% colnames(data)) {
                ppe@NAMES <- as.character(data[["modification_sites"]])
        } else {
                warning("Column 'modification_sites' not found in input data. NAMES will be left empty in PhosphoExperiment object.")
                # If modification_sites is missing, limma rownames will be auto-generated by PhosR,
                # which might be less descriptive. Consider stopping or generating it.
        }

        # Determine sample groups for design matrix based on column names
        # This assumes column names are e.g., "control_1", "treat_A"
        grps <- sub("^[^_]+_", "", colnames(ppe)) # Extracts "1", "A" etc.

        # Log2 transformation and handling infinite values
        logmat <- log2(SummarizedExperiment::assay(ppe, "Quantification"))
        logmat[is.infinite(logmat)] <- NA
        SummarizedExperiment::assay(ppe, "Quantification") <- logmat

        # Filter out phosphosites with too many missing values
        ppe_filtered <- PhosR::selectGrps(ppe, grps, alpha, n = 1)

        # Imputation and Scaling
        set.seed(123) # For reproducibility of imputation
        ppe_imputed_scI <- PhosR::scImpute(ppe_filtered, beta, grps)
        ppe_imputed <- PhosR::tImpute(ppe_imputed_scI, assay = "imputed")
        ppe_imputed_scaled <- PhosR::medianScaling(ppe_imputed, scale = FALSE, assay = "imputed")

        # Differential expression analysis using limma
        design <- stats::model.matrix(~ grps - 1)

        fit <- limma::lmFit(ppe_imputed_scaled@assays@data$scaled, design)
        contrast.matrix <- limma::makeContrasts(grpstreat-grpscontrol,
                                                levels=design)
        fit2 <- limma::contrasts.fit(fit, contrast.matrix)
        fit2 <- limma::eBayes(fit2)

        norm_data <- limma::topTable(fit2, coef = "grpstreat - grpscontrol", sort.by = "none", number = Inf) %>%
                tibble::rownames_to_column(var = "phosphosites")

        # Handle no plot return
        if (plot == "no") {
                return(norm_data)
        }

        # Generate volcano plot
        if (plot == "volcano_plot") {
                y_var <- if (value == "p_value") { "P.Value" } else { "adj.P.Val" }
                plot_df <- norm_data %>%
                        dplyr::mutate(point_color = dplyr::case_when(
                                .data[[y_var]] < 0.05 & logFC < -0.2 ~ "down", # significantly down
                                .data[[y_var]] < 0.05 & logFC > 0.2 ~ "up", # significantly up
                                TRUE ~ "NS") # not significant
                        )

                # Determine y-axis variable based on 'value' parameter
                y_lab <- if (value == "p_value") { expression("-log"[10]*"P.Value") } else { expression("-log"[10]*"adj.P.Val") }

                p <- ggplot2::ggplot(data = plot_df,
                                     ggplot2::aes_string(x = "logFC", y = paste0("-log10(", y_var, ")"), col = "point_color")) +
                        ggplot2::geom_vline(xintercept = c(-0.2, 0.2),
                                            col = 'grey',
                                            linetype = 'dashed') +
                        ggplot2::geom_hline(yintercept = -log10(0.05),
                                            col = 'grey',
                                            linetype = 'dashed') +
                        ggplot2::geom_point(size = 1) +
                        ggplot2::scale_color_manual(values = c("down" = "blue", "NS" = "grey", "up" = "red"),
                                                    labels = c("Downregulated", "Not significant", "Upregulated"),
                                                    name = 'Phosphorylation State') +
                        ggplot2::labs(x = 'log2 fold change', y = y_lab) +
                        ggplot2::ggtitle('Volcano plot of Differential Phosphorylation') +
                        ggprism::theme_prism() +
                        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) # Center title

                return(list(table = norm_data, plot = p))
        }
}
