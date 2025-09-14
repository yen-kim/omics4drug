#' Perform KEGG Pathway Enrichment Analysis and Visualization
#'
#' This function performs Kyoto Encyclopedia of Genes and Genomes (KEGG) pathway
#' enrichment analysis on a provided list of gene identifiers using the
#' `clusterProfiler` package. It supports human and mouse organisms. Optionally,
#' it can generate bar plots or bubble plots to visualize the top enriched KEGG pathways.
#'
#' @param input_list A character vector of gene identifiers, ideally ENTREZ IDs,
#'   for which to perform KEGG enrichment analysis. If other key types are used,
#'   ensure they can be mapped to ENTREZ IDs by `clusterProfiler`.
#' @param organism A character string specifying the organism. Must be one of
#'   "Homo sapiens" (for human, uses KEGG code 'hsa') or "Mus musculus" (for mouse,
#'   uses KEGG code 'mmu').
#' @param plot A character string specifying the type of plot to generate.
#'   \itemize{
#'     \item "no": No plot is generated; only the KEGG enrichment table is returned.
#'     \item "bar_plot": Generates a bar plot of the top 20 enriched KEGG pathways.
#'       Bars are colored by `-log10(pvalue)`.
#'     \item "bubble_plot": Generates a bubble plot of the top 20 enriched KEGG pathways.
#'       Bubble size represents `Count` (number of genes in the pathway), and
#'       color represents `-log10(pvalue)`.
#'   }
#'
#' @returns
#' If `plot = "no"`, returns a `data.frame` of KEGG enrichment results.
#' If `plot = "bar_plot"` or `plot = "bubble_plot"`, returns a `list` containing:
#' \itemize{
#'   \item `KEGG_data`: A `data.frame` of KEGG enrichment results.
#'   \item `plot`: A `ggplot` object representing the combined plot of top KEGG pathways.
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Example for Homo sapiens using ENTREZIDs
#' # A gene list (replace with your actual ENTREZ IDs)
#' human_entrez_genes <- c("597", "1000", "1001", "10010", "10011", "10012",
#'                         "10013", "10014", "10015", "10016", "10017", "10018")
#'
#' # Get KEGG data only
#' kegg_results_human_table <- get_KEGG(input_list = human_entrez_genes,
#'                                           organism = "Homo sapiens",
#'                                           plot = "no")
#' print(head(kegg_results_human_table))
#'
#' # Get KEGG data and bar plot
#' kegg_results_human_bar_plot <- get_KEGG(input_list = human_entrez_genes,
#'                                              organism = "Homo sapiens",
#'                                              plot = "bar_plot")
#' print(head(kegg_results_human_bar_plot$KEGG_data))
#' kegg_results_human_bar_plot$plot # This will display the plot
#'
#' # Example for Mus musculus using ENTREZIDs
#' # A gene list (replace with your actual ENTREZ IDs)
#' mouse_entrez_genes <- c("12536", "12543", "12551", "12574", "12610", "12613",
#'                         "12614", "12615", "12618", "12620", "12621", "12623")
#'
#' # Get KEGG data and bubble plot for mouse
#' kegg_results_mouse_bubble_plot <- get_KEGG(input_list = mouse_entrez_genes,
#'                                                  organism = "Mus musculus",
#'                                                  plot = "bubble_plot")
#' print(head(kegg_results_mouse_bubble_plot$KEGG_data))
#' kegg_results_mouse_bubble_plot$plot # This will display the plot
#' }
get_KEGG <- function(input_list,
                          organism = c("Homo sapiens", "Mus musculus"),
                          plot = c("no","bar_plot", "bubble_plot")) {
        organism <- match.arg(organism)
        plot <- match.arg(plot)

        # Determine KEGG organism code
        kegg_organism_code <- NULL
        if (organism == "Homo sapiens") {
                kegg_organism_code <- "hsa"
        } else if (organism == "Mus musculus") {
                kegg_organism_code <- "mmu"
        }

        # Perform KEGG enrichment analysis
        # Note: enrichKEGG typically expects ENTREZID for 'gene' and uses KEGG organism code for 'organism'
        KEGG_data_results <- clusterProfiler::enrichKEGG(gene = input_list,
                                                         organism = kegg_organism_code,
                                                         pvalueCutoff = 0.05)

        # Convert results to data frame
        KEGG_data <- as.data.frame(KEGG_data_results)

        # Handle no plot return
        if (plot == "no") {
                return(KEGG_data)
        }

        # Prepare data for plotting
        # Filter for top 20 terms, sorted by pvalue
        plot_data_KEGG <- KEGG_data %>%
                dplyr::arrange(pvalue) %>% # Arrange by pvalue for 'top' terms
                utils::head(20)

        # Define the plotting function
        # Corrected 'category' to 'Description' for KEGG pathways
        p <- ggplot2::ggplot(plot_data_KEGG, ggplot2::aes(x = FoldEnrichment, y = stats::reorder(Description, -pvalue))) +
                ggplot2::labs(x = "FoldEnrichment", y = "KEGG Pathway",
                              title = "KEGG Pathway Enrichment Analysis") +
                viridis::scale_color_viridis(option = "G", name = "-log10(pvalue)") +
                ggprism::theme_prism() +
                ggplot2::theme(
                        axis.text.x = ggplot2::element_text(size = 8),
                        axis.text.y = ggplot2::element_text(size = 9),
                        axis.title = ggplot2::element_text(size = 12, face = "bold"),
                        plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold"),
                        legend.position = "right" # Show legend for color
                )

        if (plot == "bar_plot") {
                p <- p + ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = -log10(pvalue)), color = "black", linewidth = 0.2) +
                        viridis::scale_fill_viridis(option = "G", name = "-log10(pvalue)") # Use fill for bar plots
        } else if (plot == "bubble_plot") {
                p <- p + ggplot2::geom_point(ggplot2::aes(color = -log10(pvalue)), alpha = 0.8)
        }

        return(list(KEGG_data = KEGG_data, plot = p))
}






