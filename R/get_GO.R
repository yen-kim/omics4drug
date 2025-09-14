#' Perform Gene Ontology (GO) Enrichment Analysis and Visualization
#'
#' This function performs Gene Ontology (GO) enrichment analysis on a provided
#' list of gene identifiers using the `clusterProfiler` package. It supports
#' human and mouse organisms and various key types for input genes. Optionally,
#' it can generate bar plots or bubble plots to visualize the top enriched GO
#' terms across Biological Process (BP), Cellular Component (CC), and Molecular
#' Function (MF) ontologies.
#'
#' @param input_list A character vector of gene identifiers (e.g., ENTREZIDs, Symbols)
#'   for which to perform GO enrichment analysis.
#' @param organism A character string specifying the organism. Must be one of
#'   "Homo sapiens" (for human) or "Mus musculus" (for mouse).
#' @param input_keytype A character string specifying the type of gene identifiers
#'   provided in `input_list`. This must be a valid `keyType` supported by
#'   `clusterProfiler::enrichGO` and the corresponding `OrgDb` package.
#'   Common options include "ENTREZID", "SYMBOL", "ENSEMBL", "UNIPROT", etc.
#'   Refer to `keytypes(OrgDb_package_name)` for a full list (e.g., `keytypes(org.Hs.eg.db)`
#'   for human or `keytypes(org.Mm.eg.db)` for mouse).
#' @param plot A character string specifying the type of plot to generate.
#'   \itemize{
#'     \item "no": No plot is generated; only the GO enrichment table is returned.
#'     \item "bar_plot": Generates separate bar plots for the top 20 enriched
#'       terms in each GO ontology (BP, CC, MF). The bars are colored by `-log10(pvalue)`.
#'     \item "bubble_plot": Generates separate bubble plots for the top 20 enriched
#'       terms in each GO ontology (BP, CC, MF). Bubble size represents `Count`
#'       (number of genes in the term), and color represents `-log10(pvalue)`.
#'   }
#'
#' @returns If `plot = "no"`, returns a `data.frame` of GO enrichment results.
#' If `plot = "bar_plot"` or `plot = "bubble_plot"`, returns a `list` containing:
#' \itemize{
#'   \item `GO_data`: A `data.frame` of GO enrichment results.
#'   \item `plot`: A `patchwork` object combining the plots for BP, CC, and MF ontologies.
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Example for Homo sapiens using ENTREZID
#' # A gene list (replace with your actual data)
#' human_genes <- c("100", "1000", "10000", "10001", "10002", "10003",
#'                  "10004", "10005", "10006", "10007", "10008", "10009",
#'                  "1001", "10010", "10011", "10012")
#'
#' # Get GO data only
#' go_results_human_table <- get_GO(input_list = human_genes,
#'                                       organism = "Homo sapiens",
#'                                       input_keytype = "ENTREZID",
#'                                       plot = "no")
#' print(head(go_results_human_table))
#'
#' # Get GO data and bar plots
#' go_results_human_bar_plot <- get_GO(input_list = human_genes,
#'                                          organism = "Homo sapiens",
#'                                          input_keytype = "ENTREZID",
#'                                          plot = "bar_plot")
#' print(head(go_results_human_bar_plot$GO_data))
#' go_results_human_bar_plot$plot # This will display the combined plot
#'
#' # Get GO data and bubble plots
#' go_results_human_bubble_plot <- get_GO(input_list = human_genes,
#'                                             organism = "Homo sapiens",
#'                                             input_keytype = "ENTREZID",
#'                                             plot = "bubble_plot")
#' print(head(go_results_human_bubble_plot$GO_data))
#' go_results_human_bubble_plot$plot # This will display the combined plot
#'
#' # Example for Mus musculus using SYMBOL
#' # A gene list (replace with your actual data)
#' mouse_genes <- c("Trp53", "Cdkn1a", "Mdm2", "Myc", "Fos", "Jun", "Akt1", "Pik3ca")
#'
#' # Get GO data and bubble plots for mouse
#' go_results_mouse_bubble_plot <- get_GO(input_list = mouse_genes,
#'                                              organism = "Mus musculus",
#'                                              input_keytype = "SYMBOL",
#'                                              plot = "bubble_plot")
#' print(head(go_results_mouse_bubble_plot$GO_data))
#' go_results_mouse_bubble_plot$plot
#' }
#'
get_GO <- function(input_list,
                        organism = c("Homo sapiens", "Mus musculus"),
                        input_keytype = c("ACCNUM", "ALIAS", "ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS",
                                          "ENTREZID", "ENZYME", "EVIDENCE", "EVIDENCEALL", "GENENAME",
                                          "GENETYPE", "GO", "GOALL", "IPI", "MGI",
                                          "ONTOLOGY", "ONTOLOGYALL", "PATH", "PFAM", "PMID",
                                          "PROSITE", "REFSEQ", "SYMBOL", "UNIPROT"),
                        plot = c("no","bar_plot", "bubble_plot")) {
        organism <- match.arg(organism)
        plot <- match.arg(plot)
        input_keytype <- match.arg(input_keytype) # Ensure keytype is validated
        OrgDb_package <- NULL
        if (organism == "Homo sapiens") {
                OrgDb_package <- "org.Hs.eg.db"
        } else if (organism == "Mus musculus") {
                OrgDb_package <- "org.Mm.eg.db"
        }

        # Perform GO enrichment analysis, passing the package name as a string
        GO_data_results <- clusterProfiler::enrichGO(gene = input_list,
                                                     OrgDb = OrgDb_package,
                                                     keyType = input_keytype,
                                                     ont = "ALL",
                                                     pAdjustMethod = "BH",
                                                     qvalueCutoff = 0.05,
                                                     readable = TRUE)
        # Convert results to data frame
        GO_data <- as.data.frame(GO_data_results)

        # Handle no plot return
        if (plot == "no") {
                return(GO_data)
        }

        # Prepare data for plotting
        # Filter for top 20 terms per ontology, sorted by pvalue
        plot_data_BP <- GO_data %>%
                dplyr::filter(ONTOLOGY == "BP") %>%
                dplyr::mutate(GO_pathway = paste0("BP: ", Description)) %>%
                dplyr::arrange(pvalue) %>% # Arrange by pvalue for 'top' terms
                utils::head(20)

        plot_data_CC <- GO_data %>%
                dplyr::filter(ONTOLOGY == "CC") %>%
                dplyr::mutate(GO_pathway = paste0("CC: ", Description)) %>%
                dplyr::arrange(pvalue) %>%
                utils::head(20)

        plot_data_MF <- GO_data %>%
                dplyr::filter(ONTOLOGY == "MF") %>%
                dplyr::mutate(GO_pathway = paste0("MF: ", Description)) %>%
                dplyr::arrange(pvalue) %>%
                utils::head(20)

        # Base plot elements for shared aesthetics
        base_plot <- function(df, plot_type) {
                p <- ggplot2::ggplot(df, ggplot2::aes(x = -log10(pvalue), y = stats::reorder(GO_pathway, -pvalue))) + # Reorder by -pvalue for consistency
                        ggplot2::labs(x = "-log10(p-value)", y = "", title = "") +
                        viridis::scale_color_viridis(option = "G", name = "-log10(pvalue)") + # Add name for color legend
                        ggprism::theme_prism() +
                        ggplot2::theme(
                                axis.text.x = ggplot2::element_text(size = 8),
                                axis.text.y = ggplot2::element_text(size = 9),
                                axis.title = ggplot2::element_text(size = 12, face = "bold"),
                                plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
                                legend.position = "right" # Show legend for color
                        )

                if (plot_type == "bar_plot") {
                        p <- p + ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = -log10(pvalue)), color = "black", linewidth = 0.2) +
                                viridis::scale_fill_viridis(option = "G", name = "-log10(pvalue)") # Use fill for bar plots
                } else if (plot_type == "bubble_plot") {
                        p <- p + ggplot2::geom_point(ggplot2::aes(color = -log10(pvalue)), alpha = 0.8) +
                                ggplot2::scale_size_continuous() # Adjust bubble size range
                }
                return(p)
        }

        # Generate plots for each ontology
        plot_BP <- base_plot(plot_data_BP, plot) + ggplot2::ggtitle("Biological Process (BP)")
        plot_CC <- base_plot(plot_data_CC, plot) + ggplot2::ggtitle("Cellular Component (CC)")
        plot_MF <- base_plot(plot_data_MF, plot) + ggplot2::ggtitle("Molecular Function (MF)")

        # Combine plots using patchwork
        combined_plot <- patchwork::wrap_plots(plot_BP, plot_CC, plot_MF, ncol = 1) +
                patchwork::plot_annotation(title = "Gene Ontology Enrichment Analysis",
                                           theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 18, face = "bold")))

        return(list(GO_data = GO_data, plot = combined_plot))
}
