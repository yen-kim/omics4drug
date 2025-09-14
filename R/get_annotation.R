#' Map Gene Identifiers
#'
#' @param data A data frame containing the gene identifiers to be mapped.
#' @param input_id_column A character string specifying the name of the column
#'   in `data` that contains the gene identifiers to be mapped.
#' @param input_id_type A character string specifying the type of the input
#'   gene identifiers (e.g., "ENSEMBL", "SYMBOL", "REFSEQ"). This corresponds
#'   to the `keytype` argument in `AnnotationDbi::mapIds`.
#'   Common `keytype` values can be found by running `keytypes(org.Hs.eg.db)`
#'   for homo sapiens or `keytypes(org.Mm.eg.db)` for mus musculus.
#' @param output_id_type A character string specifying the type of the input
#'   gene identifiers (e.g., "ENSEMBL", "SYMBOL", "REFSEQ"). This corresponds
#'   to the `keytype` argument in `AnnotationDbi::mapIds`.
#'   Common `keytype` values can be found by running `keytypes(org.Hs.eg.db)`
#'   for homo sapiens or `keytypes(org.Mm.eg.db)` for mus musculus.
#' @param organism  character string, either "Homo sapiens" or "Mus musculus",
#'   specifying the organism for which the mapping should be performed.
#'   Defaults to "Homo sapiens".
#'
#' @returns A data frame identical to the input `data` frame, but with an
#'   additional column named `entrez_id` containing the mapped Entrez Gene IDs.
#'   If a mapping is not found for an ID, `NA` will be returned in the `entrez_id` column.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # --- Example for Homo sapiens ---
#' # Create a sample data frame
#' my_data_human <- data.frame(
#'   GeneSymbol = c("TP53", "BRCA1", "MYC", "GAPDH", "UNKNOWN_GENE"),
#'   Expression = c(10, 15, 20, 50, 5),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Map Gene Symbols to Entrez IDs for human
#' mapped_data_human <- get_annotation(
#'   data = my_data_human,
#'   input_id_column = "GeneSymbol",
#'   input_id_type = "SYMBOL",
#'   output_id_type = "ENTREZID",
#'   organism = "Homo sapiens"
#' )
#' print(mapped_data_human)
#'
#' # --- Example for Mus musculus ---
#' # Create another sample data frame
#' my_data_mouse <- data.frame(
#'   EnsemblID = c("ENSMUSG00000020717", "ENSMUSG00000026774", "ENSMUSG00000000001"),
#'   FoldChange = c(1.2, -0.8, 2.5),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Map Ensembl IDs to Entrez IDs for mouse
#' mapped_data_mouse <- get_annotation(
#'   data = my_data_mouse,
#'   input_id_column = "EnsemblID",
#'   input_id_type = "ENSEMBL",
#'   output_id_type = "ENTREZID",
#'   organism = "Mus musculus"
#' )
#' print(mapped_data_mouse)
#'
#' # Example with a non-existent input column
#' tryCatch({
#'   get_annotation(my_data_human, "NonExistentColumn", "SYMBOL", "ENTREZID")
#' }, error = function(e) {
#'   message("Caught expected error: ", e$message)
#' })
#' }
get_annotation <- function(data,
                           input_id_column,
                           input_id_type,
                           output_id_type = "ENTREZID", # Default to ENTREZID as per common use case
                           organism = c("Homo sapiens", "Mus musculus")) {

        # Validate organism input
        organism <- match.arg(organism)

        # Check if input_id_column exists in the data frame
        if (!input_id_column %in% names(data)) {
                stop(paste0("Error: The specified input_id_column '", input_id_column, "' does not exist in the data frame."))
        }

        # Determine which organism database to use
        if (organism == "Homo sapiens") {
               db_object <- org.Hs.eg.db::org.Hs.eg.db
        } else if (organism == "Mus musculus") {
                db_object <- org.Mm.eg.db::org.Mm.eg.db
        }

        # Perform the mapping
        entrez <- AnnotationDbi::mapIds(
                        x = db_object,
                        keys = data[[input_id_column]],
                        column = output_id_type,
                        keytype = input_id_type,
                        multiVals = 'first'
        )

        # Add the mapped IDs as a new column to the data frame
        data <- data %>%
                dplyr::mutate(entrez_id = entrez)

        return(data)
}
