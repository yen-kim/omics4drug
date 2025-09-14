#' Get Kinase Inhibitors from a Database
#'
#' This function filters a provided drug-kinase interaction database to find
#' drugs that inhibit a list of specified kinases in a given organism.
#' It can return either the filtered data frame or an interactive `DT` data table.
#'
#' @param kinase_list A character vector of kinase UniProt accessions or names.
#'   These are the kinases you want to find inhibitors for.
#' @param drug_database A data frame of known drug-kinase interactions.
#'   This table must contain columns for `uniprot_accessions` or
#'   a similar kinase identifier, and `Drug`.
#' @param return_datatable A logical value. If `TRUE`, the function returns a
#'   list containing the raw data frame and an interactive `DT::datatable`.
#'   If `FALSE`, it returns only the data frame.
#'
#' @returns A data frame or a list. If `return_datatable` is `FALSE`, returns
#'   a data frame with the filtered inhibitors. If `TRUE`, returns a list
#'   containing both the data frame and the interactive data table.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a database for demonstration
#' database <- data.frame(
#'   Drug = c("Imatinib", "Erlotinib", "Gefitinib", "Dasatinib", "Rapamycin"),
#'   Target = c("ABL1", "EGFR", "EGFR", "SRC", "MTOR"),
#'   uniprot_accessions = c("P00519", "P00533", "P00533", "P12931", "P42345"),
#'   Organism = c("Homo sapiens", "Homo sapiens", "Homo sapiens",
#'                "Homo sapiens", "Mus musculus")
#' )
#'
#' # Example 1: Get the data frame of inhibitors for Homo sapiens kinases
#' my_kinases <- c("P00519", "P00533", "P12931")
#' inhibitors_df <- get_inhibitor (kinase_list = my_kinases,
#'                                 drug_database = database,
#'                                 return_datatable = FALSE
#'                                 )
#' print(inhibitors_df)
#'
#' # Example 2: Get both the data frame and the interactive data table
#' my_kinases_mouse <- c("P42345")
#' inhibitors_dt <- get_inhibitor(kinase_list = my_kinases_mouse,
#'                                drug_database = database,
#'                                return_datatable = TRUE
#'                               )
#'
#' # Access the list elements
#' print(inhibitors_dt$raw_data)
#' inhibitors_dt$datatable # This will display the interactive table
#' }
get_inhibitor <- function(kinase_list,
                                  drug_database,
                                  return_datatable = FALSE) {

        # Validate inputs

        if (!is.character(kinase_list) || length(kinase_list) == 0) {
                stop("`kinase_list` must be a non-empty character vector.")
        }
        if (!is.data.frame(drug_database)) {
                stop("`drug_database` must be a data frame.")
        }
        if (!all(c("uniprot_accessions") %in% colnames(drug_database))) {
                stop("`drug_database` must contain 'uniprot_accessions' columns.")
        }

        # Filter the database
        inhibitor_list <- drug_database %>%
                dplyr::filter(uniprot_accessions %in% kinase_list)

        if (return_datatable) {
                drug_table <- DT::datatable(inhibitor_list,
                                            options = list(scrollX = TRUE, scrollY = "600px", autoWidth = TRUE))
                return(list(raw_data = inhibitor_list, datatable = drug_table))
        } else {
                return(inhibitor_list)
        }
}
