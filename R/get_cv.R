#' Calculate Coefficient of Variation (CV) and Optionally Draw Plot
#'
#' This function calculates the Coefficient of Variation (CV) for each row
#' (e.g., each protein or phosphosite) in a numerical data frame. It can
#' optionally generate a box plot, violin plot, histogram or density plot of the calculated CVs.
#' The CV is calculated as (standard deviation / mean) * 100.
#'
#' @param data A data frame containing numerical values for which CV should be
#'   calculated row-wise. Missing values (`NA`) are handled by removing them
#'   from the calculation of mean and standard deviation for each row.
#' @param name A character string specifying the name of the new column that
#'   will store the calculated Coefficient of Variation values. Default is "cv".
#' @param plot A character string indicating whether to generate a plot and
#'   which type.
#'   \itemize{
#'     \item "no": No plot is generated; only the data frame with the added CV column is returned.
#'     \item "box_plot": Generates a box plot of the calculated CV values.
#'     \item "violin_plot": Generates a violin plot of the calculated CV values.
#'     \item "histogram": Generates a histogram of the calculated CV values.
#'     \item "density_plot": Generates a density plot of the calculated CV values.
#'   }
#'
#' @returns #' If `plot = "no"`, returns the original `data` data frame with an additional
#' column (named as specified by `name`) containing the row-wise Coefficient
#' of Variation.
#' If `plot = "box_plot"` or `plot = "violin_plot"`or `plot = "histogram"` or
#' `plot = "violin_plot"` returns a `list` containing:
#' \itemize{
#'   \item `data`: The original `data` data frame with the added CV column.
#'   \item `plot`: A `ggplot` object of the generated box plot, violin plot,
#' histogram or density plot.
#' }
#' @details
#' The function calculates the CV for each row independently. It's particularly
#' useful for quality control in omics data, where CV can indicate technical
#' variability across replicates for a given feature.
#'
#' **Missing Values (`NA`)**: `NA` values in a row are ignored when calculating
#' the mean and standard deviation for that row's CV. If a row contains only
#' `NA` values, its CV will be `NaN` (Not a Number).
#'
#' **Plotting**:
#' For plotting, all calculated CV values are aggregated into a single distribution.
#' The plots do not categorize CV by groups unless further manipulation is done
#' outside this function. The x-axis for plots is set to "All Samples" to reflect
#' this global distribution.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' sample_data <- data.frame(
#'   Feature = paste0("F", 1:10),
#'   Replicate1 = rnorm(10, mean = 100, sd = 10),
#'   Replicate2 = rnorm(10, mean = 105, sd = 12),
#'   Replicate3 = rnorm(10, mean = 95, sd = 8),
#'   Replicate4 = c(rnorm(9, mean = 100, sd = 10), NA) # Add an NA for demonstration
#' )
#'
#' # Case 1: Calculate CV and return data only
#' cv_results_df <- get_cv(data = sample_data[, -1], # Exclude 'Feature' column for CV calculation
#'                         name = "CV_percent",
#'                         plot = "no")
#' print(head(cv_results_df))
#'
#' # Case 2: Calculate CV and generate a box plot
#' cv_results_boxplot <- get_cv(data = sample_data[, -1],
#'                              name = "CV_value",
#'                              plot = "box_plot")
#' print(head(cv_results_boxplot$data))
#' cv_results_boxplot$plot # Display the box plot
#'
#' # Case 3: Calculate CV and generate a violin plot
#' cv_results_violinplot <- get_cv(data = sample_data[, -1],
#'                                name = "CV_value",
#'                                plot = "violin_plot")
#' print(head(cv_results_violinplot$data))
#' cv_results_violinplot$plot # Display the violin plot
#'
#' # Case 4: Calculate CV and generate a histogram
#' cv_results_histogram <- get_cv(data = sample_data[, -1],
#'                                 name = "CV_value",
#'                                 plot = "histogram")
#' print(head(cv_results_histogram$data))
#' cv_results_histogram$plot # Display the histogram
#'
#' # Case 5: Calculate CV and generate a density_plot
#' cv_results_density_plot <- get_cv(data = sample_data[, -1],
#'                                name = "CV_value",
#'                                plot = "density_plot")
#' print(head(cv_results_density_plot$data))
#' cv_results_density_plot$plot # Display the density plot
#'
#' # Example with NA-only row
#' na_data <- data.frame(A = c(1, NA, 3), B = c(2, NA, 4))
#' na_cv <- get_cv(na_data, plot = "no")
#' print(na_cv)

#'
get_cv <- function(data,
                   name = "cv",
                   plot = c("no","box_plot", "violin_plot", "histogram", "density_plot")){

        # Validate 'plot' argument
        plot <- match.arg(plot)

        # Check for required packages
        if (!requireNamespace("dplyr", quietly = TRUE)) {
                stop("The 'dplyr' package is required but not installed. Please install it.")
        }

        if (plot %in% c("box_plot", "violin_plot")) {
                if (!requireNamespace("ggplot2", quietly = TRUE)) {
                        stop("The 'ggplot2' package is required for plotting but not installed. Please install it.")
                }
                if (!requireNamespace("ggprism", quietly = TRUE)) {
                        stop("The 'ggprism' package is required for plotting but not installed. Please install it.")
                }
        }

        # Calculate CV for each row
        # Using rowwise and summarise is more direct for row-wise operations
        output <- data %>%
                dplyr::rowwise() %>% # Operate row-by-row
                dplyr::mutate("{name}" := { # Calculate CV for the current row
                        current_row_values <- unlist(dplyr::cur_data()) # Get all values in the current row
                        s_d <- stats::sd(current_row_values, na.rm = TRUE)
                        m_ean <- mean(current_row_values, na.rm = TRUE)
                        if (is.na(m_ean) || m_ean == 0) { # Handle cases where mean is NA or 0 to avoid NaN/Inf CV
                                NA_real_
                        } else {
                                (s_d / m_ean) * 100
                        }
                }) %>%
                dplyr::ungroup() # Always ungroup after rowwise operations

        # Handle no plot return
        if (plot == "no") {
                return(output)
        }

        # Prepare data for plotting (only the CV column)
        plot_df <- output %>%
                dplyr::select(tidyselect::all_of(name)) %>% # Select only the CV column
                dplyr::rename(CV = tidyselect::all_of(name)) # Rename it to 'CV' for plotting convenience

        # Generate plot based on type
        p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = "All Samples", y = CV)) + # X-axis indicates overall distribution
                ggplot2::labs(x = "", y = "Coefficient of Variation (%)",
                              title = paste0("Distribution of ", name, " Values")) +
                ggprism::theme_prism() +
                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), # Center plot title
                               axis.text.x = ggplot2::element_blank(), # No specific groups on x-axis
                               axis.ticks.x = ggplot2::element_blank()) # No x-axis ticks


        if (plot == "box_plot") {
                p <- p + ggplot2::geom_boxplot(fill = "steelblue", color = "black", alpha = 0.7)
        } else if (plot == "violin_plot") {
                p <- p + ggplot2::geom_violin(fill = "lightgreen", color = "black", alpha = 0.7) +
                        ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) # Add narrow boxplot inside violin
        } else if (plot == "histogram") {
                # Histogram uses a single continuous variable on the x-axis
                p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = CV)) +
                        ggplot2::geom_histogram(binwidth = 5, fill = "#00AFBB", color = "black", alpha = 0.7) +
                        ggplot2::labs(x = "Coefficient of Variation (%)", y = "Frequency",
                                      title = paste0("Distribution of ", name, " Values")) +
                        ggprism::theme_prism() +
                        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        } else if (plot == "density_plot") {
                p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = CV)) +
                        ggplot2::geom_density(fill = "#00AFBB", color = "black", alpha = 0.7) +
                        ggplot2::labs(x = "Coefficient of Variation (%)", y = "Density",
                                      title = paste0("Distribution of ", name, " Values")) +
                        ggprism::theme_prism() +
                        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        }
        return(list(data = output, plot = p))
}
