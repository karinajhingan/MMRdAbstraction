#' @importFrom ggplot2 ggplot aes_string geom_bar labs theme_minimal
#' @importFrom biostatUtil wrap_ggkm
#' @importFrom dplyr filter select distinct mutate %>% any_of all_of
#' @importFrom tidyr drop_na
#' @importFrom purrr map
#' @importFrom gtsummary tbl_uvregression tbl_strata bold_p bold_labels style_pvalue
#' @importFrom rlang expr sym parse_expr :=
#' @importFrom utils globalVariables
utils::globalVariables(c("patient_id", "km_args_list"))

#' Create a Bar Plot for a Categorical Variable
#'
#' @param data A data frame containing the variable to plot.
#' @param x_var A string name of the categorical variable to plot on the x-axis.
#' @param plot_title A string for the plot title and x-axis label.
#'
#' @return A ggplot2 bar plot object.
#' @export
create_bar_plot <- function(data, x_var, plot_title) {
  ggplot(data, aes_string(x = x_var)) +
    geom_bar(stat = "bin", width = 0.7, fill = "steelblue") +
    labs(x = plot_title, y = "Count", title = plot_title) +
    theme_minimal()
}

#' Count Patients with Specific Marker Loss Patterns
#'
#' @param data A data frame containing marker columns and patient IDs.
#' @param loss_levels A vector of values considered as "loss" for markers.
#' @param must_have_loss A character vector of marker column names that should not show loss.
#' @param must_not_have_loss A character vector of marker column names that should not not show loss.
#'
#' @return An integer count of distinct patients matching the criteria.
#' @export
count_loss <- function(data, loss_levels, must_have_loss, must_not_have_loss) {
  data %>%
    filter(
      # Include markers that must show loss
      !!!lapply(must_have_loss, function(marker) {
        rlang::expr(!!rlang::sym(marker) %in% loss_levels)
      }),
      # Exclude markers that should not show loss
      !!!lapply(must_not_have_loss, function(marker) {
        rlang::expr(is.na(!!rlang::sym(marker)) | !(!!rlang::sym(marker) %in% loss_levels))
      })
    ) %>%
    select(patient_id) %>%
    distinct() %>%
    nrow()
}

#' Mark Patients with Multi-Marker Loss
#'
#' @param data A data frame containing marker columns.
#' @param markers A character vector of marker column names to check for loss.
#' @param loss_levels A vector of values considered as "loss".
#' @param tag_name A string name for the new logical column indicating multi-marker loss.
#'
#' @return A data frame with a new logical column indicating whether all markers show loss.
#' @export
tag_multi_loss <- function(data, markers, loss_levels, tag_name) {
  data %>%
    drop_na(all_of(markers)) %>%
    mutate(
      !!tag_name := ifelse(
        Reduce(`&`, lapply(markers, function(marker) {
          get(marker) %in% loss_levels
        })),
        TRUE,
        FALSE
      )
    )
}

#' Generate and Print Kaplan-Meier Curves
#'
#' @param args A named list of arguments to pass to `doKMPlots`.
#' @param caps A string title to use for the plot.
#'
#' @return A printed ggplot object with KM curves.
#' @export
plot_km_curves <- function(args, caps) {
  p <- km_args_list %>%
    map(~ exec("doKMPlots", !!!c(.x, args))) %>%
    wrap_ggkm(title = caps, tag_levels = "A")
  print(p)
}
