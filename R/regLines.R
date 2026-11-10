

#' Automatisierte Regressionsanalyse nach Gruppen
#'
#' Diese Funktion führt eine lineare Regression für Gruppen innerhalb eines Dataframes durch
#' und erstellt eine Visualisierung mitsamt einer Tabelle der wichtigsten Kennwerte der
#' linearen Regression
#'
#' @param df Ein Dataframe, der die zu analysierenden Daten enthält.
#' @param reg_group Die Spalte im Dataframe, die die Gruppen für die linearen Regressionen angibt.
#' @param point_group Die Spalte im Dataframe, die die Gruppen für die Punktfarben in der Visualisierung angibt.
#' @param y Die abhängige Variable für die Regression.
#' @param x Die unabhängige Variable für die Regression.
#' @param colors Ein benutzerdefiniertes Farbschema für die Gruppen in form eines Vektors mit c("Treatment" = "Farbe", usw)
#' @param variable_name Optionaler Name für die abhängige Variable in der Titelbeschriftung.
#' @param predictor_name Optionaler Name für die unabhängige Variable in der Titelbeschriftung.
#' @param y_label Optionales Label für die y-Achse.
#' @param x_label Optionales Label für die x-Achse.
#'
#' @return Eine Liste mit dem Plot und separat noch einmal der Modelloutput
#' \itemize{
#'   \item `regression_plot`: Ein ggplot-Objekt mit der Visualisierung der linearen Regressionen und einer Tabelle.
#'   \item `models`: Ein Dataframe mit den Modellen, deren Gleichungen, R²-Werten und p-Werten.
#' }
#'
#' @examples
#' # Beispiel-Datenframe erstellen
#' df <- data.frame(
#'   reg_group = rep(c("A", "B"), each = 10),
#'   point_group = rep(c("X", "Y"), times = 10),
#'   y = rnorm(20, mean = 50, sd = 10),
#'   x = rnorm(20, mean = 10, sd = 3)
#' )
#' colors <- c("A" = "blue", "B" = "red")
#' regLines(df, reg_group = "reg_group", point_group = "point_group", y = "y", x = "x", colors = colors)
#'
#' @export
regLines <- function(df, reg_group, point_group, y, x,
                     colors, variable_name = y, predictor_name = x, y_label = y, x_label = x){

  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(ggplotify)
  library(cowplot)

  significance_label <- function(p) {
    if (p < 0.001) return("***")
    else if (p < 0.01) return("**")
    else if (p < 0.05) return("*")
    else if (p < 0.1) return(".")
    else return("")
  }

  models <- df %>%
    group_by(!!sym(reg_group)) %>%
    do(model = lm(reformulate(x, y), data = .)) %>%
    mutate(
      eq = paste0("y = ", round(coef(model)[2], 2), "x + ", round(coef(model)[1], 2)),
      r2 = paste0("R² = ", round(summary(model)$r.squared, 2)),
      p_val_raw = summary(model)$coefficients[2, 4],
      p_val = paste0("p = ", formatC(p_val_raw, format = "e", digits = 2), " ", significance_label(p_val_raw))
    )

  plot_1 <- ggplot(data = df, aes(x = x, y = y)) +
    geom_point(aes(x = .data[[x]], y = .data[[y]], color = .data[[point_group]]))+
    geom_smooth(aes(x = .data[[x]], y = .data[[y]], color = .data[[reg_group]]), method = "lm", se = FALSE)+


    #geom_point(aes(col = point_group)) +
    #geom_smooth(aes(col = reg_group), method = "lm", se = F) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    labs(title = paste("Linear Regression", variable_name, "\nagainst", predictor_name),
         y = y_label,
         x = x_label) +
    theme(
      text = element_text(family = "serif", size = 25),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.border = element_rect(color = "darkgrey", fill = NA, size = 1),
      axis.line = element_line(color = "black"),
      plot.caption = element_text(hjust= 0.5, size = 15),
      plot.margin = margin(5, 5, 5, 5)
    )


  table_data <- models %>% select(reg_group, eq, r2, p_val)
  colnames(table_data) <- c("Mineral", "Equation", "Coefficient of determination", "p-value")


  table_grob <- gridExtra::tableGrob(table_data)


  regression_plot <- plot_grid(plot_1, table_grob, ncol = 1, rel_heights = c(5, 1))

  return(list(regression_plot = regression_plot, models = models))


}
