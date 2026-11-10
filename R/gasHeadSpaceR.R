#' Convert GC raw data from headspace gas measurements to Concentrations in ppm
#'
#' Filters gas measurement data, fits linear calibration models for carbon
#' dioxide (CO2), nitrous oxide (N2O), and methane (CH4) using standard samples 
#' (P1, P2,and P3, and uses the fitted models to calculate concentrations in ppm.
#' Calibration plots are also generated.
#'
#' The standard concentration vectors must contain three values in the order
#' P1, P2, P3. Runs listed in \code{remove} are excluded before calibration
#' and concentration calculation. By default, no runs are excluded.
#'
#' @param data A data frame containing the raw data from our gas GC. It must include
#'   the columns \code{Sample}, \code{Run}, \code{Methane_Area},
#'   \code{Carbon_Dioxide_Height}, and \code{Nitrous_Oxide_Area}.
#' @param soll_ppm_CO2 Numeric vector of known CO2 concentrations in ppm for
#'   standards P1, P2, and P3. Defaults to \code{c(265, 1433, 2197)}.
#' @param soll_ppm_N2O Numeric vector of known N2O concentrations in ppm for
#'   standards P1, P2, and P3. Defaults to \code{c(576, 1099, 2332)}.
#' @param soll_ppm_methan Numeric vector of known methane concentrations in ppm
#'   for standards P1, P2, and P3. Defaults to \code{c(0.92, 2.45, 11.16)}.
#' @param remove Values in the \code{Run} column to exclude. Defaults to
#'   \code{c()}, which excludes no runs.
#'
#' @return A list containing:
#'   \describe{
#'     \item{data}{A data frame with the selected input columns and calculated
#'     concentrations \code{CO2_ppm}, \code{N2O_ppm}, and \code{methan_ppm}.}
#'     \item{plot}{Calibration plots for CO2, N2O, and methane, combined into
#'     one plot if the plots are assembled using the \pkg{patchwork} package.}
#'   }
#'
#' @export
#'
#' @examples
#' under construction


gasHeadSpaceR <- function(data, soll_ppm_CO2 = c(265, 1433, 2197), soll_ppm_N2O = c(576, 1099, 2332), soll_ppm_methan = c(0.92, 2.45, 11.16), remove = c() ){
  
  df_gas1 <- data %>%
    dplyr::select(
      Sample, Run, Methane_Area,
      Carbon_Dioxide_Height, Nitrous_Oxide_Area
    ) %>%
    dplyr::filter(!(Run %in% remove))
  
  
  
  # standard curve ----
  
  standards <- c("P1", "P2", "P3")
  
  df_prüf <- df_gas1 %>%
    dplyr::filter(Sample %in% standards) %>%
    dplyr::mutate(
      soll_ppm_CO2 = .env$soll_ppm_CO2[
        match(.data$Sample, .env$standards)
      ],
      soll_ppm_N2O = .env$soll_ppm_N2O[
        match(.data$Sample, .env$standards)
      ],
      soll_ppm_methan = .env$soll_ppm_methan[
        match(.data$Sample, .env$standards)
      ]
    )
  
  ## CO2 ----
  
  cal_CO2 <- lm( df_prüf$Carbon_Dioxide_Height ~ df_prüf$soll_ppm_CO2)
  summary(cal_CO2)
  rsqCO2 <- round (summary(cal_CO2)$r.squared, 4)
  intercept_CO2 <- coef(cal_CO2)[1]
  slope_CO2 <- coef(cal_CO2)[2]
  
  plot_std_CO2 <- ggplot(df_prüf, aes(x = soll_ppm_CO2, y =  Carbon_Dioxide_Height)) +
    geom_abline(intercept = coef(cal_CO2)[1], slope = coef(cal_CO2)[2], color = "darkblue", linewidth = 2) +
    geom_point(color = "darkgrey", shape = 16, size = 3) +
    geom_text(aes(label = Run), hjust = 2, size = 2)+
    geom_segment(aes(xend = soll_ppm_CO2, yend = fitted(cal_CO2)), color = "darkred", size = 1.5) +
    labs(
      title = paste0("Kalibriergerade CO2"),
      subtitle = bquote(R^2 == .(rsqCO2) ~ "," ~ "slope" == .(slope_CO2)~ "," ~ "intercept" == .(intercept_CO2)),
      x = "CO2-Konzentration [ppm]",
      y = "Peak Höhe CO2"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.border = element_rect(color = "darkgrey", fill = NA, size = 0.5),
      axis.line = element_line(color = "black")
    )
  
  ## N2O ----
  
  cal_N2O <- lm( df_prüf$Nitrous_Oxide_Area ~ df_prüf$soll_ppm_N2O)
  summary(cal_N2O)
  rsqN2O <- round (summary(cal_N2O)$r.squared, 4)
  intercept_N2O <- coef(cal_N2O)[1]
  slope_N2O <- coef(cal_N2O)[2]
  
  plot_std_N2O <- ggplot(df_prüf, aes(x = soll_ppm_N2O, y =  Nitrous_Oxide_Area)) +
    geom_abline(intercept = coef(cal_N2O)[1], slope = coef(cal_N2O)[2], color = "darkblue", linewidth = 2) +
    geom_point(color = "darkgrey", shape = 16, size = 3) +
    geom_text(aes(label = Run), hjust = 2, size = 2)+
    geom_segment(aes(xend = soll_ppm_N2O, yend = fitted(cal_N2O)), color = "darkred", size = 1.5) +
    labs(
      title = paste0("Kalibriergerade N2O"),
      subtitle = bquote(R^2 == .(rsqN2O) ~ "," ~ "slope" == .(slope_N2O)~ "," ~ "intercept" == .(intercept_N2O)),
      x = "N2O-Konzentration [ppm]",
      y = "Peak Höhe N2O"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.border = element_rect(color = "darkgrey", fill = NA, size = 0.5),
      axis.line = element_line(color = "black")
    )
  
  
  ## methan ----
  
  cal_methan <- lm( df_prüf$Methane_Area ~ df_prüf$soll_ppm_methan)
  summary(cal_methan)
  rsqmethan <- round (summary(cal_methan)$r.squared, 4)
  intercept_methan <- coef(cal_methan)[1]
  slope_methan <- coef(cal_methan)[2]
  
  plot_std_CH4 <- ggplot(df_prüf, aes(x = soll_ppm_methan, y =  Methane_Area)) +
    geom_abline(intercept = coef(cal_methan)[1], slope = coef(cal_methan)[2], color = "darkblue", linewidth = 2) +
    geom_point(color = "darkgrey", shape = 16, size = 3) +
    geom_text(aes(label = Run), hjust = 2, size = 2)+
    geom_segment(aes(xend = soll_ppm_methan, yend = fitted(cal_methan)), color = "darkred", size = 1.5) +
    labs(
      title = paste0("Kalibriergerade methan"),
      subtitle = bquote(R^2 == .(rsqmethan) ~ "," ~ "slope" == .(slope_methan)~ "," ~ "intercept" == .(intercept_methan)),
      x = "methan-Konzentration [ppm]",
      y = "Peak Höhe methan"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.border = element_rect(color = "darkgrey", fill = NA, size = 0.5),
      axis.line = element_line(color = "black")
    )
  
  
  # values to ppm ----
  
  df_gas2 <- df_gas1 %>%
    mutate(
      CO2_ppm = (Carbon_Dioxide_Height - intercept_CO2) / slope_CO2, 
      
      N2O_ppm = (Nitrous_Oxide_Area - intercept_N2O) / slope_N2O, 
      
      methan_ppm = (Methane_Area - intercept_methan) / slope_methan
    )
  
  plot_std_curves <- patchwork::wrap_plots(
    plot_std_CO2, plot_std_N2O, plot_std_CH4,
    ncol = 2
  )
  return(list(data = df_gas2, plot = plot_std_curves))
  
}
  