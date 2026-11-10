

#' Processsing Data from Chloroform fumigation extraction (CFE)
#'
#' The `cfeR` function merges data from the excel output format of our JenaTOC and
#' a table with sample infomration and processes them to calculate microbial
#' carbon (Cmic), extractable organic carbon (EOC), microbial nitrogen (Nmic),
#' extractable total nitrogen (ETN) and microbial biomass C:N ratio. Furthermore,
#' it provides a grafical visualization of the five variables including three
#' possible laysers of sample information and a table with summary statistics
#' (mean, sd). Additionally, it gives the opportunity to remove
#' identified outlaysers selectively.
#'
#'
#' @param data_jeantoc imported raw data directely from JenaTOC (stored as data frame)
#'
#' @param sample_master A data frame containing sample metadata required for
#'   calculations. Must include the following columns:
#'   \itemize{
#'     \item \code{id_jenatoc}: unique sample identifier matching Jeantoc data
#'     \item \code{TS}: dry matter content (%) of the soil sample
#'     \item \code{ew_nfum}: extraction weight for non-fumigated samples
#'     \item \code{ew_fum}: extraction weight for fumigated samples
#'   }
#' @param extr_vol_addition Numeric. Additional extraction volume in milliliters
#'   (default = 40).
#' @param dil_factor Numeric. Dilution factor applied to concentrations
#'   (default = 1, applies only if no dilution is necessary,
#'   which is the case, when 0.025 M K2SO4 was used for extraction for 13-C or 15-N).
#' @param fumcode Character. Identifier suffix for fumigated samples
#'   (default = "_fum").
#' @param nfumcode Character. Identifier suffix for non-fumigated samples
#'   (default = "_nfum").
#' @param colour Character. Column name used for coloring points in plots. This is
#' the only required sample description layer. If just one sample property should be
#' visualized please interpret colour eqal to "treatment"
#' @param shape Optional character or NA. Column name used for point shapes in plots.
#'    (default = NA).
#' @param facet Optional character or NA. Column name used for faceting plots
#'   (default = NA).
#' @param remove_outliers Vector or NA. Optional vector of sample IDs to exclude
#'   from calculations and plots (default = NA).
#'
#' @details
#' The function:
#' \itemize{
#'   \item Separates fumigated and non-fumigated samples.
#'   \item Merges with sample metadata and calculates extraction volumes.
#'   \item Computes microbial biomass indices (Cmic, Nmic), extractable pools
#'         (EOC, ETN) in µg per g TS, and their ratio.
#'
#'         here the formulas used are:
#'
#' \itemize{
#'   \item \strong{EOC (Extractable Organic Carbon)}:
#'     \deqn{EOC = \frac{c_{TC,nfum} \cdot dil\_factor \cdot extr\_vol\_nfum}{ew\_nfum \cdot TS} \cdot 100}
#'
#'   \item \strong{Cmic (Microbial Biomass Carbon)}:
#'     \deqn{Cmic = \left(\frac{c_{TC,fum} \cdot dil\_factor \cdot extr\_vol\_fum}{ew\_fum \cdot TS} - \frac{c_{TC,nfum} \cdot dil\_factor \cdot extr\_vol\_nfum}{ew\_nfum \cdot TS}\right) \cdot \frac{100}{0.45}}
#'
#'   \item \strong{ETN (Extractable Total Nitrogen)}:
#'     \deqn{ETN = \frac{c_{TN,nfum} \cdot dil\_factor \cdot extr\_vol\_nfum}{ew\_nfum \cdot TS} \cdot 100}
#'
#'   \item \strong{Nmic (Microbial Biomass Nitrogen)}:
#'     \deqn{Nmic = \left(\frac{c_{TN,fum} \cdot dil\_factor \cdot extr\_vol\_fum}{ew\_fum \cdot TS} - \frac{c_{TN,nfum} \cdot dil\_factor \cdot extr\_vol\_nfum}{ew\_nfum \cdot TS}\right) \cdot \frac{100}{0.54}}
#'
#'   \item \strong{Cmic:Nmic ratio}:
#'     \deqn{CNmic = \frac{Cmic}{Nmic}}
#' }
#'
#' where:
#' \itemize{
#'   \item `c_{TC,nfum}`: measured total C concentration in non-fumigated extracts (µg/ml)
#'   \item `c_{TC,fum}`: measured total C concentration in fumigated extracts (µg/ml)
#'   \item `c_{TN,nfum}`: measured total N concentration in non-fumigated extracts (µg/ml)
#'   \item `c_{TN,fum}`: measured total N concentration in fumigated extracts (µg/ml)
#'   \item `dil_factor`: dilution factor of the extract before measurement
#'   \item `extr_vol_nfum`: extract volume of non-fumigated sample (ml)
#'   \item `extr_vol_fum`: extract volume of fumigated sample (ml)
#'   \item `ew_nfum`: wet sample weight for non-fumigated extract (g)
#'   \item `ew_fum`: wet sample weight for fumigated extract (g)
#'   \item `TS`: dry matter content of the sample (percent)
#' }

#'
#'   \item Optionally removes specified outlier samples.
#'   \item Aggregates results by grouping variables (`colour`, `shape`, `facet`)
#'         to compute means and standard deviations.
#'   \item Produces ggplot2 visualizations with error bars for each index.
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{df_results_cfe}{Data frame with results (and some intermediate steps) added to the sample master}
#'   \item{plot_cmic}{ggplot object showing microbial biomass carbon (Cmic).}
#'   \item{plot_nmic}{ggplot object showing microbial biomass nitrogen (Nmic).}
#'   \item{plot_eoc}{ggplot object showing extractable organic carbon (EOC).}
#'   \item{plot_etn}{ggplot object showing extractable total nitrogen (ETN).}
#'   \item{plot_cnmic}{ggplot object showing the Cmic:Nmic ratio.}
#'   \item{tbl_cfe}{Data frame with mean and standard deviation of all measures grouped by colour, shape and facet
#'   (non-specified grouping variables are ignored automatically)}
#' }
#'
#' @examples
#' under construction
#'
#' @export


cfeR <- function(data_jeantoc, sample_master, extr_vol_addition = 40, dil_factor = 1, fumcode = "_fum", nfumcode = "_nfum", colour, shape = NA, facet = NA, remove_outliers = NA ) {

  library(tidyverse)

  ending <- c(fumcode, nfumcode)

  df_1 <- data_jeantoc%>%
    select("Sample_ID", "c_TC", "c_TN")%>%
    mutate(
      proc = case_when(
        grepl(paste0(fumcode), Sample_ID) ~ "fum",
        grepl(paste0(nfumcode), Sample_ID) ~ "nfum"),
      id_jenatoc = str_remove_all(Sample_ID, paste0(ending, collapse = "|")))%>%
    filter(id_jenatoc != "H2O")


  df_nfum <- df_1 %>%
    filter(proc == "nfum")%>%
    rename("c_TC_nfum" = "c_TC", "c_TN_nfum" = "c_TN")

  df_fum <- df_1 %>%
    filter(proc == "fum")%>%
    rename("c_TC_fum" = "c_TC", "c_TN_fum" = "c_TN")

  df_2 <- merge(df_nfum, df_fum, by = "id_jenatoc") %>%
    select("id_jenatoc", "c_TC_nfum", "c_TN_nfum", "c_TC_fum", "c_TN_fum") %>%
    mutate(
      c_TC_fum = as.numeric(gsub(",", ".", as.character(c_TC_fum))),
      c_TC_nfum = as.numeric(gsub(",", ".", as.character(c_TC_nfum))),
      c_TN_fum = as.numeric(gsub(",", ".", as.character(c_TN_fum))),
      c_TN_nfum = as.numeric(gsub(",", ".", as.character(c_TN_nfum)))
    )


  df_3 <- merge(sample_master, df_2, by = "id_jenatoc")

  df_calc <- df_3 %>%
    mutate(
      # Extraction volume
      WG = 1 - TS/100,
      soil_water_nfum = WG * ew_nfum,
      soil_water_fum = WG * ew_fum,
      extr_vol_nfum = extr_vol_addition + soil_water_nfum, # in ml
      extr_vol_fum = extr_vol_addition + soil_water_fum,
      dil_factor = dil_factor,
      # EOC
      eoc = (((c_TC_nfum * dil_factor * extr_vol_nfum)/(ew_nfum * TS))) * (100),
      # Cmic
      cmic = (((c_TC_fum * dil_factor* extr_vol_fum)/(ew_fum * TS)) - ((c_TC_nfum * dil_factor * extr_vol_nfum)/(ew_nfum * TS))) * (100/0.45),
      # ETN
      etn = (((c_TN_nfum * dil_factor * extr_vol_nfum)/(ew_nfum * TS))) * (100),
      # Cmic
      nmic = (((c_TN_fum * dil_factor* extr_vol_fum)/(ew_fum * TS)) - ((c_TN_nfum * dil_factor * extr_vol_nfum)/(ew_nfum * TS))) * (100/0.54),
      # cnmic
      cnmic = cmic/nmic
    )

  # outlier
  if (!all(is.na(remove_outliers))) {
    df_calc <- df_calc[!(df_calc$id %in% remove_outliers), ]
    warning(
      paste(
        "Data from the samples with the following IDs were manually removed:",
        paste(remove_outliers, collapse = ", "),
        "\nThis constitutes a data manipulation. Please ensure that this removal is properly documented in your R-script or lab notes."
      )
    )
  }


  # visualization

  if (is.na(shape)) {
    df_calc$shape_ <- "NA"

    shape = "shape_"
  }

  vars <- c(colour, shape, facet)
  vars <- vars[!is.na(vars)]   # entfernt NAs

  tbl_cfe <- df_calc %>%
    group_by(across(all_of(vars))) %>%
    summarise(
      mean_cmic = mean(cmic),
      sd_cmic   = sd(cmic),
      mean_nmic = mean(nmic),
      sd_nmic   = sd(nmic),
      mean_eoc  = mean(eoc),
      sd_eoc    = sd(eoc),
      mean_etn  = mean(etn),
      sd_etn    = sd(etn),
      mean_cnmic = mean(cnmic),
      sd_cnmic = sd(cnmic),
      .groups = "drop"
    )


  plot_cmic <- ggplot(data = df_calc,
                      aes_string(x = colour, y = "cmic", shape = shape)) +
    geom_point(data = tbl_cfe, aes_string(x = colour, y = "mean_cmic", colour = colour), size = 4)+
    geom_errorbar(data = tbl_cfe, aes_string( x = colour, ymin = "mean_cmic - sd_cmic", ymax = "mean_cmic + sd_cmic", width = 0.5, colour = colour), inherit.aes = FALSE,  linetype = "solid")+
    geom_point(size = 2, col = "black") +
    geom_text(aes(label = id), vjust = -0.5, hjust = -0.5, size = 3)+
    facet_wrap(as.formula(paste("~", facet)),
               labeller = labeller(.cols = function(x) paste(facet, "=", x))) +
    labs(title = bquote(C[mic]))+
    theme_minimal()+
    theme(
      text = element_text(family = "serif", size = 17),
      plot.title = element_text(hjust = 0, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5),
      panel.border = element_rect(color = "darkgrey", fill = NA, size = 0.5),
      axis.line = element_line(color = "black")
    )

  plot_nmic <- ggplot(data = df_calc,
                      aes_string(x = colour, y = "nmic", shape = shape)) +
    geom_point(data = tbl_cfe, aes_string(x = colour, y = "mean_nmic", colour = colour), size = 4)+
    geom_errorbar(data = tbl_cfe, aes_string( x = colour, ymin = "mean_nmic - sd_nmic", ymax = "mean_nmic + sd_nmic", width = 0.5, colour = colour), inherit.aes = FALSE,  linetype = "solid")+
    geom_point(size = 2, col = "black") +
    geom_text(aes(label = id), vjust = -0.5, hjust = -0.5, size = 3)+
    facet_wrap(as.formula(paste("~", facet)),
               labeller = labeller(.cols = function(x) paste(facet, "=", x))) +
    labs(title = bquote(N[mic]))+
    theme_minimal()+
    theme(
      text = element_text(family = "serif", size = 17),
      plot.title = element_text(hjust = 0, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5),
      panel.border = element_rect(color = "darkgrey", fill = NA, size = 0.5),
      axis.line = element_line(color = "black")
    )

  plot_eoc <- ggplot(data = df_calc,
                     aes_string(x = colour, y = "eoc", shape = shape)) +
    geom_point(data = tbl_cfe, aes_string(x = colour, y = "mean_eoc", colour = colour), size = 4)+
    geom_errorbar(data = tbl_cfe, aes_string( x = colour, ymin = "mean_eoc - sd_eoc", ymax = "mean_eoc + sd_eoc", width = 0.5, colour = colour), inherit.aes = FALSE,  linetype = "solid")+
    geom_point(size = 2, col = "black") +
    geom_text(aes(label = id), vjust = -0.5, hjust = -0.5, size = 3)+
    facet_wrap(as.formula(paste("~", facet)),
               labeller = labeller(.cols = function(x) paste(facet, "=", x))) +
    labs(title = bquote(EOC))+
    theme_minimal()+
    theme(
      text = element_text(family = "serif", size = 17),
      plot.title = element_text(hjust = 0, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5),
      panel.border = element_rect(color = "darkgrey", fill = NA, size = 0.5),
      axis.line = element_line(color = "black")
    )

  plot_etn <- ggplot(data = df_calc,
                     aes_string(x = colour, y = "etn", shape = shape)) +
    geom_point(data = tbl_cfe, aes_string(x = colour, y = "mean_etn", colour = colour), size = 4)+
    geom_errorbar(data = tbl_cfe, aes_string( x = colour, ymin = "mean_etn - sd_etn", ymax = "mean_etn + sd_etn", width = 0.5, colour = colour), inherit.aes = FALSE,  linetype = "solid")+
    geom_point(size = 2, col = "black") +
    geom_text(aes(label = id), vjust = -0.5, hjust = -0.5, size = 3)+
    facet_wrap(as.formula(paste("~", facet)),
               labeller = labeller(.cols = function(x) paste(facet, "=", x))) +
    labs(title = bquote(ETN))+
    theme_minimal()+
    theme(
      text = element_text(family = "serif", size = 17),
      plot.title = element_text(hjust = 0, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5),
      panel.border = element_rect(color = "darkgrey", fill = NA, size = 0.5),
      axis.line = element_line(color = "black")
    )

  plot_cnmic <- ggplot(data = df_calc,
                       aes_string(x = colour, y = "cnmic", shape = shape)) +
    geom_point(data = tbl_cfe, aes_string(x = colour, y = "mean_cnmic", colour = colour), size = 4)+
    geom_errorbar(data = tbl_cfe, aes_string( x = colour, ymin = "mean_cnmic - sd_cnmic", ymax = "mean_cnmic + sd_cnmic", width = 0.5, colour = colour), inherit.aes = FALSE,  linetype = "solid")+
    geom_point(size = 2, col = "black") +
    geom_text(aes(label = id), vjust = -0.5, hjust = -0.5, size = 3)+
    facet_wrap(as.formula(paste("~", facet)),
               labeller = labeller(.cols = function(x) paste(facet, "=", x))) +
    labs(title = bquote(C[mic]~ ":" ~ N[mic]))+
    theme_minimal()+
    theme(
      text = element_text(family = "serif", size = 17),
      plot.title = element_text(hjust = 0, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5),
      panel.border = element_rect(color = "darkgrey", fill = NA, size = 0.5),
      axis.line = element_line(color = "black")
    )



  return(list(df_results_cfe = df_calc, tbl_cfe = tbl_cfe, plot_cmic = plot_cmic, plot_nmic = plot_nmic, plot_eoc = plot_eoc, plot_etn = plot_etn, plot_cnmic = plot_cnmic))

}
