#' Processing qPCR Data
#'
#' This function processes qPCR and NanoDrop data, merges it with sample information,
#' calculates the number of copies per g dry soil or any converted unit, and optionally exports
#' the processed data as an Excel file. Moreover, it optionally visualizes your results.
#'
#'
#' Calculation of DNA Copies per g dry soil
#'
#' The number of copies per g dry soil (`copies_per_dw`) is calculated as follows:
#'
#' \deqn{ \text{copies\_per\_dw} = \frac{\left(\frac{\text{quantity}}{\text{dil\_factor}} \times \text{vol\_eb} \right)}{\left( \text{ew} \times \frac{\text{TS}}{100} \right)} }
#'
#' where:
#' \itemize{
#'   \item `quantity`: Measured DNA concentration from qPCR
#'   \item `dil_factor`: Dilution factor of the extract before measurement
#'   \item `vol_eb`: Volume of elution buffer in µl
#'   \item `ew`: Sample weight (wet)
#'   \item `TS`: Dry matter content of the sample in percent
#' }
#'
#'Calculation of DNA concentration per g dry soil
#'
#' \deqn{ \text{conc\_per\_dw} = \frac{\text{mean\_conc\_nanodrop}}{\text{dw}} }
#'
#' where:
#' \itemize{
#'   \item `mean_conc_nanodrop`: Average DNA concentration measured with Nanodrop (ng/µl)
#'   \item `dw`: Dry weight of the sample in g
#' }
#'
#'
#' @param data_path_raw_data File path to the raw qPCR data.
#' @param df_dna_master Data frame containing DNA sample information, including weigh-in, sample ID,
#'        sample name in qPCR, dry matter content in Prozent , etc.
#' @param df_nanodrop Data frame containing NanoDrop measurement results.
#' @param sample_name_qPCR Column name for qPCR sample identification.
#' @param id_nanodrop Column name for NanoDrop sample identification.
#' @param dil_factor Column name for dilution factor (eg. if you have diluted 10 µl of your extract with 90 µl, it must be 0.1).
#' @param ew Column name for sample weight (wet).
#' @param TS Column name for dry matter content (Prozent).
#' @param conv_factor Column name for conversion factor from dry weight soil to any unit.
#' @param group Column name for grouping variable to group measurements according to e.g. treatments when producing the plot
#' @param extraction_block Column name for block of extraction (extraction run) sample identification.
#' @param threshold_nanodrop Threshold for acceptable deviation from mean concentration.
#' @param vol_eb Volume of elution buffer in µl.
#' @param apply_conversion_factor Logical (TRUE/FALSE) indicating whether the conversion factor should be applied.
#' @param auto_export Logical (TRUE/FALSE) indicating whether the processed data should be automatically exported.
#' @param export_file Name (incl. directory) of the export file.
#' @param plot Logical (TRUE/FALSE) indicating whether you wish to get a plot of your results (this will require (`group`))
#' @param extraction_block_points Logical (TRUE/FALSE) indicating whether you wish to include the results of every measurement as individual points in your plot (this will require (`extraction_block`) )
#'
#' @return A processed data frame containing calculated DNA copy numbers per dry soil (or converted unit).
#' @export


qPCinR <- function(data_path_raw_data, df_dna_master, df_nanodrop,  sample_name_qPCR = "sample_name_qPCR",
                   id_nanodrop = "id_nanodrop" , dil_factor = "dil_factor" ,
                   ew = "ew", TS = "TS", conv_factor = "conv_factor", group = "gruppe", extraction_block = "dna_block", threshold_nanodrop,
                   vol_eb, apply_conversion_factor = FALSE, auto_export = FALSE, export_file = "qPCinR_export_file",
                   plot = FALSE, extraction_block_points = FALSE)
  {


  library(tidyverse)
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
  library(patchwork)



  ## renaming columns

  df_dna_master <- df_dna_master %>%
    rename(!!"sample_name_qPCR" := !!sym(sample_name_qPCR)) %>%
    rename(!!"ew" := !!sym(ew)) %>%
    rename(!! "TS" := !! sym(TS)) %>%
    rename(!! "conv_factor" := !! sym(conv_factor))%>%
    rename(!! "group" := !! sym(group))



  df_nanodrop <- df_nanodrop %>%
    rename(!!"id_nanodrop" := !!sym(id_nanodrop))%>%
    rename(!!"dil_factor" := !!sym(dil_factor)) %>%
    na.omit()

  # importing raw data

  df_dna_roh <- read_excel(data_path_raw_data, skip = 6) %>%
    setNames(c("well", "sample_name_qPCR", "gen", "task", "reporter",
               "quencher", "ct", "mean_ct", "sd_ct", "quantity", "mean_quantity",
               "sd_quantity", "automatic_ct_treshhold", "ct_treshold", "automatic_baseline",
               "baseline_start", "baseline_end", "tm1", "tm2", "tm3", "comments", "")) %>%
    select(sample_name_qPCR, gen, ct, quantity, automatic_ct_treshhold, ct_treshold)

  ## process NanoDrop results (and give a warning message if values are suspicious)


  tbl_nanodrop <- df_nanodrop %>%
    group_by(id_nanodrop) %>%
    summarise(
      mean_conc_nanodrop = mean(conc),
      sd_conc_nanodrop = sd(conc),
      dil_factor = first(dil_factor)
    ) %>%
    left_join(df_nanodrop, by = "id_nanodrop") %>%
    mutate(
      diff_to_mean = abs(conc - mean_conc_nanodrop),
      check_threshold = if_else(diff_to_mean > threshold_nanodrop, "TRESHOLD ÜBERSCHRITTEN!", "OK")
    ) %>%
    rename(!! "dil_factor" := !! "dil_factor.x")

  exceeded_values <- tbl_nanodrop %>%
    filter(diff_to_mean > threshold_nanodrop) %>%
    select(id_nanodrop, diff_to_mean)

  if (nrow(exceeded_values) > 0) {
    warning(paste("The following NanoDrop measurements have exceeded  the given treshhold (threshold_nanodrop) for differendes to group mean concentration:",
                  paste("id_nanodrop =", exceeded_values$id_nanodrop, "( diff =", exceeded_values$diff_to_mean, ")", collapse = ", ")))
  }


  ## merge information from different data sources

  filtered_tbl_nanodrop <- tbl_nanodrop %>%
    select("id_nanodrop", "mean_conc_nanodrop", "sd_conc_nanodrop", "dil_factor") %>%
    distinct()

  df_dna <- merge(df_dna_master, df_dna_roh, by = "sample_name_qPCR") %>%
    merge(filtered_tbl_nanodrop, by = "id_nanodrop")


  ## Check for correct dilution
  dilution_groups <- filtered_tbl_nanodrop %>%
    count(dil_factor) %>%
    rename(dilution_factor = dil_factor, number_of_samples = n)

  num_unique_dil_factor <- nrow(dilution_groups)
  num_unique_sample_names <- df_dna_master %>% pull(sample_name_qPCR) %>% unique() %>% length()

  if (num_unique_dil_factor == num_unique_sample_names) {
    message("You have diluted your samples individually")
  } else {
    message(paste("You have diluted your extracts in", num_unique_dil_factor,
                  "groups"))

    message("The following tibble reports how many samples you have diluted with the respect dilution factor:")
    print(dilution_groups)
  }

  ## Check for correct dilution confirmation
  message("Please confirm that this information is correct:")
  message("1 = How did you know? Have you observed me in the lab?!")
  message("0 = Hm... I will check my lab book and get back to you later.")

  confirm <- readline(prompt = "Enter 1 or 0: ")

  if (confirm == "1") {
    message("Great. Here we go!")
  } else if (confirm == "0") {
    stop("Gehe zurück auf Los. Ziehe nicht 200 ein. Begib dich sofort dort hin")
  } else {
    stop("Invalid input. Please enter 1 or 0.")
  }


  ## calculate number of copies per dw or any converted unit

  df_dna$copies_per_elution_vol <- ((df_dna$quantity / df_dna$dil_factor)) * vol_eb
  df_dna$dw <- (df_dna$ew * df_dna$TS / 100)
  df_dna$copies_per_dw <- df_dna$copies_per_elution_vol / df_dna$dw

  df_dna$conc_per_dw <- df_dna$mean_conc_nanodrop / df_dna$dw



  if(apply_conversion_factor){
    df_dna <- df_dna %>%
      mutate(copies_per_converted_sample =  df_dna$copies_per_dw * df_dna$conv_factor)
    df_dna <- df_dna %>%
      mutate(conc_per_converted_sample =  df_dna$conc_per_dw * df_dna$conv_factor)
  }else {df_dna <- df_dna %>%
    mutate(copies_per_converted_sample =  df_dna$copies_per_dw * 1)
  df_dna <- df_dna %>%
    mutate(conc_per_converted_sample =  df_dna$conc_per_dw * 1)}


  if (auto_export) {
    write.xlsx(df_dna, file = paste0(export_file, ".xlsx"),
               overwrite = TRUE, asTable = FALSE)
    message(paste("Your processed data were successfully exported to",
                  paste0(export_file, ".xlsx")))
  } else {
    message("Your processed data were not exported, as you have set auto_export = FALSE")
  }



  ## Plot the data if plot = TRUE
  if (plot) {

    tbl_plot <- df_dna%>%
      group_by(group) %>%
      summarize(
        mean_copies = mean(copies_per_converted_sample),
        sd_copies = sd(copies_per_converted_sample),
        mean_conc = mean(conc_per_converted_sample),
        sd_conc = sd(conc_per_converted_sample),
        dna_block = dna_block,
        copies_per_converted_sample = copies_per_converted_sample,
        conc_per_converted_sample = conc_per_converted_sample
      )
    random_colors <- sample(scales::hue_pal()(length(unique(tbl_plot$group))))

    plot_copies <- ggplot(tbl_plot, aes(x = group, y = mean_copies, fill = group)) +
      geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
      geom_errorbar(aes(ymin = mean_copies - sd_copies, ymax = mean_copies + sd_copies),
                    width = 0.2, position = position_dodge(width = 0.9)) +
      theme_minimal() +
      labs(caption = "Error bars represent standard deviation from group mean \n copies per dry soil or converted unit (if TRUE)") +
      theme(
        text = element_text(family = "serif", size = 17),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "darkgrey", fill = NA, size = 1),
        axis.line = element_line(color = "black"),
        plot.caption = element_text(hjust= 0.5, size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 12, color = random_colors)
      ) +
      scale_fill_manual(values = random_colors, guide = "none")

    if (extraction_block_points) {
      plot_copies <- plot_copies +
        geom_point(data = df_dna,
                   aes(x = group, y = copies_per_converted_sample, color = dna_block),
                   size = 3)



    }


    plot_nanodrop <- ggplot(tbl_plot, aes(x = group, y = mean_conc, fill = group)) +
      geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
      geom_errorbar(aes(ymin = mean_conc - sd_conc, ymax = mean_conc + sd_conc),
                    width = 0.2, position = position_dodge(width = 0.9)) +
      theme_minimal() +
      labs(caption = "Error bars represent standard deviation from group mean \n concentration per dry soil or converted unit (if TRUE)") +
      theme(
        text = element_text(family = "serif", size = 17),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.border = element_rect(color = "darkgrey", fill = NA, size = 1),
        axis.line = element_line(color = "black"),
        plot.caption = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 12, color = random_colors)
      ) +
      scale_fill_manual(values = random_colors, guide = "none")


    if (extraction_block_points) {
      plot_nanodrop <- plot_nanodrop +
        geom_point(data = df_dna,
                   aes(x = group, y = conc_per_converted_sample, color = dna_block),
                   size = 3)

    }


    # Extrahiere den ersten Wert aus der Spalte "Target"
    title_text <- df_dna$gen[1]

    # Setze den Titel für den kombinierten Plot
    combined_plot <- (plot_copies + plot_nanodrop) +
      plot_layout(guides = "collect") +
      plot_annotation(title = title_text)


    print(combined_plot)

    return(list(df_dna = df_dna, combined_plot = combined_plot))

  }
}
