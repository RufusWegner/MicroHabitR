#' Processing of Heat Flow Raw Data
#'
#' This function filters heat flow raw data based on a selected time interval,
#' formats the data into an R-compatible structure, and visualizes the raw data
#' to provide an initial impression. As reference point, the "ampoule removed"
#' mark in the device software in hours after any arbitrary t0 is required
#' (i.e., insertion of the ampoules).
#'
#' @param data Raw data directly exported from the TAM Air software (Excel format)
#' @param amprem Time of "Ampoule removed" in hours after the start of the experiment. The exact time
#' must be documented at which the corresponding command (remove ampoules -> next) was made in the software.
#' @param start Start of the time interval to be extracted, in hours after the beginning of the experiment (t0)
#' @param end End of the time interval to be extracted, in hours after the beginning of the experiment (t0)
#' @param trt_id Vector containing the IDs of the treatments to be considered,
#' in the order of calorimeter slots 1 to 8. A default naming convention using "t" + slot number is applied.
#'
#' @return Filtered dataframe with the columns "time" in hours after experiment start, "ambient temperature" in °C,
#' "inner temperature" in °C, "heat flow" in µW, "heat" in J (from start = 0 J) and "treatment". Additionally, a ggplot is automatically generated
#' showing the treatments with color-coded distinction.
#' @export
#'
#' @examples
#' under construction


heatSelect <- function(data, amprem, start, end, trt_id = c("ch1", "ch2", "ch3", "ch4", "ch5", "ch6", "ch7", "ch8")){


  # Parameters: data = Excel file from the calorimeter as a data frame
  # amprem = Ampoule removal time in hours after the start of the experiment (t0; i.e., insertion of the ampoules)
  ## start = Start of the selected sequence in hours after t0
  ## end = End of the selected sequence in hours after t0


  remove_mark <- subset(data, data$`Time.markers` == "Ampoule removed")

  remove.zw <- amprem * 60*60


  kalorimeter_t0 <- (remove.zw - as.numeric(remove_mark[1,1])) / (60 * 60)



  ## select and rename columns
  data <- data[, c(1,2, 3, 4, 8, 12, 16,20,24,28,32, 5, 9, 13, 17, 21, 25, 29, 33)]
  names(data) <- c("time", "ambient_temp", "inner_temp",  "heat_flow_1","heat_flow_2","heat_flow_3","heat_flow_4","heat_flow_5",
                   "heat_flow_6","heat_flow_7","heat_flow_8", "heat_1", "heat_2", "heat_3", "heat_4", "heat_5", "heat_6", "heat_7", "heat_8")


  data$time <- as.numeric(data$time)/(60*60)
  data$ambient_temp <- as.numeric(data$ambient_temp)
  data$inner_temp <- as.numeric(data$inner_temp)
  data$heat_flow_1 <- as.numeric(data$heat_flow_1)*10^6
  data$heat_flow_2 <- as.numeric(data$heat_flow_2)*10^6
  data$heat_flow_3 <- as.numeric(data$heat_flow_3)*10^6
  data$heat_flow_4 <- as.numeric(data$heat_flow_4)*10^6
  data$heat_flow_5 <- as.numeric(data$heat_flow_5)*10^6
  data$heat_flow_6 <- as.numeric(data$heat_flow_6)*10^6
  data$heat_flow_7 <- as.numeric(data$heat_flow_7)*10^6
  data$heat_flow_8 <- as.numeric(data$heat_flow_8)*10^6

  data$heat_1 <- as.numeric(data$heat_1)
  data$heat_2 <- as.numeric(data$heat_2)
  data$heat_3 <- as.numeric(data$heat_3)
  data$heat_4 <- as.numeric(data$heat_4)
  data$heat_5 <- as.numeric(data$heat_5)
  data$heat_6 <- as.numeric(data$heat_6)
  data$heat_7 <- as.numeric(data$heat_7)
  data$heat_8 <- as.numeric(data$heat_8)



  # removing NaNs
  data <- na.omit(data)

  data <- subset(data, data$time > start - kalorimeter_t0)
  data <- subset(data, data$time < end - kalorimeter_t0)

  heat_cols <- paste0("heat_", 1:8)

  data[heat_cols] <- lapply(data[heat_cols], as.numeric)


  data[heat_cols] <- lapply(
    data[heat_cols],
    function(x) x - x[1]
  )

  ch_id <- c("ch1", "ch2", "ch3", "ch4", "ch5", "ch6", "ch7", "ch8")

  index <- c(1:8)
  for (i in index) {
    var_name <- paste("df", i, sep = "_")
    assign(var_name, as.data.frame(cbind(ch_id[i], trt_id[i], data$time, data$ambient_temp, data$inner_temp,  data[[paste("heat_flow", i, sep = "_")]], data[[paste("heat", i, sep = "_")]])))
  }


  data <- rbind(df_1, df_2, df_3, df_4, df_5, df_6, df_7, df_8)
  names(data) <- c("channel", "trt", "time", "ambient_temp", "inner_temp" ,  "heat_flow", "heat")
  data$heat_flow <- as.numeric(data$heat_flow)
  data$heat <- as.numeric(data$heat)
  data$time <- as.numeric(data$time) + kalorimeter_t0
  data$ambient_temp <- as.numeric(data$ambient_temp)
  data$inner_temp <- as.numeric(data$inner_temp)
  data$trt <- factor(data$trt)

  plot <- ggplot2::ggplot(data = data, aes(x = time, y = heat_flow, color = trt)) +
    geom_point(size = 1) +
    theme_minimal() +
    labs(y = "Heat flow per Vial [µW]", x = "time [h]") +
    theme(
      text = element_text(family = "serif", size = 20),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.border = element_rect(color = "darkgrey", fill = NA, size = 0.5),
      axis.line = element_line(color = "black")
    )
  #print(plot)

  return(list(data = data, plot = plot))
}
