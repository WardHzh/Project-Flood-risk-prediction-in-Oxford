if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

if (!requireNamespace("zoo", quietly = TRUE)) {
  install.packages("zoo")
}

if (!requireNamespace("pracma", quietly = TRUE)) {
  install.packages("pracma")
}

if (!requireNamespace("lubridate", quietly = TRUE)) {
  install.packages("lubridate")
}

if (!requireNamespace("tibble", quietly = TRUE)) {
  install.packages("tibble")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

#load packages
library(dplyr)
library(tidyverse)
library(zoo)
library(pracma)
library(lubridate)
library(tibble)
library(ggplot2)

#read csv file into dataframe with "station_id" file name
retrieve_data <- function(station_id) {
  filename_base <- "River height/"
  filename <- paste0(filename_base, station_id, ".csv")
  read.csv(filename)
}

#Arrange dataframes with only needed columns, change the format of column dateTime
dataArranging <- function(station_id) {
  df <- retrieve_data(station_id)

  #datetime type, filter missing value, delete duplicated dates
  df <- df %>% mutate(date = as.POSIXct(sub("T", " ", dateTime), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
    select(date, value) %>%
    unique() %>%
    arrange(date) %>%
    filter(!is.na(value)) %>%
    distinct(date, .keep_all = TRUE)

  #complete missing dates
  startTime <- min(df$date)
  endTime <- max(df$date)
  completeTimeSeries <- seq(from = startTime, to = endTime, by = "15 min")
  complete_df <- data.frame(date = completeTimeSeries)
  merged_df <- merge(complete_df, df, by = "date", all.x = TRUE) %>%
    arrange(date) %>%
    unique() %>%
    mutate(station_id=station_id)
  return(merged_df)
}


#interpolate NA with methods with linear and LOESS regression.
PCHIPinterpolate <- function(x, y) {
  naIndices <- which(is.na(y))
  nonNaIndices <- which(!is.na(y))

  yInterpolated <- interp1(as.numeric(x[nonNaIndices]),
                           y[nonNaIndices],
                           as.numeric(x),
                           method = "linear")

  y[naIndices] <- yInterpolated[naIndices]
  return(y)
}

LOESSinterpolate <- function(x,y, windowSize = 193) {
  df <- data.frame(date = x, value = y)
  naIndices <- which(is.na(df$value))

  #only interpolate for the datapoints where within the centred window,
  #more than half are non-missing values
  fillNaLoess <- function(index, data, dates, windowSize) {
    start <- max(1, index - floor(windowSize / 2))
    end <- min(length(data), index + floor(windowSize / 2))
    windowData <- data[start:end]
    windowDates <- dates[start:end]

    if (sum(!is.na(windowData)) < floor(windowSize / 2)) {
      return(NA)
    }

    loessFit <- loess(windowData ~ as.numeric(windowDates), span = 0.9, na.action = na.exclude)
    return(predict(loessFit, newdata = as.numeric(dates[index])))
  }

  for (i in naIndices) {
    df$value[i] <- fillNaLoess(i, df$value, df$date, windowSize)
  }

  return(df$value)
}

interpolateNAs <- function(df) {
  df$newValue <- LOESSinterpolate(df$date, df$value)
  df$newValue <- PCHIPinterpolate(df$date, df$newValue)
  df <- df %>% select(date, value, newValue, station_id)
  return(df)
}


#Sub-functions to be used in function "outlierDetection", calculating the
#differences mentioned in section 3.1
avgDifferenceMid <- function(vec) {
  index <- floor((length(vec) + 1) / 2)
  medianPositions <- which(vec < (median(vec) + 0.015) & vec > (median(vec) - 0.015))
  d <- abs(medianPositions - index)
  distance <- min(d)

  if (distance == 0) {
    return(0)
  } else {
    ratio <- abs((vec[index] - vec[medianPositions]) / d)
    return(max(ratio))
  }
}

avgDifferenceLeft <- function(vec) {
  index <- 1
  medianPositions <- which(vec < (median(vec) + 0.015) & vec > (median(vec) - 0.015))
  d <- abs(medianPositions - index)
  distance <- min(d)

  if (distance == 0) {
    return(0)
  } else {
    ratio <- abs((vec[index] - vec[medianPositions]) / d)
    return(max(ratio))
  }
}

avgDifferenceRight <- function(vec) {
  index <- length(vec)
  medianPositions <- which(vec < (median(vec) + 0.015) & vec > (median(vec) - 0.015))
  d <- abs(medianPositions - index)
  distance <- min(d)

  if (distance == 0) {
    return(0)
  } else {
    ratio <- abs((vec[index] - vec[medianPositions]) / d)
    return(max(ratio))
  }
}

#function "outlierDetection", above subfunctions appear in rollapply.
#it's a function to generate column "Difference", removed detected outliers whose
#"Difference" value exceeds 0.03
outlierDetection <- function(df, windowSize = 193) {
  left_df <- df %>%
    slice(1:floor(3*(windowSize - 1)/2)) %>%
    mutate(Difference = rollapply(newValue, width = windowSize, FUN = avgDifferenceLeft, fill = NA, align = "left")) %>%
    slice(1:floor((windowSize - 1)/2))

  middle_df <- df %>%
    mutate(Difference = rollapply(newValue, width = windowSize, FUN = avgDifferenceMid, fill = NA, align = "center"))  %>%
    slice((floor((windowSize + 1)/2)) : (nrow(df) - floor((windowSize - 1)/2)))

  right_df <- df %>%
    slice((nrow(df) + floor((5- 3 * windowSize) / 2)) : nrow(df)) %>%
    mutate(Difference = rollapply(newValue, width = windowSize, FUN = avgDifferenceRight, fill = NA, align = "right"))
  right_df <- right_df %>% slice((nrow(right_df) - floor((windowSize - 2) / 2)) : nrow(right_df))
  merged_df <- bind_rows(left_df, middle_df, right_df)

  threshold <- 0.03
  merged_df <- merged_df %>% mutate(newValue = if_else(Difference > threshold, NA_real_, value))
  merged_df$height <- PCHIPinterpolate(merged_df$date, merged_df$newValue)
  merged_df <- merged_df %>% select(date, value, newValue, height, station_id)

  return(merged_df)
}

#function to: detect outliers, set to NA, re-interpolate, save as new file
riverThames <- c("7027", "7034", "7035", "7036", "7037", "7038", "7046", "7048", "7049", "7055", "7056", "7057", "7072", "7405", "7073", "7083", "7085", "7086")
riverOck <- c("7402", "7078", "7077")
riverCherwell <- c("7059", "1414", "7063", "7064", "7066", "7071")
streamsBrooks <- c("7039", "7001", "7075", "7074", "7076", "7003")
riverStations <- c(riverThames, riverOck, riverCherwell, streamsBrooks)
riverStationsList <- tibble::tibble(station_id = c(riverThames, riverOck, riverCherwell, streamsBrooks))

#transfer timescales: hourly (4); dayly (96); weekly(672).
timeScaleTransfer <- function(df, timeScale) {
  df <- df %>% select(date, height, station_id) %>%
    mutate(movingAverage = rollapply(height, width = timeScale+1, FUN = mean, fill = NA, align = "center", partial = FALSE))

  startTime <- min(df$date)
  endTime <- max(df$date)
  timeGap <- 15 * timeScale
  timeGap_str <- paste(timeGap, "mins")
  completeTimeSeries <- seq(from = startTime, to = endTime, by = timeGap_str)

  complete_df <- data.frame(date = completeTimeSeries)
  merged_df <- inner_join(df, complete_df, by = "date") %>%
    arrange(date) %>% drop_na() %>% select(date, movingAverage, station_id) %>%
    rename(height = movingAverage)

  return(merged_df)
}
