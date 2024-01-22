library(tidyverse)
library(tidyr)
library(lubridate)
library(survival)
library(ggfortify)
library(ranger)

RmExtra<-function(df){
  df$Date.Anesthesia.Start<- gsub("\n", "", df$Date.Anesthesia.Start)
  df$Date.of.Graft.failure<- gsub("\n", "", df$Date.of.Graft.failure)
  df$Date.Anesthesia.Start<- gsub(" ", "", df$Date.Anesthesia.Start)
  df$Date.of.Graft.failure<- gsub(" ", "", df$Date.of.Graft.failure)

  df$censored=1
  df<-df %>%
    filter(!is.na(Date.Anesthesia.Start)) %>%
    filter(!is.na(ClassI...Mismatch.All))
  return(df)
}

CensorGraftTime<-function(df){
  df$censored=1
  today_date <- ymd(Sys.Date())
  df <- df %>%
    mutate(censored = ifelse(is.na(Date.of.Graft.failure) | Date.of.Graft.failure == "", 0, censored),
           Date.of.Graft.failure = ifelse(is.na(Date.of.Graft.failure) | Date.of.Graft.failure == "", today_date, Date.of.Graft.failure))
  return(df)
}

MFtransform<-function(df){
  # Filter out rows where 'Sex' is neither 'Male' nor 'Female'
  df <- df %>%
    filter(Sex %in% c("Male", "Female"))

  # Substitute 'Male' with 0 and 'Female' with 2
  df <- df %>%
    mutate(Sex = ifelse(Sex == "Male", 0, 2))

  return(df)
}

plot_FacetHisto <- function(df, output) {
  # Convert the dataframe to long format for numeric columns
  long_df <- gather(df, key = "variable", value = "value",
                    which(sapply(df, is.numeric)))

  # Add 'TCMR' to the long format dataframe
  long_df$TCMR <- df$TCMR[rep(1:nrow(df),
                                                                   times = length(which(sapply(df, is.numeric))))]

  # Create grouped histograms
  p_grouped <- ggplot(long_df, aes(x = value, fill = TCMR)) +
    geom_histogram(position = "stack", binwidth = 1) +  # Adjust binwidth as needed
    facet_wrap(~ variable, scales = "free") +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    labs(x = "Value", y = "Frequency", title = "Grouped Histograms by Likely Immune Rejection Status")

  print(p_grouped)
  ggsave(output, width = 10, height = 10, dpi = 300)
}
plot_DensityLine <- function(df, output) {
  # Load necessary libraries
  library(ggplot2)
  library(tidyr)

  # Convert the dataframe to long format for numeric columns
  long_df <- gather(df, key = "variable", value = "value", which(sapply(df, is.numeric)))

  # Add 'TCMR' to the long format dataframe
  long_df$TCMR <- df$TCMR[rep(1:nrow(df), times = length(which(sapply(df, is.numeric))))]

  # Number of unique categories plus one for 'Total'
  num_categories <- length(unique(df$TCMR)) + 1

  # Create density line graphs including total density
  p_grouped <- ggplot(long_df, aes(x = value, color = TCMR)) +
    geom_density() +  # Density for each category
    geom_density(data = long_df, aes(x = value, color = "Total"), linetype = "dashed") +  # Total density
    scale_color_brewer(palette = "Set1", name = "Category") +
    guides(color = guide_legend(override.aes = list(linetype = c(rep("solid", num_categories - 1), "dashed")))) +
    facet_wrap(~ variable, scales = "free") +
    theme_minimal() +
    labs(x = "Value", y = "Density", title = "Density Line Graphs by Likely Immune Rejection Status")

  print(p_grouped)
  ggsave(output, width = 10, height = 10, dpi = 300)
}

plot_CountLine <- function(df, cat_var, num_bins = 10, output) {
  if(!cat_var %in% names(df)) {
    stop("Categorical variable not found in the dataframe.")
  }

  # Ensure num_bins is numeric and positive
  if(!is.numeric(num_bins) || num_bins <= 0) {
    stop("num_bins must be a positive numeric value.")
  }

  # Convert the dataframe to long format for numeric columns
  long_df <- gather(df, key = "variable", value = "value", which(sapply(df, is.numeric)))

  # Check if long_df$value is numeric
  if(!is.numeric(long_df$value)) {
    stop("Values to be binned must be numeric.")
  }

  # Add the categorical variable to the long format dataframe
  long_df[cat_var] <- df[[cat_var]][rep(1:nrow(df), times = length(which(sapply(df, is.numeric))))]

  # Bin values into fewer bins
  long_df$value_bin <- cut(long_df$value, breaks = num_bins, labels = FALSE)

  # Calculate counts
  count_df <- long_df %>%
    group_by(variable, value_bin, .data[[cat_var]]) %>%
    summarize(count = n(), .groups = 'drop')

  # Total counts
  total_counts <- count_df %>%
    group_by(variable, value_bin) %>%
    summarize(count = sum(count), .groups = 'drop')

  # Create count-based line graphs including total count
  p_grouped <- ggplot() +
    geom_line(data = count_df, aes(x = value_bin, y = count, color = .data[[cat_var]], group = .data[[cat_var]]), size = 1) +
    geom_line(data = total_counts, aes(x = value_bin, y = count, color = "Total", group = 1), linetype = "dashed", size = 1) +
    scale_color_brewer(palette = "Set1", name = "Category") +
    facet_wrap(~ variable, scales = "free") +
    theme_minimal() +
    labs(x = "Binned Value", y = "Count", title = paste("Count Line Graphs by", cat_var)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels

  print(p_grouped)
  ggsave(output, width = 12, height = 10, dpi = 300) # Adjust the plot size
}

AddVals<-function(df){
  df <- df %>%
    mutate(AMR = ifelse(AMR == "", 'Unavailable', AMR),
           TCMR = ifelse(TCMR != "", TCMR, 'Unavailable'))
  return(df)
}

AddTimes <- function(df){
  df$Date.Anesthesia.Start<-mdy(df$Date.Anesthesia.Start)
  df$Date.of.Graft.failure<-mdy(df$Date.of.Graft.failure)
  df$Time<-df$Date.of.Graft.failure-df$Date.Anesthesia.Start
  return(df)
}

MakeCategorical <- function(continuous_var, number_of_bins) {
  # Check if the input is numeric
  if(!is.numeric(continuous_var)) {
    stop("The variable must be numeric.")
  }

  # Check if number_of_bins is a positive integer
  if(!is.numeric(number_of_bins) || number_of_bins <= 0 || number_of_bins != round(number_of_bins)) {
    stop("Number of bins must be a positive integer.")
  }

  # Use cut to create categories
  categorical_var <- cut(continuous_var, breaks = number_of_bins, include.lowest = TRUE, labels = FALSE)

  return(categorical_var)
}

plot_SurvivalCurve <- function(df, output) {
  km <- with(df, Surv(Time, censored))
  km_fit <- survfit(km ~ TotalMismatch_cat, data = df)
  autoplot(km_fit)
}

plot_cox <- function(df, output) {
  # Fit the Cox proportional hazards model
  cox <- coxph(Surv(Time, censored) ~ A.Mismatch.No_cat
    + B.Mismatch.No_cat
    + DQ..Mismatch.No_cat + Age.at.Transplant_cat
    + DRB345.Mismatch.No_cat
    + C.Mismatch.No_cat
    + Sex, data = df)

  cox_fit <- survfit(cox)
  autoplot(cox_fit)
}

plot_aareg <- function(df, output) {
  # Fit the Cox proportional hazards model
  aareg <- aareg(Surv(Time, censored) ~ A.Mismatch.No_cat
    + B.Mismatch.No_cat
    + DQ..Mismatch.No_cat + Age.at.Transplant_cat
    + DRB345.Mismatch.No_cat
    + C.Mismatch.No_cat
    + DP.Mismatch.No_cat
    + DRB1.Mismatch.No_cat
    + Sex, data = df)
  autoplot(aareg)

}

perform_LogRankTests <- function(df, time_var, status_var, covariates) {
  results <- list()
  for(covar in covariates) {
    surv_obj <- Surv(df[[time_var]], df[[status_var]])
    test_result <- survdiff(surv_obj ~ df[[covar]])
    results[[covar]] <- test_result
  }
  return(results)
}