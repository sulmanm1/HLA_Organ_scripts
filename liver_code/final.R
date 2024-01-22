library(tidyverse)
library(lubridate)
library(survival)
library(ggfortify)
library(ranger)
library(anytime)
library(icenReg)
library(survminer)

load_data <- function(){
  df <- read.csv('data/data.csv')
  df$Age.at.Transplant <- as.numeric(df$Age.at.Transplant)
  df <- df %>% filter(!Age.at.Transplant == "Failing PCR x2 - Try nanopore when available" | !is.na(Age.at.Transplant))
  df <- df %>% filter(!is.na(A.Mismatch.No.))
  df$FullName<-paste(df$First.Name, df$Last.Name)
  df$FullName
  return(df)
}

add_time <- function(df){
  df$Date.Anesthesia.Start<-as.Date(df$Date.Anesthesia.Start, format="%m/%d/%Y")
  df$status<-ifelse(is.na(df$Date.of.Graft.failure), 0, 1)
  df$Date.of.Graft.failure<-ifelse(is.na(df$Date.of.Graft.failure), "01/15/2024", df$Date.of.Graft.failure)
  df$Date.of.Graft.failure<-as.Date(df$Date.of.Graft.failure, format="%m/%d/%Y")
  df$Time<-df$Date.of.Graft.failure-df$Date.Anesthesia.Start
  df<-df[!duplicated(df$FullName),]
  return(df)
}

df<-load_data()
df<-add_time(df)
df<- df %>% filter(Time > 0)


colnames<-c("A.Mismatch.No.", "B.Mismatch.No.", "C.Mismatch.No.", "ClassI...Mismatch.All", "DRB1.Mismatch.No.","DRB345.Mismatch.No.", "DP.Mismatch.No.","DQ..Mismatch.No.", "ClassII...Mismatch.All")
filtered_df <- df[, colnames]
df_long <- pivot_longer(filtered_df, cols = everything(), names_to = "Variable", values_to = "Values")
ggplot(df_long, aes(sample = Values)) +
  geom_qq() +
  geom_qq_line(color = "red") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal() +
  xlab("Theoretical Quantiles") +
  ylab("Sample Quantiles")

#df<-df%>%filter(Time>365)
surv_obj<-Surv(df$Time, df$status)
survfit_obj<-survfit(surv_obj~1, data = df)
ggsurvplot(survfit_obj, fun = 'cumhaz')


coxmodel<-coxph(Surv(Time, status)~ A.Mismatch.No.
                  + B.Mismatch.No.
                  + C.Mismatch.No.
                  + DQ..Mismatch.No.
                  + DRB1.Mismatch.No.
                  + DRB345.Mismatch.No.
                  + DP.Mismatch.No.
                  + Age.at.Transplant, data=df)

ggcoxzph(cox.zph(coxmodel))
survfit_data1 <- survfit_obj$surv
survfit_data2 <- survfit(coxmodel)$surv
aaregstuff
# Create time points
time_points <- survfit_obj$time

# Create a dataframe for plotting
plot_data <- data.frame(
  Time = rep(time_points, 2), # Repeat time points for both survival curves
  Survival_Probability = c(survfit_data1, survfit_data2),
  Model = rep(c("Survfit", "coxmodel"), each = length(time_points))
)

# Create a custom survival plot using ggplot2
ggplot(plot_data, aes(x = Time, y = Survival_Probability, color = Model)) +
  geom_step() +
  labs(
    title = "Survival Curves Comparison",
    x = "Time",
    y = "Survival Probability",
    color = "Model"
  ) +
  theme_minimal()


aaregmodel<-aareg(formula=Surv(Time, status)~ A.Mismatch.No.
                + B.Mismatch.No.
                + C.Mismatch.No.
                + DQ..Mismatch.No.
                + DRB1.Mismatch.No.
                + DRB345.Mismatch.No.
                + DP.Mismatch.No.
                + Age.at.Transplant, data=df)

aaregmodel
autoplot(aaregmodel)

modelaft<-function(df, var){
  types=c("logistic", "weibull", "lognormal", "gaussian", "loglogistic", "exponential")
  AIC.data<-data.frame(
    model=c(),
    AIC_val=c()
  )
  for (model in types){
    print(model)
    aftmodel <- survreg(Surv(Time, status) ~ var + Age.at.Transplant, data = df, dist = model)
    model.data<-data.frame(
      model=c(model),
      AIC_val=c(AIC(aftmodel))
    )
    print(summary(aftmodel))
    AIC.data<-rbind(AIC.data, model.data)
  }
  return(AIC.data)
}

modelaft(df, var)

