library(icenReg)
library(survminer)
library(tidyverse)
library(anytime)
library(ggfortify)

load_data<-function(){
  df<-read.csv("data/data_clean.csv")
  df$Time.upper.G1 <- ifelse(df$Time.lower.G1 == df$Time.upper.G1, Inf, df$Time.upper.G1)
  df$Time.upper.G2 <- ifelse(df$Time.lower.G2 == df$Time.upper.G2, Inf, df$Time.upper.G2)
  df$Time.upper.G3 <- ifelse(df$Time.lower.G3 == df$Time.upper.G3, Inf, df$Time.upper.G3)
  df$Transplant.Date<-anydate(df$Transplant.Date)
  df$Date.of.Graft.failure<-anydate(df$Date.of.Graft.failure)
  df<-df %>% filter(!df$Pre.DSAs=="Yes")
  df$time2fail<-ifelse(is.na(df$Date.of.Graft.failure), anydate("2024-01-15")-df$Transplant.Date, df$Date.of.Graft.failure-df$Transplant.Date)
  df$status<-ifelse(is.na(df$Date.of.Graft.failure), 0, 1)
  return(df)
}

create_survival_object<-function(df, lower, upper){
  surv_obj <- with(df, Surv(lower, upper, type = "interval2"))
}

#load and transform data
df<-load_data()

plot(survfit(Surv(time2fail, status) ~ 1, data = df), xlab = "Time (days)", ylab = "Survival Probability", main = "Kaplan-Meier Estimate")

cox<-coxph(Surv(time2fail, status) ~ A.Mismatch.No.
      + B.Mismatch.No.
      + C.Mismatch.No.
      + DRB1.Mismatch.No.
      + DRB345.Mismatch.No.
      + DP.Mismatch.No.
      + DQ.Mismatch.No., data = df)

autoplot(survfit(Surv(time2fail, status) ~ 1, data = df))

coxA<-coxph(Surv(time2fail, status) ~ A.Mismatch.No., data=df)
coxB<-coxph(Surv(time2fail, status) ~ B.Mismatch.No., data=df)
coxC<-coxph(Surv(time2fail, status) ~ C.Mismatch.No., data=df)
coxDRB1<-coxph(Surv(time2fail, status) ~ DRB1.Mismatch.No., data=df)
coxDRB345<-coxph(Surv(time2fail, status) ~ DRB345.Mismatch.No., data=df)
coxDP<-coxph(Surv(time2fail, status) ~ DP.Mismatch.No., data=df)
coxDQ<-coxph(Surv(time2fail, status) ~ DQ.Mismatch.No., data=df)
coxc1<-coxph(Surv(time2fail, status) ~ ClassII....Mismatch..All, data=df)
coxc2<-coxph(Surv(time2fail, status) ~ ClassI....Mismatch..All, data=df)

cox.zph(coxc2)
cox.zph(coxc1)
cox.zph(coxA)
cox.zph(coxB)
cox.zph(coxC)
cox.zph(coxDP)
cox.zph(coxDRB1)
cox.zph(coxDRB345)
cox.zph(coxDQ)

coxDP<-coxph(Surv(time2fail, status) ~ tt(DP.Mismatch.No.), data=df)
DP<-aareg(formula=Surv(time2fail, status)~DP.Mismatch.No., data=df)
autoplot(DP)

median_var1 <- median(df$DQ.Mismatch.No., na.rm = TRUE)
df$DQ_mismatch <- ifelse(df$DQ.Mismatch.No. > median_var1, "Top 50%", "Bottom 50%")
km_fit <- survfit(Surv(time2fail, status) ~ DQ_mismatch, data = df)
ggsurvplot(km_fit, 
           data = df, 
           xlab = "Time to Failure", 
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curves",
           palette = c("#00AFBB", "#FC4E07"))
log_rank_test <- survdiff(Surv(time2fail, status) ~ DQ_mismatch, data = df)
log_rank_test
#create survival object

surv_obj2<-create_survival_object(df, df$Time.lower.G2, df$Time.upper.G2)
surv_obj<-create_survival_object(df, df$Time.lower.G1, df$Time.upper.G1)
survfit1<-survfit(surv_obj~1, data=df)
survfit2<-survfit(surv_obj2~1, data=df)
ggsurvplot(survfit1, title = "Survival Grade 1")
ggsurvplot(survfit2, title = "Survival Grade 2")

#univariate coxph
icA<-ic_sp(surv_obj ~ A.Mismatch.No., data = df, model = 'ph', bs_samples = 500)
icB<-ic_sp(surv_obj ~ B.Mismatch.No., data = df, model = 'ph', bs_samples = 500)
icC<-ic_sp(surv_obj ~ C.Mismatch.No., data = df, model = 'ph', bs_samples = 500)
icDRB1<-ic_sp(surv_obj ~ DRB1.Mismatch.No., data = df, model = 'ph', bs_samples = 500)
icDRB345<-ic_sp(surv_obj ~ DRB345.Mismatch.No., data = df, model = 'ph', bs_samples = 500)
icDP<-ic_sp(surv_obj ~ DP.Mismatch.No., data = df, model = 'ph', bs_samples = 500)
icDQ<-ic_sp(surv_obj ~ DQ.Mismatch.No., data = df, model = 'ph', bs_samples = 500)
icC2<-ic_sp(surv_obj ~ ClassII....Mismatch..All, data = df, model = 'ph', bs_samples = 500)
icC1<-ic_sp(surv_obj ~ ClassI....Mismatch..All, data = df, model = 'ph', bs_samples = 500)
diag_covar(object=icA, model='po')


#list of covariates
ic_sp(surv_obj ~ A.Mismatch.No.
  + B.Mismatch.No.
  + C.Mismatch.No.
  + DRB1.Mismatch.No.
  + DRB345.Mismatch.No.
  + DP.Mismatch.No.
  + DQ.Mismatch.No., data = df, model = 'ph', bs_samples = 500)

summary(coxph_fit)

