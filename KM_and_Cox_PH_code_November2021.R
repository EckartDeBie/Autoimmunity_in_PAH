#Cox-PH and Kaplan Meier code autoimmunity in PAH paper
#Original code produced by Stefan Gräf and Chris Wallace
#Edited/adjusted by Chris Wallace and Eckart de Bie on 8 and 10 November 2021

options(stringsAsFactors=FALSE)
library("survival")
library("survminer")
library("survMisc")
library("lubridate")
library("ggsci")
library("magrittr")
library("tidyverse")
library("openxlsx")

analysis.dir="../"

## Cohort Study start date
start_of_study <- ymd("2013-06-01")

## Census date
date_of_census <- ymd("2021-01-31")

oc_data <- read_excel("/mortality PAH adjusted for age and sex/oc_data_2021-03-20_corrected.xlsx")
head(oc_data)
table(oc_data$cohort_group) # all cases
table(oc_data$cluster) # clusters 1, 2, 3

#quick check if clusters correspond
definitive.clusterings.3k.PAM.labels.3March2021 <- read.csv("~/My publications/Data new analyses March 2021/definitive clusterings 3k PAM labels 3March2021.csv", sep=";")
check_df <- merge(definitive.clusterings.3k.PAM.labels.3March2021, oc_data, by.x = "Row.names", by.y = "id_cohort")
same_labels <- check_df$definitive_labels == check_df$cluster
summary(same_labels)
#everything is true so the labels are correct!

#run the KM code, but first check length in cohort
ggplot(oc_data, aes(color=as.factor(event))) +
  geom_segment(aes(x=ttlt, xend=ttrc,
                   y=as.numeric(rownames(oc_data)),
                   yend=as.numeric(rownames(oc_data))), lwd=0.3) +
  geom_point(aes(ttv, as.numeric(rownames(oc_data))), shape="+", colour = "black") +
  scale_colour_nejm(name="Event") +
  xlab(paste0("Days (0 = start of study (", start_of_study, "))")) +
  ylab("Cohort patients") +
  geom_vline(xintercept=c(-10*365.25,0), linetype="dashed", color="darkgrey")
ggsave(paste0(analysis.dir, "oc_data_time_in_study_", Sys.getenv("OC_REL"), "_corrected.pdf"))


## Kaplan-Meier plots
plotKM <- function(.df=tbl, out.dir=analysis.dir, suffix=name) {
  
  ggplot(.df, aes(color=as.factor(event))) +
    geom_segment(aes(x=ttlt, xend=ttrc,
                     y=as.numeric(rownames(.df)),
                     yend=as.numeric(rownames(.df))), lwd=0.3) +
    geom_point(aes(ttv, as.numeric(rownames(.df))), shape="+", colour = "black") +
    scale_colour_nejm(name="Event") +
    xlab(paste0("Days (0 = start of study (", start_of_study, "))")) +
    ylab("Cohort patients") +
    geom_vline(xintercept=c(-10*365.25,0), linetype="dashed", color="darkgrey")
  ggsave(paste0(out.dir, "oc_data_time_in_study_", Sys.getenv("OC_REL"), "_corrected", suffix, ".pdf"))
  
  fit <- with(.df,
              survfit(Surv(tte, event) ~ cluster))
  summary(fit)$table
  km <- ggsurvplot(fit, .df,
                   conf.int=TRUE,
                   test.for.trend=FALSE,
                   pval=TRUE,
                   pval.method = TRUE,
                   ##surv.median.line="hv",
                   risk.table=TRUE,
                   risk.table.height=3,
                   xscale="d_y",
                   palette="nejm",
                   legend = c(.95,.9),
                   legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3"),
                   xlab="Time [years after diagnosis]")
  return(km)
  ## ggsave(paste0(out.dir, "survival_analysis_", Sys.getenv("OC_REL"), suffix, ".pdf"), km$plot, width=15, height=8)
  ## ggsave(paste0(out.dir, "survival_analysis_", Sys.getenv("OC_REL"), suffix, "_table.pdf"), km$table, width=15.35, height=2)
  
}

## all left truncated right censored (ltrc)
tbl <- filter(oc_data, diagnosis_verified %in% c("HPAH", "IPAH", "PVOD", "PCH"))
print(tbl) #, width=Inf)
select(tbl, sex, cluster) %>% table()
select(tbl, diagnosis_verified, cluster) %>% table()
basename <- "_ltrc_PAH-PVOD-PCH"
p=plotKM(tbl, analysis.dir, basename)
p$plot
p$table

tbl$cluster = as.factor(tbl$cluster)
m0=coxph(Surv(tte, event) ~ cluster, data=tbl)
m1=coxph(Surv(tte, event) ~ cluster + sex, data=tbl)
m2=coxph(Surv(tte, event) ~ sex + age_diagnosis + cluster , data=tbl)

survdiff(Surv(tte, event) ~ cluster, data=tbl)

## Pairwise survdiff

sobj <- with(tbl, Surv(tte, event))
names(fit)

res <- pairwise_survdiff(Surv(tte, event) ~ cluster,
                         data = tbl, p.adjust.method="BH")
res


for (i in seq(10, 20, 5)) {
  print(i)
  .tbl <- filter(tbl, tdxse<=i*365.25)
  print(dim(.tbl))
  name <- paste0(basename, "_tdxse_le", sprintf("%02d", i), "yrs")
  plotKM(.tbl, analysis.dir, name)
  
  res <- pairwise_survdiff(Surv(tte, event) ~ cluster,
                           data = .tbl, p.adjust.method="BH")
  print(res)
}


## Fit a Cox proportional hazards model
## unadjusted
fit.coxph1 <- coxph(Surv(tte, event) ~ cluster,
                    data = tbl)
ggforest(fit.coxph1, data = tbl)
ggsave("coxph_unadjusted.pdf", width=8, height=2)

## add treatment effects
tbl <- mutate(tbl, age_group = ifelse(age_diagnosis >= 50, "old", "young"))
treat = readRDS("C:/Users/Gebruiker/Documents/My publications/revisions/treatment_per_patient_PAH_cohort_November2021.RDS")
head(as.data.frame(treat))
head(tbl)
tbl2=merge(tbl, treat[,c("id_cohort.x","type_of_therapy")], by.x="id_cohort", by.y="id_cohort.x")
stopifnot(nrow(tbl)==nrow(tbl2))
head(as.data.frame(tbl2))

## adjusted, with age grouped
fit.coxph2 <- coxph(Surv(tte, event) ~  sex + age_group + type_of_therapy + cluster,
                    data = tbl2)
anova(fit.coxph2)
ggforest(fit.coxph2, data = tbl2)

## adjusted with age as continuous (preferable)
fit.coxph2.a <- coxph(Surv(tte, event) ~  sex + age_diagnosis + type_of_therapy + cluster,
                      data = tbl2)
anova(fit.coxph2.a)
summary(fit.coxph2.a)
ggforest(fit.coxph2.a, data = tbl2)
ggsave("coxph_adjusted.pdf", width=8, height=8)

#this all runs fine, now just adjust the labels so that they correspond to the labels used in the paper
tbl3 <- tbl2
tbl3$type_of_therapy <- ifelse(tbl3$type_of_therapy == "tripple", "triple", tbl3$type_of_therapy)
tbl3$type_of_therapy <- ifelse(tbl3$type_of_therapy == "no_medication_recorded", "no medication recorded", tbl3$type_of_therapy)
#also name the clusters like they're named in the paper
tbl3$cluster <- ifelse(tbl3$cluster == "1", "High autoantibody cluster", tbl3$cluster)
tbl3$cluster <- ifelse(tbl3$cluster == "2", "Low autoantibody cluster", tbl3$cluster)
tbl3$cluster <- ifelse(tbl3$cluster == "3", "Intermediate autoantibody cluster", tbl3$cluster)

#Now check if labels are nice in the Cox-PH
fit.coxph3 <- coxph(Surv(tte, event) ~  sex + age_diagnosis + type_of_therapy + cluster,
                    data = tbl3)
anova(fit.coxph3)
summary(fit.coxph3)
ggforest(fit.coxph3, data = tbl3)
ggsave("coxph_adjusted_labels_November_2021.pdf", width=11, height=8)
