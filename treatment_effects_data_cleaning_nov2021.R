#EdB started here on 4-11-2021
#data cleaning of treatment to assess treatment effects
library(tidyverse)
autoimmunity_dataframe <- read_excel("TBR PAH cohort study/Origional data with minor transformations/Original data/autoimmunity.dataframe.xlsx")

#first step = rename the colnames with first row
colnames(autoimmunity_dataframe) <- autoimmunity_dataframe[1,]
autoimmunity_dataframe <- autoimmunity_dataframe[-c(1),]

treatment_df <- autoimmunity_dataframe %>% select(starts_with("dh_") | ends_with("_base"), c("vasoresponder", "vasoresponder2", "initial_therapy", "id_cohort.x"))

#problem: dates are imported weirdly
library(lubridate)
#see if it's easy to fix.....
treatment_df$dh_pde5_date_start <- as.numeric(treatment_df$dh_pde5_date_start)
treatment_df$dh_pde5_date_start <- as.Date(treatment_df$dh_pde5_date_start, origin="1899-12-30")

#so this has worked!
#apply this for the entire df
treatment_df[,c(4,6,8,10,12)] <- lapply(treatment_df[,c(4,6,8,10,12)], as.numeric)
dates_normal <- function(x) {as.Date(x, origin = "1899-12-30")}
treatment_df$dh_era_date_start <- dates_normal(treatment_df$dh_era_date_start)
treatment_df[,c(4,6,8,10,12)] <- lapply(treatment_df[,c(4,6,8,10,12)], dates_normal)

#all the dates are fixed now!

#make the following categories:
#1) single therapy
#2) combination therapy
#3) triple therapy
#4) IV therapy 

#first get a col of IV therapy
summary(as.factor(treatment_df$dh_pa_drug))
iv_meds <- c("Epoprostenol iv", "Iloprost iv", "Treprostinil iv")
treatment_df$IV_treatment <- ifelse(treatment_df$dh_pa_drug %in% iv_meds, "IV", "no")

#now get numbers for therapy
#note to self: if there is a NA for the recorded drug, but the recorded drug
#does have a start date --> we can safely assume this drug was being used.
#first for ERA
summary(as.factor(treatment_df$dh_era_drug))
treatment_df$era1 <- ifelse(!is.na(treatment_df$dh_era_drug), 1, 0)
#also add if a date was given for start of drug, but drug isn't noted as being given
treatment_df$era1 <- ifelse(!is.na(treatment_df$dh_era_date_start), 1, treatment_df$era1)

#now for PDE5
treatment_df$pde1 <- ifelse(!is.na(treatment_df$dh_pde5_drug), 1, 0)
#also add if a date was given for start of drug, but drug isn't noted as being given
treatment_df$pde1 <- ifelse(!is.na(treatment_df$dh_pde5_date_start), 1, treatment_df$pde1)

#now for PA drugs
treatment_df$pa1 <- ifelse(!is.na(treatment_df$dh_pa_drug), 1, 0)
treatment_df$pa1 <- ifelse(!is.na(treatment_df$dh_pa_date_start), 1, treatment_df$pa1)

#now for PRA drugs
summary(as.factor(treatment_df$dh_pra_drug))
treatment_df$pra1 <- ifelse(!is.na(treatment_df$dh_pra_drug), 1, 0)
treatment_df$pra1 <- ifelse(!is.na(treatment_df$dh_pra_date_start), 1, treatment_df$pra1)

#now sgc drugs
summary(as.factor(treatment_df$dh_sgc_drug))
treatment_df$sgc1 <- ifelse(!is.na(treatment_df$dh_sgc_drug), 1, 0)
treatment_df$sgc1 <- ifelse(!is.na(treatment_df$dh_sgc_date_start), 1, treatment_df$sgc1)

#now at drugs
summary(as.factor(treatment_df$dh_at_drug))
treatment_df$at1 <- ifelse(!is.na(treatment_df$dh_at_drug), 1, 0)
treatment_df$at1 <- ifelse(!is.na(treatment_df$dh_at_date_start), 1, treatment_df$at1)

#also include vasoresponders (used vasoresponder2 as less NAs) if no treatment is recorded. 
#If people have no recorded treatment but are vasoresponders, we can assume they were on CCBs
treatment_df$at1 <- ifelse(treatment_df$vasoresponder2 %in% c("vasoresponder"), 1, treatment_df$at1)
#(this gave an eror when using ==, fixed it by using %in%)

#now get al col with the total therapy #
treatment_df$total_therapy <- treatment_df$era1 + treatment_df$pa1 + treatment_df$pra1 + treatment_df$sgc1 + treatment_df$at1 + treatment_df$pde1
#================================================================
#now also check if the recorded initial therapy (in case other missing data) corresponds with the calculation above
summary(as.factor(treatment_df$initial_therapy))
#intial therapy overlaps nicely, only need to check if patients with initial combination therapy aren't 
#also counted as being on single therapy
treatment_df %>% filter(initial_therapy == "combination therapy") %>% summary(as.factor(total_therapy))
#nothing is classed as <2! --> so that's okay
#==================================
#get a single column with all the variables
treatment_df$tripple_therapy <- ifelse(treatment_df$total_therapy >= 3, "tripple", "no")
treatment_df$combination_therapy <- ifelse(treatment_df$total_therapy == 2, "combination", "no")
treatment_df$single_therapy <- ifelse(treatment_df$total_therapy == 1, "single", "no")

#now generate a col that has all the data for the comparison
treatment_df$type_of_therapy <- treatment_df$IV_treatment
treatment_df$type_of_therapy <- ifelse(treatment_df$type_of_therapy == "no", treatment_df$tripple_therapy, treatment_df$type_of_therapy)
treatment_df$type_of_therapy <- ifelse(treatment_df$type_of_therapy == "no", treatment_df$combination_therapy, treatment_df$type_of_therapy)
treatment_df$type_of_therapy <- ifelse(treatment_df$type_of_therapy == "no", treatment_df$single_therapy, treatment_df$type_of_therapy)


#people with a ''no'' for treatment have no recoded treatment --> note this correctly
treatment_df$type_of_therapy <- ifelse(treatment_df$type_of_therapy == "no", "no_medication_recorded", treatment_df$type_of_therapy)
summary(as.factor(treatment_df$type_of_therapy))

#get a nice dataframe to share
saveRDS(treatment_df, file = "treatment_per_patient_PAH_cohort_November2021.RDS")
df <- treatment_df[,c(22:30, 34)]

write.csv2(df, file="treatment_per_patient_inPAHCohortNov2021.csv")

#this is ready to include in the Cox-PH!

#brief check with Chi-square
df$type_of_therapy <- as.factor(df$type_of_therapy)
definitive.clusterings.3k.PAM.labels.3March2021 <- read.csv("~/My publications/Data new analyses March 2021/definitive clusterings 3k PAM labels 3March2021.csv", sep=";")
names(definitive.clusterings.3k.PAM.labels.3March2021)[[1]] <- "id_cohort.x"
df <- merge(df, definitive.clusterings.3k.PAM.labels.3March2021, by="id_cohort.x")
df$definitive_labels <- as.factor(df$definitive_labels)
chisq.test(df$definitive_labels, df$type_of_therapy)

library(Publish)
abc <- univariateTable(definitive_labels ~ type_of_therapy, data = df, show.totals = TRUE, column.percent = TRUE)
abc <- summary(abc)
write.csv2(abc, file="cluster_differences_per_treatment_type.csv")
