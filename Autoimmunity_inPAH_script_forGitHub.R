#Autoimmunity in PAH paper code for GitHub
##reproducible code for AAB analysis
## original by Eckart de Bie
##02/03/2020

## rerun and edited by Chris Wallace 02/03/2021
## some changes in style to favour data.table due to familiarity or to automate some manual steps
## otherwise should follow to original total_script_Eck.R
##code made ready for publication on GitHub on 01/05/21 By Eckart de Bie
#(by removing some comments & minor analyses not reported in the paper)


library(readxl)
library(data.table)
library(magrittr)
library(MatchIt)
library(tidyverse)
set.seed(42) # reproducility, not sure if matchit uses random


#First step match the cases and controls


data_dir="~/My publications/data re-run code" # where all the data are stored
data.Eoin.v1.4March <- read.csv2(file.path(data_dir,"data_Eoin_v1_4March.csv"), stringsAsFactors=FALSE)
autoimmune_profiling_2nd_batch_5_missing_samples_removed <-
  read_excel(file.path(data_dir,"autoimmune_profiling_2nd_batch__5_missing_samples_removed.xlsx")) %>%
  as.data.table()
dob_data_v1_20may <- read.csv2(file.path(data_dir,"dob_data_v1_20may.csv")) %>% as.data.table()
dob_data_v1_20may$DOB <- as.Date(dob_data_v1_20may$DOB, format = "%Y-%m-%d")

#===================================================

##the visit_date variable is the correct date of sampling for the cases and this
##is the one that should be matched on! the age at this time can be calculated
##with the DOB data and the autoimmune_profiling file

#===================================================

#match data from Eoin with our own data 2:1
data_Eoin <- data.Eoin.v1.4March
propensity_demographics <- as.data.table(autoimmune_profiling_2nd_batch_5_missing_samples_removed)
propensity_demographics[,visit_date:=as.Date(visit_date, format = "%Y-%m-%d")]

propensity_demographics$sample <- 1
propensity_Eoin <- as.data.table(data_Eoin[c(2,22,23)])

#generate numeric data for gender
table(propensity_demographics$sex,exclude=NULL) # no missing, all male / female
propensity_demographics[,S:=ifelse(sex=="male",1,0)]


table(propensity_Eoin$SEX,exclude=NULL) # male, female, and NA
propensity_Eoin[SEX=="M",S:=1]
propensity_Eoin[SEX=="F",S:=0]

##=============================================
##for this bit the new df needs to be cleaned where the ''age'' variable is the
##age at sampling and all other cols apart from the ones below can be removed
##the cases need a 1 for the sample_group variable if this one does not exist
##yet
names(propensity_demographics)


propensity_demographics[,diag_date:=as.Date(cfh_date_of_diagnosis,"YYYY-MM-DD")]
head(propensity_demographics[,.(AGE=age_diagnosis, diag_date, visit_date,
                                diff=as.numeric(visit_date - diag_date))])
propensity_demographics[,AGE:=age_diagnosis + as.numeric(visit_date-diag_date)/365.25]
prop_PAH_for_match <- propensity_demographics[,.(AGE, S,
                                                 sample=pah_subject_id_bridge_validated,
                                                 sample_group=1)]
summary(prop_PAH_for_match$AGE)
propensity_demographics[is.na(cfh_date_of_diagnosis)]
missing_ids=propensity_demographics[is.na(cfh_date_of_diagnosis)]$pah_subject_id_cohort_validated
dob_data_v1_20may[Patient %in% missing_ids]
propensity_demographics=merge(propensity_demographics, dob_data_v1_20may[,.(Patient,DOB)],
                              by.x="pah_subject_id_cohort_validated",
                              by.y="Patient",
                              all.x=TRUE)
propensity_demographics[pah_subject_id_cohort_validated %in% missing_ids]
propensity_demographics[pah_subject_id_cohort_validated %in% missing_ids, AGE:=(visit_date - DOB)/365.25]
propensity_demographics[,DOB:=NULL]
prop_PAH_for_match <- propensity_demographics[,.(AGE, S,
                                                 sample=pah_subject_id_cohort_validated,
                                                 sample_group=1)]
##====================================================
#give Eoins pts a number
propensity_Eoin$sample_group <- 0

#combine both df's for MatchIt
df_match <- rbind(prop_PAH_for_match[,.(sample,sample_group,AGE,S)],
                  propensity_Eoin[,.(sample,sample_group,AGE,S)])
lapply(df_match, is.numeric)
## PresenceNA <- lapply(df_match, is.na)

#the formula indicates there is NA data
#will omit NA data (now 4477 obs and 4 vars)
df_match_no_na <- na.omit(df_match)
table(df_match$sample_group)
table(df_match_no_na$sample_group) # lost some controls (expected) but no cases, so ok
df_match[sample_group==1 & is.na(sample)]


match_age_sex_group <- MatchIt::matchit(as.factor(sample_group) ~ as.factor(S) + AGE,
                                        data = df_match_no_na,
                                        method = "nearest",
                                        ratio= 2)

match_age_sex_group %>% summary()

matrix_match <- as.data.frame(match_age_sex_group$match.matrix)

table_of_match <- match_age_sex_group %>% match.data() %>% #as.tibble() %>%
  as.data.frame() %>% as.data.table()


Selected_samples_Eoin = table_of_match[sample_group==0,]$sample
head(Selected_samples_Eoin)

#now merge data so only select relevant samples from Eoins data

healthy_age_sex_matched_ctrls <- data_Eoin[data_Eoin$sample %in% Selected_samples_Eoin,]
healthy_age_sex_matched_ctrls <- healthy_age_sex_matched_ctrls[,-c(2,3)]
head(healthy_age_sex_matched_ctrls)
dim(healthy_age_sex_matched_ctrls)

#create a csv file
write.csv2(healthy_age_sex_matched_ctrls, file = file.path(data_dir,"Healthy_crtls_propensity_matched_by_age_and_sex_v2_3March2021.csv"))
write.csv2(Selected_samples_Eoin, file=file.path(data_dir,"list of healthy controls from Eoin v2 3March2021.csv"))

#now normalise the data from Eoin as the plates do not have normalised AAB levels
#========================
(load(file.path(data_dir,"prop.matching.RData"))) # Aab_P1

##per plate there was a row which required normalisations (everything needed to
##be devided by a certain value for a positive control)
#========================
##normalisation of the data from Eoin by the positive controls
plate1_E <- as.data.frame(Aab_P1[["1"]])
plate1_E_normalised <- rbind(plate1_E, plate1_E[rep(1, 90), ])
plate1_E_pos_crtls <- plate1_E_normalised[c(1,92:180),]
normalised_PE1 <- plate1_E[2:91,]/plate1_E_pos_crtls
## #this works --> this is the right approach!

f=function(x)
  x[-1,]/x[rep(1,nrow(x)-1),]
norm2=f(plate1_E)
sum(normalised_PE1-norm2)
identical(normalised_PE1,norm2)

## this gives same answer, so use it
Aab_norm=lapply(Aab_P1, f)

#now merge all normalised plates into one large df which can be filtered later on
lapply(Aab_norm,colnames) %>% unique() # not all the same - same Aab, but some abbreviated some not
for(i in 2:length(Aab_norm))
  colnames(Aab_norm[[i]])=colnames(Aab_norm[[1]])
total_normalised_healthy_control_df <- do.call("rbind",Aab_norm) %>% as.data.frame()

#now only select cases from this df that have actually been selected
## healt_crtls_match <- list.of.healthy.controls.from.Eoin.v1.5March
total_normalised_healthy_control_df$sample <- sub("X","",rownames(total_normalised_healthy_control_df))
select_health_ctrls=total_normalised_healthy_control_df[match(Selected_samples_Eoin,total_normalised_healthy_control_df$sample),]
table(Selected_samples_Eoin %in% total_normalised_healthy_control_df$sample)
dim(select_health_ctrls)
length(Selected_samples_Eoin)
sum(is.na(select_health_ctrls))

## ##create a csv file
write.csv2(select_health_ctrls,
           file = file.path(data_dir, "normalised FI-BKG levels of healthy controls Eoin v3_3March2021.csv"))

#===========================================
#now I had to log transform these variables together with the cases
#===========================================
#log normalisation of all auto-antibodies
normalised.antibodies.per.patient <- read.csv2(file.path(data_dir,"normalised_antibodies_per_patient.csv"))
normalised.antibodies.for.the.healthy.controls <- select_health_ctrls

#clean the df so that it matches
pts_abs <- normalised.antibodies.per.patient %>% as.data.table()
crtl_abs <- normalised.antibodies.for.the.healthy.controls %>% as.data.table()

#generate a file for pts with all the info
## pts_abs <- merge(pts_abs, outlier_free_df, by="Row.names")
dim(pts_abs)
pts_abs <- merge(pts_abs, propensity_demographics,
                 by.x="Row.names",
                 by.y="pah_subject_id_cohort_validated")
setnames(pts_abs,"Row.names","sample")
dim(pts_abs)
head(pts_abs)

pts_abs$PAH <- "yes"

##total file with all info for patients is done --> now need to generate a file for crtls
colnames(pts_abs)
colnames(crtl_abs)
crtl_abs$PAH <- "no"
crtl_abs=merge(crtl_abs, propensity_Eoin[,.(sample,S,AGE)], by.x="sample", by.y="sample")

## make colnames match
library(stringdist)
D=outer(colnames(pts_abs),colnames(crtl_abs),stringdist)
heatmap(D)
D.pts=apply(D,1,which.min)
D.ctrl=apply(D,2,which.min)
options(width=120)
data.frame(pts=colnames(pts_abs),
           crtl=colnames(crtl_abs)[D.pts],
           d=apply(D,1,min)) # 5:8,9:24 are good matches of aab names, Human.IgG missing from crtl
idx=c(5:8,10:24)
options(width=120)
cbind(colnames(crtl_abs)[D.pts[idx]], colnames(pts_abs)[idx])

colnames(crtl_abs)[D.pts[idx]] <- colnames(pts_abs)[idx]
crtl_abs[,sex:=ifelse(S==0,"female","male")]

##now the datasets are similar and can be merged
i=intersect(colnames(pts_abs),colnames(crtl_abs))
i
total_abs <- rbind(pts_abs[,i,with=FALSE], crtl_abs[,i,with=FALSE])

#export to csv to save this file
write.csv2(total_abs, file = "autoAbs of ctrls and pts merge v1 3March2021.csv")

#make sure all values are positive + 1 before log-transforming
#========================================================

aab_idx=2:20
offset=min(total_abs[,aab_idx,with=FALSE])
for(aab in names(total_abs)[aab_idx]) {
  ## if(any(total_abs[[aab]] < 1)) {
  ## offset=1-min(total_abs[[aab]])
  total_abs[[aab]] <- total_abs[[aab]] + 1 - offset
  ## }
  ##now it is possible to log-transform the complete dataset
  total_abs[[aab]] <- log(total_abs[[aab]])
}

#save as a csv file
write.csv2(total_abs, file="log-transformed autoAbs of ctrls and pts v1 3March2021.csv")

#=====================================
#now the clustering can start
#=====================================
#1) get elbow and silhoutte plots
#plots of shilouette & elbow clusters
#import the log-transformed AABs (pts & crtls)
log_AAB <- copy(total_abs)
##rownames(log_AAB) <- log_AAB$Row.names

#do hc for all these
library(factoextra)
aab_idx=2:20
mat=as.matrix(log_AAB[,aab_idx,with=FALSE])

#first get euclidean distances
m=as.matrix(log_AAB[,aab_idx,with=FALSE])
dist_hc <- dist(m, method = "euclidean")
## print(dist_hc)

#do the hc
#default method is complete linkage (will use this as well, as eclust for HC is also run on default settings)
library(dendextend)
fit <- hclust(dist_hc)
d=as.dendrogram(fit)
plot(d)
cols=ifelse(log_AAB$PAH=="yes",3,2)
colored_bars(matrix(cols,ncol=1),dend=d)
## some clustering apparent

dendogram_h=numeric(13)
for(i in 2:13) dendogram_h[i] <- fit$height[i-1]
plot(13:1, dendogram_h, type = "b", xlab = "number of clusters", ylab = "dendogram height", main = "Scree plot of HC")

#repeat this for total graph
dendogram_2=numeric(1419)
for(i in 2:1419) dendogram_2[i] <- fit$height[i-1]
plot(1419:1, dendogram_2, type = "b", xlab = "number of clusters", ylab = "dendogram height", main = "scree_plot_hc v2")

## install.packages("cluster")
library(cluster)
library(dplyr)
library(factoextra)
#retry it now for a silhouette plot
fviz_nbclust(m, cluster::pam, method="silhouette")+theme_classic() + ggtitle("PAM")

#now do it for an elbow plot
fviz_nbclust(m, cluster::pam, method="wss")+theme_classic() + ggtitle("PAM") + geom_vline(xintercept=3, lty=2,col="#66666666")

#now repeat the elbow plot for hclust
fviz_nbclust(m, hcut, method="wss")+theme_classic() + ggtitle("hclust")
fviz_nbclust(m, hcut, method="silhouette")+theme_classic() + ggtitle("hclust")

library(ggplot2)


#definitive clustering analysis 25 May
## library(tidyverse)
## library(dplyr)
library(Publish)
## library(factoextra)

#import the dataset 
imput_df <- m
rownames(imput_df) <- log_AAB$sample
#do the clustering
definitive_clustering <- eclust(imput_df,FUNcluster="pam", k=3,hc_metric = "euclidean")
#get clustering labels
definitive_labels <- definitive_clustering$clustering
definitive_labels <- as.data.frame(definitive_labels)

definitive_labels$Row.names <- rownames(definitive_labels)
#write csv
write.csv2(definitive_labels, file="definitive clusterings 3k PAM labels 3March2021.csv")

#get new df ready to analyse differences
log_AAB$definitive_labels=as.factor(definitive_labels[ log_AAB$sample, "definitive_labels" ])

#=================
#17 Aug --> quick tests which groups are comparable in PAH status

## CWcomment: this is the correct way to test for a difference in PAH proportions between groups

tt=with(log_AAB,table(PAH,definitive_labels))
print(tt)
chisq.test(tt)

log_AAB[,y:=ifelse(PAH=="yes",1,0)]
mod=glm(y ~ definitive_labels -1, data=as.data.frame(log_AAB), family="binomial")
## group 2 has lower frequency of PAH than others
summary(mod)
exp(coef(mod))
## but groups 1 and 3 also differ slightly
l <- cbind(1, 0, -1)
aod::wald.test(b = coef(mod), Sigma = vcov(mod), L = l)

log_AAB[,.(N=.N,
           mn_age=mean(AGE),med_age=median(AGE),
           prop_female=mean(sex=="female"),prop_PAH=mean(PAH=="yes")),by="definitive_labels"]

## df1 <- df_aabs %>% dplyr::filter(definitive_labels == 1)
## df2 <- df_aabs %>% dplyr::filter(definitive_labels == 2)
## df3 <- df_aabs %>% dplyr::filter(definitive_labels == 3)

## df1$PAH <- as.factor(df1$PAH)
## df2$PAH <- as.factor(df2$PAH)
## df3$PAH <- as.factor(df3$PAH)

## c12 <- rbind(df1, df2)
## chisq.test(c12$PAH, c12$definitive_labels)

## c13 <- rbind(df1, df3)
## chisq.test(c13$PAH, c13$definitive_labels)

## c23 <- rbind(df2, df3)
## chisq.test(c23$PAH, c23$definitive_labels)

## CWcomment: don't do bonferroni here.  The correct way is the chisq.test(tt) above, which tests for any difference between groups.  then we can detect which are different conditional on a difference being there
## #bonferroni correct
## pah_prev <- c(2.798e-07, 0.001162, 0.003452)
## p.adjust(pah_prev, method="bonferroni")

#Everything between the following lines isn't reported in the paper but still relevant to know
#=======================================================
#do kruskal wallis test for AABs (loop adapted from: https://stackoverflow.com/questions/50194609/r-kruskal-wallis-test-in-loop-over-specified-columns-in-data-frame)
results_kruskal <- list()
for(i in colnames(log_AAB)[aab_idx]){
  results_kruskal[[i]] <- kruskal.test(formula(paste(i, "~ definitive_labels")), data = log_AAB)
}
results_kruskal
#this works!
## or collate p values
pvals_kruskal <- sapply(colnames(log_AAB)[aab_idx], function(i) {
  kruskal.test(formula(paste(i, "~ definitive_labels")), data = log_AAB)$p.value
})
names(pvals_kruskal)=colnames(log_AAB)[aab_idx]
pvals_kruskal
#this works!

#generate a univariate table
table_definitive_difference_per_cluster <- univariateTable(definitive_labels ~ Cardiolipin_norm + Centromere.Protein.B..CENP.B._norm + H2a.F2a2..and.H4.F2a1._norm + Histone.type.IIA_norm + Jo.1_norm + La.SS.B.Antigen_norm + Mi.2b_norm + Myeloperoxidase_norm  + Proteinase.3_norm + Pyruvate.Dehydrogenase.from.porcine.hear_norm + RNP.complex_norm + Ro.SS.A..Antigen_norm + SCL.70.Antigen_norm + Scl.34..Fibrillarin.his.tag._norm + Smith.antigen..Sm._norm + Thyroglobulin_norm + Thyroid.peroxidase_norm + Transglutaminase_norm + U1.snRNP.68_norm + sex + PAH + AGE, data = log_AAB, show.totals = TRUE, column.percent = TRUE)

cluster_diff <- summary(table_definitive_difference_per_cluster)

is.numeric(log_AAB$AGE)
anova_age_meas <- aov(AGE ~ definitive_labels, data=log_AAB)
summary(anova_age_meas)

cluster_diff <- as.data.frame(cluster_diff)

#problem: if pval is very small, R treats it as character (due to < sign), this needs to be numeric
#this poses a problem for multiple testing correction
#calculate pvals for these by hand and add these to the list
#after that multiple testing correction can be properly applied

cluster_diff$v <- pvals_kruskal[match(cluster_diff$Variable, names(pvals_kruskal))]
cluster_diff[24,8] <- 0.0112

chisq.test(log_AAB$PAH, log_AAB$definitive_labels)
cluster_diff[23,8] <- 9.603e-10

chisq.test(log_AAB$sex, log_AAB$definitive_labels)
cluster_diff[21,8] <- 5.106e-06
cluster_diff$q <- p.adjust(cluster_diff$v, method = "fdr")

#write as csv
write.csv2(cluster_diff, file="differences per PAH cluster 3k v6 DEF 3March2021.csv")

#=======================================================================================================================

## HERE

#now generate clusters for the clinical data
#first import the clinical file
list.files(data_dir)
clinical_numeric_file <- read.csv2(file.path(data_dir,"NA_33_df_with_outliers_removed.csv"))
#note to self: it's called ith outliers removed because variables with loads of NAs have been removed from this df
head(log_AAB,2)
head(clinical_numeric_file,2)
dim(clinical_numeric_file)
dim(log_AAB)
all(clinical_numeric_file$Row.names %in% log_AAB$sample)

clin_df <- merge(log_AAB[,.(sample,definitive_labels,AGE,sex)],clinical_numeric_file, by.x="sample", by.y ="Row.names")
is.factor(clin_df$definitive_labels)
is.factor(clin_df$sex)

#generate the univariate table
univariate_num_clin <- univariateTable(definitive_labels ~ Cardiolipin_norm + Centromere.Protein.B..CENP.B._norm + H2a.F2a2..and.H4.F2a1._norm + Histone.type.IIA_norm + Jo.1_norm + La.SS.B.Antigen_norm + Mi.2b_norm + Myeloperoxidase_norm  + Proteinase.3_norm + Pyruvate.Dehydrogenase.from.porcine.hear_norm + RNP.complex_norm + Ro.SS.A..Antigen_norm + SCL.70.Antigen_norm + Scl.34..Fibrillarin.his.tag._norm + Smith.antigen..Sm._norm + Thyroglobulin_norm + Thyroid.peroxidase_norm + Transglutaminase_norm + U1.snRNP.68_norm + cbt_haem_platelets_x10e9pl + cbt_haem_wbc_x10e9pl + cbt_haem_haematocrit_ratio + cbt_haem_hb_gpl + cbt_inflammation_crp_mgpl + cbt_liver_alkaline_phos_iupl + cbt_liver_albumin_gpl + cbt_liver_alt_iupl + cbt_liver_bilirubin_umolpl + cbt_liver_total_prot_gpl + cbt_renal_sodium_mmolpl +  cbt_renal_potassium_mmolpl + cbt_renal_urea_mmolpl + cbt_renal_creatinine_umolpl + cbt_thyr_tsh_mupl + cfe_rest_spo2 + cfe_bp_systolic + cfe_bp_diastolic + cfe_heart_rate + ep_1_distance_meters + ep_1_oxygen_saturation_pre + ep_1_oxygen_saturation_post + hb_heart_rate + hb_sp_syst + hb_sp_dias + hb_sp_mean + hb_pawp_m + hb_pap_s + hb_pap_d + hb_pap_m + hb_rap_m + hb_cardiac_output_value_1 + hb_sa_o2 + hb_sv_o2 + hb_pvr_calc + lf_fev1_liters + lf_fev1_pc + lf_fvc_liters + lf_fvc_pc + lf_tlc_liters + lf_kco_mmol + lf_kco_pc + lf_va_liters + fev1_fvc + pcwp_adj + pvr + ci + pp + ca + rc + hr + lf_fev1_fvc_ratio + hb_sv_value + hb_pac_value + age_diagnosis + bs_height + bs_weight + bs_bsa + bs_bmi + AGE, data=clin_df, column.percent = TRUE, show.totals = TRUE)
univariate_num_clin <- summary(univariate_num_clin)
univariate_num_clin
univariate_num_clin2 <- as.data.frame(univariate_num_clin)

#do a kruskal test for the AABs

names(clin_df)
aab_idx=66:87
## or collate p values
pvals_kruskal <- sapply(colnames(clin_df)[aab_idx], function(i) {
  kruskal.test(formula(paste(i, "~ definitive_labels")), data = clin_df)$p.value
})
names(pvals_kruskal)=colnames(clin_df)[aab_idx]
pvals_kruskal


#do an ANOVA for the normal clinical variables

colnames(clin_df)
clin_idx=c(3,6:64) # again, not sure if I have the right variables - 4:61 is one
# fewer than 6:64, but hopefully all are included

pvals_anova_clin = sapply(colnames(clin_df)[clin_idx], function(x) {
  message(x)
  f=as.formula(paste0(x, " ~ definitive_labels"))
  summ=summary(aov(f, data = na.omit(clin_df[,c(x,"definitive_labels"),with=FALSE])))
  summ[[1]]$"Pr(>F)"[1]
})
names(pvals_anova_clin) = colnames(clin_df)[clin_idx]
## histogram shows nice enrichment of smaller p values - otherwise this plot would be flat if all nulls were true
hist(pvals_anova_clin)

#now the adjusted p-values need to be included as this package does not do ANOVA properly
univariate_num_clin2$v <- pvals_kruskal[match(univariate_num_clin2$Variable, names(pvals_kruskal))]

univariate_num_clin_new <- as.data.frame(univariate_num_clin)
univariate_num_clin_new$v <- pvals_anova_clin[match(univariate_num_clin2$Variable, names(pvals_anova_clin))]

#combine these files
u2 <- univariate_num_clin_new[c(20:138),]
u3 <- univariate_num_clin2[c(1:19),]

#combine these files
univariate_numeric <- bind_rows(u2, u3)
univariate_numeric <- univariate_numeric[,-7]
names(univariate_numeric)[[7]] <- "p-value"

#calculate separate pvals for TSH as this - even despite the large n - clearly is not normal
anov1 <- aov(log(cbt_thyr_tsh_mupl)~ definitive_labels, data = na.omit(clin_df[,.(cbt_thyr_tsh_mupl,definitive_labels)]))
summary(anov1)

lm_new_tsh <- lm(log(cbt_thyr_tsh_mupl)~ definitive_labels, data = na.omit(clin_df[,.(cbt_thyr_tsh_mupl,definitive_labels)]))
anova(lm_new_tsh)

#integrate this pval in the list for FDR
univariate_numeric[29,7] <- 0.06623


write.csv2(univariate_numeric, file=file.path(data_dir,"ANOVA and Kruskal-Wallis results clin vars DEF cluster v2_3March2021.csv"))
#===============================================================
#EdB note to self: only run this after complete FDR so that I am sure
#which variables are significant and which ones aren't
#===============================================================

#significant differences were further assessed with an ANCOVA
## CWcomment, I can't run this - what is lm_df?  I can't find diagnosis_verified in any df to date
#EdB Comment: sorry!lm_df just is the autoimmune dataframe that I renamed for simplicity;
#I will change it so that it runs properly
autoimmunity.dataframe_v1 <- fread(file.path(data_dir,"autoimmunity.dataframe_v1.csv"), dec = ",")
lm_df <- autoimmunity.dataframe_v1
label123 <- clin_df[,c(1,2)]
names(lm_df)[[4]] <- "sample"
lm_df$diagnosis_verified.x <- as.factor(lm_df$diagnosis_verified.x)
lm_df$sex <- as.factor(lm_df$sex)
lm_df <- merge(lm_df, label123, by="sample")

lm_df$definitive_labels <- as.factor(lm_df$definitive_labels)

lm_pvr2 <- lm(pvr ~ age_diagnosis + sex + diagnosis_verified.x + bs_bmi + definitive_labels, data=lm_df)

anova_pvr2 <- anova(lm_pvr2)
anova_pvr2
#there's still a significant association, although age sex and PAH subtype also influence PVR

#repeat for PVR-calc
lm_pvrc2 <- lm(hb_pvr_calc ~ age_diagnosis + sex + diagnosis_verified.x + bs_bmi + definitive_labels, data=lm_df)

anova_pvrc2 <- anova(lm_pvrc2)
anova_pvrc2
## #similar result as for PVR

## #repeat for CO
lm_co2 <- lm(hb_cardiac_output_value_1 ~ age_diagnosis + sex + diagnosis_verified.x + bs_bmi + definitive_labels, data=lm_df)

anova_co2 <- anova(lm_co2)
anova_co2

## #CO still is significant

## #now finaly do spo2 post 6mwd
lm_spo2 <- lm(ep_1_oxygen_saturation_post ~ age_diagnosis + sex + diagnosis_verified.x + bs_bmi + definitive_labels, data=lm_df)

anova_spo2 <- anova(lm_spo2)
anova_spo2

#=========================================
#I selected the pvals for the clusters and did Bonferroni correction via p.adjust for these
#performed in the console so not saved in the script
#==========================================


#now analyse the differences for the character data
#import the character dataset
autoimmunity.dataframe_v1 <- fread(file.path(data_dir,"autoimmunity.dataframe_v1.csv"))

df_character <- autoimmunity.dataframe_v1
df_character[df_character==""]<-NA
## colnames(df_character) <- df_character[1,]
## df_character <- df_character[-c(1),]
df_character$Row.names <- df_character$id_cohort.x

#sometimes a variables says not done --> change this to NA
df_character <- df_character %>% naniar::replace_with_na_all(condition =  ~.x == "not done")
#do the same if a variable is unkonw
df_character <- df_character %>% naniar::replace_with_na_all(condition =  ~.x == "unknown")
df_character <- df_character %>% naniar::replace_with_na_all(condition =  ~.x == "UNK")
df_character <- df_character %>% naniar::replace_with_na_all(condition =  ~.x == "not-seen")
df_character <- df_character %>% naniar::replace_with_na_all(condition =  ~.x == "not-recorded")
df_character$ep_1_supplemental_oxygen <- gsub('%', 'Lmin', df_character$ep_1_supplemental_oxygen)
df_character$hb_supplemental_oxygen <- gsub('%', 'Lmin', df_character$hb_supplemental_oxygen)


##the df is clean now --> merge with labels
df_character <- merge(df_character, definitive_labels, by="Row.names")

#the variables of interest all are factors --> change entire df to factors for analysis
df_character <- lapply(df_character, as.factor)

lapply(df_character, is.factor)

df_character <- as.data.frame(df_character)

#now generate the univariate table
univar_DEF_char <- univariateTable(definitive_labels ~ sex + diagnosis + cbt_aab_ana + cbt_aab_anticardiolipin_ab + cbt_aab_antidsdna_ab + cbt_aab_antiscl70_ab + cbt_aab_anticentromere_ab + cbt_aab_antirho_ab + cbt_aab_antiena_ab + cbt_aab_anca + cbt_sero_hep_b_sero + cfe_jugular_venous_pressure_select + cfe_digital_clubbing + reveal_risk + cfe_spider_naevi + cfe_hf_ankle_swelling + cfe_hf_ascites + cfh_symptoms_onset_infection + cfh_syncope + cfh_haemoptysis_hospitalization + cfh_raynauds + ec_percicardial_effusion + ec_right_ventricle + el_rhythm + el_rbbb + el_dominant_r_wave + el_ecg + ep_1_supplemental_oxygen + functional_class + hb_supplemental_oxygen + hv_vasodilator_responder + img_vq_pe + img_emphysema + img_fibrosis + img_thromboembolic_disease + img_portal_venous_flow + lf_sleep_study_diagnosis + vasoresponder + vasoresponder2 + obstructive_pfts + autoimmunity + LtoRshunt + highsbp + highdbp + tg + hdl + img_fibrosis_category + img_emphysema_category + events10 + status_d_ltx + smoking_history + smoking_current + smoking_habit + dh_pde5_drug + dh_era_drug + dh_pa_drug + who1_2 + initial_therapy + comorbid_copd + comorbid_emphysema + comorbid_asthma + comorbid_OSA + comorbid_CAD + comorbid_CVA + comorbid_PAD + comorbid_HTN + comorbid_PE + comorbid_cirrhosis + comorbid_hepatitis + comorbid_AVM + comorbid_bleed + comorbid_epistaxis + comorbid_diabetes + comorbid_hypothyroidism + comorbid_HHT + comorbid_CKD + comorbid_GERD + comorbid_DM1 + comorbid_DM2 + comorbid_asplenia + comorbid_heterotaxy + comorbid_ca + comorbid_sjogren + comorbid_ankylosing.spondylitis + comorbid_pfo + comorbid_overlap_syndrome + comorbid_polymyalgia_rheumatica + FHx + drug_exposure + amphetamines + dasatinib + anorexigen + sex + ethnic_category + BMI_category + svo2 + case, data = df_character, column.percent = TRUE, show.totals = TRUE, compare.groups = TRUE)
univar_DEF_char <- summary(univar_DEF_char)

#only 1 pval has the < sign --> manually adjust this
chisq.test(df_character$definitive_labels, df_character$comorbid_hypothyroidism)
#Pearson's Chi-squared test
#data:  df_character$definitive_labels and df_character$comorbid_hypothyroidism
#X-squared = 24.621, df = 2, p-value = 4.503e-06

univar_DEF_char[231,7] <- 4.503e-06


univar_DEF_char$`p-value` <- as.numeric(univar_DEF_char$`p-value`)
class(univar_DEF_char$`p-value`)


#the publish packages automatically does either Chi-square (>2 groups) 

univariate_character <- as.data.frame(univar_DEF_char)


u1 <- bind_rows(univariate_character, univariate_numeric)
#EdB on 9/3/21: so u1 here is the total df for which only the BMPR comparisons need to be added for the final FDR


#also analysed the effect of clusters on (BMPR2) mutations
#total mutation analysis pts Cohort AABs
#28Sept 2020
## library(tidyverse)
library(Publish)
#BMPR2 mutations 31August

#========================================
#EdB Stopped here on 5/3/21
#=======================================

bmpr2_df <- readRDS(file.path(data_dir,"BMPR2_genotype_patients.rds")) %>% as.data.frame() %>% as.data.table()
table(bmpr2_df$Project)

bmpr2_df <- bmpr2_df[Project=="PAH" & ! (geno %in% c("KDR.het.mis", "EIF2AK4.het.ptv", "EIF2AK4.het.mis", "SMAD4.het.del:SMAD4.het.del:SMAD4.het.del:SMAD4.het.del:SMAD4.het.del", "EIF2AK4.het.ptv:KDR.het.mis"))]

#check which patients are in the cohort
label_wgs <- autoimmunity.dataframe_v1[,c(1:4)]
bmpr2_df <- merge(bmpr2_df, label_wgs, by="WGS.ID")
dim(bmpr2_df)


#this can be analysed further
bmpr_per_clust <- bmpr2_df
bmpr_per_clust <- merge(label_wgs, bmpr_per_clust, by = "id_cohort.x", all=TRUE)
names(bmpr_per_clust)[[1]] <- "Row.names"
bmpr_per_clust <- merge(bmpr_per_clust, definitive_labels, by="Row.names")

#now do some quick stats for mutations per cluster


#Not all patients have WGS data --> exclude these as they can generate a bias, so only compare patients who actaully have had their genes analysed
bmpr_per_clust <- bmpr_per_clust %>% filter (!is.na(WGS.ID.x))

#problem: empty cells aren't considered to be NA
#function from: https://stackoverflow.com/questions/24172111/change-the-blank-cells-to-na
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}
bmpr_per_clust <- bmpr_per_clust %>% mutate_each(funs(empty_as_na)) 

#verified this function --> it worked
#now filter
bmpr_per_clust <- bmpr_per_clust %>% filter (!is.na(WGS.ID.x))
#399 patients are kept in which is as expected

#clean this df
bmpr_per_clust <- bmpr_per_clust[,c("Row.names", "definitive_labels", "Project", "geno")]

#if any NAs are present in this dataset, we know there are no relevant PAH mutations (as these patients were all sequenced and no variants were ''flagged'' as PAH mutations)
bmpr_per_clust[is.na(bmpr_per_clust)] <- "no_mutation"

bmpr_per_clust$definitive_labels <- as.factor(bmpr_per_clust$definitive_labels)
bmpr_per_clust$Project <- as.factor(bmpr_per_clust$Project)
bmpr_per_clust$geno <- as.factor(bmpr_per_clust$geno)

#now compare mutations between groups
chisq.test(bmpr_per_clust$definitive_labels, bmpr_per_clust$Project)
#Pearson's Chi-squared test
#data:  bmpr_per_clust$definitive_labels and bmpr_per_clust$Project
#X-squared = 8.5206, df = 2, p-value = 0.01412

univBMPR1 <- univariateTable(definitive_labels ~ Project, data = bmpr_per_clust, column.percent = TRUE, show.totals = TRUE, compare.groups = TRUE)
univBMPR1 <- summary(univBMPR1)


#now compare if only BMPR2 mutation frequency differs
library(stringr)
bmpr_per_clust$geno <-as.character(bmpr_per_clust$geno)
list_bmpr <- bmpr_per_clust %>% filter(str_detect(geno, "BMPR"))
list_bmpr$bmpr2_mut <- "yes"
list_bmpr2 <- bmpr_per_clust %>% filter(!str_detect(geno, "BMPR"))
list_bmpr2$bmpr2_mut <- "no"

list_bmpr <- bind_rows(list_bmpr, list_bmpr2)
list_bmpr$bmpr2_mut <- as.factor(list_bmpr$bmpr2_mut)
chisq.test(list_bmpr$definitive_labels, list_bmpr$bmpr2_mut)
#Pearson's Chi-squared test
#data:  list_bmpr$definitive_labels and list_bmpr$bmpr2_mut
#X-squared = 4.3304, df = 2, p-value = 0.1147

univBMPR2 <- univariateTable(definitive_labels ~ bmpr2_mut, data = list_bmpr, column.percent = TRUE, show.totals = TRUE, compare.groups = TRUE)
univBMPR2 <- summary(univBMPR2)

bmpr2 <- bind_rows(univBMPR1, univBMPR2)
bmpr2$`p-value` <- as.numeric(bmpr2$`p-value`)

univariate_def <- bind_rows(u1, bmpr2)

#now get an FDR
univariate_def$q_value <- p.adjust(univariate_def$`p-value`, method = "fdr")
write.csv2(univariate_def, file="pvals_all_comparisons_pts_clusters_9March.csv")

## get summary's per patient cluster
bmpr2_1 <- bmpr_per_clust %>% dplyr::filter(definitive_labels == 1)
bmpr2_2 <- bmpr_per_clust %>% dplyr::filter(definitive_labels == 2)
bmpr2_3 <- bmpr_per_clust %>% dplyr::filter(definitive_labels == 3)

summary(bmpr2_1$Project)
summary(bmpr2_1$geno)

summary(bmpr2_2$Project)
summary(bmpr2_2$geno)

summary(bmpr2_3$Project)
summary(bmpr2_3$geno)


#==========================================================
#now get a measure of differences in positivity between cases and controls
#============================================================


#====================================================================
#2 September --> change threshold to 0.75Q + 2IQR 

#calculate positivity for the specific dataset of AAB pts and crtls
#import the dataset

#now measure for positivity is 0.75Q + 2IQR
#if higher than this --> positivit

##calculate this for each AAB in the crtl pop
log_aabs=copy(log_AAB)
calc_threshold=function(x) {
  a1 <- quantile(x, 0.75)
  a2 <- 2*(IQR(x))
  a1+a2
}
colnames(log_aabs)
cols=colnames(log_aabs)[2:20]
thresholds=lapply(cols, function(x) calc_threshold(log_aabs[PAH=="no"][[x]]))
names(thresholds)=cols

#================================================================
#now assess positivity

pos_log <- copy(log_aabs)
for(x in cols)
  pos_log[[x]] = as.factor(ifelse(pos_log[[x]] > thresholds[[x]], "yes", "no"))

newnames=c(
  "cardiolipin_positivity"="Cardiolipin_norm",
  "cenpb_positivity"="Centromere.Protein.B..CENP.B._norm",
  "h2ah4_positivity"="H2a.F2a2..and.H4.F2a1._norm",
  "HistoneIIa_positivity"="Histone.type.IIA_norm",
  "jo1_positivity"="Jo.1_norm",
  "mi2b_positivity"="Mi.2b_norm",
  "lassb_positivity"="La.SS.B.Antigen_norm",
  "mpo_positivity"="Myeloperoxidase_norm",
  "prot3_positivity"="Proteinase.3_norm",
  "pdh_positivity"="Pyruvate.Dehydrogenase.from.porcine.hear_norm",
  "rnpcompl_positivity"="RNP.complex_norm",
  "rossa_positivity"="Ro.SS.A..Antigen_norm",
  "scl70_positivity"="SCL.70.Antigen_norm",
  "scl34_positivity"="Scl.34..Fibrillarin.his.tag._norm",
  "smith_positivity"="Smith.antigen..Sm._norm",
  "tgb_positivity"="Thyroglobulin_norm",
  "tpo_positivity"="Thyroid.peroxidase_norm",
  "transglut_positivity"="Transglutaminase_norm",
  "u1snRNP_positivity"="U1.snRNP.68_norm")
setnames(pos_log, newnames, names(newnames))

names(pos_log)

#publish automatically runs chi-square tests for these groups
univar_pos <- Publish::univariateTable(PAH ~ cardiolipin_positivity + cenpb_positivity + h2ah4_positivity + HistoneIIa_positivity + jo1_positivity + mi2b_positivity + lassb_positivity + mpo_positivity + prot3_positivity + pdh_positivity + rnpcompl_positivity + rossa_positivity + scl70_positivity + scl34_positivity + smith_positivity + tgb_positivity + tpo_positivity + transglut_positivity + u1snRNP_positivity, data=pos_log, show.totals=TRUE, column.percent=TRUE)

univar_pos$p.values %<>% p.adjust(., method="fdr")
univar_pos <- summary(univar_pos)


#save this file
write.csv2(univar_pos, file="positivity_differences_aabs_1sept.csv")

#=========================
#get a measure of positive vs negative pt

#change no in 0 and yes in 1 for this
#can't get this to work quickly..... --
write.csv2(pos_log, file=file.path(data_dir,"num_pos2_3March2021.csv"))


num_pos <- copy(pos_log)

class(num_pos$cardiolipin_positivity)

#now determine positivity
names(num_pos)
num_pos$sum <- rowSums(num_pos[,2:20,with=FALSE]=="yes")
#get row for positivity
num_pos$positivity <- ifelse(num_pos$sum >= 1, "positive", "negative")

#now calculate differences
num_pos$positivity <- as.factor(num_pos$positivity)
chisq.test(num_pos$PAH,num_pos$positivity)

pts_num <- num_pos %>% filter(PAH == "yes")
crtls_num <- num_pos %>% filter(PAH == "no")

summary(pts_num$positivity)
summary(crtls_num$positivity)
#====================================
#now assess the effects of the clusters

## cluster_pos <- merge(num_pos, definitive.clusterings.3k.PAM.labels.25may, by= "Row.names")
cluster_pos=copy(num_pos)
cluster_pos$definitive_labels <- definitive_labels[cluster_pos$sample,"definitive_label"]

#check the class before doing a univariate analysis
lapply(cluster_pos, class)
#need to make some adjustments
#cluster_pos[,2:25] <- lapply(cluster_pos[,2:25], as.factor)

univar_pos <- Publish::univariateTable(definitive_labels ~ cardiolipin_positivity + cenpb_positivity + h2ah4_positivity + HistoneIIa_positivity + jo1_positivity + mi2b_positivity + lassb_positivity + mpo_positivity + prot3_positivity + pdh_positivity + rnpcompl_positivity + rossa_positivity + scl70_positivity + scl34_positivity + smith_positivity + tgb_positivity + tpo_positivity + transglut_positivity + u1snRNP_positivity + positivity, data=cluster_pos, show.totals=TRUE, column.percent=TRUE)
univar_pos$p.values %<>% p.adjust(., method="fdr")
unvar_pos2 <- summary(univar_pos)

write.csv2(unvar_pos2, file="positivity_per_cluster2Sept2020.csv")

#=========================
#now repeat this for clusters with solely patients
pts_pos <- cluster_pos %>% dplyr::filter(PAH == "yes")
pts_pos_nw <- pts_pos
saveRDS(pts_pos_nw, file="AAB_positivity_per_pt_cluster.rds")
AAB_positivity_per_pt_cluster <- pts_pos_nw
definitive.clusterings.3k.PAM.labels.3March2021 <- definitive_labels
names(definitive.clusterings.3k.PAM.labels.3March2021)[[2]] <- "sample"
AAB_positivity_per_pt_cluster <- merge(AAB_positivity_per_pt_cluster, definitive.clusterings.3k.PAM.labels.3March2021, by="sample")
AAB_positivity_per_pt_cluster$definitive_labels <- as.factor(AAB_positivity_per_pt_cluster$definitive_labels)

univar_pos_nw <- Publish::univariateTable(definitive_labels ~ cardiolipin_positivity + cenpb_positivity + h2ah4_positivity + HistoneIIa_positivity + jo1_positivity + mi2b_positivity + lassb_positivity + mpo_positivity + prot3_positivity + pdh_positivity + rnpcompl_positivity + rossa_positivity + scl70_positivity + scl34_positivity + smith_positivity + tgb_positivity + tpo_positivity + transglut_positivity + u1snRNP_positivity + positivity, data= AAB_positivity_per_pt_cluster, show.totals=TRUE, column.percent=TRUE)
univar_pos_nw$p.values %<>% p.adjust(., method="fdr")
univar_pos_nw <- summary(univar_pos_nw)
write.csv2(univar_pos_nw, file="positivity_per_clusterPTS2Sept2020.csv")



#===============================================================
#now plot figures
#first for positivity
AAB_pos_percent_cases_vs_crtls_9March <- read_excel(file.path(data_dir,"AAB_pos_percent_cases_vs_crtls_9March.xlsx"))
clust <- AAB_pos_percent_cases_vs_crtls_9March
clust <- as.data.frame(clust)
rownames(clust) <- clust$Antibody
clust <- clust[,-c(1)]

AAB_pos_percent_clusters_9March <- read_excel(file.path(data_dir,"AAB_pos_percent_clusters_9March.xlsx"))
pah <- AAB_pos_percent_clusters_9March
pah <- as.data.frame(pah)
rownames(pah) <- pah$Antibody
pah <- pah[,-c(1)]

library(pheatmap)

ph1 <- pheatmap(clust, main = "Autoantibody positivity between clusters of solely patients", cluster_rows = F, cluster_cols = F)
ph2 <- pheatmap(pah, main = "Autoantibody positivity between cases and controls", cluster_rows = F, cluster_cols = F)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(ph1, "heatmap_PAH_clusters_AAB_positivity.pdf")
save_pheatmap_pdf(ph2, "heatmap_PAH_vs_crtls_AAB_positivity.pdf")

#===============================================================
#visualise the distributions of the cases and controls
list.files()
data=readRDS(file.path(data_dir,"data_complete_for_paper_aabs_figure4.rds"))

head(data)
class(data)

log_aabs=copy(log_AAB)
calc_threshold=function(x) {
  a1 <- quantile(x, 0.75)
  a2 <- 2*(IQR(x))
  a1+a2
}
colnames(log_aabs)
cols=colnames(log_aabs)[2:20]
thresholds=lapply(cols, function(x) calc_threshold(log_aabs[PAH=="no"][[x]]))
names(thresholds)=cols

colnames(data[,c(25:43)])
names(thresholds) <- colnames(data[,c(25:43)])

log_plot <- log_aabs

log_plot$cardiolipin_threshold <- thresholds$cardiolipin_threshold
log_plot$cenp_b_threshold <- thresholds$cenp_b_threshold
log_plot$h2h4_threshold <- thresholds$h2h4_threshold
log_plot$histIIA_threshold <- thresholds$histIIA_threshold
log_plot$jo1_threshold <- thresholds$jo1_threshold
log_plot$lassb_threshold <- thresholds$lassb_threshold
log_plot$mi2b_threshold <- thresholds$mi2b_threshold
log_plot$mpo_threshold <- thresholds$mpo_threshold
log_plot$proteinase3_threshold <- thresholds$proteinase3_threshold
log_plot$pdh_threshold <- thresholds$pdh_threshold
log_plot$rnp_complex_treshold <- thresholds$rnp_complex_treshold
log_plot$rossa_threshold <- thresholds$rossa_threshold
log_plot$SCL70_threshold <- thresholds$SCL70_threshold
log_plot$scl34_threshold <- thresholds$scl34_threshold
log_plot$smith_antigen_threshold <- thresholds$smith_antigen_threshold
log_plot$tgb_threshold <- thresholds$tgb_threshold
log_plot$TPO_threshold <- thresholds$TPO_threshold
log_plot$transglutaminase_threshold <- thresholds$transglutaminase_threshold
log_plot$u1snrnp_threshold <- thresholds$u1snrnp_threshold

#EdB on 9/3/21 --> further adapt Chris' code

data=log_plot
head(data)
class(data)
library(data.table)
pattern.names=function (..., cols = character(0L))
{
  data.table:::patterns(...,cols=colnames(data))
  p = unlist(list(...), use.names = any(nzchar(names(...))))
  if (!is.character(p))
    stop("Input patterns must be of type character.")
  lapply(p, grep, cols,value=TRUE)
}
nm=pattern.names(list(threshold="_threshold$",value="_norm$"),cols=colnames(data)) # manually checked order matches

setnames(data,"rnp_complex_treshold","rnp_complex_threshold")
names(data)
m=melt(as.data.table(data),
       measure.vars=patterns(threshold="_threshold$",value="_norm$"))
m[,variable:=sub("_norm","",nm$value)[variable]]
head(m)
renamers=c("PDH"= "Pyruvate.Dehydrogenase.from.porcine.hear",
           "CENP-B"= "Centromere.Protein.B..CENP.B.")
for(i in seq_along(renamers))
  m[variable==renamers[i], variable:=names(renamers)[i]]
m[,variable:=factor(variable,levels=rev(sort(unique(m$variable))))] # rev makes the y axis order top to bottom
m[,x:=as.numeric(variable)]

m[,variable:=factor(variable)]
m[,row:=as.numeric(variable) < max(as.numeric(variable))/2]
m[,x:=as.numeric(variable)]

## store median in hc
mmed=m[PAH=="no",.(med0=median(value)),by="variable"]
head(mmed)
nrow(mmed)
m=merge(m,mmed,by="variable")

## center and scale
## use this for untransformed
m[,value_cs:=(value-med0)/(threshold-med0)]
## use this for log transformed
## m[,value_cs:=(log(value)-log(med0))/(log(threshold)-log(med0))]
m[,threshold_cs:=1]

## set 25 as max limit
lim=20
m[,truncated:=value_cs > lim | value_cs < -lim]
m[truncated==TRUE,value_cs:=sign(value_cs) * lim]


mx=min(m[row==FALSE]$x)
m[row==FALSE,x:=x-mx+1]

## summarise by positivity
ms=m[,.(pos=100*mean(value > threshold)),by=c("variable","PAH")]
ms[,PAH:=factor(PAH)]
ms[,x:=as.numeric(variable) + as.numeric(PAH)/5-0.4]
head(ms)

library(ggplot2)
library(ggbeeswarm)
library(cowplot); theme_set(theme_cowplot(font_size=10))


## alternative split panel
library(ggrastr) # add geom_boxplot_jitter
m[,group:=ifelse(PAH=="no", NA, PAH)]
m[,group:=factor(as.character(group),levels=c("1","3","2"))]
levels(m$group)=c("high","intermediate", "low")
m[variable=="RNP.complex",.(RNP=median(value_cs)),by="group"]
m[variable!="RNP.complex",.(notRNP=median(value_cs)),by="group"]

## identify outliers to jitter
find_outliers <- function(y, coef = 1.5) {
  stats <- as.numeric(quantile(y, c(0.25,0.75)))
  iqr <- diff(stats)
  list(lower=stats[1] - coef*iqr, upper=stats[2]+coef*iqr)
}
m[,c("lower","upper"):=find_outliers(value_cs),by=c("variable","PAH")]
m[,c("lower_group","upper_group"):=find_outliers(value_cs),by=c("variable","group")]

def1 <- m
m <- def1
m1 <- m %>% filter(variable == "U1.snRNP.68")
m2 <- m %>% filter(variable == "Thyroid.peroxidase")
m3 <- m %>% filter(variable == "Smith.antigen..Sm.")
m4 <- m %>% filter(variable == "Scl.34..Fibrillarin.his.tag.")
m5 <- m %>% filter(variable == "RNP.complex")
m6 <- m %>% filter(variable == "PDH")
m7 <- m %>% filter(variable == "Mi.2b")
m8 <- m %>% filter(variable == "Jo.1")
m9 <- m %>% filter(variable == "H2a.F2a2..and.H4.F2a1.")
m10 <- m %>% filter(variable == "Cardiolipin")
m11 <- m %>% filter(variable == "Transglutaminase")
m12 <- m %>% filter(variable == "Thyroglobulin")
m13 <- m %>% filter(variable == "SCL.70.Antigen")
m14 <- m %>% filter(variable == "Ro.SS.A..Antigen")
m15 <- m %>% filter(variable == "Proteinase.3")
m16 <- m %>% filter(variable == "Myeloperoxidase")
m17 <- m %>% filter(variable == "La.SS.B.Antigen")
m18 <- m %>% filter(variable == "Histone.type.IIA")
m19 <- m %>% filter(variable == "CENP-B")



b1=ggplot(m1, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m1[!is.na(PAH) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  ggtitle(m1$variable) +
  theme(legend.position="top")

b2=ggplot(m2, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m2[!is.na(PAH) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m2$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b3=ggplot(m3, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m3[!is.na(PAH) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m3$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b4=ggplot(m4, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m4[!is.na(PAH) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m4$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b5=ggplot(m5, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m5[!is.na(PAH) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m5$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b6=ggplot(m6, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m6[!is.na(PAH) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m6$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b7=ggplot(m7, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m7[!is.na(PAH) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m7$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b8=ggplot(m8, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m8[!is.na(PAH) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m8$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b9=ggplot(m9, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m9[!is.na(PAH) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m9$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b10=ggplot(m10, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m10[!is.na(PAH) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m10$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b11=ggplot(m11, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m11[!is.na(PAH) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m11$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b12=ggplot(m12, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m12[!is.na(PAH) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m12$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b13=ggplot(m13, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m13[!is.na(PAH) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m13$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b14=ggplot(m14, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m14[!is.na(PAH) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m14$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b15=ggplot(m15, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m15[!is.na(PAH) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m15$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b16=ggplot(m16, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m16[!is.na(PAH) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m16$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b17=ggplot(m17, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m17[!is.na(PAH) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m17$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b18=ggplot(m18, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m18[!is.na(PAH) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m18$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

b19=ggplot(m19, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=PAH),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=PAH),
             position=position_jitterdodge(jitter.width=0.1),
             data=m19[!is.na(PAH) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m19$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("Healthy controls", "PAH patients"), values=c("darkgreen", "black")) +
  theme(legend.position="top")

m <- def1
m <- m %>% filter(PAH == "yes")

m$group <- m$definitive_labels

m1 <- m %>% filter(variable == "U1.snRNP.68")
m2 <- m %>% filter(variable == "Thyroid.peroxidase")
m3 <- m %>% filter(variable == "Smith.antigen..Sm.")
m4 <- m %>% filter(variable == "Scl.34..Fibrillarin.his.tag.")
m5 <- m %>% filter(variable == "RNP.complex")
m6 <- m %>% filter(variable == "PDH")
m7 <- m %>% filter(variable == "Mi.2b")
m8 <- m %>% filter(variable == "Jo.1")
m9 <- m %>% filter(variable == "H2a.F2a2..and.H4.F2a1.")
m10 <- m %>% filter(variable == "Cardiolipin")
m11 <- m %>% filter(variable == "Transglutaminase")
m12 <- m %>% filter(variable == "Thyroglobulin")
m13 <- m %>% filter(variable == "SCL.70.Antigen")
m14 <- m %>% filter(variable == "Ro.SS.A..Antigen")
m15 <- m %>% filter(variable == "Proteinase.3")
m16 <- m %>% filter(variable == "Myeloperoxidase")
m17 <- m %>% filter(variable == "La.SS.B.Antigen")
m18 <- m %>% filter(variable == "Histone.type.IIA")
m19 <- m %>% filter(variable == "CENP-B")


p1=ggplot(m1, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m1[!is.na(group) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m1$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p2=ggplot(m2, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m2[!is.na(group) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m2$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p3=ggplot(m3, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m3[!is.na(group) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m3$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p4=ggplot(m4, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m4[!is.na(group) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m4$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p5=ggplot(m5, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m5[!is.na(group) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m5$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p6=ggplot(m6, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m6[!is.na(group) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m6$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p7=ggplot(m7, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m7[!is.na(group) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m7$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")


p8=ggplot(m8, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m8[!is.na(group) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m8$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p9=ggplot(m9, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m9[!is.na(group) &
                       (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m9$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p10=ggplot(m10, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m10[!is.na(group) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m10$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p11=ggplot(m11, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m11[!is.na(group) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m11$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red","blue", "orange")) +
  theme(legend.position="top")

p12=ggplot(m12, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m12[!is.na(group) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m12$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p13=ggplot(m13, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m13[!is.na(group) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m13$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p14=ggplot(m14, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m14[!is.na(group) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m14$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p15=ggplot(m15, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m15[!is.na(group) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m15$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p16=ggplot(m16, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m16[!is.na(group) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m16$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p17=ggplot(m17, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m17[!is.na(group) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m17$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p18=ggplot(m18, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m18[!is.na(group) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m18$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", labels = c("cluster 1", "cluster 2", "cluster 3"), values=c("red", "blue", "orange")) +
  theme(legend.position="top")

p19=ggplot(m19, aes(x=variable)) +
  geom_boxplot(aes(y=value_cs,col=group),outlier.shape=NA) +
  geom_point(aes(y=value_cs,col=group),
             position=position_jitterdodge(jitter.width=0.1),
             data=m19[!is.na(group) &
                        (value_cs < lower_group | value_cs > upper_group)],alpha=0.4,size=0.6) +
  geom_hline(yintercept=c(1),col="black",linetype="dashed",size=0.4) +
  ylab("Scaled autoantibody level") +
  ggtitle(m19$variable) +
  theme_bw() +
  scale_colour_manual(name = "groups", c("cluster 1", "cluster 2", "cluster 3"),  values=c("red", "blue", "orange")) +
  theme(legend.position="top")


btw_tick <- 1
library(ggpubr)
c1<- ggarrange(p1 + scale_y_continuous(breaks = seq(0,6, btw_tick))+ylim(-0.5, 6), b1 + scale_y_continuous(breaks = seq(0,6, btw_tick))+ylim(-0.5, 6), align = "v")
c2<- ggarrange(p2, b2)
c3<- ggarrange(p3, b3)
c4<- ggarrange(p4, b4)
c5<- ggarrange(p5, b5)
c6<- ggarrange(p6, b6)
c7<- ggarrange(p7 + scale_y_continuous(breaks = seq(0,4, btw_tick))+ylim(-0.5, 4), b7 + scale_y_continuous(breaks = seq(0,4, btw_tick)) +ylim(-0.5, 4), align = "v")
c8<- ggarrange(p8, b8)
c9<- ggarrange(p9, b9)
c10<- ggarrange(p10, b10)
c11<- ggarrange(p11, b11)
c12<- ggarrange(p12, b12)
c13<- ggarrange(p13, b13)
c14<- ggarrange(p14, b14)
c15<- ggarrange(p15, b15)
c16<- ggarrange(p16, b16)
c17<- ggarrange(p17, b17)
c18<- ggarrange(p18, b18)
c19<- ggarrange(p19, b19)

ggsave(file="positivity_plot1.pdf",plot=c1,dpi=300)
ggsave(file="positivity_plot2.pdf",plot=c2,dpi=300)
ggsave(file="positivity_plot3.pdf",plot=c3,dpi=300)
ggsave(file="positivity_plot4.pdf",plot=c4,dpi=300)
ggsave(file="positivity_plot5.pdf",plot=c5,dpi=300)
ggsave(file="positivity_plot6.pdf",plot=c6,dpi=300)
ggsave(file="positivity_plot7.pdf",plot=c7,dpi=300)
ggsave(file="positivity_plot8.pdf",plot=c8,dpi=300)
ggsave(file="positivity_plot9.pdf",plot=c9,dpi=300)
ggsave(file="positivity_plot10.pdf",plot=c10,dpi=300)
ggsave(file="positivity_plot11.pdf",plot=c11,dpi=300)
ggsave(file="positivity_plot12.pdf",plot=c12,dpi=300)
ggsave(file="positivity_plot13.pdf",plot=c13,dpi=300)
ggsave(file="positivity_plot14.pdf",plot=c14,dpi=300)
ggsave(file="positivity_plot15.pdf",plot=c15,dpi=300)
ggsave(file="positivity_plot16.pdf",plot=c16,dpi=300)
ggsave(file="positivity_plot17.pdf",plot=c17,dpi=300)
ggsave(file="positivity_plot18.pdf",plot=c18,dpi=300)
ggsave(file="positivity_plot19.pdf",plot=c19,dpi=300)


