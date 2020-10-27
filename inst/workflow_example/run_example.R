library(data.table)
library(dplyr)
library(ukbpheno)
library(ggplot2)
library(ggforce)
library(tableone)
library(survminer)
library("MatchIt")


#############################################################################
# specify path to UKB data and to the package
#############################################################################
# folder where the files are downloaded and decrypted
pheno_dir<-".../data/ukb12345/"
fukbtab <- paste(pheno_dir,"ukb12345.tab",sep="")
fhtml<-paste(pheno_dir,"ukb12345.html",sep="")
fhesin<-paste(pheno_dir,"hesin.txt",sep="")
fhesin_diag<-paste(pheno_dir,"hesin_diag.txt",sep="")
fhesin_oper<-paste(pheno_dir,"hesin_oper.txt",sep="")
fgp_clinical<-paste(pheno_dir,"gp_clinical.txt",sep="")
fgp_scripts<-paste(pheno_dir,"gp_scripts.txt",sep="")
fdeath_portal<-paste(pheno_dir,"death.txt",sep="")
fdeath_cause_portal<-paste(pheno_dir,"death_cause.txt",sep="")

#####################################
#  Withdrawal list includes participants who have dropped out of the study.
#  This list is sent to the principal investigator and delegates periodically
####################################

f_withdrawal<-paste(pheno_dir,"w12345_20210809.csv",sep="")


#############################################################################
# load the example definition table and data setting included in the package
#############################################################################
extdata_dir<-paste0(system.file("extdata", package="ukbpheno"),"/")
fdefinitions <- paste0(extdata_dir,"definitions_cardiometabolic_traits.tsv")
fdata_setting <- paste0(extdata_dir,"data.settings.tsv")

# ##########################################
# read setting
############################################
dfData.settings <-fread(fdata_setting)
####################################################
# # read definitions.
###################################################


# ##################################################################################################################################################
# #  *******************this is an optional block -> for wiki?****************************
# ##################################################################################################################################################
# # to process the definition table, we need all available codes
# # currently 2 coding systems exist for gp_clinical
# cf_read3<-dfData.settings[dfData.settings$classification=="CTV3",]$code_map
# cf_read2<-dfData.settings[dfData.settings$classification=="READ2",]$code_map
# # corresponding column names from the .txt file
# get_all_exsiting_codes(fgp_clinical,c("read_2","read_3"),c(paste0(extdata_dir,cf_read2),paste0(extdata_dir,cf_read3)))
#
# #  3 gp_script
# cf_dmd<-dfData.settings[dfData.settings$classification=="DMD",]$code_map
# # bnf.england bnf.scotland are separated because they appear to be coded slightly differently according to the official documentation...hence different post-processing may be needed
# cf_bnf<-unique(dfData.settings[dfData.settings$classification=="BNF",]$code_map)
# cf_read2d<-dfData.settings[dfData.settings$classification=="READ2_drugs",]$code_map
# get_all_exsiting_codes(fgp_scripts,c("read_2","bnf_code","dmd_code"),c(paste0(extdata_dir,cf_read2d),paste0(extdata_dir,cf_bnf),paste0(extdata_dir,cf_dmd)))
#
# # get_all_exsiting_codes(fhesin_diag,c('diag_icd9','diag_icd10'),c(paste0(pheno_dir,'hesin_icd9.code'),paste0(pheno_dir,'hesin_icd10.code')))
# # get_all_exsiting_codes(fdeath_cause_portal,c('cause_icd10'),c(paste0(pheno_dir,'death_icd10.code')))
# # get_all_exsiting_codes(fhesin_oper,c('oper3','oper4'),c(paste0(pheno_dir,'hesin_oper3.code'),paste0(pheno_dir,'hesin_oper4.code')))
# ##################################################################################################################################################



# data.table::fwrite(list(missing_codes),paste(pheno_dir,"ukbphenodata_feb2021.Rdata_notInData.code"))
dfDefinitions_processed_expanded<-read_definition_table(fdefinitions,fdata_setting,extdata_dir)

# ##################################################
# Prepare UKB data:
####################################################

# harmonize data without definition table, default fields including fields from nurse interview are taken
# lst.harmonized.data<-harmonize_ukb_data(f.ukbtab = fukbtab,f.html = fhtml,f.gp_clinical = fgp_clinical,f.gp_scripts = fgp_scripts,f.hesin = fhesin,f.hesin_diag = fhesin_diag,f.hesin_oper =fhesin_oper,f.death_portal = fdeath_portal,f.death_cause_portal = fdeath_cause_portal )

# with definition table, fields required in definition table are also included
# rm(lst.harmonized.data)
lst.harmonized.data<-harmonize_ukb_data(f.ukbtab = fukbtab,f.html = fhtml,dfDefinitions=dfDefinitions_processed_expanded,f.gp_clinical = fgp_clinical,f.gp_scripts = fgp_scripts,f.hesin = fhesin,f.hesin_diag = fhesin_diag,f.hesin_oper =fhesin_oper,f.death_portal = fdeath_portal,f.death_cause_portal = fdeath_cause_portal,f.withdrawal_list = f_withdrawal,allow_missing_fields = TRUE)

# ##################################################
# Get case/control with an example phenotype:
####################################################
# we need 1) definition of the target trait, 2)harmonized data tables, 3) data.setting dataframe and # 4) target set of individuals either specified via df_reference_date or vct.identifiers


# 1) definition of the target trait
trait<-"DmRxT2"
dfDefinitions_processed_expanded[dfDefinitions_processed_expanded$TRAIT==trait,]$DESCRIPTION
# [1] "Atrial fibrillation/Atrial flutter"

# 2)harmonized data table -> lst.harmonized.data
# 3) data.setting dataframe -> dfData.settings
# 4.1) target set of individuals specified via df_reference_date
# df_reference_date is a data.table with identifiers in first column and corresponding reference dates in second column
# Time to disease will be calculated from these specified reference dates
df_reference_dt_v0<-lst.harmonized.data$dfukb[,c("identifier","f.53.0.0")]
# remove the individuals who have withdrawn from the study
df_reference_dt_v0<-df_reference_dt_v0[! identifier  %in% df_withdrawal$V1]

# individuals with DmT2 based on self-reported codes , medication and linkage
lst.DmRxT2.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT==trait), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)

View(lst.DmRxT2.case_control$all_event_dt.Include_in_cases.summary)
View(lst.DmRxT2.case_control$df.casecontrol)

# get a disease time line to see the relative contribution of different data sources over time
DmRxT2_timeline<-plot_disease_timeline_by_source(definition=dfDefinitions_processed_expanded%>%filter(TRAIT==trait),lst.harmonized.data$lst.data,dfData.settings,df_reference_dt_v0$identifier)
DmRxT2_timeline
# get a UpSet plot to see the relative contribution of different data sources at baseline
upset_plot<-make_upsetplot(definition=dfDefinitions_processed_expanded%>%filter(TRAIT==trait),lst.harmonized.data$lst.data,dfData.settings,df.reference.dates = df_reference_dt_v0)
upset_plot
# extract all hospital admission records
all_DmRxT2_evnt<-lst.DmRxT2.case_control$all_event_dt.Include_in_cases
DmRxT2_hesin_rec<-all_DmRxT2_evnt[grepl ("hesin" ,all_DmRxT2_evnt$.id)]
#  get some descriptive statistics on the records on a code level
hesin_stats<-get_stats_for_events(DmRxT2_hesin_rec)
hesin_stats$stats.codes.summary.p

#  get some summary statistics on the records on individual level
DmRxT2_rec_cnt<-DmRxT2_hesin_rec[,.(count=.N),by=c("identifier")]
max(DmRxT2_rec_cnt$count)
median(DmRxT2_rec_cnt$count)
mean(DmRxT2_rec_cnt$count)
quantile(DmRxT2_rec_cnt$count)

# visualization count with barplot
ggplot2::ggplot(DmRxT2_rec_cnt, ggplot2::aes(x=count)) +
  ggplot2::geom_bar(fill="#0073C2FF" ) +  ggplot2::xlab("Number fo secondary care record per person") +
  ggplot2::ylab("Frequency") + #theme with white background
  ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size=22),panel.grid.minor =ggplot2::element_blank(),panel.grid.major =ggplot2::element_blank()) +  ggforce::facet_zoom(xlim = c(0, 50))

# plot individual time line
plot_individual_timeline(df.data.settings = dfData.settings,lst.data=lst.harmonized.data$lst.data,ind_identifier = 1111111)



################################
# refine the diagnoses
###############################
# identify individuals with specific DmT2 codes
lst.DmT2.case_control<-get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="DmT2"), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)
# identify individuals with specific DmT1 codes
lst.DmT1.case_control<-get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="DmT1"), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)
# identify individuals with DmG
lst.DmG.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="DmG"), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)
# identify individuals with metformin use , 1140884600 contributes 15,085
lst.RxMet.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="RxMet"), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)
# identify indviduals with diagnosis codes related diabetes excluding medication use
lst.Dm.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Dm"), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)
# identify young onset self reported diabetes (european origin)
lst.SrDmYEw.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="SrDmYEw"), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)
# identify young onset self reported diabetes (carribean african origin)
lst.SrDmYSaCa.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="SrDmYSaCa"), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)
#identify use of insulin/oral diabetic medication other than metformin
lst.RxDmNoMet.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="RxDmNoMet"), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)
#identify individuals that are on metformin but no self-reported diabetes nor use of insulin/oral diabetic medication other than metformin
RxMet_DmUnlikely<-setdiff(lst.RxMet.case_control$df.casecontrol[Hx==2]$identifier,union(lst.Dm.case_control$df.casecontrol[Hx==2]$identifier,lst.RxDmNoMet.case_control$df.casecontrol[Hx==2]$identifier))

# identify individuals with self-report insulin <12 months post-diagnosis
lst.RxDmInsFirstYear.case_control<-get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="RxDmInsFirstYear"), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)
#identify individuals with insulin
lst.RxDmIns.case_control<-get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="RxDmIns"), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)

# cross examine various diagnoses, for example to get young onset diabetes who do not have evidence of  other types of diabetes (on insulin, on insulin within 1 year of diagnosis, specific type 1 diabetes codes, specific gestational diabetes codes)
# individuals of young onset diabetes
ind_young_onset<- union(lst.SrDmYSaCa.case_control$df.casecontrol[Any==2]$identifier,lst.SrDmYEw.case_control$df.casecontrol[Any==2]$identifier)

# individuals with evidence of other types of diabetes reported
ind_RxInsFirstYear_DmT1_DmG<- union(union(lst.RxDmInsFirstYear.case_control$df.casecontrol[Any==2]$identifier,lst.DmT1.case_control$df.casecontrol[Hx==2]$identifier),lst.DmG.case_control$df.casecontrol[Hx==2]$identifier)

# young onset but no DM type 1/ gestational diabetes specific codes nor self report of insulin within first year of diagnosis
inds_young_onset_potential_DmT2 <-setdiff(ind_young_onset,ind_RxInsFirstYear_DmT1_DmG)

# further filter out individuals who are on the metformin but likely NOT diabetes
inds_young_onset_potential_DmT2 <-setdiff(inds_young_onset_potential_DmT2,RxMet_DmUnlikely)

# check if these individuals have specific DmT2 codes
# and inspect the data on those individuals who don't  i.e. the uncertain cases
lst.DmRxT2.case_control$all_event_dt.Include_in_cases[identifier %in% setdiff(inds_young_onset_potential_DmT2,lst.DmT2.case_control$df.casecontrol[Hx==2]$identifier)]

###########################################################################
# Part 2 generate phenotype in batch and make a clinical characteristic table
##########################################################################
# read the definitions table
fdefinitions <- paste(extdata_dir,"definitions_cardiometabolic_traits.tsv",sep="")
dfDefinitions_processed_expanded<-read_defnition_table(fdefinitions,fdata_setting,extdata_dir)

# extract clinical variables from the main dataset using read_ukb_tabdata()
# we need the metadata (.html) file for read_ukb_tabdata()
dfhtml <- read_ukb_metadata(fhtml)
# rename the identifier column in the metadata
dfhtml[which(dfhtml$field.tab=="f.eid"),]$field.tab<-"identifier"

# age at assessment centre visit, sex, BMI , HbA1c, glucose,insulin within 1 year of diagnosis,UK Biobank assessment center location,visit date
baseline_fields<-c(21003,31,21001,30750,30740,2986,54,53)
# extract these variables from main dataset
dfukb_baseline <- read_ukb_tabdata(fukbtab,dfhtml,fields_to_keep = baseline_fields)
gc()
# the target disease traits we will generate in batch
diseases<-c("Af","Cad","DmT2","Hcm","Hf","HtRx","HyperLipRx")

# make an output folder to store the result
out_folder<-paste0(pheno_dir,"output/")
if(!dir.exists(file.path(out_folder))){
 dir.create(file.path(out_folder))
}

df_withdrawal<-fread(f_withdrawal)
dfukb_baseline_pheno<-dfukb_baseline[! identifier %in% df_withdrawal$V1]
# loop through the traits, including family history of related diseases and the diabetes medication use
for (disease in c(diseases,"HxDm","HxHrt","HxHt","RxDmOr","RxDmIns")){
  print(disease)
  lst.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT==disease), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)
  # add the trait to the  col names
  colnames(lst.case_control$df.casecontrol) <- paste(disease,"0",colnames(lst.case_control$df.casecontrol),  sep = "_")
  # except for participant identifier
  names(lst.case_control$df.casecontrol)[names(lst.case_control$df.casecontrol) == paste(disease,"0","identifier",  sep = "_")]<-"identifier"
  # merge these columns with dfukb_baseline_pheno
  dfukb_baseline_pheno<-merge(dfukb_baseline_pheno,lst.case_control$df.casecontrol,by="identifier",all.x = TRUE,all.y = FALSE)
}


##################################
# Baseline Characteristics table
#################################

# keep only the variables needed for the table
dfukb_baseline_pheno_fortable1<-dfukb_baseline_pheno[,c('identifier',"DmT2_0_Hx","f.21003.0.0","f.21001.0.0","f.30740.0.0","f.30750.0.0","DmT2_0_first_diagnosis_days","f.31.0.0","HxDm_0_Any","HxHrt_0_Any","HxHt_0_Any","HtRx_0_Hx","HyperLipRx_0_Hx","Af_0_Hx","Hcm_0_Hx","Hf_0_Hx","RxDmOr_0_Hx","RxDmIns_0_Hx","f.2986.0.0"),with=FALSE]
# -ve first diagnosis day indicates history while +ve indicates follow-up cases
dfukb_baseline_pheno_fortable1$DmT2_0_first_diagnosis_years<-(-1*dfukb_baseline_pheno_fortable1$DmT2_0_first_diagnosis_days)/365.25

# rename for readability
colnames(dfukb_baseline_pheno_fortable1)<-c("identifier","Type 2 diabetes","Age","BMI","Glucose","HbA1c","Days since type 2 diabetes diagnosis","Sex","Family history of diabetes","Family history of heart disease","Family history of hypertension","Hypertension","Hyperlipidemia","Atrial fibrillation","Hypertrophic cardiomyopathy","Heart failure","Oral diabetes medication","Insulin","Insulin within 1 year of diagnosis","Years since type 2 diabetes diagnosis")
# below the parameters for CreateTableOne
# the full variable list
vars<-c("Age","BMI","Glucose","HbA1c","Years since type 2 diabetes diagnosis","Sex","Family history of diabetes","Family history of heart disease","Family history of hypertension","Hypertension","Hyperlipidemia","Atrial fibrillation","Hypertrophic cardiomyopathy","Heart failure","Oral diabetes medication","Insulin","Insulin within 1 year of diagnosis")
# the categorical variables on the clinical characteristics table
factorVars<-setdiff(vars,c("Age","BMI","Glucose","HbA1c","Years since type 2 diabetes diagnosis"))

## Create the clinical characteristic table stratified by type 2 diabetes
tableOne <- CreateTableOne(vars = vars, strata = "Type 2 diabetes", data = dfukb_baseline_pheno_fortable1, factorVars = factorVars)
hist(dfukb_baseline_pheno_fortable1$`Years since type 2 diabetes diagnosis`)
tableOne
tab1Mat <- print(tableOne, quote = FALSE, noSpaces = TRUE, printToggle = FALSE,nonnormal =c("Glucose","HbA1c","Years since type 2 diabetes diagnosis") )
## Save the table to a CSV file
write.csv(tab1Mat, file =paste0(out_folder,"BaselineTable.csv"))



##################################
# Survival analysis
##################################
# get death dates from data
deathdt<-unique(lst.harmonized.data$lst.data$tte.death.icd10.primary[,.(identifier,eventdate)])
# rename the column
colnames(deathdt)<-c("identifier","deathdt")
# merge
dfukb_baseline_pheno<-merge(dfukb_baseline_pheno,deathdt,by="identifier",all.x=TRUE,all.y = FALSE)
# HESIN censoring date are different by regions
# retrieve this info using UK Biobank assessment center location attended by the participants
england<-c("10003","11001","11002","11007","11008","11009","11010","11011","11012","11013","11014","11016","11017","11018","11019","11020","11021", "11006")
scotland<-c("11004","11005")
wales<-c("11003","11022","11023")
# corresponding censoring dates
dfukb_baseline_pheno[dfukb_baseline_pheno$f.54.0.0 %in% england,"censordateHES"]<-as.Date("2021-03-31")
dfukb_baseline_pheno[dfukb_baseline_pheno$f.54.0.0 %in% scotland,"censordateHES"]<-as.Date("2021-03-31")
dfukb_baseline_pheno[dfukb_baseline_pheno$f.54.0.0 %in% wales,"censordateHES"]<-as.Date("2018-02-28")

# time-to-event/observed time is determined at earliest of date of event, date of death and censoring date of HESIN data (last follow up)
# this is already calculated for those who have events
range(dfukb_baseline_pheno[Hf_0_Fu==2,Hf_0_Fu_days])
# ppl wihtout event but died before HESIN censoring date (end of follow up)
dfukb_baseline_pheno[Hf_0_Fu==1 & !is.na(deathdt) & deathdt-censordateHES<=0,Hf_0_Fu_days:=deathdt-as.Date(f.53.0.0)]
# people censored at last fu
# ppl wihtout event but died after censoring date (HESIN),
dfukb_baseline_pheno[Hf_0_Fu==1 &!is.na(deathdt)& deathdt-censordateHES>0 ,Hf_0_Fu_days:=censordateHES-as.Date(f.53.0.0)]
# ppl wihtout event and alive by censoring date
dfukb_baseline_pheno[Hf_0_Fu==1 &is.na(deathdt) ,Hf_0_Fu_days:=censordateHES-as.Date(f.53.0.0)]

#Estimate risk of new onset heart failure by presence/absence of type 2 diabetes at baseline
fit<-survival::survfit(survival::Surv(Hf_0_Fu_days/365.25,Hf_0_Fu ) ~ DmT2_0_Hx, data = dfukb_baseline_pheno[DmT2_0_Hx>0])
# summary(fit)
# Make Kaplan-Meier plot
ggsurvplot(fit, data =  dfukb_baseline_pheno[DmT2_0_Hx>0], size = 0.8,
           break.time.by=2,
           xlab = "Follow up (years)",
           censor.size=2,
           palette = c("#072A6C", "#FF8400"),
           conf.int = TRUE,          # Add confidence interval
           pval = TRUE,              # Add p-value
           risk.table = TRUE,        # Add risk table
           risk.table.col = "strata",# Risk table color by groups
           legend.labs = c("No type 2 diabetes at baseline","Type 2 diabetes at baseline"),
           risk.table.height = 0.2)

########################################
# 1:2 case control matching with MatchIt
########################################
# Remove individuals with either missing or excluded phenotype for target phenotype (type 2 diabetes at baseline)
df_to_matchit<-dfukb_baseline_pheno[!is.na(DmT2_0_Hx) & DmT2_0_Hx>0]

# Pick three covariates age at assessment center visit, sex and BMI for matching
df_to_matchit<-na.omit(df_to_matchit[,.(identifier,DmT2_0_Hx,f.21003.0.0,f.31.0.0,f.21001.0.0)])
# Format the data for the matchit function
# Control/case: 1/2 to 0/1
df_to_matchit$DmT2_0_Hx<-df_to_matchit$DmT2_0_Hx-1
# Name the rows
rownames(df_to_matchit)<-df_to_matchit$identifier
colnames(df_to_matchit)<-c("identifier","Type 2 diabetes","Age","Sex","BMI")
# Run matchit
m.dm2<-matchit(`Type 2 diabetes`~Age + Sex+BMI,data=df_to_matchit,ratio=2)
#Check result
summary(m.dm2)
# Each row in the match.matrix shows identifier of one case with 2 matched controls
m.dm2$match.matrix
