library(data.table)
library(dplyr)
library(ukbpheno)
library(ggplot2)
library(ggforce)
library(tableone)

#############################################################################
# specify path to UKB data and to the package
#############################################################################

repo_dir<-".../repos/ukbpheno/"
pheno_dir<-".../data/ukb12345/"


fukbtab <- paste(pheno_dir,"ukb47316.tab",sep="")
fhtml<-paste(pheno_dir,"ukb47316.html",sep="")
fhesin<-paste(pheno_dir,"hesin.txt",sep="")
fhesin_diag<-paste(pheno_dir,"hesin_diag.txt",sep="")
fhesin_oper<-paste(pheno_dir,"hesin_oper.txt",sep="")
fgp_clinical<-paste(pheno_dir,"gp_clinical.txt",sep="")
fgp_scripts<-paste(pheno_dir,"gp_scripts.txt",sep="")
fdeath_portal<-paste(pheno_dir,"death.txt",sep="")
fdeath_cause_portal<-paste(pheno_dir,"death_cause.txt",sep="")


#############################################################################
# load the example definition table and data setting included in the package
#############################################################################
data_dir<-paste(repo_dir,"inst/extdata/",sep="")
fdefinitions <- paste(data_dir,"definitions_cardiometabolic_traits.tsv",sep="")
fdata_setting <- paste(data_dir,"data.settings.tsv",sep="")

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
# get_all_exsiting_codes(fgp_clinical,c("read_2","read_3"),c(paste0(data_dir,cf_read2),paste0(data_dir,cf_read3)))
#
# #  3 gp_script
# cf_dmd<-dfData.settings[dfData.settings$classification=="DMD",]$code_map
# # bnf.england bnf.scotland are separated because they appear to be coded slightly differently according to the official documentation...hence different post-processing may be needed
# cf_bnf<-unique(dfData.settings[dfData.settings$classification=="BNF",]$code_map)
# cf_read2d<-dfData.settings[dfData.settings$classification=="READ2_drugs",]$code_map
# get_all_exsiting_codes(fgp_scripts,c("read_2","bnf_code","dmd_code"),c(paste0(data_dir,cf_read2d),paste0(data_dir,cf_bnf),paste0(data_dir,cf_dmd)))
#
# # get_all_exsiting_codes(fhesin_diag,c('diag_icd9','diag_icd10'),c(paste0(pheno_dir,'hesin_icd9.code'),paste0(pheno_dir,'hesin_icd10.code')))
# # get_all_exsiting_codes(fdeath_cause_portal,c('cause_icd10'),c(paste0(pheno_dir,'death_icd10.code')))
# # get_all_exsiting_codes(fhesin_oper,c('oper3','oper4'),c(paste0(pheno_dir,'hesin_oper3.code'),paste0(pheno_dir,'hesin_oper4.code')))
# ##################################################################################################################################################



# data.table::fwrite(list(missing_codes),paste(pheno_dir,"ukbphenodata_feb2021.Rdata_notInData.code"))
dfDefinitions_processed_expanded<-read_defnition_table(fdefinitions,fdata_setting,data_dir)

# ##################################################
# Prepare UKB data:
####################################################

# harmonize data without definition table, default fields including fields from nurse interview are taken
# lst.harmonized.data<-harmonize_ukb_data(f.ukbtab = fukbtab,f.html = fhtml,f.gp_clinical = fgp_clinical,f.gp_scripts = fgp_scripts,f.hesin = fhesin,f.hesin_diag = fhesin_diag,f.hesin_oper =fhesin_oper,f.death_portal = fdeath_portal,f.death_cause_portal = fdeath_cause_portal )

# with definition table, fields required in definition table are also included
# rm(lst.harmonized.data)
lst.harmonized.data<-harmonize_ukb_data(f.ukbtab = fukbtab,f.html = fhtml,dfDefinitions=dfDefinitions_processed_expanded,f.gp_clinical = fgp_clinical,f.gp_scripts = fgp_scripts,f.hesin = fhesin,f.hesin_diag = fhesin_diag,f.hesin_oper =fhesin_oper,f.death_portal = fdeath_portal,f.death_cause_portal = fdeath_cause_portal,allow_missing_fields = TRUE)

#####################################
# read withdrawal list and remove them from the data
####################################

f_particip_withdraw<-paste(pheno_dir,"w12345_20210809.csv",sep="")
df_withdrawal<-fread(f_particip_withdraw)


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
DmRxT2_timeline<-plot_disease_timeline_by_source(definition=dfDefinitions_processed_expanded%>%filter(TRAIT==trait),lst.harmonized.data$lst.data,dfData.settings,lst.harmonized.data$vct.identifiers)
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
inds_young_onset_probable_DmT2 <-setdiff(ind_young_onset,ind_RxInsFirstYear_DmT1_DmG)

# further filter out individuals who are on the metformin but likely NOT diabetes
inds_young_onset_probable_DmT2 <-setdiff(inds_young_onset_probable_DmT2,RxMet_DmUnlikely)

# check if these individuals have specific DmT2 codes
# and inspect the data on those individuals who don't  i.e. the uncertain cases
lst.DmRxT2.case_control$all_event_dt.Include_in_cases[identifier %in% setdiff(inds_young_onset_probable_DmT2,lst.DmT2.case_control$df.casecontrol[Hx==2]$identifier)]


###########################################################################
# Part 2 generate phenotype in batch and make a clinical characteristic table
##########################################################################

# extract clinical variables from the main dataset using read_ukb_tabdata()
# we need the metadata (.html) file for read_ukb_tabdata()
dfhtml <- read_ukb_metadata(fhtml)
# rename the identifier column in the metadata
dfhtml[which(dfhtml$field.tab=="f.eid"),]$field.tab<-"identifier"

# age at assessment centre visit, sex, BMI , HbA1c, glucose,insulin within 1 year of diagnosis
baseline_fields<-c(21003,31,21001,30750,30740,2986)
# extract these variables from main dataset
dfukb_baseline <- read_ukb_tabdata(fukbtab,dfhtml,fields_to_keep = baseline_fields)
gc()

# read the definitions table
fdefinitions <- paste(data_dir,"definitions_cardiometabolic_traits.tsv",sep="")
# the target disease traits we will generate in batch
diseases<-c("Af","Cad","DmRxT2","Hcm","Hf","HtRx","HyperLipRx")

# make an output folder to store the result
out_folder<-paste0(pheno_dir,"output/")
if(!dir.exists(file.path(out_folder))){
 dir.create(file.path(out_folder))
}


dfukb_baseline_pheno<-dfukb_baseline
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

dfukb_baseline_pheno$DmRxT2_0_first_diagnosis_years<-(-1*dfukb_baseline_pheno$DmRxT2_0_first_diagnosis_days)/365.25

# keep only the variables needed for the table
dfukb_baseline_pheno<-dfukb_baseline_pheno[,c('identifier',"DmRxT2_0_Hx","f.21003.0.0","f.21001.0.0","f.30740.0.0","f.30750.0.0","DmRxT2_0_first_diagnosis_years","f.31.0.0","HxDm_0_Any","HxHrt_0_Any","HxHt_0_Any","HtRx_0_Hx","HyperLipRx_0_Hx","Af_0_Hx","Hcm_0_Hx","Hf_0_Hx","RxDmOr_0_Hx","RxDmIns_0_Hx","f.2986.0.0"),with=FALSE]
# rename for readability
colnames(dfukb_baseline_pheno)<-c("identifier","Type 2 diabetes","Age","BMI","Glucose","HbA1c","Years of diabetes","Sex","Family history of diabetes","Family history of heart disease","Family history of hypertension","Hypertension","Hyperlipidemia","Atrial fibrillation","Hypertropic cardiomyopathy","Heart failure","Oral diabetes medication","Insulin","Insulin within 1 year of diagnosis")
# below the parameters for CreateTableOne
# the full variable list
vars<-c("Age","BMI","Glucose","HbA1c","Years of diabetes","Sex","Family history of diabetes","Family history of heart disease","Family history of hypertension","Hypertension","Hyperlipidemia","Atrial fibrillation","Hypertropic cardiomyopathy","Heart failure","Oral diabetes medication","Insulin","Insulin within 1 year of diagnosis")
# the categorical variables on the clinical characteristics table
factorVars<-setdiff(vars,c("Age","BMI","Glucose","HbA1c","Years of diabetes"))

## Create the clinical characteristic table stratified by type 2 diabetes
tableOne <- CreateTableOne(vars = vars, strata = "Type 2 diabetes", data = dfukb_baseline_pheno, factorVars = factorVars)
hist(dfukb_baseline_pheno$`Years of diabetes`)
tableOne
tab1Mat <- print(tableOne, quote = FALSE, noSpaces = TRUE, printToggle = FALSE,nonnormal =c("Glucose","HbA1c","Years of diabetes") )
## Save the table to a CSV file
write.csv(tab1Mat, file =paste0(out_folder,"BaselineTable.csv"))



