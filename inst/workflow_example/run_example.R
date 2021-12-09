library(data.table)
library(dplyr)
library(ukbpheno)

#############################################################################
# specify path to UKB data and to the package
#############################################################################

repo_dir=".../repos/ukbpheno/"
pheno_dir=".../data/ukb12345/"


fukbtab = paste(pheno_dir,"ukb47316.tab",sep="")
fhtml = paste(pheno_dir,"ukb47316.html",sep="")
fhesin=paste(pheno_dir,"hesin.txt",sep="")
fhesin_diag=paste(pheno_dir,"hesin_diag.txt",sep="")
fhesin_oper=paste(pheno_dir,"hesin_oper.txt",sep="")
fgp_clinical =paste(pheno_dir,"gp_clinical.txt",sep="")
fgp_scripts =paste(pheno_dir,"gp_scripts.txt",sep="")
fdeath_portal=paste(pheno_dir,"death.txt",sep="")
fdeath_cause_portal=paste(pheno_dir,"death_cause.txt",sep="")


#############################################################################
# load the example definition table and data setting included in the package
#############################################################################
data_dir<-paste(repo_dir,"inst/extdata/",sep="")
fdefinitions = paste(data_dir,"definitions_CadDcmHcmAfHf_CTV_READ2.tsv",sep="")
fdefinitions = paste(data_dir,"definitions_v3.1.tsv",sep="")

fdata_setting = paste(data_dir,"data.settings.tsv",sep="")

# ##########################################
# read setting
############################################
dfData.settings <-fread(fdata_setting)
####################################################
# # read definitions.
###################################################


##################################################################################################################################################
#  *******************this is an optional block -> for wiki?****************************
##################################################################################################################################################
# to process the definition table, we need all available codes
# currently 2 coding systems exist for gp_clinical
cf_read3<-dfData.settings[dfData.settings$classification=="CTV3",]$code_map
cf_read2<-dfData.settings[dfData.settings$classification=="READ2",]$code_map
# corresponding column names from the .txt file
get_all_exsiting_codes(fgp_clinical,c("read_2","read_3"),c(paste0(data_dir,cf_read2),paste0(data_dir,cf_read3)))

#  3 gp_script
cf_dmd<-dfData.settings[dfData.settings$classification=="DMD",]$code_map
# bnf.england bnf.scotland are separated because they appear to be coded slightly differently according to the official documentation...hence different post-processing may be needed
cf_bnf<-unique(dfData.settings[dfData.settings$classification=="BNF",]$code_map)
cf_read2d<-dfData.settings[dfData.settings$classification=="READ2_drugs",]$code_map
get_all_exsiting_codes(fgp_scripts,c("read_2","bnf_code","dmd_code"),c(paste0(data_dir,cf_read2d),paste0(data_dir,cf_bnf),paste0(data_dir,cf_dmd)))

# get_all_exsiting_codes(fhesin_diag,c('diag_icd9','diag_icd10'),c(paste0(pheno_dir,'hesin_icd9.code'),paste0(pheno_dir,'hesin_icd10.code')))
# get_all_exsiting_codes(fdeath_cause_portal,c('cause_icd10'),c(paste0(pheno_dir,'death_icd10.code')))
# get_all_exsiting_codes(fhesin_oper,c('oper3','oper4'),c(paste0(pheno_dir,'hesin_oper3.code'),paste0(pheno_dir,'hesin_oper4.code')))
##################################################################################################################################################



# data.table::fwrite(list(missing_codes),paste(pheno_dir,"ukbphenodata_feb2021.Rdata_notInData.code"))
dfDefinitions_processed_expanded<-read_defnition_table(fdefinitions,fdata_setting,data_dir)


# ##################################################
# Prepare UKB data:
####################################################


# harmonize data without definition table, default fields including fields from nurse interview are taken
lst.harmonized.data<-harmonize_ukb_data(f.ukbtab = fukbtab,f.html = fhtml,f.gp_clinical = fgp_clinical,f.gp_scripts = fgp_scripts,f.hesin = fhesin,f.hesin_diag = fhesin_diag,f.hesin_oper =fhesin_oper,f.death_portal = fdeath_portal,f.death_cause_portal = fdeath_cause_portal )

# with definition table, fields required in definition table are also included
lst.harmonized.data<-harmonize_ukb_data(f.ukbtab = fukbtab,f.html = fhtml,dfDefinitions=dfDefinitions_processed_expanded,f.gp_clinical = fgp_clinical,f.gp_scripts = fgp_scripts,f.hesin = fhesin,f.hesin_diag = fhesin_diag,f.hesin_oper =fhesin_oper,f.death_portal = fdeath_portal,f.death_cause_portal = fdeath_cause_portal )


# ##################################################
# Get case/control with an example phenotype:
####################################################
# we need 1) definition of the target trait, 2)harmonized data tables, 3) data.setting dataframe and # 4) target set of individuals either specified via df_reference_date or vct.identifiers
# 1) definition of the target trait
trait<-"DmRxT2"
# dfDefinitions_processed_expanded[dfDefinitions_processed_expanded$TRAIT==trait,]$DESCRIPTION
# [1] "Atrial fibrillation/Atrial flutter"

# 2)harmonized data table -> lst.harmonized.data
# 3) data.setting dataframe -> dfData.settings
# 4.1) target set of individuals specified via df_reference_date
# df_reference_date is a data.table with identifiers in first column and corresponding reference dates in second column
# Time to disease will be calculated from these specified reference dates
df_reference_dt_v0<-lst.harmonized.data$dfukb[,c("identifier","f.53.0.0")]
df_reference_dt_v0$f.53.0.0<-as.Date(df_reference_dt_v0$f.53.0.0)
lst.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT==trait), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)


df_reference_dt_v2<-lst.harmonized.data$dfukb[,c("identifier","f.53.2.0")]
lst.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT==trait), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v2)
# TODO CHECK numbers of cases not a sum from the messages , can be confusing


# 4.2) target set of individuals specified via vct.identifiers i.e. no specified reference dates
# The dates of the first event will be taken as reference date for the calculation of time to disease
lst.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT==trait), lst.harmonized.data$lst.data,dfData.settings, vct.identifiers=lst.harmonized.data$vct.identifiers)

# it is a good practice to check if all the identifiers are in the main dataset, as non-existing identifier will be classified as a control (no events)
my.curated.identifiers<-df_reference_dt_v2[1:20]$identifier
my.curated.identifiers[1]<-10000011
all(my.curated.identifiers %in% lst.harmonized.data$vct.identifiers) #FALSE

all_event_dt <- get_all_events(definition=dfDefinitions_processed_expanded%>%filter(TRAIT==trait)%>%filter(Definitions=="Include_in_cases"),lst.data=lst.harmonized.data$lst.data,df.data.settings=dfData.settings)   #Af


plot_individual_timeline(df.data.settings = dfData.settings,lst.data=lst.harmonized.data$lst.data,ind_identifier = 1111111)

DmRxT2_timeline<-plot_disease_timeline_by_source(definition=dfDefinitions_processed_expanded%>%filter(TRAIT==trait),lst.harmonized.data$lst.data,dfData.settings,lst.harmonized.data$vct.identifiers)
DmRxT2_timeline


# TODO check the Hx part is not working
DmRxT2_status_by_source<-get_case_status_by_source(definition=dfDefinitions_processed_expanded%>%filter(TRAIT==trait),lst.harmonized.data$lst.data,df.data.settings,df.reference.dates = df_reference_dt_v0)
df_reference_dt_today$f.53.0.0<-Sys.Date()
DmRxT2_today_status_by_source<-get_case_status_by_source(definition=dfDefinitions_processed_expanded%>%filter(TRAIT==trait),lst.harmonized.data$lst.data,df.data.settings,df.reference.dates = df_reference_dt_today)
