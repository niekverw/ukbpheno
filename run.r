# TODO: Add output for death (primary and primary+secondary) in get_incidence_prevalence() )
#library(ukbtools)  # chmod -R u+w /usr/local/Cellar/r
library(matrixStats)
library(fasttime)
# library(disk.frame) # not really working like I thought.. 
library(data.table)
library(dplyr)
library(tictoc)
require(XML)
require(stringr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(pheatmap)
library(ggdendro)

if (Sys.getenv("USER")=="niek"){
  repo_dir="/Users/niek/repos/ukbpheno/"
  pheno_dir="/Volumes/data/ukb/"
}else if (Sys.getenv("USER")=="mw") {
  pheno_dir="/home/mw/Analyses/Ukbpheno_data/"
  repo_dir="/home/mw/Repos/ukbpheno/"
}

source(paste(repo_dir,"convert_nurseinterview_to_episodedata.r",sep=""))
source(paste(repo_dir,"ProcessdfDefinitions.R",sep=""))
source(paste(repo_dir,"read-data-ukb-functions.R",sep=""))
source(paste(repo_dir,"query-event-tables.r",sep=""))
source(paste(repo_dir,"process_event_tables.r",sep=""))
source(paste(repo_dir,"convert_touchscreen_to_episodedata.r",sep=""))


### set paths. 
fukbtab = paste(pheno_dir,"ukb41823.tab",sep="") 
#fukbtab = paste(pheno_dir,"ukb41823.tab.head10k",sep="") # header only for testing.
# fukbtab = paste(pheno_dir,"ukb38326.tab.head300",sep="") # small subset only for testing.
fhtml = paste(pheno_dir,"ukb41823.html",sep="")
fhesin=paste(pheno_dir,"hesin.txt",sep="")
fhesin_diag=paste(pheno_dir,"hesin_diag.txt",sep="")
fhesin_oper=paste(pheno_dir,"hesin_oper.txt",sep="")
fgp_clinical =paste(pheno_dir,"gp_clinical.txt",sep="")
# fgp_scripts =paste(pheno_dir,"gp_scripts.txt",sep="")  
fdefinitions = paste(repo_dir,"definitions.tsv",sep="")
fdeath_portal=paste(pheno_dir,"death.txt",sep="")
fdeath_cause_portal=paste(pheno_dir,"death_cause.txt",sep="")
# output file:
fukbphenodata <- paste(pheno_dir,"ukbphenodata.Rdata",sep="") #where to store final object

##########################################
# read definitions. 
dfDefinitions <- fread(fdefinitions, colClasses = 'character', data.table = FALSE)
dfDefinitions_processed <- ProcessDfDefinitions(dfDefinitions)
########################################## 
# Prepare UKB data: 
# read ukb data from ukbconv (.html + .tab) 
tic("converting data")
# ukb's .tab meta data
dfhtml <- read_ukb_metadata(fhtml) # 9 columns in dfhtml: "field.number","field.count","field.showcase","field.html","field.tab","field.description","col.type","col.name","fread_column_type"
# ukb's .tab file; extract only relevant fields 
dfDefinitions_ukb_fields <- get_allvarnames(dfDefinitions_processed) ## regarding downstream functions, what happens if å field is not present, will we get an error? 
dfukb <- read_ukb_tabdata(fukbtab,dfhtml,fields_to_keep = dfDefinitions_ukb_fields$all_ukb_fields) # 439.52 sec
print(format(object.size(dfukb), units = "Gb"))
# converting data
lst.data <- list()
# touscreen, which uses information from the dfDefinitions_processed (event==2: only the first occurence is an event)
lst.data$ts <- convert_touchscreen_to_episodedata(dfukb,ts_conditions = dfDefinitions_processed$TS)
# self reported data  (event==2: only the first occurence is an event)
lst.data$tte.sr.20001 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20001",field_sr_date = "20006",qc_treshold_year = 10) # cancer
lst.data$tte.sr.20002 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20002",field_sr_date = "20008",qc_treshold_year = 10) # non cancer
lst.data$tte.sr.20004 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20004",field_sr_date = "20010",qc_treshold_year = 10) # operation
# medication (event==0, since no age of diagnosis)
lst.data$sr.20003 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20003",field_sr_date = "",qc_treshold_year = 10) 
# death registry from tab file (event==1: every occurence is a real event)
lst.data$tte.death.icd10.primary <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40001",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10) # death
lst.data$tte.death.icd10.secondary <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40002",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10) # death
# death registry from data portal , same data as the main dataset but more up to date, refer document DeathLinkage
# this is merged with the death from tab file for completeness, but it is the same data now. 
lst.data_dth<-read_death_data(fdeath_portal,fdeath_cause_portal)
lst.data$tte.death.icd10.primary <-union(lst.data_dth$primary,lst.data$tte.death.icd10.primary)
lst.data$tte.death.icd10.secondary <-union(lst.data_dth$secondary,lst.data$tte.death.icd10.secondary)
rm(lst.data_dth)
# hesin  (event==1)
lst.data <- append(lst.data,read_hesin_data(fhesin ,fhesin_diag ,fhesin_oper )) #tte.hes.primary + tte.hes.secondary
# primary care, gp  (event==1)
lst.data <- append(lst.data,read_gp_clinical_data(fgp=fgp_clinical )) # 462.085 

#attr(lst.data,all_identifiers <- as.character(dfukb$f.eid))

# setkey on code. 
toc() #  1111.306 sec elapsed, 18min. 

# load lst.data.settings
lst.data.settings <- data.frame(fread("
    datasource classification  datatype  expand_codes  diagnosis ignore.case
    tte.sr.20002  f.20002 numeric 0 1 FALSE
    tte.sr.20001  f.20001 numeric 0 1 FALSE
    tte.sr.20004  f.20004 numeric 0 1 FALSE
    sr.20003  f.20003 numeric 0 1 FALSE
    tte.death.icd10.primary ICD10 character 1 1 TRUE
    tte.death.icd10.secondary ICD10 character 1 2 TRUE
    tte.hesin.oper3.primary OPCS3 character 1 1 TRUE
    tte.hesin.oper3.secondary OPCS3 character 1 2 TRUE
    tte.hesin.oper4.primary OPCS4 character 1 1 TRUE
    tte.hesin.oper4.secondary OPCS4 character 1 2 TRUE
    tte.hesin.icd10.primary ICD10 character 1 1 TRUE
    tte.hesin.icd10.secondary ICD10 character 1 2 TRUE
    tte.hesin.icd9.primary  ICD9  character 1 1 TRUE
    tte.hesin.icd9.secondary  ICD9  character 1 2 TRUE
    tte.gpclincal.read2 READ  character 0 2 FALSE
    tte.gpclincal.read3 CTV3  character 0 2 FALSE
    tte.gpscript.dmd.england  DMD character 0 2 FALSE
    tte.gpscript.bnf.england  BNF character 0 2 FALSE
    tte.gpscript.bnf.scotland BNF character 0 2 FALSE
    tte.gpscript.read2.wales  READ  character 1 2 FALSE
    ts  TS  character 0 1 TRUE"))

# generate meta data dynamically,  returns a list with the number of rows per code based on lst.data.settings
lst.counts <- get_lst_counts(lst.data,lst.data.settings = lst.data.settings)

lst.identifiers <- as.character(dfukb$f.eid)
# save 
save(dfhtml,dfukb,lst.data,lst.data.settings,lst.counts,file=fukbphenodata)
# load(fukbphenodata)

##########################################
## analyse 1 single definition
##########################################
# expand the definitions based on the data that is loaded
dfDefinitions_processed_expanded <- expand_dfDefinitions_processed(dfDefinitions_processed,lst.data.settings=lst.data.settings,lst.counts = lst.counts)
#all collapsed to 1 datatable
all_event_dt <- get_all_events(dfDefinitions_processed_expanded[14,],lst.data,lst.data.settings)   #MI
# all_event_dt <- get_all_events(dfDefinitions_processed_expanded[8,],lst.data)  #DmT2
# all_event_dt <- get_all_events(dfDefinitions_processed_expanded[17,],lst.data)  #Ht #all collapsed to 1 datatable
# all_event_dt <- get_all_events(dfDefinitions_processed_expanded[9,],lst.data)  #DmT2
all_event_dt.stats <- get_stats_for_events(all_event_dt) #should generate several plots
all_event_dt.stats$stats.codes.summary.p
all_event_dt.stats$stats.class.cooccur.p
all_event_dt.stats$stats.codes.cooccur.filtered.p.dendro
all_event_dt.stats$stats.codes.cooccur.filtered.p.heat

# get incidence/prevalence from baseline
all_event_dt.summary <- get_incidence_prevalence(all_event_dt = all_event_dt,lst.data.settings, reference_date = setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
View(all_event_dt.summary %>% filter(is.na(Hx) & is.na(Fu)))

# get occcurence  from first event and recurrence of primary events: 
all_event_dt.summary <- get_incidence_prevalence(all_event_dt = all_event_dt,lst.data.settings,reference_date = NULL,window_fu_days_mask = 15)
hist(all_event_dt.summary %>% pull(Fu_days) %>% as.numeric,breaks=100)
hist(all_event_dt.summary %>% filter(Fu_days<100) %>% pull(Fu_days) %>% as.numeric,breaks=100)
hist(all_event_dt.summary$reference_date,breaks=200) ## <-  actual first date of diagnosis
print(format(object.size(lst.data), units = "Mb")) #"2014.1 Mb"

##########################################
############# DEFINE CASE & CONTROLS for 1 definition, Tryout; please check on consistency, features it should have and if everything behaves ok... it can get quite confusing with the reference date, 
##########################################
### get only get cases (with exclusions)
TEST <- get_cases(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Cad"), lst.data,lst.data.settings, reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))

### get everything; population, cases (with exclusions), and controls (with exclusions)
TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Nicm"), lst.data,lst.data.settings, reference_date=NULL)
TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Nicm"), lst.data,lst.data.settings,  reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="HfInCad"), lst.data,lst.data.settings, reference_date=NULL)
TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="HfInCad"), lst.data,lst.data.settings, reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Cad"), lst.data,lst.data.settings, reference_date=NULL)
TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Cad"), lst.data,lst.data.settings, reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Cad"), lst.data,lst.data.settings, reference_date=setNames(as.Date(as.character(dfukb$f.53.1.0),format="%Y-%m-%d"),dfukb$f.eid))
## ~50% has family history of heart disease, that is a alot !?
TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="HxHrt"), lst.data,lst.data.settings, reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))



gc()