
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
source(paste(repo_dir,"read-data.R",sep=""))
source(paste(repo_dir,"query-event-tables.r",sep=""))
source(paste(repo_dir,"process_event_tables.r",sep=""))

fukbtab = paste(pheno_dir,"ukb41823.tab",sep="") 
fukbtab = paste(pheno_dir,"ukb41823.tab.head10k",sep="") # header only for testing.
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

fukbphenodata <- paste(pheno_dir,"ukbphenodata.Rdata",sep="") #where to store final object

##########################################
# read definitions. 
dfDefinitions <- fread(fdefinitions, colClasses = 'character', data.table = FALSE)
dfDefinitions_processed <- ProcessDfDefinitions(dfDefinitions)
dfDefinitions_ukb_fields <- get_allvarnames(dfDefinitions_processed)
########################################## 
# > loading data
# read ukb data from ukbconv (.html + .tab)
# 9 columns in dfhtml: "field.number","field.count","field.showcase","field.html","field.tab","field.description","col.type","col.name","fread_column_type"
dfhtml <- read_ukb_metadata(fhtml)
# TODO: unit of operation - phenotype such that read only required fields from the ukbxxxxx.tab file 
# TODO: create a loop : for each definition in definitions ,  read ukbdata + fetch all required info + generate the variable 
dfukb <- read_ukb_tabdata(fukbtab,dfhtml,fields_to_keep = dfDefinitions_ukb_fields$all_ukb_fields)
print(format(object.size(dfukb), units = "Gb"))

#################################################################################### \
# Should we put the following in a function?
tic("converting data")
lst.data <- list()
##########################################
# SELF REPORT (TIME TO EVENT DATA); data is unique (eid,code) contains first events, but also some without date, so can't say if its event. 
lst.data$tte.sr.20002 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20002",field_sr_date = "20008",qc_treshold_year = 10) # non cancer
lst.data$tte.sr.20001 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20001",field_sr_date = "20006",qc_treshold_year = 10) # cancer
lst.data$tte.sr.20004 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20004",field_sr_date = "20010",qc_treshold_year = 10) # operation
##########################################
# MEDICATION (data is not unique for eid, code)
lst.data$sr.20003 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20003",field_sr_date = "",qc_treshold_year = 10) 
##########################################
# DEATH REGISTRY REPORT
lst.data$tte.death.icd10.primary <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40001",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10) # death
## secondary death, use only 1 date... 
lst.data$tte.death.icd10.secondary <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40002",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10) # death
# DEATH from data portal , same data as the main dataset but more up to date, refer document DeathLinkage
lst.data_dth<-read_death_data(fdeath_portal,fdeath_cause_portal)
# merge the records   dplyr union
lst.data$tte.death.icd10.primary <-union(lst.data_dth$primary,lst.data$tte.death.icd10.primary)
lst.data$tte.death.icd10.secondary <-union(lst.data_dth$secondary,lst.data$tte.death.icd10.secondary)
setkey(lst.data$tte.death.icd10.primary,'code')
setkey(lst.data$tte.death.icd10.secondary,'code')
rm(lst.data_dth)
##########################################
# HESIN (data is not unique for eid, code; could be used for reevents) contains duration 
# add 8 lists from HES data primary/secondary x oper3/oper4/icd9/icd10 , with columns eid,eventdate,epidur,<diag>,event
lst.data <- append(lst.data,read_hesin_data(fhesin ,fhesin_diag ,fhesin_oper )) #tte.hes.primary + tte.hes.secondary
# GP # add 2 lists with read2 /read3
lst.data <- append(lst.data,read_gp_clinical_data(fgp=fgp_clinical ))
#TS
lst.data$ts <- convert_touchscreen_to_episodedata(dfukb,ts_conditions = dfDefinitions_processed$TS)

lst.data<-lapply(lst.data,function(x) {setkey(x,code) }) # double check that everything has the same setkey. 
##########################################

lst.data.settings <- data.frame(fread("datasource	classification	datatype	expand_codes	diagnosis	ignore.case
tte.sr.20002	n_20002	numeric	0	1	FALSE
tte.sr.20001	n_20001	numeric	0	1	FALSE
tte.sr.20004	n_20004	numeric	0	1	FALSE
sr.20003	n_20003	numeric	0	1	FALSE
tte.death.icd10.primary	ICD10	character	1	1	TRUE
tte.death.icd10.secondary	ICD10	character	1	2	TRUE
tte.hesin.oper3.primary	OPCS3	character	1	1	TRUE
tte.hesin.oper3.secondary	OPCS3	character	1	2	TRUE
tte.hesin.oper4.primary	OPCS4	character	1	1	TRUE
tte.hesin.oper4.secondary	OPCS4	character	1	2	TRUE
tte.hesin.icd10.primary	ICD10	character	1	1	TRUE
tte.hesin.icd10.secondary	ICD10	character	1	2	TRUE
tte.hesin.icd9.primary	ICD9	character	1	1	TRUE
tte.hesin.icd9.secondary	ICD9	character	1	2	TRUE
tte.gpclincal.read2	READ	character	0	2	FALSE
tte.gpclincal.read3	CTV3	character	0	2	FALSE
tte.gpscript.dmd.england	DMD	character	0	2	FALSE
tte.gpscript.bnf.england	BNF	character	0	2	FALSE
tte.gpscript.bnf.scotland	BNF	character	0	2	FALSE
tte.gpscript.read2.wales	READ	character	1	2	FALSE
ts	TS	character	0	1	TRUE"))
toc()

##########################################
# generate meta data dynamically,  returns a list with the number of rows per code based on default_datatable_defCol_pair()
lst.counts <- get_lst_counts(lst.data,datatable_defCol_pair = default_datatable_defCol_pair() )
# to dataset make more lean: retain identifier , visitdates and additional fields needed for definitions besides default (as cols in file)
#dfukb<- dfukb[,dfhtml[dfhtml$field.showcase %in% c("eid", "53",ukb_fields$nondefault_ukb_fields),]$field.tab,with=FALSE]
# save 
save(dfhtml,dfukb,lst.data,lst.data.settings,lst.counts,file=fukbphenodata)
# load(fukbphenodata)
##########################################

##########################################
## analyse 1 definition
##########################################
# dfDefinitions_processed_expanded[12,]$TRAIT #"Cad"
# expand the definitions as exact match much faster 
dfDefinitions_processed_expanded <- expand_dfDefinitions_processed(dfDefinitions_processed,datatable_defCol_pair=default_datatable_defCol_pair(),lst.counts = lst.counts)

all_event_dt <- get_all_events(dfDefinitions_processed_expanded[12,],lst.data)  #all collapsed to 1 datatable
dev.off()
rm(all_event_dt.stats)
all_event_dt.stats <- get_stats_for_events(all_event_dt) #should generate several plots
all_event_dt.stats$stats.codes.summary.p
all_event_dt.stats$stats.class.cooccur.p
all_event_dt.stats$stats.codes.cooccur.filtered.p.dendro
all_event_dt.stats$stats.codes.cooccur.filtered.p.heat

all_event_dt.summary <- get_incidence_prevalence(all_event_dt = all_event_dt,reference_date = setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
View(all_event_dt.summary %>% filter(is.na(Hx) & is.na(Fu)))

### without referencec datet, study recurrencce of primary events: 
all_event_dt.summary <- get_incidence_prevalence(all_event_dt = all_event_dt,reference_date = NULL,window_fu_days_mask = 15)
hist(all_event_dt.summary %>% pull(Fu_days) %>% as.numeric,breaks=100)
hist(all_event_dt.summary %>% filter(Fu_days<100) %>% pull(Fu_days) %>% as.numeric,breaks=100)
hist(all_event_dt.summary$reference_date,breaks=200)
print(format(object.size(lst.data), units = "Mb")) #"2014.1 Mb"

all_event_dt.summary

