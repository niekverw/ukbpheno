
#library(ukbtools)  # chmod -R u+w /usr/local/Cellar/r
library(matrixStats)
library(fasttime)
# library(disk.frame) # not really working like I thought.. 
library(data.table)
library(dplyr)
library(tictoc)
require(XML)
require(stringr)

if (Sys.getenv("USER")=="niek"){
  repo_dir="/Users/niek/repos/ukbpheno/"
  pheno_dir="/Volumes/data/ukb/"
}else if (Sys.getenv("USER")=="ming") {
  pheno_dir="/home/ming/UKB/Ukbpheno_data/"
  repo_dir="/home/ming/Repos/ukbpheno/"
}

source(paste(repo_dir,"convert_nurseinterview_to_episodedata.r",sep=""))
source(paste(repo_dir,"ProcessdfDefinitions.R",sep=""))
source(paste(repo_dir,"read-data.R",sep=""))
source(paste(repo_dir,"query-event-tables.r",sep=""))

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

# ### These aren't used, are they? 
# fsr_coding=paste(repo_dir,"data/20003.coding4.tsv",sep="")
# fcncr_coding=paste(repo_dir,"data/20001.coding3.tsv",sep="")
# fnoncncr_coding=paste(repo_dir,"data/20002.coding6.tsv",sep="")
# foper_coding=paste(repo_dir,"data/20004.coding5.tsv",sep="")

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
lst <- list()
##########################################
# SELF REPORT (TIME TO EVENT DATA); data is unique (eid,code) contains first events, but also some without date, so can't say if its event. 
lst$tte.sr.20002 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20002",field_sr_date = "20008",qc_treshold_year = 10) # non cancer
lst$tte.sr.20001 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20001",field_sr_date = "20006",qc_treshold_year = 10) # cancer
lst$tte.sr.20004 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20004",field_sr_date = "20010",qc_treshold_year = 10) # operation
##########################################
# MEDICATION (data is not unique for eid, code)
lst$sr.20003 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20003",field_sr_date = "",qc_treshold_year = 10) 
##########################################
# DEATH REGISTRY REPORT
lst$tte.death.icd10.primary <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40001",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10) # death
## secondary death, use only 1 date... 
lst$tte.death.icd10.secondary <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40002",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10) # death
# DEATH from data portal , same data as the main dataset but more up to date, refer document DeathLinkage
lst_dth<-read_death_data(fdeath_portal,fdeath_cause_portal)
# merge the records   dplyr union
lst$tte.death.icd10.primary <-union(lst_dth$primary,lst$tte.death.icd10.primary)
lst$tte.death.icd10.secondary <-union(lst_dth$secondary,lst$tte.death.icd10.secondary)
setkey(lst$tte.death.icd10.primary,'code')
setkey(lst$tte.death.icd10.secondary,'code')
rm(lst_dth)
##########################################
# HESIN (data is not unique for eid, code; could be used for reevents) contains duration 
# add 8 lists from HES data primary/secondary x oper3/oper4/icd9/icd10 , with columns eid,eventdate,epidur,<diag>,event
lst <- append(lst,read_hesin_data(fhesin ,fhesin_diag ,fhesin_oper )) #tte.hes.primary + tte.hes.secondary
# GP # add 2 lists with read2 /read3
lst <- append(lst,read_gp_clinical_data(fgp=fgp_clinical ))
lst<-lapply(lst,function(x) {setkey(x,code) }) # double check that everything has the same setkey. 
##########################################
toc()

##########################################
# generate meta data dynamically,  returns a list with the number of rows per code based on default_datatable_defCol_pair()
lst.counts <- get_lst_counts(lst,datatable_defCol_pair = default_datatable_defCol_pair() )
##########################################
# View(lst.counts$ICD10)
# to dataset make more lean: retain identifier , visitdates and additional fields needed for definitions besides default (as cols in file)
#dfukb<- dfukb[,dfhtml[dfhtml$field.showcase %in% c("eid", "53",ukb_fields$nondefault_ukb_fields),]$field.tab,with=FALSE]
# save 
save(dfhtml,dfukb,lst,lst.counts,file=fukbphenodata)
##########################################
#load(fukbphenodata)
##### I changed get_all_events() code a little bit with this  principle:
# ### 1) get a vector with definitions to use as input. 
# definitions <- dfDefinitions_processed[10,]
# ### 2)   expand ICD codes  using the counts, now incorporated in the function expand_dfDefinitions_processed()
# lookupquery <- strsplit(definitions[["ICD10"]],split = ",")[[1]] 
# lookuptarget <- lst.counts[['ICD10']][,get("code")] 
# lookupquery <- grep( paste(sep="","^",lookupquery, collapse='|'),lookuptarget,ignore.case = T,value = T)
# ### 3) then use datatable lookup:
# lst[["tte.hesin.icd10.primary"]][.(lookupquery) ]


##########################################
## test: 
dfDefinitions_processed_expanded <- expand_dfDefinitions_processed(dfDefinitions_processed,datatable_defCol_pair=default_datatable_defCol_pair(),lst.counts = lst.counts)
test <- get_all_events(dfDefinitions_processed_expanded[9,],lst) #list of 11 dfs 
# some dfs are NA on every column, not sure why that is.. 


# print(format(object.size(lst), units = "Mb")) #"2014.1 Mb"




