
#library(ukbtools)  # chmod -R u+w /usr/local/Cellar/r
library(matrixStats)
library(fasttime)
# library(disk.frame) # not really working like I thought.. 
library(data.table)
library(dplyr)
library(tictoc)
require(XML)
require(stringr)
require(disk.frame)

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
fsr_coding=paste(repo_dir,"data/20003.coding4.tsv",sep="")
fcncr_coding=paste(repo_dir,"data/20001.coding3.tsv",sep="")
fnoncncr_coding=paste(repo_dir,"data/20002.coding6.tsv",sep="")
foper_coding=paste(repo_dir,"data/20004.coding5.tsv",sep="")
fdeath_portal=paste(pheno_dir,"death.txt",sep="")
fdeath_cause_portal=paste(pheno_dir,"death_cause.txt",sep="")




fukbphenodata <- paste(pheno_dir,"ukbphenodata.Rdata",sep="") #where to store final object






# write.table(dfCodesheetREAD_SR.Coding, file = "CodesheetREAD_SR.Coding", append = FALSE, quote = TRUE, sep = "\t",
#             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
#             col.names = TRUE, qmethod = c("escape", "double"),
#             fileEncoding = "")
#  med field 20003 code table 
dfCodesheetREAD_SR.Coding<-fread(fsr_coding)



# read definitions. 
dfDefinitions <- fread(fdefinitions, colClasses = 'character', data.table = FALSE)
dfDefinitions_processed <- ProcessDfDefinitions(dfDefinitions)



ukb_fields <- get_allvarnames(dfDefinitions_processed)



# read ukb data from ukbconv (.html + .tab)
# 9 columns in dfhtml: "field.number","field.count","field.showcase","field.html","field.tab","field.description","col.type","col.name","fread_column_type"
dfhtml <- read_ukb_metadata(fhtml)
# TODO: unit of operation - phenotype such that read only required fields from the ukbxxxxx.tab file 
# TODO: create a loop : for each definition in definitions ,  read ukbdata + fetch all required info + generate the variable 
dfukb <- read_ukb_data(fukbtab,dfhtml,fields_to_keep = ukb_fields$all_ukb_fields)
print(format(object.size(dfukb), units = "Gb"))
ncol(dfukb)

tic("converting data")
lst <- list()
# SELF REPORT (TIME TO EVENT DATA); data is unique (eid,code) contains first events, but also some without date, so can't say if its event. 
lst$tte.sr.20002 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20002",field_sr_date = "20008",qc_treshold_year = 10) # non cancer
lst$tte.sr.20001 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20001",field_sr_date = "20006",qc_treshold_year = 10) # cancer
lst$tte.sr.20004 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20004",field_sr_date = "20010",qc_treshold_year = 10) # operation
# MEDICATION (data is not unique for eid, code)
lst$sr.20003 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20003",field_sr_date = "",qc_treshold_year = 10) 
# DEATH REGISTRY REPORT
lst$tte.death.icd10.primary <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40001",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10) # death
## secondary death, use only 1 date... 
lst$tte.death.icd10.secondary <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40002",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10) # death
# DEATH from data portal , same data as the main dataset but more up to date, refer document DeathLinkage
lst_dth<-read_death_data(fdeath_portal,fdeath_cause_portal)
# merge the records   dplyr union
lst$tte.death.icd10.primary <-union(lst_dth$primary,lst$tte.death.icd10.primary)
lst$tte.death.icd10.secondary <-union(lst_dth$secondary,lst$tte.death.icd10.secondary)
rm(lst_dth)
# HESIN (data is not unique for eid, code; could be used for reevents) contains duration 
# add 8 lists from HES data primary/secondary x oper3/oper4/icd9/icd10 , with columns eid,eventdate,epidur,<diag>,event
lst <- append(lst,read_hesin_data(fhesin ,fhesin_diag ,fhesin_oper )) #tte.hes.primary + tte.hes.secondary
# GP # add 2 lists with read2 /read3
lst <- append(lst,read_gp_clinical_data(fgp=fgp_clinical ))
toc()

# generate meta data dynamically,  returns a list with the number of rows per code based on default_datatable_defCol_pair()
lst.counts <- get_lst_counts(lst)
# View(lst.counts$ICD10)
# to dataset make more lean: retain identifier , visitdates and additional fields needed for definitions besides default (as cols in file)
#dfukb<- dfukb[,dfhtml[dfhtml$field.showcase %in% c("eid", "53",ukb_fields$nondefault_ukb_fields),]$field.tab,with=FALSE]
# save 
save(dfhtml,dfukb,lst,lst.counts,file=fukbphenodata)

##### test:
### get a vector with definitions to use as input. 
Vctdef <- dfDefinitions_processed[10,c("TRAIT","DESCRIPTION", unique(default_datatable_defCol_pair()))]

### ICD10 expand first using counts, then use datatable lookup:
lookupquery <- strsplit(Vctdef[["ICD10"]],split = ",")[[1]] 
lookuptarget <- lst.counts[['ICD10']][,get("code")] 
lookupquery <- grep( paste(sep="","^",lookupquery, collapse='|'),lookuptarget,ignore.case = T,value = T)
lst[["tte.hesin.icd10.secondary"]][.(lookupquery) ] 

### note that self reported data is numeric: maybe we can store this in default_datatable_defCol_pair
lookupquery <- as.numeric(strsplit(Vctdef[["n_20002"]],split = ",")[[1]])
lst[["tte.sr.20002"]][.(lookupquery)] ## lookup multiple codes super fast. 


# # check object sizes
# print(format(object.size(lst), units = "Mb")) #"2014.1 Mb"

# # TOUCSCHREEN Self reported, unstructured. leave this unstructured? 
# # "6150=1[3627]" # (field == code (age))
# # "3581≥0[3581]" #Menopause, only age. 
# # 6153_=5 # anticonception
# # 20110=4 (Bowel cancer),20107=4, 20111=4 # family history. 
# lst$df_6150 <- read_ukb_data(f, fields_to_keep = c("6150","3627:numeric")) 
# 





