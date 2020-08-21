
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

# fukbtab = paste(pheno_dir,"ukb38326.tab",sep="") 
fukbtab = paste(pheno_dir,"ukb41823.tab.head10k",sep="") # header only for testing.
# fukbtab = paste(pheno_dir,"ukb38326.tab.head300",sep="") # small subset only for testing.


fhtml = paste(pheno_dir,"ukb41823.html",sep="")
fhesin=paste(pheno_dir,"hesin.txt_",sep="")
fhesin_diag=paste(pheno_dir,"hesin_diag.txt_",sep="")
fhesin_oper=paste(pheno_dir,"hesin_oper.txt_",sep="")
fgp_clinical =paste(pheno_dir,"gp_clinical.txt_",sep="")
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


# GP 
# add 2 lists with read2 /read3
lst <- append(lst,read_gp_clinical_data(fgp=fgp_clinical ))
toc()

# meta data
# returns a new variable with the number of rows per code
lst.counts <- lapply(lst, function(x) x[, .N, by=.(code)] )

lst.counts$icd10 <- sumcounts(list(tte.hesin.icd10.primary=lst.counts$tte.hesin.icd10.primary,
                                   tte.hesin.icd10.secondary=lst.counts$tte.hesin.icd10.secondary,
                                   tte.death.icd10.primary=lst.counts$tte.death.icd10.primary,
                                   tte.death.icd10.secondary=lst.counts$tte.death.icd10.secondary)) # sumcounts function  in read-data.R
lst.counts$icd9 <- sumcounts(list(tte.hesin.icd9.primary=lst.counts$tte.hesin.icd9.primary,
                                  tte.hesin.icd9.secondary=lst.counts$tte.hesin.icd9.secondary)) # sumcounts function in read-data.R
lst.counts$oper3 <- sumcounts(list(tte.hesin.oper3.primary=lst.counts$tte.hesin.oper3.primary,
                                   tte.hesin.oper3.secondary=lst.counts$tte.hesin.oper3.secondary)) # sumcounts function in read-data.R
lst.counts$oper4 <- sumcounts(list(tte.hesin.oper4.primary=lst.counts$tte.hesin.oper4.primary,
                                   tte.hesin.oper4.secondary=lst.counts$tte.hesin.oper4.secondary)) # sumcounts function in read-data.R

# View(lst.counts$icd10)

# filter dfukb. 
# retain identifier , visitdates and additional fields needed for definitions besides default (as cols in file)
dfukb<- dfukb[,dfhtml[dfhtml$field.showcase %in% c("eid", "53",ukb_fields$nondefault_ukb_fields),]$field.tab,with=FALSE]


# save 
save(dfhtml,dfukb,lst,lst.counts,file=fukbphenodata)

View(lst$tte.death.icd10.primary)


# codes.ts <- dfDefinitions_processed[1,]$TS

View(dfDefinitions_processed)
# 
# 

codes.icd10 <- dfDefinitions_processed$ICD10CODES[8]
codes.icd10.expanded <- unique(lst.counts[['icd10']]$code[grep(paste(sep="","^",strsplit(codes.icd10,",")[[1]], collapse='|'),lst.counts$icd10$code)])

#grepoper <-unlist( mclapply(  , function(col) grep(paste(sep="","^",VctCodes, collapse='|'), col, ignore.case=FALSE),mc.cores =detectCores()/2 ) ) ## PARALLEL of the above.


 

# TODO GP SCRIPT
# TODO TOUCSCHREEN? 

# extract_case(dfDefinitions[1,],lst,dfukb[,c("f.eid", "f.53.0.0")] ) # extract case data. _per source (HESIN, NURSE-TOUCHSCREEN, GP, DEATH) and everything together? 
# make venn diagram// stats 

# filter dfukb to reduce memory? dont need these fields anymore? 
# keep <- !names(dfukb) %in% dfhtml[dfhtml$field.showcase %in% c("20001","20002","20004","20003","20006","20008","20010","53","40000","40001","40002"),]$field.tab
# dfukb <- dfukb[,..keep]

# # check object sizes
# print(format(object.size(lst), units = "Mb")) #"2014.1 Mb"

# 
# View(dfhes[(dfhes$level.x != dfhes$level.y) & !is.na(dfhes$level.x) & !is.na(dfhes$level.y),])
# View(dfukb[,grepl("4000", names(dfukb)),with=FALSE ])
# 
# # TOUCSCHREEN Self reported, unstructured. leave this unstructured? 
# # "6150=1[3627]" # (field == code (age))
# # "3581≥0[3581]" #Menopause, only age. 
# # 6153_=5 # anticonception
# # 20110=4 (Bowel cancer),20107=4, 20111=4 # family history. 
# lst$df_6150 <- read_ukb_data(f, fields_to_keep = c("6150","3627:numeric")) 
# 





