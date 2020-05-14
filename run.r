
#library(ukbtools)  # chmod -R u+w /usr/local/Cellar/r
library(matrixStats)
library(fasttime)
# library(disk.frame) # not really working like I thought.. 
library(data.table)
library(dplyr)
library(tictoc)
require(XML)

source("/Users/niek/repos/ukbpheno/convert_nurseinterview_to_episodedata.r")
source("/Users/niek/repos/ukbpheno/ProcessdfDefinitions.R")
source("/Users/niek/repos/ukbpheno/read-data.R")

fukbtab = "/Volumes/data/ukb/ukb38326.tab"
#fukbtab = "/Volumes/data/ukb/ukb38326.tab.head" # header only for testing. 
fhtml = "/Volumes/data/ukb/ukb38326.html"
fhesin="/Volumes/data/ukb/hesin.txt"
fhesin_diag="/Volumes/data/ukb/hesin_diag.txt"
fhesin_oper="/Volumes/data/ukb/hesin_oper.txt"
fgp_clinical = "/Volumes/data/ukb/gp_clinical.txt"
fdefinitions = "/Users/niek/repos/ukbpheno/definitions.tsv"

# read definitions. 
dfDefinitions <- fread(fdefinitions, colClasses = 'character', data.table = FALSE)
dfDefinitions_processed <- ProcessDfDefinitions(dfDefinitions)
fields_to_keep <- get_allvarnames(dfDefinitions_processed)

# read ukb data from ukbconv (.html + .tab)
dfhtml <- read_ukb_metadata(fhtml)
dfukb <- read_ukb_data(fukbtab,dfhtml,fields_to_keep = fields_to_keep)

tic("converting data")
lst <- list()
# SELF REPORT (TIME TO EVENT DATA); data is unique (eid,code) contains first events, but also some without date, so can't say if its event. 
lst$tte.sr_20002 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20002",field_sr_date = "20008",qc_treshold_year = 10) # non cancer
lst$tte.sr_20001 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20001",field_sr_date = "20006",qc_treshold_year = 10) # cancer
lst$tte.sr_20004 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20004",field_sr_date = "20010",qc_treshold_year = 10) # operation
# MEDICATION (data is not unique for eid, code)
lst$sr_medication <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20003",field_sr_date = "",qc_treshold_year = 10) 
# DEATH RECORDS
lst$tte.sr_40001 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40001",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10) # operation
## secondary death, use only 1 date... 
lst$tte.sr_40002 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40002",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10) # operation
# HESIN (data is not unique for eid, code; could be used for reevents) contains duration 
toc();tic("convert HESIN data")
lst <- append(lst,read_hesin_data(fhesin ,fhesin_diag ,fhesin_oper )) #tte.hes.primary + tte.hes.secondary
# GP 
lst <- append(lst,read_gp_clinical_data(fgp=fgp_clinical ))
toc()

# filter dfukb to reduce memory? dont need these fields anymore? 
keep <- !names(dfukb) %in% dfhtml[dfhtml$field.showcase %in% c("20001","20002","20004","20003","20006","20008","20010","53","40000","40001","40002"),]$field.tab
dfukb <- dfukb[,..keep]

# # check object sizes
print(format(object.size(lst), units = "Mb")) #"2014.1 Mb"

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





