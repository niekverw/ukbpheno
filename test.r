
#library(disk.frame)
#library(ukbtools)  # chmod -R u+w /usr/local/Cellar/r

source("/repos/ukbpheno/convert_nurseinterview_to_episodedata.r")
source("/repos/ukbpheno/ProcessdfDefinitions.R")
source("/repos/ukbpheno/read-data.R")

fukb = "/Volumes/data/ukb/ukb38326.tab"
fukb = "/Volumes/data/ukb/ukb38326.tab.head"
fhtml = "/Volumes/data/ukb/ukb38326.html"
fdefinitions = "/repos/ukbpheno/definitions.tsv"

dfDefinitions <- fread(fdefinitions, colClasses = 'character', data.table = FALSE)
dfDefinitions_processed <- ProcessDfDefinitions(dfDefinitions)
fields_to_keep <- get_allvarnames(dfDefinitions_processed)

dfhtml <- read_ukb_metadata(fhtml)
dfukb <- read_ukb_data(fukb,dfhtml,fields_to_keep = fields_to_keep)


# install.packages('unixtools', repos = 'http://www.rforge.net/'); unixtools::set.tempdir("/new/tmp/path") 

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

print(format(object.size(lst$tte.sr_20002), units = "Mb"))

# HESIN (data is not unique for eid, code; could be used for reevents) contains duration 

# GP (data is not unique for eid, code)
View(dfukb[,grepl("4000", names(dfukb)),with=FALSE ])

# TOUCSCHREEN Self reported, unstructured. leave this unstructured? 
# "6150=1[3627]" # (field == code (age))
# "3581≥0[3581]" #Menopause, only age. 
# 6153_=5 # anticonception
# 20110=4 (Bowel cancer),20107=4, 20111=4 # family history. 
lst$df_6150 <- read_ukb_data(f, fields_to_keep = c("6150","3627:numeric")) 


