
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

if (Sys.getenv("USER")=="niek"){
  repo_dir="/Users/niek/repos/ukbpheno/"
  pheno_dir="/Volumes/data/ukb/"
}else if (Sys.getenv("USER")=="ming") {
  pheno_dir="/home/mw/Analyses/Ukbpheno_data/"
  repo_dir="/home/mw/Repos/ukbpheno/"
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
lst.data<-lapply(lst.data,function(x) {setkey(x,code) }) # double check that everything has the same setkey. 
##########################################
toc()

##########################################
# generate meta data dynamically,  returns a list with the number of rows per code based on default_datatable_defCol_pair()
lst.counts <- get_lst_counts(lst.data,datatable_defCol_pair = default_datatable_defCol_pair() )

##########################################
##########################################
#### We should make it into one data-object, which will make things more dynamic.
### ukb.data.object <- list(data=lst.data,lst.counts=lst.counts,settings=default_datatable_defCol_pair())
### lst.counts <- NULL
### lst.data <- NULL

# to dataset make more lean: retain identifier , visitdates and additional fields needed for definitions besides default (as cols in file)
#dfukb<- dfukb[,dfhtml[dfhtml$field.showcase %in% c("eid", "53",ukb_fields$nondefault_ukb_fields),]$field.tab,with=FALSE]
# save 
save(dfhtml,dfukb,lst.data,lst.counts,file=fukbphenodata)
# load(fukbphenodata)
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
# lst.data[["tte.hesin.icd10.primary"]][.(lookupquery) ]


##########################################
## analyse 1 definition
##########################################
dfDefinitions_processed_expanded <- expand_dfDefinitions_processed(dfDefinitions_processed,datatable_defCol_pair=default_datatable_defCol_pair(),lst.counts = lst.counts)
all_event_dt <- get_all_events(dfDefinitions_processed_expanded[12,],lst.data) #list of 11 dfs 
all_event_dt[, .(count = .N, sum.event = sum(event,na.rm = T),sum.epidur= sum(epidur,na.rm = T),median.epidur= median(epidur,na.rm = T),max.epidur= max(epidur,na.rm=T)), by = f.eid]


# show stats on codes
stats.codes <- all_event_dt[, .(count=.N,sum.event = sum(event,na.rm = T),sum.epidur= sum(epidur,na.rm = T),median.epidur= median(epidur,na.rm = T),max.epidur= max(epidur,na.rm=T) ), keyby=list(f.eid,classification, code)]
stats.codes <- stats.codes %>% dplyr::group_by(classification, code) %>% summarise(count=n() )
stats.codes <- stats.codes %>% arrange(count)
stats.codes$rank <- 1:nrow(stats.codes)
p1 <- ggplot(stats.codes, aes(rank, count,label=code,color=classification)) + geom_point() + ylim(-((max(stats.codes$count))/3),NA)  + geom_text_repel(size =3,segment.size=0.5)
p2 <- ggplot(stats.codes, aes(rank, count,label=code,color=classification)) + scale_y_continuous(trans='log2') + geom_point()   + geom_text_repel(size =3,segment.size=0.5)
ggarrange(p1,p2,nrow = 1, ncol = 2,common.legend = TRUE)

# show co-occurences of classifications 
stats.coocurrence <- table(all_event_dt[,c("f.eid" ,"classification")])
stats.coocurrence[stats.coocurrence>0] <-1
mat <- crossprod(as.matrix(stats.coocurrence))     
mat <- floor(t(mat * 100 / diag(mat)))                 # calculate the percentage
pheatmap::pheatmap(mat,display_numbers=mat,cluster_cols = F,cluster_rows = F )

# show co-occurences of codes
stats.coocurrence <- table(all_event_dt[,c("f.eid" ,"code")])
stats.coocurrence[stats.coocurrence>0] <-1
mat <- crossprod(as.matrix(stats.coocurrence))     
mat <- floor(t(mat * 100 / diag(mat)))                 # calculate the percentage
pheatmap::pheatmap(mat,fontsize = 6,cluster_cols = F,cluster_rows = F)


reshape(d,direction = 'wide',v.names='classification',idvar='f.eid',timevar='value')
### get prevalence (move this.. )
get_incidence_prevalence <- function(all_event_dt,reference_date,include_secondary_recurrence=TRUE,datatable_defCol_pair=default_datatable_defCol_pair() ) {
  # event==0, event cannot be used forr age - of - diagnosiis or new events.
  # event==1, event can be used for any type of future events, as it is based on ICD10 type of data 
  # event==2, event can be used for any type of future events and age of diagnoses but only if there is no evidence in history, as it is self reported.
  
  # define window arround reference_dates? -- e.g. if you're diagnosed with XX, count follow-up events only after 20 days of no-events. 
  # define window to score the diagnosis on refereence date.
  # reference_dates # should be a vector of dates, with the identifier as name. 
  #reference_date = setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid)
  
  ### Future: 1) count any new event==1 code (icd10 etc) 
  ###         2) count any primary event==1 if Hx==1 (recurrence) if 'include_secondary_recurrence'==TRUE
  ###         3) countany new event==2 (self report) age of diagnosis in future
  
  
  # 4139206
  df.ref <- data.frame(reference_date)
  df.ref$pid <- rownames(df.ref)
  
  df <- merge(all_event_dt,df.ref,by.x = 'f.eid' , by.y='pid',all.x = T) %>% arrange(eventdate) %>% as.data.table()
  df <- df %>% filter(!is.na(reference_date)) # comment out if missing f.eids is fixed. 
  
  df$days <- df$eventdate - df$reference_date
  setkey(df,days) # i don't know why, but setkey was alreaday on f.eid and cannot refresh..
  setkey(df,f.eid)
  
  if(include_secondary_recurrence){
    sources_recurrence_events <- datatable_defCol_pair %>% filter(diagnosis==1) %>% pull(datasource)
  } else {
    sources_recurrence_events <- datatable_defCol_pair$datasource
  } 
  ### History
  dfHx <- df[days<=0]
  Hxd <- dfHx[event>0][, .(Hx_days=min(days,na.rm=T) ), keyby=list(f.eid)]
  dfHx[,Hx:=1]
  
  ### Future
  dfFu <- df[  ((days>(0) & event==1) & (!f.eid %in% dfHx$Hx)) |
               (days>(0) & event==1 & (f.eid %in% dfHx$Hx) & .id %in% sources_recurrence_events ) |
               (days>(0) & event==2  & (!f.eid %in% dfHx$Hx)) ]
  dfFu[,Fu:=1] #unique(dfFu$f.eid)
  Fud <- dfFu[,.(Fu_days= min(days,na.rm=T) ), keyby=list(f.eid)]
  
  ### age of diagnosis
  #system.time({ df %>% filter(event>0) %>% group_by(f.eid) %>% summarise(first_diagnosis_days=min(days)) }) # <- slow.. 
  #system.time({ df[df$event>0][,.(first_diagnosis_days=min(days,na.rm=T) ), by=f.eid] }) # <- fast..
  #df[df$event>0][df[, .I[which.max(days)], by=f.eid]$V1] # <- aanother way.. 
  first_diagnosis_days <- df[df$event>0][,.(first_diagnosis_days=min(days,na.rm=T) ), by=f.eid] # 
  # View(Reduce(function(...) merge(..., all = TRUE,by='f.eid'), list(df,
  #                                                                   Hxd,
  #                                                                   Fud,
  #                                                                   unique(dfHx[,c("f.eid","Hx")]),
  #                                                                   unique(dfFu[,c("f.eid","Fu")]),
  #                                                                   first_diagnosis_days)))

  df <- Reduce(function(...) merge(..., all = TRUE,by='f.eid'), list( Hxd,
                                                                Fud,
                                                                unique(dfHx[,c("f.eid","Hx")]),
                                                                unique(dfFu[,c("f.eid","Fu")]),
                                                                first_diagnosis_days))
  df[,Any:=1]
  df
  return(df)
}


df <- get_incidence_prevalence(all_event_dt = all_event_dt,reference_date = setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))




  
# print(format(object.size(lst.data), units = "Mb")) #"2014.1 Mb"




