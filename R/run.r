# # TODO: Add output for death (primary and primary+secondary) in get_incidence_prevalence() )
# #library(ukbtools)  # chmod -R u+w /usr/local/Cellar/r
# library(matrixStats)
# library(fasttime)
# # library(disk.frame) # not really working like I thought..
# library(data.table)
# library(dplyr)
# library(tictoc)
# require(XML)
# require(stringr)
# library(ggplot2)
# library(ggrepel)
# library(ggpubr)
# library(pheatmap)
# library(ggdendro)
# library(readxl)
#
# if (Sys.getenv("USER")=="niek"){
#   repo_dir="/Users/niek/repos/ukbpheno/R/"
#   pheno_dir="/Volumes/data/ukb/"
# }else if (Sys.getenv("USER")=="mw") {
#   pheno_dir="/home/mw/Analyses/Ukbpheno_data/"
#   repo_dir="/home/mw/Repos/ukbpheno/R/"
# }
#
# source(paste(repo_dir,"convert_nurseinterview_to_episodedata.r",sep=""))
# source(paste(repo_dir,"ProcessdfDefinitions.R",sep=""))
# source(paste(repo_dir,"read-data-ukb-functions.R",sep=""))
# source(paste(repo_dir,"query-event-tables.r",sep=""))
# source(paste(repo_dir,"process_event_tables.r",sep=""))
# source(paste(repo_dir,"convert_touchscreen_to_episodedata.r",sep=""))
# source(paste(repo_dir,"plot_timeline.R",sep=""))
#
#
# ### set paths.
# fukbtab = paste(pheno_dir,"ukb41823.tab",sep="")
# #fukbtab = paste(pheno_dir,"ukb41823.tab.head10k",sep="") # header only for testing.
# # fukbtab = paste(pheno_dir,"ukb38326.tab.head300",sep="") # small subset only for testing.
# fhtml = paste(pheno_dir,"ukb41823.html",sep="")
# fhesin=paste(pheno_dir,"hesin.txt",sep="")
# fhesin_diag=paste(pheno_dir,"hesin_diag.txt",sep="")
# fhesin_oper=paste(pheno_dir,"hesin_oper.txt",sep="")
# fgp_clinical =paste(pheno_dir,"gp_clinical.txt",sep="")
# # fgp_scripts =paste(pheno_dir,"gp_scripts.txt",sep="")
# fdeath_portal=paste(pheno_dir,"death.txt",sep="")
# fdeath_cause_portal=paste(pheno_dir,"death_cause.txt",sep="")
# # output file:
# fukbphenodata <- paste(pheno_dir,"ukbphenodata.Rdata",sep="") #where to store final object
# fdefinitions = paste(repo_dir,"definitions.tsv",sep="")
# fdata_setting = paste(repo_dir,"data.settings.tsv",sep="")
# # code_map_dir<-paste(repo_dir,"data/",sep="")
# ##########################################
# # read definitions.
# dfDefinitions <- fread(fdefinitions, colClasses = 'character', data.table = FALSE)
# dfDefinitions_processed <- ProcessDfDefinitions(dfDefinitions)
#
#
#
# ##########################################
# # Prepare UKB data:
# # read ukb data from ukbconv (.html + .tab)
# tic("converting data")
# # ukb's .tab meta data
# dfhtml <- read_ukb_metadata(fhtml) # 9 columns in dfhtml: "field.number","field.count","field.showcase","field.html","field.tab","field.description","col.type","col.name","fread_column_type"
# # ukb's .tab file; extract only relevant fields
# dfDefinitions_ukb_fields <- get_allvarnames(dfDefinitions_processed,dfhtml) ## regarding downstream functions, what happens if å field is not present, will we get an error? YES a check now is performed in get_allvarnames()
#
#
# dfukb <- read_ukb_tabdata(fukbtab,dfhtml,fields_to_keep = dfDefinitions_ukb_fields$all_ukb_fields) # 439.52 sec
# print(format(object.size(dfukb), units = "Gb"))
# # load lst.data.settings
# lst.data.settings <-fread(fdata_setting)
# lst.identifiers <- as.character(dfukb$f.eid)
#
# # converting data
# lst.data <- list()
# # touscreen, which uses information from the dfDefinitions_processed (event==2: only the first occurence is an event)
# lst.data$ts <- convert_touchscreen_to_episodedata(dfukb,ts_conditions = dfDefinitions_processed$TS)
# # self reported data  (event==2: only the first occurence is an event)
# lst.data$tte.sr.20001 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20001",field_sr_date = "20006",qc_treshold_year = 10) # cancer
# lst.data$tte.sr.20002 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20002",field_sr_date = "20008",qc_treshold_year = 10) # non cancer
# lst.data$tte.sr.20004 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20004",field_sr_date = "20010",qc_treshold_year = 10) # operation
# # medication (event==0, since no age of diagnosis)
# lst.data$sr.20003 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20003",field_sr_date = "",qc_treshold_year = 10)
# # death registry from tab file (event==1: every occurence is a real event)
# # lst.data$tte.death.icd10.primary <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40001",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10,event_code=1) # death
# # lst.data$tte.death.icd10.secondary <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40002",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10,event_code=1) # death
# # death registry from data portal , same data as the main dataset but more up to date, refer document DeathLinkage
# # this is merged with the death from tab file for completeness, but it is the same data now.
# lst.data_dth <- read_death_data(fdeath_portal,fdeath_cause_portal)
#
#
#
# # take only the death record from portal for now
# lst.data$tte.death.icd10.primary <-lst.data_dth$primary
# lst.data$tte.death.icd10.secondary<-lst.data_dth$secondary
# rm(lst.data_dth)
# View(lst.data$tte.death.icd10.primary)
# # hesin  (event==1)
# lst.data <- append(lst.data,read_hesin_data(fhesin ,fhesin_diag ,fhesin_oper )) #tte.hes.primary + tte.hes.secondary
# # primary care, gp  (event==1)
#
# lst.data <- append(lst.data,read_gp_clinical_data(fgp=fgp_clinical,min_instance = lst.data.settings[lst.data.settings$classification=="CTV3",minimum_instance])) # 462.085 - > #822.694 sec with the instance filter!
#
#
#
#
# # make sure eeverything is in the right format:
# lst.data <- lapply(lst.data,function(x) {setkey(x,code) })
# lst.data <- lapply(lst.data,function(x) {x[, ('f.eid') := lapply(.SD, as.character), .SDcols = 'f.eid'] })
# lst.data <- lapply(lst.data,function(x) {x[,'eventdate'] <-  round(x$eventdate);return(x) })
#
# toc() #  1111.306 sec elapsed, 18min.
#
#
# # generate meta data dynamically,  returns a list with the number of rows per code based on lst.data.settings
# # lst.counts <- get_lst_counts(lst.data,lst.data.settings = lst.data.settings)
#
# # save
# save(dfhtml,dfukb,lst.data,lst.data.settings,lst.counts,file=fukbphenodata)
# # load(fukbphenodata)
#
#
# tic()
# test_expandDef<-expand_dfDefinitions_processed2(dfDefinitions_processed,lst.data.settings)
# toc()
# #4.949 sec elapsed with grep read2_drugs     # 2.699  taking only codes that exists in data
# #1.231 sec elapsed without
# # orginal codes in definition table : f1,f2,x006f  (insulin)
# length (unlist(strsplit(test_expandDef$READ2_drugs[4], ","))) #267 from 2 codes f1,f2 # 190taking only codes that exists in data
#
#
# ##########################################
# ## analyse 1 single definition
# ##########################################
# # expand the definitions based on the data that is loaded
# dfDefinitions_processed_expanded <- expand_dfDefinitions_processed(dfDefinitions_processed,lst.data.settings=lst.data.settings,lst.counts = lst.counts)
# dfDefinitions_processed_expanded <-expand_dfDefinitions_processed2(dfDefinitions_processed,lst.data.settings)
# #all collapsed to 1 datatable
# all_event_dt <- get_all_events(dfDefinitions_processed_expanded[14,],lst.data,lst.data.settings)   #MI
#
#
# # TODO refine get_stats_for_events
# all_event_dt.stats <- get_stats_for_events(all_event_dt) #should generate several plots
# View(all_event_dt.stats$stats.codes.summary.table)
# all_event_dt.stats$stats.codes.summary.p
# all_event_dt.stats$stats.class.cooccur.p
# all_event_dt.stats$stats.codes.cooccur.filtered.p.dendro
# all_event_dt.stats$stats.codes.cooccur.filtered.p.heat
#
# # get incidence/prevalence from baseline
# all_event_dt.summary <- get_incidence_prevalence(all_event_dt = all_event_dt,lst.data.settings, reference_date = setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
# View(all_event_dt.summary %>% filter(is.na(Hx) & is.na(Fu)))
#
# #########################################################################################
# # TODO: Add output for death (primary and primary+secondary) in get_incidence_prevalence() )
# all_event_dt.summary <- get_incidence_prevalence(all_event_dt = all_event_dt,lst.data.settings, reference_date = NULL)
# head(all_event_dt.summary)
# # death_event_dt.summary<- get_survival_data(dfDefinitions_processed_expanded[14,],lst.data,lst.data.settings) #1559
# # death_event_dt.summary<- get_survival_data(dfDefinitions_processed_expanded[14,],lst.data,lst.data.settings,window_days_mask = 5)
#
# #########################################################################################
#
# # get occcurence  from first event and recurrence of primary events:
# all_event_dt.summary <- get_incidence_prevalence(all_event_dt = all_event_dt,lst.data.settings,reference_date = NULL,window_fu_days_mask = 15)
# hist(all_event_dt.summary %>% pull(Fu_days) %>% as.numeric,breaks=100)
# hist(all_event_dt.summary %>% filter(Fu_days<100) %>% pull(Fu_days) %>% as.numeric,breaks=100)
# hist(all_event_dt.summary$reference_date,breaks=200) ## <-  actual first date of diagnosis
# print(format(object.size(lst.data), units = "Mb")) #"2014.1 Mb"
#
# ##########################################
# ############# DEFINE CASE & CONTROLS for 1 definition, Tryout; please check on consistency, features it should have and if everything behaves ok... it can get quite confusing with the reference date,
# ##########################################
# ### get only get cases (with exclusions)
# TEST <- get_cases(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Nicm"), lst.data,lst.data.settings, reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
#
#
# ### get everything; population, cases (with exclusions), and controls (with exclusions)
# TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Nicm"), lst.data,lst.data.settings, reference_date=NULL,lst.identifiers)
# TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Nicm"), lst.data,lst.data.settings,  reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
# TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="HfInCad"), lst.data,lst.data.settings, reference_date=NULL,lst.identifiers)
#
#
#
#
#
# TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="HfInCad"), lst.data,lst.data.settings, reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
# TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Cad"), lst.data,lst.data.settings, reference_date=NULL,lst.identifiers)
# TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Cad"), lst.data,lst.data.settings, reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
# TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="Cad"), lst.data,lst.data.settings, reference_date=setNames(as.Date(as.character(dfukb$f.53.1.0),format="%Y-%m-%d"),dfukb$f.eid))
# ## ~50% has family history of heart disease, that is a alot !?
# TEST <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT=="HxHrt"), lst.data,lst.data.settings, reference_date=setNames(as.Date(as.character(dfukb$f.53.0.0),format="%Y-%m-%d"),dfukb$f.eid))
#
# ### get first date of diagnosis by combining reference_date +  first_diagnosis_days (should triple check for some specific participants if dates etc are correct, appropriate NA's etc. . found one error alreday)
# hist(TEST$df.casecontrol %>% filter(Any==2) %>% mutate(date_of_first_diag=reference_date+first_diagnosis_days) %>% pull(date_of_first_diag),breaks=200 )
#
#
#
#
# ###########individual timeline########################
# lst.data.f.eid<-lapply(lst.data,function(x) {x[, ('f.eid') := lapply(.SD, as.numeric), .SDcols = 'f.eid'] }) # set eid to numeric
# lst.data.f.eid<-lapply(lst.data,function(x) {setkey(x,f.eid) }) # double check that everything has the same setkey.
#
# plot_individual_timeline(lst.data.settings,NULL,lst.data.f.eid,identifier="1234567")
#
