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
fdata_setting = paste(data_dir,"data.settings.tsv",sep="")

# ##########################################
# read setting
############################################
rm(dt.data.settings)
dt.data.settings <-fread(fdata_setting)
class(dt.data.settings)
####################################################
# # read definitions.
###################################################
dfDefinitions <- fread(fdefinitions, colClasses = 'character', data.table = FALSE)
dfDefinitions_processed <- ProcessDfDefinitions(dfDefinitions)
# Next we check the codes in the definition table against the code maps
# The code maps are either downloaded from UK biobank showcase (codingxx.tsv) [search corresponding Data-Coding in showcase] for which has been checked to include all codes present in data # or to be created from the current data at hand (.code)

# currently 2 coding systems exist for gp_clinical
cf_read3<-dt.data.settings[dt.data.settings$classification=="CTV3",]$code_map
cf_read2<-dt.data.settings[dt.data.settings$classification=="READ2",]$code_map
# corresponding column names from the .txt file
get_all_exsiting_codes(fgp_clinical,c("read_2","read_3"),c(paste0(data_dir,cf_read2),paste0(data_dir,cf_read3)))

#  3 gp_script
cf_dmd<-dt.data.settings[dt.data.settings$classification=="DMD",]$code_map
# bnf.england bnf.scotland are separated because they appear to be coded slightly differently according to the official documentation...hence different post-processing may be needed
cf_bnf<-unique(dt.data.settings[dt.data.settings$classification=="BNF",]$code_map)
cf_read2d<-dt.data.settings[dt.data.settings$classification=="READ2_drugs",]$code_map
get_all_exsiting_codes(fgp_scripts,c("read_2","bnf_code","dmd_code"),c(paste0(data_dir,cf_read2d),paste0(data_dir,cf_bnf),paste0(data_dir,cf_dmd)))




# With the code maps, we check the codes in the definition table
#
missing_codes<-check_dfDefinitions_codes(dfDefinitions_processed,dt.data.settings,code_map_dir=data_dir,F)
data.table::fwrite(list(missing_codes),paste(pheno_dir,"ukbphenodata_feb2021.Rdata_notInData.code"))

# WHAT IT MEANS: Codes that are not present will be removed in the expand_dfDefinitions...() as the they will be removed
# in theory these extra codes will not cause crashes of the pipeline unless the whole line is empty (no codes from any source)
dfDefinitions_processed_expanded <-expand_dfDefinitions_processed2(dfDefinitions_processed,dt.data.settings,code_map_dir =data_dir )

########################################################################
########################################################################
# TODO remove the hierchical part of the expand_dfDefinitions_processed2 () and grep everything?
# reason not to do this : ICD 9/10 codes are used HESIN (portal),DEATH (portal)and CANCER (ukb.tab column) -> to gather all available codes in the data means we need to first get lst.data which is the same as dfDefinitions_processed()!

########################################################################
########################################################################








# ##################################################
# # Prepare UKB data:
####################################################

#TODO wrap the below(?) session to this function
generate_harmonized_long_format_datatable <- function(f.ukbtab=NULL,f.html=NULL,fhesin=NULL,....){



  # # read ukb data from ukbconv (.html + .tab)
  message("Start data harmonization")


  if (!is.null(f.html)){
    message("Read metadata file ")
  # # ukb's .tab meta data
  dfhtml <- read_ukb_metadata(fhtml) # 9 columns in dfhtml: "field.number","field.count","field.showcase","field.html","field.tab","field.description","col.type","col.name","fread_column_type"`
  }else{
    message("Required metadata file (.html) is missing, please generate it with ukbconv. Exit!")
    return()
  }

  message("Verify all required fields from the definitions are present in the main dataset (.tab)")

  # # ukb's .tab file; extract only relevant fields
  # we then use this meta-data to check if all fields we need are inside the .tab file
  # the function outputs a list of fields that are required
  dfDefinitions_ukb_fields <- get_allvarnames(dfDefinitions_processed,dfhtml)
  if (is.null(dfDefinitions_ukb_fields)){
    message("Please ensure all required fields are present in the main dataset, abort!")
    return()
  }else{
  message("Read .tab file and retrieve required fields.")
  tictoc::tic()
  # Next is to extract those columns required from the .tab file
  dfukb <- read_ukb_tabdata(fukbtab,dfhtml,fields_to_keep = dfDefinitions_ukb_fields$all_ukb_fields)
  tictoc::toc("time elapsed processing .tab file")
  }
  vct.identifiers <- as.character(dfukb$f.eid)
  gc()


  ################################################################################
  # converting data
  ################################################################################
  lst.data <- list()
  ################################################################################
  # touchscreen (event==2: only the first occurence is an event)
  ################################################################################
  message("Convert touchscreen data")
  lst.data$ts <- convert_touchscreen_to_episodedata(dfukb,ts_conditions = dfDefinitions_processed$TS) #154.964sec
  ################################################################################
  # self reported data  (event==2: only the first occurence is an event)
  ################################################################################
  message("Convert self-reported data")
  message(" ->Self-reported cancer codes (field 20001)")
  # cancer
  lst.data$tte.sr.20001 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20001",field_sr_date = "20006",qc_threshold_year = 10)
  message(" ->Self-reported non-cancer codes(field 20002)")
  # non cancer
  lst.data$tte.sr.20002 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20002",field_sr_date = "20008",qc_threshold_year = 10)
  message(" ->Self-reported operation codes(field 20004)")
  # operation
  lst.data$tte.sr.20004 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20004",field_sr_date = "20010",qc_threshold_year = 10)
  message(" ->Self-reported medication codes(field 20003)")
  # medication (event==0, since no age of diagnosis)
  lst.data$sr.20003 <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "20003",field_sr_date = "",qc_threshold_year = 10)



  ################################################################################
  # cancer registry
  ################################################################################
  message("Convert cancer registry data(field 40006/40013)")
  # qc_threshold set to zeros, because cancer recurrence is fairly common (?)
  lst.data$tte.cancer.icd10 <- convert_cancerregister_to_episodedata(dfukb,field_diagnosis = "40006",field_date = "40005",field_date_type="date",qc_threshold_year = 0,codetype = "character")
  lst.data$tte.cancer.icd9 <- convert_cancerregister_to_episodedata(dfukb,field_diagnosis = "40013",field_date = "40005",field_date_type="date",qc_threshold_year = 0,codetype ="character")

  ################################################################################
  # death record
  ################################################################################
  # death registry from tab file (event==1: every occurence is a real event)
  # lst.data$tte.death.icd10.primary <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40001",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10,event_code=1) # death
  # lst.data$tte.death.icd10.secondary <- convert_nurseinterview_to_episodedata(dfukb,field_sr_diagnosis = "40002",field_sr_date = "40000",field_sr_date_type="date",qc_treshold_year = 10,event_code=1) # death
  # death registry from data portal , same data as the main dataset but more up to date, refer document DeathLinkage
  # this is merged with the death from tab file for completeness, but it is the same data now.
  message("Convert death registry data(data portal)")
  lst.data_dth <- read_death_data(fdeath_portal,fdeath_cause_portal)
  # take only the death record from portal for now
  lst.data$tte.death.icd10.primary <-lst.data_dth$primary
  lst.data$tte.death.icd10.secondary<-lst.data_dth$secondary
  rm(lst.data_dth)
  # View(lst.data$tte.death.icd10.primary)

  ################################################################################
  # hesin  (event==1)
  ################################################################################
  tictoc::tic("Convert Hospital Inpatient data(data portal)")

  lst.data <- append(lst.data,read_hesin_data(fhesin ,fhesin_diag ,fhesin_oper )) #tte.hes.primary + tte.hes.secondary

  tictoc::toc()

  ################################################################################
  # primary care, gp  (event==1)
  ################################################################################
  tictoc::tic("Convert GP diagnosis data(data portal")
  #individuals in gp clinical: 229436 ~55.6 millions entries
  lst.data <- append(lst.data,read_gp_clinical_data(fgp=fgp_clinical)) # 311.128 sec, also thousands!
  tictoc::toc()

  gc()
  tictoc::tic("Convert GP prescription data(data portal")
  # without prescriptions #individuals in gp script: 222104 ~57.8 million entries
  # this one does medication
  lst.data <- append(lst.data,read_gp_script_data(fgp=fgp_scripts)) # 2798.576 sec
  gc()
  tictoc::toc()

  # make sure eeverything is in the right format:
  lst.data <- lapply(lst.data,function(x) {setkey(x,code) })
  lst.data <- lapply(lst.data,function(x) {x[, ('f.eid') := lapply(.SD, as.character), .SDcols = 'f.eid'] })
  lst.data <- lapply(lst.data,function(x) {x[,'eventdate'] <-  round(x$eventdate);return(x) })


  return(lst.data,vct.identifiers)

}






# fukbphenodata <- paste(pheno_dir,"ukbphenodata_dataPortal_june2021_defv3.1.Rdata",sep="") #where to store final object



# tictoc::toc() #  1111.306 sec elapsed, 18min.

# here it is stored as .Rdata
save(dfhtml,dfukb,lst.data,dt.data.settings,file=fukbphenodata)



#####################################################################
# if the steps above were performed already simply load the R data

load(fukbphenodata)


#######################################################################
# make outcome files
######################################################################


  visitdatefields = names(dfukb)[grepl(paste0("[^0-9]","53","[^0-9]"), names(dfukb))]
  out_folder=paste(pheno_dir,"output/2021-02-10_HESIN-Jun21.def3.1/",sep="")

  if(!dir.exists(file.path(out_folder))){
    dir.create(file.path(out_folder))
  }
  out_visit_folders=vector()
  for (visit in c(0,2,3)){
  out_visit_folder=paste(out_folder,"visit",visit,"/",sep="")
  out_visit_folders<-c(out_visit_folders,out_visit_folder)
  if(!dir.exists(file.path(out_visit_folder))){
    dir.create(file.path(out_visit_folder))
  }
  }
  # dfwTs<-dfDefinitions_processed_expanded[!is.na(dfDefinitions_processed_expanded$TS),]
 gc()
  i=1


  tictoc::tic()
  msgcon <- file(paste(pheno_dir,"creatPhenotype-msg_Jun2021.def3.1.txt",sep=''), open = "a")
  sink(paste(pheno_dir,"creatPhenotype-out_Jun2021.def3.1.txt",sep=''), type = "output", append = TRUE, split = TRUE)
  sink(msgcon, type = "message")


  # for (trait in dfwTs[dfwTs$Definitions=="Include_in_cases","TRAIT"]){
  # loop over the trait
  for (trait in dfDefinitions_processed_expanded[dfDefinitions_processed_expanded$Definitions=="Include_in_cases","TRAIT"]) {
    message(glue::glue("Identify case and control for : {unique(dfDefinitions_processed_expanded[dfDefinitions_processed_expanded$TRAIT==trait,]$DESCRIPTION)}"))


    if (trait %in% dfDefinitions_processed_expanded[dfDefinitions_processed_expanded$Definitions=="Study_population","TRAIT"]){

      # !!!!!!if there is study population reference date will not be considered! HENCE dealt separately

      message(glue::glue("Reference date for **{unique(dfDefinitions_processed_expanded[dfDefinitions_processed_expanded$TRAIT==trait,]$DESCRIPTION)} ** will be first event of trait in Study Population"))
      lst.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT==trait), lst.data,dt.data.settings,reference_date = NULL,vct.identifiers)
      # keep a subset to reduce file size
      lst.case_control$df.casecontrol<-lst.case_control$df.casecontrol[,c("f.eid","survival_days",  "Death_any", "Hx_days","Fu_days","Hx","Fu","Any")]
      # rename col names
      colnames(lst.case_control$df.casecontrol) <- paste(trait,colnames(lst.case_control$df.casecontrol), sep = "_")
      names(lst.case_control$df.casecontrol)[names(lst.case_control$df.casecontrol) == paste( trait,"f.eid", sep = "_")]<-"f_eid"

      # this is for mergeall in stata which give type mismatch error if ID is not numeric
      lst.case_control$df.casecontrol$f_eid<-as.integer(lst.case_control$df.casecontrol$f_eid)
      # for stata which doesn't accept period in var name
      colnames(lst.case_control$df.casecontrol)<-gsub("\\.","_",colnames(lst.case_control$df.casecontrol))
      # write to file
      save_path<-paste(out_folder,paste("ukb41823",trait,"dta",sep="."),sep="")
      haven::write_dta(lst.case_control$df.casecontrol,path=save_path)

    }

    else{

      # query by visits baseline and imaging
      for (visit in c(0,2,3)){
        message(glue::glue("********************Query visit {visit}**********************"))
        # visit=0
        visitdatefield = visitdatefields[visit+1]
        # trait="RxChol"
        lst.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT==trait), lst.data,dt.data.settings, reference_date=setNames(as.Date(as.character(dfukb[[visitdatefield]]),format="%Y-%m-%d"),dfukb$f.eid))

        if (stringr::str_detect(trait,"Death")){
          # keep only death related columns to avoid confusion

          lst.case_control$df.casecontrol<-lst.case_control$df.casecontrol[,c("f.eid", "reference_date","survival_days","Death_primary", "Death_any", "Ref")]
        }else{

          # keep a subset to reduce file size
          lst.case_control$df.casecontrol<-lst.case_control$df.casecontrol[,c("f.eid","survival_days","Death_primary",  "Death_any", "Hx_days","Fu_days","Hx","Fu","Any")]

        }


        # rename col names

        colnames(lst.case_control$df.casecontrol) <- paste(trait,visit,colnames(lst.case_control$df.casecontrol),  sep = "_")
        names(lst.case_control$df.casecontrol)[names(lst.case_control$df.casecontrol) == paste(trait,visit,"f.eid",  sep = "_")]<-"f_eid"

        # this is for mergeall in stata which give type mismatch error if ID is not numeric
        lst.case_control$df.casecontrol$f_eid<-as.integer(lst.case_control$df.casecontrol$f_eid)

        # for stata which doesn't accept period in var name
        colnames(lst.case_control$df.casecontrol)<-gsub("\\.","_",colnames(lst.case_control$df.casecontrol))
        # write to file
        save_path<-paste(out_folder,"visit",visit,paste("/ukb41823",trait,visit,"dta",sep="."),sep="")
        haven::write_dta(lst.case_control$df.casecontrol,path=save_path)

           }

    }
        print(paste("*******Finish phenotype #",paste(i,':',trait,"****************",sep=" "),sep=""))
        i<-i+1

    }
# case_TS<-lst.case_control$all_event_dt.Include_in_cases[lst.case_control$all_event_dt.Include_in_cases$.id=="TS",]
# unique(lst.case_control$all_event_dt.Include_in_cases$.id)

  sink(NULL, type = "message")
  sink(NULL, type = "output")
  gc()
  message(glue::glue("DONE!Processed {i} phenotypes."))
  tictoc::toc()

# trait="Cad"
# lst.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT==trait), lst.data,dt.data.settings,reference_date = NULL,vct.identifiers)
# View(lst.case_control$df.casecontrol)
# colnames(lst.case_control$df.casecontrol)
# rm(trait)
# rm(lst.case_control)
