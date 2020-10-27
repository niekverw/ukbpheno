#!/usr/bin/env Rscript

# idea from
# https://termbrowser.nhs.uk/?perspective=full&conceptId1=404684003&edition=uk-edition&release=v20200610&server=https://termbrowser.nhs.uk/sct-browser-api/snomed&langRefset=999001261000000100,999000691000001104
# snomedbrowser.com

# require(optparse)
# require(data.table)
# require(readxl)
# require(dplyr)
# require(stringr)
# require(shiny)
# require(DT)
# require(ukbpheno)

required_pckg = c(
  "optparse",
  "data.table",
  "readxl",
  "dplyr",
  "stringr",
  "shiny",
  "DT",
  "devtools",
  "ukbpheno"
)
new.pckg <-
  required_pckg[!(required_pckg %in% installed.packages()[, "Package"])]
if (length(new.pckg) > 0) {
  print("*********The following required packages are not installed*********")
  print(new.pckg)
  print("Trying to install the required packages")
  for (pckg in new.pckg) {
    install.packages(pckg)
    if (require(pckg)) {
      print("successfully installed and loaded")
      print(pckg)
    } else {
      print(pckg)
      stop("fail to install required packag
e, exit.")
    }
  }
} else{
  for (pckg in required_pckg[!(required_pckg %in% new.pckg)]) {
    print(pckg)
    suppressPackageStartupMessages(library(pckg, character.only = TRUE))

  }
}



# source("ProcessdfDefinitions.R")
# library(readxl)
# library(dplyr)
# library(data.table)
# library(stringr)

# library(optparse)
# library(data.table)
# library(ukbpheno)

option_list <- list(
  optparse::make_option(
    "--fcoding_xls",
    type = "character",
    default = NULL,
    help = "Filepath to clinical coding classification systems and maps"
  ),
  optparse::make_option(
    "--f_med_readSR",
    type = "character",
    default = "../extdata/dfCodesheetREAD_SR.Coding.RData",
    help = "Filepath to code map beteen READ medication codes and self report medication codes in RData format(included in package)"
  ),
  optparse::make_option(
    "--fcoding_icd10",
    type = "character",
    default = "../extdata/ICD10.coding19.tsv",
    help = "Filepath to ICD10 coding classification system"
  ),
  optparse::make_option(
    "--fcoding_icd9",
    type = "character",
    default = "../extdata/ICD9.coding87.tsv",
    help = "Filepath to ICD9 coding classification system"
  ),
  optparse::make_option(
    "--fcoding_opcs4",
    type = "character",
    default = "../extdata/OPCS4.coding240.tsv",
    help = "Filepath to OPCS4 coding classification system"
  ),
  optparse::make_option(
    "--fcoding_20003",
    type = "character",
    default = "../extdata/20003.coding4.tsv",
    help = "Filepath to code map for field 20003"
  ),
  optparse::make_option(
    "--fbnf",
    type = "character",
    default = "../extdata/read2_bnf",
    help = "Filepath to code map between read2 and bnf"
  ),
  optparse::make_option(
    "--fdmd",
    type = "character",
    default = "../extdata/read2_dmd",
    help = "Filepath to code map between read2 and DMD"
  )

)


opt_parser <- optparse::OptionParser(option_list = option_list)

opt <- optparse::parse_args(opt_parser)

if (is.null(opt$fcoding_xls)) {
  print_help(opt_parser)
  stop("Please provide the code maps (.xlsx). Refer resource 592 in Showcase.",
       call. = FALSE)
}
if (is.null(opt$f_med_readSR)) {
  print_help(opt_parser)
  stop("Please check the path to code map for self-report & READ medications",
       call. = FALSE)

}




load_data <-
  function(fcoding.xls,
           fmed.readSR,
           ficd10,
           ficd9,
           fopcs4,
           f20003,fread2_bnf,fread2_dmd) {
    #### READ V2
    dfCodesheet.read_v2_read_ctv3 <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "read_v2_read_ctv3"))[, c(2, 7)]
    colnames(dfCodesheet.read_v2_read_ctv3) <- c("read_code", "CTV3")
    dfCodesheet.read_v2_icd9 <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "read_v2_icd9"))[, c(1, 2)]
    dfCodesheet.read_v2_icd10 <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "read_v2_icd10"))[, c(1, 2)]
    dfCodesheet.read_v2_opcs4 <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "read_v2_opcs4"))[, c(1, 2)]

    dfCodesheet.read_v2_lkp <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "read_v2_lkp"))
    dfCodesheet.read_v2_lkp <-
      as.data.frame(dfCodesheet.read_v2_lkp %>% dplyr::arrange(read_code, term_code))
    #dfCodesheet.read_v2_lkp <- dfCodesheet.read_v2_lkp[dfCodesheet.read_v2_lkp$term_code==0,]
    dfCodesheet.read_v2_bnf <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "read_v2_drugs_bnf"))[, c(1, 2)]
    # According to documentation: "BNF codes in the TPP extract follow the format 00.00.00.00.00. However, the coding structure does not
    # always map to codes provided by the NHSBSA. The first six digits of the code typically relate to BNF
    # chapter, section and paragraph in the NHSBSA code lists, although this is not consistent. Digits 7 and 8
    # do not appear to correspond to subparagraphs in the NHSBSA codes, and digits 9 and 10 are always
    # coded as 00."
    # It is therefore necessary to preserve the "." inbetween to know this is TPP BNF code??
    # dfCodesheet.read_v2_bnf$bnf_code<-gsub("\\.","",dfCodesheet.read_v2_bnf$bnf_code)
    colnames(dfCodesheet.read_v2_bnf)<-c("read_code", "bnf")
    dfCodesheet.read_v2_bnf<-dfCodesheet.read_v2_bnf[dfCodesheet.read_v2_bnf$read_code!=""& dfCodesheet.read_v2_bnf$bnf!="",]
    dfCodesheet.read_v2_bnf<-dfCodesheet.read_v2_bnf[!is.na(dfCodesheet.read_v2_bnf$read_code) & !is.na(dfCodesheet.read_v2_bnf$bnf),]

    dfCodesheet.read_v2_dmd <-
    as.data.frame(data.table::fread(fread2_dmd,sep = "\t",colClasses = c("character","character","character")))[,c(1,2)]


    colnames(dfCodesheet.read_v2_dmd)<-c("read_code", "dmd")
    dfCodesheet.read_v2_dmd<-dfCodesheet.read_v2_dmd[dfCodesheet.read_v2_dmd$read_code!=""& dfCodesheet.read_v2_dmd$dmd!="",]
    dfCodesheet.read_v2_dmd<-dfCodesheet.read_v2_dmd[dfCodesheet.read_v2_dmd$dmd!="0",]
    dfCodesheet.read_v2_dmd<-dfCodesheet.read_v2_dmd[!is.na(dfCodesheet.read_v2_dmd$read_code)& !is.na(dfCodesheet.read_v2_dmd$dmd),]


    dfCodesheet.read_v2_drugs_lkp <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "read_v2_drugs_lkp"))
    dfCodesheet.read_v2_lkp <-
      rbind(dfCodesheet.read_v2_lkp[, c(1, 3)], dfCodesheet.read_v2_drugs_lkp[, 1:2])
    dfCodesheet.read_v2_lkp <-
      as.data.frame(
        dfCodesheet.read_v2_lkp %>% dplyr::group_by(read_code) %>% dplyr::summarize(text = stringr::str_c(term_description, collapse = "/"))
      )
    #dfCodesheet.read_v2_lkp$read_code <- gsub("\\.","", dfCodesheet.read_v2_lkp$read_code)
    dfCodesheet.read_v2_lkp <-
      dfCodesheet.read_v2_lkp %>% unique()  %>% dplyr::arrange(read_code)
    dfCodesheet.read_v2_lkp$text <-
      stringr::str_replace_all(dfCodesheet.read_v2_lkp$text, "[^/[:^punct:]]", "")

    # meds; read to ukb code
    message("read dfCodesheetREAD_SR.Coding.RData")
    load(fmed.readSR)
    colnames(dfCodesheetREAD_SR.Coding) <-
      c("UKB.Coding", "text", "read_code", "term.id")
    dfCodesheet.read_v2_UKBmeds <- dfCodesheetREAD_SR.Coding[, c(3, 1)]
    dfCodesheet.read_v2_lkp <-
      rbind(dfCodesheet.read_v2_lkp, dfCodesheetREAD_SR.Coding[!dfCodesheetREAD_SR.Coding$READ.CODE %in% dfCodesheet.read_v2_lkp$read_code, c("read_code", "text")])

    dfCodesheet.READ <-
      merge(
        dfCodesheet.read_v2_read_ctv3,
        dfCodesheet.read_v2_icd9,
        by = "read_code",
        all = T
      )
    dfCodesheet.READ <-
      merge(dfCodesheet.READ,
            dfCodesheet.read_v2_icd10,
            by = "read_code",
            all = T)
    dfCodesheet.READ <-
      merge(dfCodesheet.READ,
            dfCodesheet.read_v2_opcs4,
            by = "read_code",
            all = T)
    dfCodesheet.READ <-
      merge(dfCodesheet.READ,
            dfCodesheet.read_v2_UKBmeds,
            by = "read_code",
            all = T)
    dfCodesheet.READ <-
      merge(dfCodesheet.READ,
            dfCodesheet.read_v2_bnf,
            by = "read_code",
            all = T)
    dfCodesheet.READ <-
      merge(dfCodesheet.READ,
            dfCodesheet.read_v2_dmd,
            by = "read_code",
            all = T)
    dfCodesheet.READ <-
      merge(dfCodesheet.READ,
            dfCodesheet.read_v2_lkp,
            by = "read_code",
            all = T)

    dfCodesheet.READ <- unique(dfCodesheet.READ)
    colnames(dfCodesheet.READ) <-
      c("READ", "CTV3", "ICD10", "ICD9", "OPCS4", "n_20003","BNF","DMD","text")
    dfCodesheet.READ$source = "READ"
    dfCodesheet.READ <- data.table::data.table(dfCodesheet.READ)
    dfCodesheet.READ$text <-
      stringr::str_replace_all(dfCodesheet.READ$text, "[^/[:^punct:]]", "")

    data.table::setkey(dfCodesheet.READ, "READ")
    ##############################
    ##### READ V3
    message("read ctv3 to read2")

    dfCodesheet.read_ctv3_readv2 <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "read_ctv3_read_v2"))[, c(1, 5)]
    dfCodesheet.read_ctv3_readv2 <-
      dfCodesheet.read_ctv3_readv2[!dfCodesheet.read_ctv3_readv2$READV2_CODE %in% "_NONE", ]
    colnames(dfCodesheet.read_ctv3_readv2) <- c("read_code", "READ")

    dfCodesheet.read_ctv3_icd9 <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "read_ctv3_icd9"))[, c(1, 2)]
    dfCodesheet.read_ctv3_icd10 <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "read_ctv3_icd10"))[, c(1, 2)]
    dfCodesheet.read_ctv3_opcs4 <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "read_ctv3_opcs4"))[, c(1, 2)]
    dfCodesheet.read_ctv3_lkp <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "read_ctv3_lkp"))[, c(1, 2)]
    dfCodesheet.read_ctv3_lkp <-
      as.data.frame(
        dfCodesheet.read_ctv3_lkp %>% dplyr::group_by(read_code) %>%  dplyr::summarize(text = stringr::str_c(term_description, collapse = "/"))
      )

    dfCodesheet.CTV3 <-
      merge(
        dfCodesheet.read_ctv3_readv2,
        dfCodesheet.read_ctv3_icd10,
        by = "read_code",
        all = T
      )
    dfCodesheet.CTV3 <-
      merge(dfCodesheet.CTV3,
            dfCodesheet.read_ctv3_icd9,
            by = "read_code",
            all = T)
    dfCodesheet.CTV3 <-
      merge(dfCodesheet.CTV3,
            dfCodesheet.read_ctv3_opcs4,
            by = "read_code",
            all = T)
    dfCodesheet.CTV3 <-
      merge(dfCodesheet.CTV3,
            dfCodesheet.read_ctv3_lkp,
            by = "read_code",
            all = T)
    dfCodesheet.CTV3 <- unique(dfCodesheet.CTV3)
    colnames(dfCodesheet.CTV3) <-
      c("CTV3", "READ", "ICD10", "ICD9", "OPCS4", "text")
    dfCodesheet.CTV3$n_20003 <- NA
    dfCodesheet.CTV3$BNF <- NA
    dfCodesheet.CTV3$DMD <- NA

    dfCodesheet.CTV3$source = "CTV3"
    dfCodesheet.CTV3$text <-
      stringr::str_replace_all(dfCodesheet.CTV3$text, "[^/[:^punct:]]", "")

    dfCodesheet.CTV3 <- data.table::data.table(dfCodesheet.CTV3)
    data.table::setkey(dfCodesheet.CTV3, "CTV3")

    ####################BNF
    message("BNF")
    dfCodesheet.bnf_readv2 <-as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "read_v2_drugs_bnf"))[, c(2,1)]

                                               colnames(dfCodesheet.bnf_readv2)<-c( "bnf","read_code")
    # see above for not removing the dot!
    # dfCodesheet.bnf_readv2$bnf<-gsub("\\.","",dfCodesheet.bnf_readv2$bnf)
    dfCodesheet.bnf_readv2<-dfCodesheet.bnf_readv2[dfCodesheet.bnf_readv2$bnf!="",]
    dfCodesheet.bnf_readv2<-dfCodesheet.bnf_readv2[ !is.na(dfCodesheet.bnf_readv2$bnf),]

    dfCodesheet.bnf_lookup<-as.data.frame(data.table::fread(fread2_bnf,sep = "\t"))[,c(2,3)]
    colnames(dfCodesheet.bnf_lookup)<-c( "bnf","text")
    # dfCodesheet.bnf_lookup$bnf<-gsub("\\.","",dfCodesheet.bnf_lookup$bnf)

    # dfCodesheet.bnf_lookup[dfCodesheet.bnf_lookup$bnf=="",]

    dfCodesheet.bnf_lookup<-dfCodesheet.bnf_lookup[dfCodesheet.bnf_lookup$bnf!="",]
    dfCodesheet.bnf_lookup<-dfCodesheet.bnf_lookup[ !is.na(dfCodesheet.bnf_lookup$bnf),]

    dfCodesheet.BNF <-
      merge(
        dfCodesheet.bnf_readv2,
        dfCodesheet.bnf_lookup,
        by = "bnf",
        all = T
      )

    dfCodesheet.BNF <- unique(dfCodesheet.BNF)
    colnames(dfCodesheet.BNF) <-
      c("BNF", "READ","text")
    dfCodesheet.BNF$n_20003 <- NA
    dfCodesheet.BNF$ICD10 <- NA
    dfCodesheet.BNF$ICD9 <- NA
    dfCodesheet.BNF$OPCS4 <- NA
    dfCodesheet.BNF$CTV3 <- NA
    dfCodesheet.BNF$DMD <- NA
    dfCodesheet.BNF$source = "BNF"
    dfCodesheet.BNF$text <-
      stringr::str_replace_all(dfCodesheet.BNF$text, "[^/[:^punct:]]", "")

    dfCodesheet.BNF <- data.table::data.table(dfCodesheet.BNF)
    data.table::setkey(dfCodesheet.BNF, "BNF")





    #######################DMD



    dfCodesheet.read_v2_dmd <-
      as.data.frame(data.table::fread(fread2_dmd))
    colnames(dfCodesheet.read_v2_dmd)<-c("read_code", "dmd")
    dfCodesheet.read_v2_dmd<-dfCodesheet.read_v2_dmd[dfCodesheet.read_v2_dmd$read_code!=""& dfCodesheet.read_v2_dmd$dmd!="",]
    dfCodesheet.read_v2_dmd<-dfCodesheet.read_v2_dmd[!is.na(dfCodesheet.read_v2_dmd$read_code)& !is.na(dfCodesheet.read_v2_dmd$dmd),]


    dfCodesheet.DMD<-as.data.frame(data.table::fread(fread2_dmd,sep = "\t",colClasses = c("character","character","character")))
    colnames(dfCodesheet.DMD)<-c( "READ","DMD","text")
    # dfCodesheet.bnf_lookup$bnf<-gsub("\\.","",dfCodesheet.bnf_lookup$bnf)
    dfCodesheet.DMD$DMD<-as.character(dfCodesheet.DMD$DMD)
    # dfCodesheet.bnf_lookup[dfCodesheet.bnf_lookup$bnf=="",]
    # clean out non-useful mapping
    dfCodesheet.DMD<-dfCodesheet.DMD[dfCodesheet.DMD$DMD!="",]
    dfCodesheet.DMD<-dfCodesheet.DMD[dfCodesheet.DMD$DMD!="0",]

    dfCodesheet.DMD<-dfCodesheet.DMD[ !is.na(dfCodesheet.DMD$DMD),]

    dfCodesheet.DMD <- unique(dfCodesheet.DMD)

    dfCodesheet.DMD$n_20003 <- NA
    dfCodesheet.DMD$ICD10 <- NA
    dfCodesheet.DMD$ICD9 <- NA
    dfCodesheet.DMD$OPCS4 <- NA
    dfCodesheet.DMD$CTV3 <- as.character(NA)
    dfCodesheet.DMD$BNF <- NA
    dfCodesheet.DMD$source = "DMD"
    dfCodesheet.DMD$text <-
      stringr::str_replace_all(dfCodesheet.DMD$text, "[^/[:^punct:]]", "")

    dfCodesheet.DMD <- data.table::data.table(dfCodesheet.DMD)
    data.table::setkey(dfCodesheet.DMD, "DMD")




    #######################33


    ####### ICD10
    message("ICD10")
    dfCodesheet.icd10_icd9 <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "icd9_icd10"))[, c(3, 1)]
    colnames(dfCodesheet.icd10_icd9) <- c("ICD10", "ICD9")
    dfCodesheet.icd10_lkp <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "icd10_lkp"))[, c(2, 5)] # not complete (e.g. X*)
    names(dfCodesheet.icd10_lkp) <- c("ICD10", "text")
    dfCodesheet.ICD10.coding19 <-
      data.frame(data.table::fread(ficd10))[, 1:2] # <- contains all of the above.
    names(dfCodesheet.ICD10.coding19) <- c("ICD10", "text")

    dfCodesheet.icd10_lkp <-
      rbind(dfCodesheet.icd10_lkp, dfCodesheet.ICD10.coding19[!dfCodesheet.ICD10.coding19$ICD10 %in% dfCodesheet.icd10_lkp$ICD10, ])

    dfCodesheet.ICD10 <-
      merge(dfCodesheet.icd10_icd9,
            dfCodesheet.icd10_lkp,
            by = "ICD10",
            all = T)

    dfCodesheet.ICD10$READ <- NA
    dfCodesheet.ICD10$CTV3 <- NA
    dfCodesheet.ICD10$OPCS4 <- NA
    dfCodesheet.ICD10$n_20003 <- NA
    dfCodesheet.ICD10$BNF <- NA
    dfCodesheet.ICD10$DMD <- NA
    dfCodesheet.ICD10$source = "ICD10"
    dfCodesheet.ICD10$text <-
      stringr::str_replace_all(dfCodesheet.ICD10$text, "[^/[:^punct:]]", "")

    dfCodesheet.ICD10 <- data.table::data.table(dfCodesheet.ICD10)
    data.table::setkey(dfCodesheet.ICD10, "ICD10")

    ####### # ICD9 depscription, dfCodesheet.icd9_lkp
    message("ICD9")
    dfCodesheet.icd9_lkp <-
      as.data.frame(readxl::read_xlsx(fcoding.xls, sheet = "icd9_lkp")) # certainly not complete!
    colnames(dfCodesheet.icd9_lkp) <- c("ICD9", "text")
    dfCodesheet.ICD9.coding87 <-
      data.frame(data.table::fread(ficd9))[, 1:2] # UKB
    names(dfCodesheet.ICD9.coding87) <- c("ICD9", "text")
    dfCodesheet.icd9_lkp <-
      rbind(dfCodesheet.icd9_lkp, dfCodesheet.ICD9.coding87[!dfCodesheet.ICD9.coding87$ICD9 %in% dfCodesheet.icd9_lkp$ICD9, ])

    dfCodesheet.ICD9 <- dfCodesheet.icd9_lkp

    dfCodesheet.ICD9$ICD10 <- NA
    dfCodesheet.ICD9$READ <- NA
    dfCodesheet.ICD9$CTV3 <- NA
    dfCodesheet.ICD9$OPCS4 <- NA
    dfCodesheet.ICD9$n_20003 <- NA
    dfCodesheet.ICD9$BNF <- NA
    dfCodesheet.ICD9$DMD <- NA
    dfCodesheet.ICD9$source = "ICD9"
    dfCodesheet.ICD9$text <-
      stringr::str_replace_all(dfCodesheet.ICD9$text, "[^/[:^punct:]]", "")
    dfCodesheet.ICD9 <- data.table::data.table(dfCodesheet.ICD9)
    data.table::setkey(dfCodesheet.ICD9, "ICD9")

    ###### OPCS4
    message("OPSC4")
    dfCodesheet.OPCS4.coding240 <-
      data.frame(data.table::fread(fopcs4))[, c(1, 2)]
    names(dfCodesheet.OPCS4.coding240) <- c("OPCS4", "text")
    dfCodesheet.OPCS4 <- dfCodesheet.OPCS4.coding240
    dfCodesheet.OPCS4$ICD10 <- NA
    dfCodesheet.OPCS4$ICD9 <- NA
    dfCodesheet.OPCS4$READ <- NA
    dfCodesheet.OPCS4$CTV3 <- NA
    dfCodesheet.OPCS4$n_20003 <- NA
    dfCodesheet.OPCS4$BNF <- NA
    dfCodesheet.OPCS4$DMD <- NA
    dfCodesheet.OPCS4$source = "OPCS4"
    dfCodesheet.OPCS4$text <-
      stringr::str_replace_all(dfCodesheet.OPCS4$text, "[^/[:^punct:]]", "")
    dfCodesheet.OPCS4 <- data.table::data.table(dfCodesheet.OPCS4)
    data.table::setkey(dfCodesheet.OPCS4, "OPCS4")

    ### n_20003
    message("f.20003")
    dfCodesheet.n_20003 <-
      data.frame(data.table::fread(f20003))[, c(1, 2)]
    names(dfCodesheet.n_20003) <- c("n_20003", "text")
    dfCodesheet.n_20003$n_20003 <-
      as.character(dfCodesheet.n_20003$n_20003) ### LOOKUPS ASSUME CHARACTER, SO NEEDS TO BE CHARACTER
    dfCodesheet.n_20003$OPCS4 <- NA
    dfCodesheet.n_20003$ICD10 <- NA
    dfCodesheet.n_20003$ICD9 <- NA
    dfCodesheet.n_20003$READ <- NA
    dfCodesheet.n_20003$CTV3 <- NA
    dfCodesheet.n_20003$BNF <- NA
    dfCodesheet.n_20003$DMD <- NA
    dfCodesheet.n_20003$source = "n_20003"
    dfCodesheet.n_20003$text <-
      stringr::str_replace_all(dfCodesheet.n_20003$text, "[^/[:^punct:]]", "")
    dfCodesheet.n_20003 <- data.table::data.table(dfCodesheet.n_20003)
    data.table::setkey(dfCodesheet.n_20003, "n_20003")

    #####  FINAL OBJECT.
    LstdfCodesheets <- list(
      READ = dfCodesheet.READ[, c("READ",
                                  "CTV3",
                                  "ICD10",
                                  "ICD9",
                                  "OPCS4",
                                  "n_20003",
                                  "BNF",
                                  "DMD",
                                  "text",
                                  "source")],
      CTV3 = dfCodesheet.CTV3[, c("READ",
                                  "CTV3",
                                  "ICD10",
                                  "ICD9",
                                  "OPCS4",
                                  "n_20003",
                                  "BNF",
                                  "DMD",
                                  "text",
                                  "source")],
      ICD10 = dfCodesheet.ICD10[, c("READ",
                                    "CTV3",
                                    "ICD10",
                                    "ICD9",
                                    "OPCS4",
                                    "n_20003",
                                    "BNF",
                                    "DMD",
                                    "text",
                                    "source")],
      ICD9 = dfCodesheet.ICD9[, c("READ",
                                  "CTV3",
                                  "ICD10",
                                  "ICD9",
                                  "OPCS4",
                                  "n_20003",
                                  "BNF",
                                  "DMD",
                                  "text",
                                  "source")],
      OPCS4 = dfCodesheet.OPCS4[, c("READ",
                                    "CTV3",
                                    "ICD10",
                                    "ICD9",
                                    "OPCS4",
                                    "n_20003",
                                    "BNF",
                                    "DMD",
                                    "text",
                                    "source")],
      n_20003 = dfCodesheet.n_20003[, c("READ",
                                        "CTV3",
                                        "ICD10",
                                        "ICD9",
                                        "OPCS4",
                                        "n_20003",
                                        "BNF",
                                        "DMD",
                                        "text",
                                        "source")],
      BNF = dfCodesheet.BNF[, c("READ",
                                        "CTV3",
                                        "ICD10",
                                        "ICD9",
                                        "OPCS4",
                                        "n_20003",
                                        "BNF",
                                        "DMD",
                                        "text",
                                        "source")],
      DMD = dfCodesheet.DMD[, c("READ",
                                        "CTV3",
                                        "ICD10",
                                        "ICD9",
                                        "OPCS4",
                                        "n_20003",
                                        "BNF",
                                        "DMD",
                                        "text",
                                        "source")],
      ALL = rbind(
        dfCodesheet.DMD[, c("READ", "CTV3", "ICD10", "ICD9", "OPCS4", "n_20003", "BNF", "DMD","text", "source")],
        dfCodesheet.READ[, c("READ", "CTV3", "ICD10", "ICD9", "OPCS4", "n_20003", "BNF", "DMD","text", "source")],
        dfCodesheet.CTV3[, c("READ", "CTV3", "ICD10", "ICD9", "OPCS4", "n_20003", "BNF", "DMD","text", "source")],
        dfCodesheet.ICD10[, c("READ", "CTV3", "ICD10", "ICD9", "OPCS4", "n_20003", "BNF", "DMD","text", "source")],
        dfCodesheet.ICD9[, c("READ", "CTV3", "ICD10", "ICD9", "OPCS4", "n_20003", "BNF", "DMD","text", "source")],
        dfCodesheet.OPCS4[, c("READ", "CTV3", "ICD10", "ICD9", "OPCS4", "n_20003", "BNF", "DMD","text", "source")],
        dfCodesheet.BNF[, c("READ", "CTV3", "ICD10", "ICD9", "OPCS4", "n_20003", "BNF", "DMD","text", "source")]

      )
      ## CAN PROBABLY REMOVE ALL AND CREATE IT WHEN NEEDED.
    )

    return(LstdfCodesheets)
  }




########################################
##### FUNCTIONS to convert different codings.
########################################

convert.coding <- function(Str,
                           from.code = "READ.CODE",
                           to.code = "UKB.Coding",
                           lookuptable = dfCodesheetREAD_SR.Coding,
                           ignore.case = FALSE) {
  # Str<-"f3...,f36z.,f31"
  Str <- as.character(Str)
  Str = gsub(pattern = ".",
             replacement = "\\.",
             Str,
             fixed = T) ### INMPORTANT FOR READ CODEES!!!!!
  if (is.na(Str)) {
    return(NA)
  }
  VctStr <- unlist(strsplit(Str, ","))
  #VctRXstrings<-strsplit(df[!is.na(df$n_20003_),]$n_20003_,",")[[1]]
  c <- paste(unique(unlist(
    lapply(VctStr,  function(x)
      lookuptable[, get(to.code)] [grep(paste("^", x, sep = ""), lookuptable[, get(from.code)] , ignore.case =
                                          ignore.case)])
  )), collapse = ",")

  return(c)
}


add.description.to.codes <-
  function(Str,
           code.id = "UKB.Coding",
           description.id = "Meaning",
           description.lookuptable = dfCodesheetREAD_SR.Coding,
           ignore.case = FALSE,
           firstcodeonly = TRUE) {
    if (is.na(Str) | Str == "NA") {
      return(Str)
    }
    Str <- as.character(Str)
    if (is.na(Str)) {
      return(NA)
    }

    VctStr <- unlist(strsplit(Str, ","))



    c <- sapply(VctStr,  function(x) {
      x.d <-
        description.lookuptable[, get(description.id)] [grep(paste("^", x, sep =
                                                                     ""),
                                                             description.lookuptable[, get(code.id)] ,
                                                             ignore.case = ignore.case)]
      if (length(x.d) == 0) {
        return(paste0(x, " (NA)"))
      }
      x.d <-
        stringr::str_replace_all(x.d,  "[^/[:^punct:]]", "") # replace all symbols to not mess up downstream things.
      if (firstcodeonly == TRUE) {
        x.d <- x.d[1]
      }
      x.d <- paste0(x.d, collapse = " /")
      x.d <- paste0(x, " (", x.d, ")")
      x.d
    }, USE.NAMES = F)

    c <- paste(unique(unname(c)), collapse = ",")
    return(c)
  }

add.description.to.vectorofcodes <-
  function(vctcodes = c("G551.", "G55.."),
           code.id = "CTV3",
           description.id = "text",
           description.lookuptable = LstdfCodesheets$CTV3,
           ignore.case = FALSE,
           firstcodeonly = TRUE) {
    #c <- unlist(sapply(vctcodes,na.omit,USE.NAMES = F))

    if (length(vctcodes) == 0) {
      return(c())
    }

    c <- as.character(vctcodes)
    df_c <- description.lookuptable[.(c)] #[[code.id]]
    df_c <-
      df_c %>% dplyr::select(c = eval(code.id), text) %>% unique %>% dplyr::filter(!is.na(c)) %>% dplyr::group_by(c) %>% dplyr::mutate(text = paste0(text, collapse = "/")) %>% unique %>% dplyr::mutate(c_text =
                                                                                                                                                                                                           paste0(c, " (", text, ")"))
    df_c <- data.frame(df_c)
    #c <- paste(unique(unname(c)),collapse=",")
    return(df_c)
  }


expand_clean_codes <- function(col = df$ICD10,
                               from.code = "ALT_CODE",
                               description.id = 'DESCRIPTION',
                               lookuptable = dfCodesheet.icd10_lkp,
                               add_description = T) {
  # col=df$ICD10
  # lookuptable=LstdfCodesheets$ICD10
  # from.code="ICD10"
  # description.id='text'

  col <-
    ukbpheno::PreProcessDfDefinitions(
      df = data.frame(col = col, tmp = rep("NA", length(col))),
      VctAllColumns = c("col", "tmp"),
      VctColstoupper = NULL
    )[, 1]

  to.code = "self"
  lookuptable$self <- lookuptable[, get(from.code)]
  #c <- unlist(lapply(col, convert.coding,from.code=from.code,to.code=to.code,lookuptable=lookuptable))

  c1 <-
    paste(col, unlist(
      lapply(
        col,
        convert.coding,
        from.code = from.code,
        to.code = to.code,
        lookuptable = lookuptable
      )
    ), sep = ",")
  c2 <-
    unlist(lapply(c1, function(x) {
      x = unique(strsplit(x, ",")[[1]])
      if (length(x) == 1 &
          x[1] == "NA") {
        return("NA")
      } else{
        return(paste(x[x != "NA"], collapse = ","))
      }
    }))

  if (add_description == T) {
    c2 <- sapply(
      c2,
      add.description.to.codes,
      code.id = from.code,
      description.id = description.id,
      description.lookuptable = lookuptable,
      USE.NAMES = F
    )
  }
  return(c2)
}

lookup_list_in_df <- function(lst, df.lookup) {
  # lookup list in df, as long as names are corresponding to df headers..
  df.all <- data.frame()
  for (i in 1:length(lst)) {
    lookup = na.omit(lst[i][[1]])
    d <-
      df.lookup[df.lookup[, get(names(lst[i]))] %in% lookup , ] %>% unique()
    #if(nrow(d)>1) d$source = names(lst[i])
    df.all <- rbind(df.all, d)

  }
  df.all <- df.all %>% unique()
  return(df.all)
}

annotate_codes <- function(c, LstdfCodesheets = LstdfCodesheets) {
  if (length(c) == 0) {
    return(c())
  }
  lst_lookup_annotated <- list()
  for (i in 1:length(c)) {
    c_anno <-
      add.description.to.vectorofcodes(
        vctcodes = c[i][[1]],
        code.id = names(c[i]),
        description.id = "text",
        description.lookuptable = LstdfCodesheets[[names(c[i])]]
      )
    lst_lookup_annotated[[names(c[i])]] = c_anno
    #add.description.to.vectorofcodes(lst_lookup$ICD10,code.id = "ICD10",description.id = "text",description.lookuptable = LstdfCodesheets$ICD10)
  }
  return(lst_lookup_annotated)
}

# suggest new codes;
lookup_codes <-
  function(codes = row,
           LstdfCodesheets = LstdfCodesheets,
           expand_input = F) {
    # codes=row
    fromcodes = c("ICD10", "ICD9", "READ", "CTV3", "OPCS4", "n_20003","BNF","DMD")
    cols = c("ICD10", "ICD9", "READ", "CTV3", "OPCS4", "n_20003","BNF","DMD")
    # fromcodes = c("ICD10", "ICD9", "READ", "CTV3", "OPCS4", "n_20003")
    # cols = c("ICD10", "ICD9", "READ", "CTV3", "OPCS4", "n_20003")
    codes <-
      PreProcessDfDefinitions(df = codes,
                              VctAllColumns = cols ,
                              VctColstoupper = NULL)
    input = codes
    if (expand_input) {
      for (code in fromcodes) {
        #colname <- glue::glue("{code}CODES")
        if (is.null(codes[[code]])) {
          next
        }
        codes[[code]] <-
          expand_clean_codes(
            col = codes[[code]] ,
            from.code = code,
            description.id = 'text',
            lookuptable = LstdfCodesheets[code][[1]],
            add_description = F
          )
      }
    }

    c <-
      sapply(codes[cols], function(x)
        unique(unlist(strsplit(x, ","))))

    df_lookup <-
      lookup_list_in_df(lst = c, df.lookup = LstdfCodesheets$ALL)
    lst_lookup <-
      apply(df_lookup[, fromcodes, with = FALSE], 2, function(x)
        c(na.omit(unique(x))))
    #lst_lookup annotated:
    lst_lookup_anno <-
      annotate_codes(c = lst_lookup, LstdfCodesheets = LstdfCodesheets)
    df_lookup_anno <- df_lookup
    for (col in fromcodes) {
      #if(!is.data.frame(codes.lookup[['OPCS4']])){next}

      df_lookup_anno[[col]] <-
        lst_lookup_anno[[col]]$c_text[match(df_lookup[, get(col)] , lst_lookup_anno[[col]]$c)]

    }
    return(
      list(
        cols = cols,
        input = input,
        input_c = c,
        lst_lookup = lst_lookup,
        df_lookup = df_lookup,
        lst_lookup_anno = lst_lookup_anno,
        df_lookup_anno = df_lookup_anno
      )
    )
  }


# #### SHINY:
#
#
# library(shiny)
# library(DT)
# Define UI for app that draws a histogram ----


ui <- shiny::fluidPage(
  # App title ----
  shiny::titlePanel("UKB code explorer"),
  # Sidebar layout with input and output definitions ----
  shiny::sidebarLayout(
    # Sidebar panel for inputs ----
    shiny::sidebarPanel(
      # Input: Slider for the number of bins ----
      shiny::textAreaInput(
        inputId = "iICD10",
        label = "ICD10",
        value = "Z951 (cabg)",
        width = NULL,
        placeholder = NULL
      ),
      shiny::textAreaInput(
        inputId = "iICD9",
        label = "ICD9",
        value = "",
        width = NULL,
        placeholder = NULL
      ),
      shiny::textAreaInput(
        inputId = "iREAD",
        label = "READ",
        value = "",
        width = NULL,
        placeholder = NULL
      ),
      shiny::textAreaInput(
        inputId = "iCTV3",
        label = "CTV3",
        value = "",
        width = NULL,
        placeholder = NULL
      ),
      # shiny::textAreaInput(
      #   inputId = "iBNF",
      #   label = "BNF",
      #   value = "",
      #   width = NULL,
      #   placeholder = NULL
      # ),
      # shiny::textAreaInput(
      #   inputId = "iDMD",
      #   label = "DMD",
      #   value = "",
      #   width = NULL,
      #   placeholder = NULL
      # ),
      shiny::textAreaInput(
        inputId = "iOPCS4",
        label = "OPCS4",
        value = "K40,K41,K43,K44,K45,K46(cabg),K471(endarterectomy),K49, K50,K75 (pci)",
        width = NULL,
        placeholder = NULL
      ),
      shiny::checkboxInput(inputId = "iExpandcodes", "Expand codes, e.g. I50 -> I501,I502, etc.  ", FALSE),
      shiny::actionButton("goButton", "Go!"),
      shiny::HTML(
        "<br><br>Note; this is a tryout version for exploration - translations are not reliable and should be manually verified! <br><br>BNF/DMD codes: contains information from NHS Digital licenced under the current version of Open Government Licence <br>[https://www.nationalarchives.gov.uk/doc/open-government-licence/open-government-licence.htm]"
      )
    ),
    # Main panel for displaying outputs ----
    shiny::mainPanel(
      #verbatimTextOutput("oICD10")
      #DT::dataTableOutput("table_input")

      shiny::tabsetPanel(
        shiny::tabPanel("ICD10 ", DT::dataTableOutput("table_oICD10")),
        shiny::tabPanel("ICD9",  DT::dataTableOutput("table_oICD9")),
        shiny::tabPanel("READ",  DT::dataTableOutput("table_oREAD")),
        shiny::tabPanel("CTV3",  DT::dataTableOutput("table_oCTV3")),
        shiny::tabPanel("OPCS4",  DT::dataTableOutput("table_oOPCS4")),
        shiny::tabPanel("n_20003",  DT::dataTableOutput("table_on_20003")),
        # shiny::tabPanel("BNF",  DT::dataTableOutput("table_oBNF")),
        # shiny::tabPanel("DMD",  DT::dataTableOutput("table_oDMD")),
        shiny::tabPanel(
          "output",
          shiny::checkboxInput(inputId = "iIncludetext", "include description", FALSE),
          shiny::h4("ICD10"),
          shiny::textOutput("codes_oICD10"),
          shiny::h4("ICD9"),
          shiny::textOutput("codes_oICD9"),
          shiny::h4("READ"),
          shiny::textOutput("codes_oREAD"),
          shiny::h4("CTV3"),
          shiny::textOutput("codes_oCTV3"),
          shiny::h4("OPCS4"),
          shiny::textOutput("codes_oOPCS4"),
          shiny::h4("n_20003"),
          shiny::textOutput("codes_on_20003")
          # shiny::h4("BNF"),
          # shiny::textOutput("codes_oBNF"),
          # shiny::h4("DMD"),
          # shiny::textOutput("codes_oDMD")
        ),
        shiny::tabPanel("lookuptable",  DT::dataTableOutput("table_oraw"))



      )
    )
  )
)


server <- function(input, output) {
  values_lookup <-
    shiny::reactiveValues(codes.lookup.shinyready = NULL)

  shiny::observeEvent(input$goButton, {
    shiny::showModal(modalDialog("please wait" , easyClose = FALSE, footer =
                                   NULL))
    row <- data.frame(
      ICD10 = gsub("\\|", ",", input$iICD10),
      ICD9 = gsub("\\|", ",", input$iICD9),
      READ = gsub("\\|", ",", input$iREAD),
      CTV3 = gsub("\\|", ",", input$iCTV3),
      # BNF = gsub("\\|", ",", input$iBNF),
      # DMD = gsub("\\|", ",", input$iDMD),
      OPCS4 = gsub("\\|", ",", input$iOPCS4)
    )


    codes.lookup <-
      lookup_codes(
        codes = row,
        LstdfCodesheets = LstdfCodesheets,
        expand_input = input$iExpandcodes
      )


    convert_lookup_to_df <- function(codes, input_c = row_exp$ICD10) {
      df <-
        data.frame(
          codes = codes$c,
          text = codes$text,
          c_text = codes$c_text,
          input = codes$c %in% input_c
        )

      if (nrow(df) == 0) {
        df <- data.frame(
          codes = "",
          text = "",
          c_text = "",
          input = FALSE
        )
      }
      return(df)
    }
    values_lookup$codes.lookup.shinyready <<-
      lapply(codes.lookup$cols, function(x)
        convert_lookup_to_df(
          codes = codes.lookup$lst_lookup_anno[[x]],
          input_c = codes.lookup$input_c[[x]]
        ))
    names(values_lookup$codes.lookup.shinyready) <<-
      codes.lookup$cols

    dtoptions = list(
      pageLength = 25,
      info = FALSE,
      lengthMenu = list(c(25, 50, 100, 200,-1), c("25", "50", "100", "200", "All"))
    )

    output$table_oICD10 = DT::renderDataTable({
      values_lookup$codes.lookup.shinyready$ICD10[, c(1, 2, 4)]
    }, filter = "top", options = dtoptions,
    selection = list(
      mode = 'multiple',
      selected = which(values_lookup$codes.lookup.shinyready$ICD10$input)
    ))

    output$table_oICD9 = DT::renderDataTable({
      values_lookup$codes.lookup.shinyready$ICD9[, c(1, 2, 4)]
    }, filter = "top", options = dtoptions,
    selection = list(
      mode = 'multiple',
      selected = which(values_lookup$codes.lookup.shinyready$ICD9$input)
    ))


    output$table_oREAD = DT::renderDataTable({
      values_lookup$codes.lookup.shinyready$READ[, c(1, 2, 4)]
    }, filter = "top", options = dtoptions,
    selection = list(
      mode = 'multiple',
      selected = which(values_lookup$codes.lookup.shinyready$READ$input)
    ))

    output$table_oCTV3 = DT::renderDataTable({
      values_lookup$codes.lookup.shinyready$CTV3[, c(1, 2, 4)]
    }, filter = "top", options = dtoptions,
    selection = list(
      mode = 'multiple',
      selected = which(values_lookup$codes.lookup.shinyready$CTV3$input)
    ))

    # output$table_oBNF = DT::renderDataTable({
    #   values_lookup$codes.lookup.shinyready$BNF[, c(1, 2, 4)]
    # }, filter = "top", options = dtoptions,
    # selection = list(
    #   mode = 'multiple',
    #   selected = which(values_lookup$codes.lookup.shinyready$BNF$input)
    # ))
    #
    # output$table_oDMD = DT::renderDataTable({
    #   values_lookup$codes.lookup.shinyready$DMD[, c(1, 2, 4)]
    # }, filter = "top", options = dtoptions,
    # selection = list(
    #   mode = 'multiple',
    #   selected = which(values_lookup$codes.lookup.shinyready$CTV3$input)
    # ))

    output$table_oOPCS4 = DT::renderDataTable({
      values_lookup$codes.lookup.shinyready$OPCS4[, c(1, 2, 4)]
    }, filter = "top", options = dtoptions,
    selection = list(
      mode = 'multiple',
      selected = which(values_lookup$codes.lookup.shinyready$OPCS4$input)
    ))

    output$table_on_20003 = DT::renderDataTable({
      values_lookup$codes.lookup.shinyready$n_20003[, c(1, 2, 4)]
    }, filter = "top", options = dtoptions,
    selection = list(
      mode = 'multiple',
      selected = which(values_lookup$codes.lookup.shinyready$n_20003$input)
    ))

    # data.frame(x="123",t="asd")
    removeModal()
  })

  output$codes_oICD10 <- shiny::renderPrint({
    s = input$table_oICD10_rows_selected
    if (input$iIncludetext) {
      paste(values_lookup$codes.lookup.shinyready$ICD10[s, ]$c_text,
            collapse = ", ")
    } else{
      paste(values_lookup$codes.lookup.shinyready$ICD10[s, ]$codes,
            collapse = ", ")
    }

  })
  output$codes_oICD9 <- shiny::renderPrint({
    s = input$table_oICD9_rows_selected
    if (input$iIncludetext) {
      paste(values_lookup$codes.lookup.shinyready$ICD9[s, ]$c_text,
            collapse = ", ")
    } else {
      paste(values_lookup$codes.lookup.shinyready$ICD9[s, ]$codes,
            collapse = ", ")
    }

  })
  output$codes_oREAD <- shiny::renderPrint({
    s = input$table_oREAD_rows_selected
    if (input$iIncludetext) {
      paste(values_lookup$codes.lookup.shinyready$READ[s, ]$c_text,
            collapse = ", ")
    } else{
      paste(values_lookup$codes.lookup.shinyready$READ[s, ]$codes,
            collapse = ", ")
    }
  })
  output$codes_oCTV3 <- shiny::renderPrint({
    s = input$table_oCTV3_rows_selected
    if (input$iIncludetext) {
      paste(values_lookup$codes.lookup.shinyready$CTV3[s, ]$c_text,
            collapse = ", ")
    } else {
      paste(values_lookup$codes.lookup.shinyready$CTV3[s, ]$codes,
            collapse = ", ")
    }
  })
  output$codes_oOPCS4 <- shiny::renderPrint({
    s = input$table_oOPCS4_rows_selected
    if (input$iIncludetext) {
      paste(values_lookup$codes.lookup.shinyready$OPCS4[s, ]$c_text,
            collapse = ", ")
    } else {
      paste(values_lookup$codes.lookup.shinyready$OPCS4[s, ]$codes,
            collapse = ", ")
    }
  })
  output$codes_on_20003 <- shiny::renderPrint({
    s = input$table_on_20003_rows_selected
    if (input$iIncludetext) {
      paste(values_lookup$codes.lookup.shinyready$n_20003[s, ]$c_text,
            collapse = ", ")
    } else {
      paste(values_lookup$codes.lookup.shinyready$n_20003[s, ]$codes,
            collapse = ", ")
    }
  })

  # output$codes_oBNF <- shiny::renderPrint({
  #   s = input$table_oBNF_rows_selected
  #   if (input$iIncludetext) {
  #     paste(values_lookup$codes.lookup.shinyready$BNF[s, ]$c_text,
  #           collapse = ", ")
  #   } else {
  #     paste(values_lookup$codes.lookup.shinyready$BNF[s, ]$codes,
  #           collapse = ", ")
  #   }
  # })
  # output$codes_oDMD <- shiny::renderPrint({
  #   s = input$table_oDMD_rows_selected
  #   if (input$iIncludetext) {
  #     paste(values_lookup$codes.lookup.shinyready$BNF[s, ]$c_text,
  #           collapse = ", ")
  #   } else {
  #     paste(values_lookup$codes.lookup.shinyready$BNF[s, ]$codes,
  #           collapse = ", ")
  #   }
  # })

  #
  output$table_oraw = DT::renderDataTable({
    LstdfCodesheets$ALL
  }, filter = "top")
  #output$oICD10 <- renderText({ })
}


# # # codes.lookup$lst_lookup_anno$ICD10$c_text[match(codes.lookup$df_lookup[,get("ICD10")] ,codes.lookup$lst_lookup_anno[["ICD10"]]$c)]
# # #######################################
# # #### LOAD DATA.
# # #######################################
# # source("ProcessdfDefinitions.R")
# fcoding_xls="../inst/extdata/all_lkps_maps_v3.xlsx"
# if(!exists("LstdfCodesheets")){ LstdfCodesheets <- load_data(fcoding.xls=opt$fcoding_xls,fmed.readSR = opt$f_med_readSR) }
# #
# gc()
# #######################################
# #### STAND ALONE EXAMPLE FOR 1 ROW.
# #######################################
# # ICD10, icd9, read2,ct3,
# dfDefinitions_file="../inst/extdata/definitions_cardiometabolic_traits.tsv"
# df = data.frame(fread(dfDefinitions_file))
# df <- ProcessDfDefinitions(df=df,fill_dependencies = F)
# df <- df[,c("TRAIT","DESCRIPTION", "ICD10","ICD9","READ2","CTV3","OPCS4")]
# names(df) <- c("TRAIT","DESCRIPTION", "ICD10","ICD9","READ","CTV3","OPCS4") # READ2 > READ
# # # DO IIT FOR 1 ROW suggest codes for 1 selected row.
# irow=13 #12 #3#12
# row <- df[irow,]
# #row$OPCS4 <- "K02"
# codes.lookup <- lookup_codes(codes = row,LstdfCodesheets=LstdfCodesheets,expand_input=T)
#
# codes.lookup$df_lookup

# #######################################
# #### LOAD DATA.
# #######################################
# fcoding_xls="../inst/extdata/all_lkps_maps_v3.xlsx"
print("Prepare the code maps...")
if (!exists("LstdfCodesheets")) {
  LstdfCodesheets <-
    load_data(
      fcoding.xls = opt$fcoding_xls,
      fmed.readSR = opt$f_med_readSR,
      ficd10 = opt$fcoding_icd10,
      ficd9 = opt$fcoding_icd9,
      fopcs4 = opt$fcoding_opcs4,
      f20003 = opt$fcoding_20003,fread2_bnf =opt$fbnf,fread2_dmd = opt$fdmd
    )
}
#
gc()

shiny::shinyApp(ui, server)

