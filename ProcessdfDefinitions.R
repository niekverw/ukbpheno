#' @export
ConvertFactorsToStringReplaceNAInDf<-function (df) {
  df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE) ## CHANGE Factors to strings; everything is now a string.
  df[df==""]  <- NA ### REPLACE EMPTY WITH NA.
  return(df)
}

#' @export
pasteRemoveNA <- function(..., sep = " ", collapse = NULL, na.rm = F) {
  if (na.rm == F)
    paste(..., sep = sep, collapse = collapse)
  else
    if (na.rm == T) {
      paste.na <- function(x, sep) {
        x <- gsub("^\\s+|\\s+$", "", x)
        ret <- paste(na.omit(x), collapse = sep)
        is.na(ret) <- ret == ""
        return(ret)
      }
      df <- data.frame(..., stringsAsFactors = F)
      ret <- apply(df, 1, FUN = function(x) paste.na(x, sep))

      if (is.null(collapse))
        ret
      else {
        paste.na(ret, sep = collapse)
      }
    }
}

#' @export
CheckDuplicateTRAITS<-function(df){
  if(length(unique(duplicated(df["TRAIT"])))>1){stop("TRAIT column contains duplicate ID's")}
}



#' @export
PreProcessDfDefinitions<-function(df,VctAllColumns,VctColstoupper=NULL ){ # c("ICD10","ICD9","OPCS4","OPCS3")
## df<-dfDefinitions
  # check if nrows==1
  checkr=0
  checkc=0
  if(nrow(df)==1){df<-rbind(df,df);checkr=1}
  if(ncol(df)==1 & length(VctAllColumns)==1){df<-cbind(df,df);checkc=1; names(df)<-c(VctAllColumns,"tmp"); VctAllColumns=c(VctAllColumns,"tmp") }

  ## for the names: remove everything between dots (R converts symbols to dots "(,.-)/" etc )
  names(df) <- gsub( " *\\..*?\\. *", "", names(df) )
  ## add missing columns
  df[, VctAllColumns[!VctAllColumns %in% colnames(df)]] <- NA
  ## remove everything between brackets
  df[,VctAllColumns]<- data.frame(apply(df[,VctAllColumns], 2, function(y) gsub( " *\\(.*?\\) *", "", y)) )
  ### remove dots(.): names(df)
  VctColumnsRemoveDots=c("ICD10","ICD9","OPCS4","OPCS3")
  VctColumnsRemoveDots <-VctColumnsRemoveDots[VctColumnsRemoveDots %in% colnames(df)]
  df[,VctColumnsRemoveDots]<- data.frame(apply(df[,VctColumnsRemoveDots],2,function(x) gsub(".", "", x, fixed = TRUE)))
  ### remove spaces:
  df[,VctAllColumns]<- data.frame(apply(df[,VctAllColumns],2,function(x) gsub(" ", "", x, fixed = TRUE)))
  ### remove trailing commas:
  trim.commas <- function (x) gsub("(?<=[\\,])\\,*|^\\,+|\\,+$", "", x, perl=TRUE)
  df[,VctAllColumns]<- data.frame(apply(df[,VctAllColumns],2,function(x) trim.commas(x)))

  ### replace one equal sign to logical equal if needed
  VctCustomFields="TS"
  df[,VctCustomFields]<-gsub("\\b[=]+\\b","==",df[,VctCustomFields],perl=TRUE)
  df[,VctCustomFields]<-gsub("\\b[≥]\\b",">=",df[,VctCustomFields],perl=TRUE)
  df[,VctCustomFields]<-gsub("\\b[≤]\\b","<=",df[,VctCustomFields],perl=TRUE)

  if(length(VctColstoupper)==1){
    df[,VctColstoupper] <- toupper(df[,VctColstoupper])
  }
  if(length(VctColstoupper)>1){
    df[,VctColstoupper] <- apply(df[,VctColstoupper],2,toupper)
  }

  
  df<-ConvertFactorsToStringReplaceNAInDf(df) #### CONVERT FACTOR TO STRING

  if(checkr==1){df<-df[1,]}
  if(checkc==1){df<-df[,1]}

  return(df)
}

#' @export
FillInSRdefinitions<-function(df,Var="SR",cols=c("n_20001","n_20002","n_20004") ) {
  # 20001 cancer code, self reported; 20002 non-cancer illness, self reported; 20004 operation code
  ## fill in SR
  df[,Var]<-as.character(df[,Var])
  try(df[is.na( as.character(df[,Var])) ,][,Var] <- "",silent = T)

 # df [ is.na(as.character( df[,Var] )) %in% "NA"  ,]

  cols<- cols[cols %in% colnames(df)] ### check if available.
  for(col in cols){
    #print(col)
    Columnmatches<-gsub ( ",","|",paste(col,unlist(df[col]),sep="=" ))
    Columnmatches [grepl("NA",Columnmatches)]<-"" ### remove NAs..
    df[,Var]<-  paste(df[,Var],Columnmatches ,sep="," )
  }
  print(df[,Var])
  ### trim commas:
  trim.commas <- function (x) gsub("(?<=[\\,])\\,*|^\\,+|\\,+$", "", x, perl=TRUE)
  df[,Var]<-trim.commas(df[,Var])
  print(df[,Var])

  df<-ConvertFactorsToStringReplaceNAInDf(df)

  return(df)
}

#' @export
CovertMednamesToUkbcoding<- function(StrRx){
  #StrRx<-"phenformin,metformin,buformin,glibenclamide,chlorpropamide,tolbutamide,glibornuride,tolazamide,carbutamide,glipizide,gliquidone,gliclazide,metahexamide,glisoxepide,glimepiride,acetohexamide,glymidine,acarbose,miglitol,voglibose,troglitazone,rosiglitazone,pioglitazone,sitagliptin,vildagliptin,saxagliptin,alogliptin,linagliptin,gemigliptin,repaglinide,nateglinide,exenatide,pramlintide,benfluorex,liraglutide,mitiglinide,dapagliflozin,lixisenatide,canagliflozin,empagliflozin,albiglutide,dulaglutide"
  StrRx<-as.character(StrRx)
  if(is.na(StrRx)) { return(NA)}
  VctRXstrings<-unlist(strsplit(StrRx,","))
  #VctRXstrings<-strsplit(df[!is.na(df$n_20003),]$n_20003,",")[[1]]
  # StrRxCodes<-paste(unique(unlist(lapply(VctRXstrings,  function(x) dfCodesheetREAD_SR.Coding[,"UKB.Coding"] [ grep(x,dfCodesheetREAD_SR.Coding[,"Meaning"] ,ignore.case=TRUE )]  ))),collapse=",")
  # StrRxCodes<-paste(unique(unlist(lapply(VctRXstrings,  function(x) dfCodesheetREAD_SR.Coding[grep(x, Meaning,ignore.case=TRUE)][,"UKB.Coding"]    ))),collapse = ",")
  # StrRxCodes <- paste(unique(unlist(lapply(VctRXstrings,   function(x) dfCodesheetREAD_SR.Coding[match(x, dfCodesheetREAD_SR.Coding$Meaning), "UKB.Coding"]  ))),collapse = ",")
  # for each medication, only attempt to match string that cannot be converted to numeric, numeric assumed to be in ukb coding already
  StrRxCodes <- paste(unique(unlist(lapply(VctRXstrings,function(x) { if(is.na(as.numeric(x))) dfCodesheetREAD_SR.Coding[match(x, dfCodesheetREAD_SR.Coding$Meaning), "UKB.Coding"]  else x }))),collapse = ",")
  # check if codes exist
  print(StrRxCodes)
  StrRxCodes_c<-unlist(strsplit(StrRxCodes,","))
  StrRxCodes_inCode<- StrRxCodes_c[StrRxCodes_c %in% dfCodesheetREAD_SR.Coding$UKB.Coding  ]
  StrRxCodes_notInCode<- StrRxCodes_c[! StrRxCodes_c %in% dfCodesheetREAD_SR.Coding$UKB.Coding  ]
  print("The following codes not found in UKB coding:")
  print(StrRxCodes_notInCode)
  StrRxCodes<-paste(StrRxCodes_inCode,collapse = ",")
  
  # inCodeBool<-unlist(lapply(unlist(strsplit(StrRxCodes,",")),function(x) x %in% dfCodesheetREAD_SR.Coding$UKB.Coding)))
  return(StrRxCodes)
}

#' @export
CovertReadcodesToSelfReportedUkbCoding<- function(StrRx){
 # StrRx<-"f3,f4,ft"
  StrRx<-as.character(StrRx)
  if(is.na(StrRx)) { return(NA)}
  VctRXstrings<-unlist(strsplit(StrRx,","))
  #VctRXstrings<-strsplit(df[!is.na(df$n_20003),]$n_20003,",")[[1]]
  StrRxCodes<-paste(unique(unlist(lapply(VctRXstrings,  function(x) dfCodesheetREAD_SR.Coding[,"UKB.Coding"] [ grep(paste("^", x,sep=""),dfCodesheetREAD_SR.Coding[,"READ.CODE"] ,ignore.case=TRUE )]  ))),collapse=",")
  return(StrRxCodes)
}


#### TODO:
#' @export
ReduceRedundancyDf<- function(df){ ### NOT really nessesary
  return(df)
}



# DfDefinitions<-read.table("/Users/niekverw/Downloads/ex",sep="\t",header=T)
#columns<-c("ICD10","ICD9","OPCS4","OPCS3","TOUCHSCREEN","TS_AGE_DIAG_COLNAME","SELFREPORTED","MEDICATION","LAB")
# print(dfDefinitions)
#dfDefinitionstmp2<-ProcessDfDefinitions(dfDefinitions,columns)

#' ProcessDfDefinitions
#'
#' Process definitions, for input
#'
#' @param dfDefinitions df
#' @param VctAllColumns Vct
#' @keywords ExtractVarsFromMasterSet CreateUKBiobankPhentoypes ProcessDfDefinitions
#' @return None
#'
#' @examples
#' #
#' #This function processes an excel file with definitions and is automtically performed in CreateUKBiobankPhentoypes().
#' #It can be usefull to run this function as a check prior to running CreateUKBiobankPhentoypes.
#' #
#' #VctAllColumns contains all column names of interest, so that it can ignore everything else.
#' #20001, 20002 and 20004 go into SR
#' #READ and 20003 is parsed into RX
#' #
#' #
#' #
#' VctAllColumns<-  c("TS", "SR", "TS_RX", "SR_RX", "LAB", "ICD10", "ICD9", "OPCS4","OPCS3", "TS_AGE_DIAG_COLNAME", "READ","CTV3","BNF","DMD", "n_20001",    "n_20002", "n_20003", "n_20004", "DEPENDENCY")
#' ProcessDfDefinitions(dfDefinitions,VctAllColumns)
#'
#' @export
ProcessDfDefinitions<-function(df,
                               VctAllColumns=c("TS",
                                               "ICD10", "ICD9", "OPCS4","OPCS3",
                                               "READ", "CTV3",
                                               "BNF","DMD",
                                               "n_20001",    "n_20002", "n_20003", "n_20004",
                                               "DEPENDENCY"),
                               VctColstoupper=c("ICD10","ICD9","OPCS4","OPCS3"),
                               fill_dependencies=T){
  
  #df<- dfDefinitions  #  df<- dfDefinitions2
  # VctAllColumns<-  c("TS", "SR", "TS_RX", "SR_RX", "LAB", "ICD10", "ICD9", "OPCS4","OPCS3", "TS_AGE_DIAG_COLNAME", "READ","CTV3", "n_20001",    "n_20002", "n_20003", "n_20004", "DEPENDENCY")

  #if(nrow(df)==1 ) {stop("please have more than 1 phenotype definition.")} ## check if excel file has more than 1 row.
  df <- data.frame(df)
  names(df) <- sub(pattern = "CODES",replacement = "",names(df) )
    names(df) <- sub(pattern = "_$",replacement = "",names(df) ) # "n_20002_" --> "n_20002"

  df <- PreProcessDfDefinitions(df,VctAllColumns,VctColstoupper=VctColstoupper)
  if(any(!VctAllColumns %in% names(df))) print(paste("WARNING missing columns:", paste(VctAllColumns[!VctAllColumns %in% names(df)],collapse=", ")))

  #################################
  CheckDuplicateTRAITS(df) # check duplicateids.
  ### df = excel matrix. 1 hij loopt een voor een over elke rij heen,
  # 2 zoekt per dependency in die rij de bijpassende rijen voor elke dependency  erbij en plakt die naast elkaar (inclusief de dependencies van de dependencies).
  # 4) dan delete hij de depenencies die hij ingevuld heeft # dit op repeat tot dat er geen dependencies meer zijn en alles is ingevuld.
  if(fill_dependencies==F){return(return(ConvertFactorsToStringReplaceNAInDf(df)))}
  #################################
  ### HELPER FUNCTION TO CROSS CHECK EVERYTHING AND LOOKUPS, SHOULD GET A SEPERATE FUNCTION OUTSIDE OF EVERYTHING.
  #################################
  # ###[unsupported] LOOKUP NAMES OF MEDICATION and put UKBIO. in RX
  # print(df$n_20003)
  # df$n_20003<-unlist(lapply( df$n_20003, CovertMednamesToUkbcoding))
  # # df$n_20003<- paste(df$n_20003, unlist(lapply( df$n_20003, CovertMednamesToUkbcoding)),sep=",")
  # print(df$n_20003)
  
  # df<-FillInSRdefinitions(df,"SR_RX",c("n_20003"))
  ### LOOKUP READ. and put UKBIO. in SR_RX

  #################################
  ### FILL SR fields with  _2000X_ 'helper' columns;
  #################################
  #df<-FillInSRdefinitions(df,"SR",c("n_20001","n_20002","n_20004"))
  #df<-FillInSRdefinitions(df,"SR_RX",c("n_20003"))
  #df[,c("n_20001","n_20002","n_20004","n_20003")] <- NA
  #VctAllColumns <- VctAllColumns[!VctAllColumns %in% c("n_20001","n_20002","n_20004","n_20003")]

  ### lookup ICD10/9/OPER and put into READ and CTV3:
  # ....? I can lookup everything in everything to make everything more complete


  repeat {
    for(i in 1:nrow(df)) {
      row <- df[i,]

      if(!is.na(row$DEPENDENCY)){
        VctDEPENDENCYs<-unlist(strsplit(row$DEPENDENCY,","))

        for (StrDEPENDENCY in VctDEPENDENCYs) {
          if(row$TRAIT == StrDEPENDENCY) {stop("Dependency is same as trait.")}
          targetrow<-df[df$TRAIT==StrDEPENDENCY,]
          if(nrow(targetrow)==0){stop(paste('Dependency: "',StrDEPENDENCY,'" not found in traits ',row$TRAIT,sep="")) }

          if( is.na(targetrow["DEPENDENCY"])){
            for(col in VctAllColumns){
              Vctcol<-unique( unlist(strsplit( c(df[i,col],df[df$TRAIT==StrDEPENDENCY,col]) ,",")) )

              df[i,col]<-pasteRemoveNA(Vctcol ,collapse=",",na.rm=T)
            }
            # remove DEPENDENCY that was just filled in:
            #df[i,"DEPENDENCY"]<-gsub(paste(StrDEPENDENCY,sep=""),"",df[i,"DEPENDENCY"],fixed=TRUE,ignore.case=FALSE)

            # remove DEPENDENCY that was just filled in:
            LstTmpDependencies<- unlist(strsplit(VctDEPENDENCYs,","))
            df[i,"DEPENDENCY"]<-paste( LstTmpDependencies [! LstTmpDependencies  %in%  StrDEPENDENCY] ,sep="",collapse = ",")

            df[i,"DEPENDENCY"]<-gsub("^,*|(?<=,),|,*$", "", df[i,"DEPENDENCY"], perl=T)
            if(df[i,"DEPENDENCY"]==""){df[i,"DEPENDENCY"]<-NA} ## if empty replace with NA
          }
        }
      }
    }
    if( length(unique(is.na(df$DEPENDENCY)))==1 ) break
  }

  return(ConvertFactorsToStringReplaceNAInDf(df))

  #write.table(df,paste(dfDefinitions_file,".processed.tsv",sep=""),quote = FALSE,row.names = FALSE,sep="\t")
}



#' @export
get_allvarnames <- function(dfDefinitions_processed){
  #  dfDefinitions_processed
  VctAllUKBVDefinitionColumns=c("TS") #set this variable to a selection of columns (dfDefinition columns) to be outputted by the _UKBV variable, default is 'VctAllUKBVDefinitionColumns=c("TS","SR","TS_RX","SR_RX","LAB")'
  # TS(Touchscreen) in definition
  TScolumns = "TS"
  defcols <- unlist(strsplit(na.omit(unname(unlist(dfDefinitions_processed[,c(TScolumns)]))),split=","))
  defcols <- unlist(strsplit(chartr("[]", "()", defcols),"[()]"))
  defcols <- gsub("(≥|=|>|<).*","",defcols)
  defcols <- gsub("≥.*","",defcols)
  defcols <- unique(defcols) ## need column classes....


  all_ukb_fields <- unique(c(defcols,default_ukb_fields()))
  all_ukb_fields<-gsub("[a-zA-Z]*?_","",all_ukb_fields)
  nondefault_ukb_fields <- all_ukb_fields[!all_ukb_fields %in% default_ukb_fields() ] 
  
  
  return(
    list(all_ukb_fields=all_ukb_fields,
         nondefault_ukb_fields=nondefault_ukb_fields,
         default_ukb_fields=default_ukb_fields())
         )
}


#' @export
convert_readv2_to_ukbmedication<-function(Vctn_20003,Vctreadcodes){
  # Vctn_20003 <- dfDefinitions$n_20003
  # Vctreadcodes <- dfDefinitions$READ

  Vctreadcodes=PreProcessDfDefinitions(data.frame(Vctreadcodes=Vctreadcodes),VctAllColumns="Vctreadcodes",VctColstoupper=F)
  Vctn_20003=PreProcessDfDefinitions(data.frame(Vctn_20003=Vctn_20003),VctAllColumns="Vctn_20003",VctColstoupper=F)

  Vctn_20003 <- paste(Vctn_20003, unlist(lapply( Vctreadcodes, CovertReadcodesToSelfReportedUkbCoding)),sep=",")
  Vctn_20003 <- unlist(lapply(Vctn_20003,function(x) {  x = unique(strsplit(x,"," )[[1]]); if(length(x)==1 & x[1] =="NA"){ return("NA")} else{ return( paste(x[x != "NA"],collapse=",") )} }))
  return(Vctn_20003)
}


#' @export
## expand definition codes based on the ccodes that are available in the raw data (using lst.counts).
## lst.data.settings includes settings that should be used, e.g. if ICD10 should be looked up case sensitive or not (incase of READ cases are important, dotts should also NOT be interpreted)
expand_dfDefinitions_processed <- function(dfDefinitions_processed,lst.data.settings,lst.counts){
  classifications <- lst.data.settings %>% filter(expand_codes==1) %>% pull (classification) %>% unique()
  for (c in classifications){
    
  for (r in 1:nrow(dfDefinitions_processed)){ # for loops just as fast as apply in this case.. 
    VctStr = unlist(strsplit(dfDefinitions_processed[r,c],",")) 
    lookuptable = lst.counts[[c]]
    ignore.case = unique(lst.data.settings[lst.data.settings$classification %in% c,'ignore.case'])[1]
    
    Str_expanded <- paste(unique(unlist(
      lapply(VctStr,  function(x)  lookuptable$code [ grep(paste("^", x,sep=""),lookuptable$code ,ignore.case=ignore.case )])
      )),collapse=",")
    dfDefinitions_processed[r,c]<- Str_expanded
    }
  }
  return(dfDefinitions_processed)
}