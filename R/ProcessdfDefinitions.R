




#' Change factors to strings in a dataframe and assign NA to empty cells
#'
#' This function takes data fields containing illness codes/time of diagnosis distributed in multiple with each row representing one individual. It processed the data and return in episodes of event for all individuals
#' @param df dataframe containing the fields
#' @return  a data.table object with all episodes
#' @keywords auxiliary
#' @export
#' @examples
#' ConvertFactorsToStringReplaceNAInDf(data.frame("a" = c(1,2,3),"p" = c("t","b",""),"y" = c(NA,"w","r")))
ConvertFactorsToStringReplaceNAInDf<-function (df) {
  df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE) ## CHANGE Factors to strings; everything is now a string.
  df[df==""]  <- NA ### REPLACE EMPTY WITH NA.
  return(df)
}


#' Custom paste function with formatting step to remove redundant spaces
#'
#' This function is used for formatting the character vectors
#' @return  formatted object
#' @keywords auxiliary
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

#' Check if there is duplicates entries of phenotypes/traits in definition tables
#'
#' This function verify traits are uniquely defined and stop execution otherwise.
#' @keywords auxiliary
#' @export
CheckDuplicateTRAITS<-function(df){
  if(length(unique(duplicated(df["TRAIT"])))>1){
    print(df[duplicated(df["TRAIT"]),])
    stop("TRAIT column contains duplicate ID's")}
}


#' Preprocess the definition table
#'
#' This function preprocess definition table - formatting and remove extra characters as well as additional description following column names (in brackets)
#' @param df definitiion table as dataframe
#' @param VctAllColumns character vector containing all column names of interest (as shown in table) where preprocessing is performed
#' @param VctColstoupper character vector containing all column names without brackets for which the characters in columns will be converted to uppercase.
#' @return  df
#' @keywords definition
#' @export
#' @examples
#' PreProcessDfDefinitions(df,c("TS(Touchscreen)","ICD10", "ICD9", "OPCS4","OPCS3","READ2","READ2_drugs", "CTV3","BNF","DMD","f.20001(sr_cancer)", "f.20002(sr_noncancer)", "f.20003(sr_med)", "f.20004(sr_oper)"),c("ICD10","ICD9","OPCS4","OPCS3"))
PreProcessDfDefinitions<-function(df,VctAllColumns,VctColstoupper=NULL ){ # c("ICD10","ICD9","OPCS4","OPCS3")
## df<-dfDefinitions
  # check if nrows==1
  checkr=0
  checkc=0
  if(nrow(df)==1){df<-rbind(df,df);checkr=1}
  if(ncol(df)==1 & length(VctAllColumns)==1){df<-cbind(df,df);checkc=1; names(df)<-c(VctAllColumns,"tmp"); VctAllColumns=c(VctAllColumns,"tmp") }

  ## for the names: remove everything between dots (R converts symbols to dots "(,.-)/" etc ) <- this is due to data.frame(check.names= TRUE) in ProcessDfDefinitions, set to FALSE
  # names(df) <- gsub( " *\\..*?\\. *", "", names(df) )

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
  if("TS" %in% names(df)){
    VctCustomFields="TS"
    df[,VctCustomFields]<-gsub("\\b[=]+\\b","==",df[,VctCustomFields],perl=TRUE)
    df[,VctCustomFields]<-gsub("\\b[≥]\\b",">=",df[,VctCustomFields],perl=TRUE)
    df[,VctCustomFields]<-gsub("\\b[≤]\\b","<=",df[,VctCustomFields],perl=TRUE)
  }
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
FillInSRdefinitions<-function(df,Var="SR",cols=c("f.20001","f.20002","f.20004") ) {
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



#' ProcessDfDefinitions
#'
#' This function processes a tsv file with definitions and is automtically performed in CreateUKBiobankPhentoypes().
#' #It can be usefull to run this function as a check prior to running CreateUKBiobankPhentoypes
#' @param df definition table as dataframe
#' @param VctAllColumns character vector containing all column names of interest (as shown in table) where preprocessing is performed
#' @param VctColstoupper character vector containing all column names without brackets for which the characters in columns will be converted to uppercase.
#' @param fill_dependencies For composite traits whether to further process dependencies in Include_definitions/Exclude_from_cases/Study_population/Exclude_from_controls columns, default: TRUE
#' @return  processed definitiion table as dataframe
#' @keywords definition
#' @export
#' @examples
#' dfDefinitions_processed <- ProcessDfDefinitions(fread("definitions.tsv", colClasses = 'character', data.table = FALSE))
ProcessDfDefinitions<-function(df,
                               VctAllColumns=c("TS(Touchscreen)",
                                               "ICD10", "ICD9", "OPCS4","OPCS3",
                                               "READ2","READ2_drugs", "CTV3",
                                               "BNF","DMD",
                                               "f.20001(sr_cancer)",    "f.20002(sr_noncancer)", "f.20003(sr_med)", "f.20004(sr_oper)"
                                               ),
                               VctColstoupper=c("ICD10","ICD9","OPCS4","OPCS3"),
                               fill_dependencies=T){

  # df<- dfDefinitions  #  df<- dfDefinitions2

  df <- data.frame(df,check.names = FALSE)
  names(df) <- sub(pattern = "CODES",replacement = "",names(df) )
  names(df) <- sub(pattern = "_$",replacement = "",names(df) ) # "n_20002_" --> "n_20002"
  # replace bracket with dot
  # names(df) <- gsub( "([^.*])\\((.*)\\)", "\\1.\\2", names(df))
  # VctAllColumns <-   gsub( "([^.*])\\((.*)\\)", "\\1.\\2", VctAllColumns)
  # remove bracket and everything in it
  names(df) <- gsub( "\\(.*\\)", "", names(df))
  VctAllColumns <-   gsub( "\\(.*\\)", "", VctAllColumns)
  df <- PreProcessDfDefinitions(df,VctAllColumns,VctColstoupper=VctColstoupper)

  if(any(!VctAllColumns %in% names(df))) print(paste("WARNING missing columns:", paste(VctAllColumns[!VctAllColumns %in% names(df)],collapse=", ")))

  #################################
  CheckDuplicateTRAITS(df) # check duplicateids.
  ### df = excel matrix. 1 hij loopt een voor een over elke rij heen,
  # 2 zoekt per dependency in die rij de bijpassende rijen voor elke dependency  erbij en plakt die naast elkaar (inclusief de dependencies van de dependencies).
  # 4) dan delete hij de depenencies die hij ingevuld heeft # dit op repeat tot dat er geen dependencies meer zijn en alles is ingevuld.
  if(fill_dependencies==F){return(return(ConvertFactorsToStringReplaceNAInDf(df)))}

  # Include_definitions : formerly DEPENDENCY for composite trait
  # the while loop rewritten to a function to be run 4 times
  # the CASE table will be concat to the main table because one will be counted if they have any of the codes ?
  df<-parseIncludeExcludeCol(df,"Include_definitions",concat_to_df = TRUE,VctAllColumns=VctAllColumns)
  Include_in_cases <- df[,c("TRAIT","DESCRIPTION",VctAllColumns)]
  Exclude_from_cases <- lookup.codes(df=df,lookupcolumn = "Exclude_from_cases",VctAllColumns)
  Study_population <- lookup.codes(df=df,lookupcolumn = "Study_population",VctAllColumns)
  Exclude_from_controls <- lookup.codes(df=df,lookupcolumn = "Exclude_from_controls",VctAllColumns)

  dfDefinition <- rbind(Include_in_cases %>% dplyr::mutate(Definitions="Include_in_cases"),
        Exclude_from_cases %>% dplyr::mutate(Definitions="Exclude_from_cases"),
        Study_population %>%dplyr::mutate(Definitions="Study_population"),
        Exclude_from_controls %>% dplyr::mutate(Definitions="Exclude_from_controls")
        )

  dfDefinition <- ConvertFactorsToStringReplaceNAInDf(dfDefinition)

  ## other ways of storing: , I don't know what's best.
  #lst.dfs <- list(Include_in_cases=Include_in_cases,Exclude_from_cases=Exclude_from_cases,Study_population=Study_population,Exclude_from_controls=Exclude_from_controls)
  #lst.dfs <- lapply(lst.dfs,function(x) ConvertFactorsToStringReplaceNAInDf(x))
  #lst.lst <- convert_dataframelist_to_lst(lst.dfs)
  return(dfDefinition)
}




#' Helper function to look up the dependencies
#'
#' This function is used for trait dependency lookup
#' @param df definition table as dataframe
#' @param VctAllColumns names of columns in which codes should be updated
#' @param lookupcolumn name of dependency column that needs lookup
#' @return  a new dataframe containing definitions as the result of the lookup
#' @keywords auxiliary
#' @export
lookup.codes <- function(df,lookupcolumn="Exclude_from_cases",VctAllColumns){
  # copy the rows that needs to be looked up
  dfInEx<-df[!(df[[lookupcolumn]] == "" |is.na(df[[lookupcolumn]])),]
  print(paste(nrow(dfInEx),"traits with dependent trait in",lookupcolumn,sep=" "))
  if (nrow(dfInEx) !=0){
  # update these columns
  dfInEx[,VctAllColumns]<-NA

    # for each row , extract the codes from the definition table
    for(i in 1:nrow(dfInEx)) {
      def=unlist(strsplit(dfInEx[i,lookupcolumn],","))
      for(col in VctAllColumns){
        Vctcol <- df[df$TRAIT %in% def,col]
        dfInEx[i,col] <- pasteRemoveNA(Vctcol ,collapse=",",na.rm=T)
      }
    }

  }
  dfInEx <- dfInEx[,c("TRAIT","DESCRIPTION",VctAllColumns)]
  return(dfInEx)
}

#' Helper function to structure the processed definition
#'
#' This function is alternative way to store and output the definition table
#' @param lst.dfs  a list of definition tables (in dataframes)
#' @return  restructured list of definition tables organized by dependency status
#' @keywords auxiliary
#' @export
convert_dataframelist_to_lst <- function(lst.dfs){
  defs <- list()
  for (t in lst.dfs$Include_in_cases$TRAIT){
    def <- list()
    def$TRAIT <- lst.dfs$Include_in_cases[lst.dfs$Include_in_cases$TRAIT %in% t,'TRAIT']
    def$DESCRIPTION <- lst.dfs$Include_in_cases[lst.dfs$Include_in_cases$TRAIT %in% t,'DESCRIPTION']
    def$Include_in_cases <- c(lst.dfs$Include_in_cases[lst.dfs$Include_in_cases$TRAIT %in% t,VctAllColumns])
    def$Exclude_from_cases <- c(lst.dfs$Exclude_from_cases[lst.dfs$Exclude_from_cases$TRAIT %in% t,VctAllColumns])
    def$Study_population <- c(lst.dfs$Study_population[lst.dfs$Study_population$TRAIT %in% t,VctAllColumns])
    def$Exclude_from_controls <- c(lst.dfs$Exclude_from_controls[lst.dfs$Exclude_from_controls$TRAIT %in% t,VctAllColumns])
    defs[t] <- list(def)
  }
  #defs$lst.dfs<- lst.dfs
  return(defs)
}

#' Helper function to parse the trait dependencies
#'
#' Each composite trait with dependencies considered a tree as chain of dependency is possible. To parse and retrieve all codes needed, this function first searches for a leaf nodes and then add all codes of nodes along the way until root is reached; this operation is done recursively until all dependent codes are filled into the df
#' @param df definition table as dataframe
#' @param InExCol name of dependency column that needs lookup
#' @param concat_to_df option to update the input df otherwise store the codes in a new df,default: TRUE
#' @param VctAllColumns names of columns in which codes should be updated
#' @return  the original dataframe or a new dataframe with codes
#' @keywords definition
#' @export
parseIncludeExcludeCol <- function (df,InExCol,concat_to_df=FALSE,VctAllColumns){
  # result dataframe with non-empty rows in the corresponding inclusion/exclusion criteria
  dfInEx<-df[!(df[[InExCol]] == "" |is.na(df[[InExCol]])),]
  print(paste(nrow(dfInEx),"traits with dependent trait in",InExCol,sep=" "))
  if (nrow(dfInEx)==0){
    return(df)
  }
  dfInEx[,VctAllColumns]<-NA

  repeat {
    for(i in 1:nrow(df)) {
      row <- df[i,]
      # maybe this not needed but somehow stuck once in a loop where !is.na(row[[InExCol]] ==TRUE && if( is.na(targetrow[InExCol])) ==FALSE (so the set to NA was never reached) ?????????
      if(df[i,InExCol]==""| is.na(df[i,InExCol])) {df[i,InExCol]<-NA}

      # row
      if(!is.na(row[[InExCol]])){
        # parse the dependent traits
        VctInEx<-unlist(strsplit(row[[InExCol]],","))
        # remove space
        VctInEx<-  gsub(" ", "", VctInEx)
        # remove brackets if applicable
        VctInEx<- gsub( " *\\(.*?\\) *", "", VctInEx)

        # for each trait in dependent traits
        for (StrInEx in VctInEx) {
          # break self-referencing loop
          if(row$TRAIT == StrInEx) {stop("Dependency trait is same as trait.Please check for circular reference of traits.")}
          #  retrieve the row for the dependent trait
          targetrow<-df[df$TRAIT==StrInEx,]
          # print(StrInEx)
          if(nrow(targetrow)==0){stop(paste('Dependent trait: "',StrInEx,'" not found in traits ',row$TRAIT,sep="")) }
          # print("is.na(targetrow[InExCol])")
          # print(is.na(targetrow[InExCol]))
          # if this target trait is a leaf in the tree  i.e itself contains no dependency, start attaching its code to its parent

          if( is.na(targetrow[InExCol])){
            # for each classification/data source, add the corresponding codes in the
            for(col in VctAllColumns){
              Vctcol<-unique( unlist(strsplit( c(df[i,col],df[df$TRAIT==StrInEx,col]) ,",")) )
              if (concat_to_df==TRUE){
                  # if to concat codes to the input df
                  df[i,col]<-pasteRemoveNA(Vctcol ,collapse=",",na.rm=T)
                } else{
                  # else put in dfInEx
                  dfInEx[i,col]<-pasteRemoveNA(Vctcol ,collapse=",",na.rm=T)
                }
            }
            # remove InExCol that was just filled in:
            #df[i,"InExCol"]<-gsub(paste(StrInclude_in_case,sep=""),"",df[i,"Include_definitions"],fixed=TRUE,ignore.case=FALSE)
            # remove InExCol trait that was just filled in
            LstTmpDependencies<- unlist(strsplit(VctInEx,","))

            df[i,InExCol]<-paste( LstTmpDependencies [! LstTmpDependencies  %in%  StrInEx] ,sep="",collapse = ",")

            df[i,InExCol]<-gsub("^,*|(?<=,),|,*$", "", df[i,InExCol], perl=T)
            if(df[i,InExCol]==""){df[i,InExCol]<-NA} ## if empty replace with NA
          }
        }
      }
    }
    if( length(unique(is.na(df[[InExCol]])))==1 ) break
  }
  if (concat_to_df==TRUE){
    return(df)
  }else{return (dfInEx)
    }
}

#' Retrieve the required data fields in the ukb dataset *ukbxxxxx.tab* as stated in the definition table
#'
#' Survey the main dataset and collect all field names that will be used as specified in definition table. These include basic data fields specified in `default_ukb_fields()`, which can be updated accordingly, and the fields in touchscreen columns. If metadata is supplied , it will check if these fields are present in the maindataset.
#' @param dfDefinitions_processed definition table as dataframe
#' @param dfhtml the metadata file for main dataset
#' @return  a list of character vectors named all_ukb_fields,nondefault_ukb_fields and default_ukb_fields
#' @keywords definition
#' @export
#' @examples
#' dfDefinitions_processed <- ProcessDfDefinitions(fread("definitions.tsv", colClasses = 'character', data.table = FALSE))
#' dfhtml <- read_ukb_metadata(fhtml)
#' dfDefinitions_ukb_fields <- get_allvarnames(dfDefinitions_processed,dfhtml)
get_allvarnames <- function(dfDefinitions_processed,dfhtml=NULL,allow_miss=TRUE){
  #  dfDefinitions_processed
  # VctAllUKBVDefinitionColumns=c("TS") #set this variable to a selection of columns (dfDefinition columns) to be outputted by the _UKBV variable, default is 'VctAllUKBVDefinitionColumns=c("TS","SR","TS_RX","SR_RX","LAB")'
  # TS(Touchscreen) in definition
  TScolumns = "TS"
  defcols <- unlist(strsplit(na.omit(unname(unlist(as.character(dfDefinitions_processed[,c(TScolumns)])))),split=","))
  defcols <- unlist(strsplit(chartr("[]", "()", defcols),"[()]"))
  defcols <- gsub("(≥|≤|=|>|<).*","",defcols)
  defcols <- gsub("≥.*","",defcols)
  defcols <- unique(defcols) ## need column classes....


  all_ukb_fields <- unique(c(defcols,default_ukb_fields()))
  all_ukb_fields<-gsub("[a-zA-Z]*?_","",all_ukb_fields)
  nondefault_ukb_fields <- all_ukb_fields[!all_ukb_fields %in% default_ukb_fields() ]

  if (!is.null(dfhtml)){
  # check for missing fields if there is metadata file
  if (all(all_ukb_fields %in%dfhtml$field.showcase)){
    message("All fields required are present in the main dataset.")

  }else{
    fields_missed<- all_ukb_fields[!all_ukb_fields %in%dfhtml$field.showcase]
    message(glue::glue("WARNING: {length(fields_missed)} fields not found in the main dataset: {glue::glue_collapse(fields_missed, sep = ',')}"))
    if (isTRUE(allow_miss)){
      # grep fields that are present for extraction
      all_ukb_fields<-all_ukb_fields[all_ukb_fields %in% dfhtml$field.showcase]
    }else{
      return(NULL)
    }

  }
  # no check for missing fields
  }else{
    message("No metadata file given, some fields may be missing.")
  }

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

#' Expand the illness/medication codes in definition table to include codes that belong to this code
#'
#' Explicitly fill in codes covered by the parent codes that are available in the raw data (using lst.counts)stated in the definition table
#' @param dfDefinitions_processed definition table as dataframe
#' @param lst.data.settings data.setting as dataframe, it includes settings that should be used, e.g. if ICD10 should be looked up case sensitive or not (incase of READ cases are important, dotts should also NOT be interpreted)
#' @param lst.counts a summary of the raw data produced by `get_lst_counts()`
#' @return  definition table as dataframe
#' @keywords definition
#' @export
expand_dfDefinitions_processed <- function(dfDefinitions_processed,lst.data.settings,lst.counts){
  classifications <- lst.data.settings %>% dplyr::filter(expand_codes==1) %>% dplyr::pull (classification) %>% unique()
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



#' Traverse hierarchical code dictionary
#'
#' Given a hierchical dictionary, identify all children codes under the input parent code.
#' @param dfcode hierarchical tree-structured code dictionary with at least columns "coding","selectable","node" and "parent".
#' @param codeVct Input code as character vector
#' @return  All selectable codes (input code inclusive) as character vector
#' @keywords auxiliary
#' @export
add_child_nodes <-function(dfcode,codeVct){
    resultCodeVct<- vector()
    currVct<-codeVct
    while (length(currVct)!=0){
      if (length(currVct)>0){
        resultCodeVct<-unique(c(resultCodeVct,currVct))}
      # get node ID of the ICD codes
      parentNodeId<-unlist(lapply(currVct,function(x) dfcode[dfcode$coding==x,]$node_id ))
     # find nodes belonging to these nodes
      childrenNodes<-unique(unlist(lapply(parentNodeId,function(x) dfcode[dfcode$parent_id==x,]$coding)))
      currVct<-childrenNodes
    }
    # remove codes not selectable
    keep<-unlist(lapply(resultCodeVct, function(x) dfcode[dfcode$coding==x,]$selectable=="Y"))
    resultCodeVct<-resultCodeVct[keep]
    return(resultCodeVct)
}

#' Expand the illness/medication codes in definition table to include codes that belong to this code
#'
#' Explicitly fill in codes covered by the parent codes stated in the definition table. This differs from `expand_dfDefinitions_processed()` in the dictionaries taken - it reads the code dictionaries stated in lst.data.setting and expand to available codes accordingly. Codes not present **will be removed**.
#' @param dfDefinitions_processed definition table as dataframe
#' @param lst.data.settings data.setting as dataframe, it includes settings that should be used, e.g. if ICD10 should be looked up case sensitive or not (incase of READ cases are important, dotts should also NOT be interpreted)
#' @param code_map_dir the directory storing the coding files
#' @return  definition table as dataframe
#' @keywords definition
#' @export
#' @example
#' expand_dfDefinitions_processed2(dfDefinitions_processed,lst.data.settings,paste0(system.file("extdata", package="ukbpheno"),"/"))
expand_dfDefinitions_processed2 <-
  function(dfDefinitions_processed,
           lst.data.settings,code_map_dir) {
    message("Expand the codes in the definition table")

    classifications <-
      lst.data.settings %>% dplyr::filter(expand_codes == 1) %>% dplyr::pull (classification) %>% unique()

    lst.codemap<-list()

    for (cls in classifications) {
      fmap=paste(code_map_dir,unique(lst.data.settings[lst.data.settings$classification==cls,]$code_map),sep="")
      #look up if case sensitive
      ignore.case = unique(lst.data.settings[lst.data.settings$classification %in% cls, 'ignore.case'])[1]
      # read file
      message(glue::glue("Read from codings for {cls} from {fmap}"))
      lst.codemap[[cls]]<-fread(fmap)
      # rename the column
      if (ncol(lst.codemap[[cls]])!=1){
        # for files downloaded from showcase
        names(lst.codemap[[cls]])[grep("^cod", names(lst.codemap[[cls]]))] <- "coding"
      }else{
        # for codes created locally / if there is only 1 column it has to be the code
        names(lst.codemap[[cls]])<-"coding"
      }
      # if case insensitive,
      # always change to upper letters as this has been done in the preprocessingg of definitiion
      if (isTRUE(ignore.case[1]$ignore.case)) {
        lst.codemap[[cls]]$coding <- toupper(lst.codemap[[cls]]$coding)
      }
      for (r in 1:nrow(dfDefinitions_processed)) {
        # for loops just as fast as apply in this case..
        VctStr = unlist(strsplit(dfDefinitions_processed[r, cls], ","))

        #if code map is hierarchical, traverse the tree to get child codes
        if ( unique(lst.data.settings[lst.data.settings$classification==cls,]$hierarchical_map)) {
          # get all nodes down from codemap
          Str_expanded <-
            paste(add_child_nodes(lst.codemap[[cls]], VctStr), collapse = ",")
          dfDefinitions_processed[r, cls] <- Str_expanded
        } else{
          # otherwise grep patterns (codes) that starts with the input code
          Str_expanded <- paste(unique(unlist(
            lapply(VctStr,  function(x){
              lst.codemap[[cls]]$coding[grep(paste("^", x, sep = ""),
                                                      lst.codemap[[cls]]$coding  ,
                                                      ignore.case = ignore.case)]}
              )


          )), collapse = ",")
          dfDefinitions_processed[r, cls] <- Str_expanded
              # print(dfDefinitions_processed[r, cls])
          # next

          ############
        }

      }
    }
    return(dfDefinitions_processed)
  }




#' Check if the illness/medication codes in definition table are found the coding files
#'
#' Given a definition table and the directory with codings, check if codes in definition table are present in the coding files.
#' @param dfDefinitions_processed definition table as dataframe
#' @param lst.data.settings data.setting as dataframe, it includes settings that should be used, e.g. if ICD10 should be looked up case sensitive or not (incase of READ cases are important, dotts should also NOT be interpreted)
#' @param code_map_dir the directory storing the coding files
#' @param check_expandable_codes_only boolean flag to indicate whether or not to check the only codes to be expanded
#' @return  character vector that are missing
#' @keywords definition
#' @export
#' @example
#' check_dfDefinitions_codes(dfDefinitions_processed,lst.data.settings,code_map_dir=paste0(system.file("extdata", package="ukbpheno"),"/"))
check_dfDefinitions_codes <-function(dfDefinitions_processed,
                                     lst.data.settings,code_map_dir="data/",check_expandable_codes_only=T){
   all_missing_codes<-list()
   #  get all codings needed in data.settings
   if (!check_expandable_codes_only){
     classifications <-
       lst.data.settings %>%  dplyr::pull (classification) %>% unique()
   }else{
     classifications <-
       lst.data.settings %>% dplyr::filter(expand_codes == 1) %>% dplyr::pull (classification) %>% unique()
   }


   lst.codemap<-list()
   for (cls in classifications) {
     # no code to compare to
     if (is.na(unique(lst.data.settings[lst.data.settings$classification==cls,]$code_map))) next
    fmap=paste(code_map_dir,unique(lst.data.settings[lst.data.settings$classification==cls,]$code_map),sep="")
    #look up if case sensitive
    ignore.case = unique(lst.data.settings[lst.data.settings$classification %in% cls, 'ignore.case'])[1]
    # read file
    message(glue::glue("Read from codings for {cls} from {fmap}"))
    lst.codemap[[cls]]<-fread(fmap)
    # rename the column
    if (ncol(lst.codemap[[cls]])!=1){
      # for files downloaded from showcase
      names(lst.codemap[[cls]])[grep("^cod", names(lst.codemap[[cls]]))] <- "coding"
    }else{
      # for codes created locally / if there is only 1 column it has to be the code
      names(lst.codemap[[cls]])<-"coding"
    }
    # if case insensitive,
    # always change to upper letters as this has been done in the preprocessingg of definitiion
    if (isTRUE(ignore.case[1]$ignore.case)){
      lst.codemap[[cls]]$coding <- toupper(lst.codemap[[cls]]$coding)
    }

    codes<-na.omit(unlist(strsplit(dfDefinitions_processed[[cls]],",")))
    # direct match
    missing_codes<-dplyr::setdiff(codes,lst.codemap[[cls]]$coding)

    # grep match
    # this is not a very ideal line as it doesn't recognize 2 different expand code options for the same classification (if needed) and will always only take the first row
    if ((! unique(lst.data.settings[lst.data.settings$classification==cls,]$hierarchical_map))&(unique(lst.data.settings[lst.data.settings$classification==cls,]$expand_codes==1) )) {
    extended_matches <- lapply(codes,  function(x) any(stringr::str_detect(lst.codemap[[cls]]$coding,stringr::regex(paste("^", x, sep = ""),ignore_case =ignore.case[1]$ignore.case ))))

    if (length(extended_matches)>0){
      missing_codes<-codes[!unlist(extended_matches)]
    }

    }

    missing_codes<-na.omit(missing_codes)

    message(glue::glue("Number of missing codes in {cls} (will be ignored):{length(unique(missing_codes))} \n******************************"))

    all_missing_codes<-append(all_missing_codes,setNames(as.data.frame(missing_codes,stringsAsFactors =F),cls))
   }
   all_missing_codes<-all_missing_codes[lapply(all_missing_codes,length)>0]
   dt_all_missing_codes <- data.table::rbindlist(
     lapply(all_missing_codes, function(x) data.table::data.table(t(unique(x)))),
     fill = TRUE
   )
   dt_all_missing_codes<-data.table::data.table(t(dt_all_missing_codes))
   colnames(dt_all_missing_codes)<-names(all_missing_codes)
   message(glue::glue("Write missing codes to file missing_code.tsv"))
   data.table::fwrite(dt_all_missing_codes,"missing_code.tsv",sep="\t")
   return(dt_all_missing_codes)
}


# View(lst.data$tte.gpclincal.read2)
# head(lst.data$tte.gpclincal.read2)
