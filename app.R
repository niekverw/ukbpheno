setwd("/Users/niek/repos/ukbpheno")

source("ProcessdfDefinitions.R")
library(readxl)
library(dplyr)
library(data.table)
library(stringr)
fcoding.xls="data/all_lkps_maps.xlsx"

#### READ V2
dfCodesheet.read_v2_read_ctv3 <- as.data.frame(read_xlsx(fcoding.xls,sheet="read_v2_read_ctv3"))[,c(2,7)]
colnames(dfCodesheet.read_v2_read_ctv3)<-c("read_code","CTV3")

dfCodesheet.read_v2_icd9 <- as.data.frame(read_xlsx(fcoding.xls,sheet="read_v2_icd9"))[,c(1,2)]
dfCodesheet.read_v2_icd10 <- as.data.frame(read_xlsx(fcoding.xls,sheet="read_v2_icd10"))[,c(1,2)]
dfCodesheet.read_v2_opcs4 <- as.data.frame(read_xlsx(fcoding.xls,sheet="read_v2_opcs4"))[,c(1,2)]

dfCodesheet.read_v2_lkp <- as.data.frame(read_xlsx(fcoding.xls,sheet="read_v2_lkp"))
dfCodesheet.read_v2_lkp <- as.data.frame(dfCodesheet.read_v2_lkp %>% arrange(read_code,term_code))
#dfCodesheet.read_v2_lkp <- dfCodesheet.read_v2_lkp[dfCodesheet.read_v2_lkp$term_code==0,]
dfCodesheet.read_v2_drugs_lkp<- as.data.frame(read_xlsx(fcoding.xls,sheet="read_v2_drugs_lkp"))
dfCodesheet.read_v2_lkp <- rbind(dfCodesheet.read_v2_lkp[,c(1,3)],dfCodesheet.read_v2_drugs_lkp[,1:2])
dfCodesheet.read_v2_lkp <- as.data.frame(dfCodesheet.read_v2_lkp %>% group_by(read_code) %>%  summarize(text = str_c(term_description, collapse = "/")))
#dfCodesheet.read_v2_lkp$read_code <- gsub("\\.","", dfCodesheet.read_v2_lkp$read_code)
dfCodesheet.read_v2_lkp <- dfCodesheet.read_v2_lkp %>% unique()  %>% arrange(read_code)
dfCodesheet.read_v2_lkp$text <- str_replace_all(dfCodesheet.read_v2_lkp$text,"[^/[:^punct:]]", "") 

# meds; read to ukb code
load("data/dfCodesheetREAD_SR.Coding.RData")
colnames(dfCodesheetREAD_SR.Coding)<- c("UKB.Coding","text","read_code","term.id")
dfCodesheet.read_v2_UKBmeds <- dfCodesheetREAD_SR.Coding[,c(3,1)]
dfCodesheet.read_v2_lkp<-rbind(dfCodesheet.read_v2_lkp,dfCodesheetREAD_SR.Coding[!dfCodesheetREAD_SR.Coding$READ.CODE %in% dfCodesheet.read_v2_lkp$read_code,c("read_code","text")])

dfCodesheet.READ <- merge(dfCodesheet.read_v2_read_ctv3,dfCodesheet.read_v2_icd9,by="read_code",all=T)
dfCodesheet.READ <- merge(dfCodesheet.READ,dfCodesheet.read_v2_icd10,by="read_code",all=T)
dfCodesheet.READ <- merge(dfCodesheet.READ,dfCodesheet.read_v2_opcs4,by="read_code",all=T)
dfCodesheet.READ <- merge(dfCodesheet.READ,dfCodesheet.read_v2_UKBmeds,by="read_code",all=T)
dfCodesheet.READ <- merge(dfCodesheet.READ,dfCodesheet.read_v2_lkp,by="read_code",all=T)
dfCodesheet.READ <- unique(dfCodesheet.READ)
colnames(dfCodesheet.READ) <- c("READ","CTV3","ICD10","ICD9","OPCS4","n_20003","text")
dfCodesheet.READ$source="READ"
dfCodesheet.READ <- data.table(dfCodesheet.READ)
dfCodesheet.READ$text <- str_replace_all(dfCodesheet.READ$text,"[^/[:^punct:]]", "") 

setkey(dfCodesheet.READ,"READ")
##############################
##### READ V3
dfCodesheet.read_ctv3_readv2 <- as.data.frame(read_xlsx(fcoding.xls,sheet="read_ctv3_read_v2"))[,c(1,5)]
dfCodesheet.read_ctv3_readv2 <- dfCodesheet.read_ctv3_readv2[!dfCodesheet.read_ctv3_readv2$READV2_CODE %in% "_NONE",]
colnames(dfCodesheet.read_ctv3_readv2)<-c("read_code","READ")
  
dfCodesheet.read_ctv3_icd9 <- as.data.frame(read_xlsx(fcoding.xls,sheet="read_ctv3_icd9"))[,c(1,2)]
dfCodesheet.read_ctv3_icd10 <- as.data.frame(read_xlsx(fcoding.xls,sheet="read_ctv3_icd10"))[,c(1,2)]
dfCodesheet.read_ctv3_opcs4 <- as.data.frame(read_xlsx(fcoding.xls,sheet="read_ctv3_opcs4"))[,c(1,2)]
dfCodesheet.read_ctv3_lkp <- as.data.frame(read_xlsx(fcoding.xls,sheet="read_ctv3_lkp"))[,c(1,2)]
dfCodesheet.read_ctv3_lkp <- as.data.frame(dfCodesheet.read_ctv3_lkp %>% group_by(read_code) %>%  summarize(text = str_c(term_description, collapse = "/")))

dfCodesheet.CTV3 <- merge(dfCodesheet.read_ctv3_readv2,dfCodesheet.read_ctv3_icd10,by="read_code",all=T)
dfCodesheet.CTV3 <- merge(dfCodesheet.CTV3,dfCodesheet.read_ctv3_icd9,by="read_code",all=T)
dfCodesheet.CTV3 <- merge(dfCodesheet.CTV3,dfCodesheet.read_ctv3_opcs4,by="read_code",all=T)
dfCodesheet.CTV3 <- merge(dfCodesheet.CTV3,dfCodesheet.read_ctv3_lkp,by="read_code",all=T)
dfCodesheet.CTV3 <- unique(dfCodesheet.CTV3)
colnames(dfCodesheet.CTV3) <- c("CTV3","READ","ICD10","ICD9","OPCS4","text")
dfCodesheet.CTV3$n_20003 <- NA
dfCodesheet.CTV3$source="CTV3"
dfCodesheet.CTV3$text <- str_replace_all(dfCodesheet.CTV3$text,"[^/[:^punct:]]", "") 

dfCodesheet.CTV3 <- data.table(dfCodesheet.CTV3)
setkey(dfCodesheet.CTV3,"CTV3")

####### ICD10
dfCodesheet.icd10_icd9<- as.data.frame(read_xlsx(fcoding.xls,sheet="icd9_icd10"))[,c(3,1)]
colnames(dfCodesheet.icd10_icd9) <- c("ICD10","ICD9")
dfCodesheet.icd10_lkp <- as.data.frame(read_xlsx(fcoding.xls,sheet="icd10_lkp"))[,c(2,5)] # not complete (e.g. X*)
names(dfCodesheet.icd10_lkp) <- c("ICD10","text")
dfCodesheet.ICD10.coding19 <- data.frame(fread("data/ICD10.coding19.tsv"))[,1:2] # <- contains all of the above.
names(dfCodesheet.ICD10.coding19) <- c("ICD10","text")

dfCodesheet.icd10_lkp <- rbind(dfCodesheet.icd10_lkp,dfCodesheet.ICD10.coding19[!dfCodesheet.ICD10.coding19$ICD10 %in% dfCodesheet.icd10_lkp$ICD10,])

dfCodesheet.ICD10 <-  merge(dfCodesheet.icd10_icd9,dfCodesheet.icd10_lkp,by="ICD10",all=T)

dfCodesheet.ICD10$READ <- NA
dfCodesheet.ICD10$CTV3 <- NA
dfCodesheet.ICD10$OPCS4 <- NA
dfCodesheet.ICD10$n_20003 <- NA
dfCodesheet.ICD10$source="ICD10"
dfCodesheet.ICD10$text <- str_replace_all(dfCodesheet.ICD10$text,"[^/[:^punct:]]", "") 

dfCodesheet.ICD10 <- data.table(dfCodesheet.ICD10)
setkey(dfCodesheet.ICD10,"ICD10")

####### # ICD9 depscription, dfCodesheet.icd9_lkp
dfCodesheet.icd9_lkp <- as.data.frame(read_xlsx(fcoding.xls,sheet="icd9_lkp")) # certainly not complete!
colnames(dfCodesheet.icd9_lkp) <- c("ICD9","text")
dfCodesheet.ICD9.coding87 <- data.frame(fread("data/ICD9.coding87.tsv"))[,1:2] # UKB
names(dfCodesheet.ICD9.coding87)<- c("ICD9","text")
dfCodesheet.icd9_lkp <- rbind(dfCodesheet.icd9_lkp,dfCodesheet.ICD9.coding87[!dfCodesheet.ICD9.coding87$ICD9 %in% dfCodesheet.icd9_lkp$ICD9,])

dfCodesheet.ICD9 <- dfCodesheet.icd9_lkp

dfCodesheet.ICD9$ICD10 <- NA
dfCodesheet.ICD9$READ <- NA
dfCodesheet.ICD9$CTV3 <- NA
dfCodesheet.ICD9$OPCS4 <- NA
dfCodesheet.ICD9$n_20003 <- NA
dfCodesheet.ICD9$source="ICD9"
dfCodesheet.ICD9$text <- str_replace_all(dfCodesheet.ICD9$text,"[^/[:^punct:]]", "") 
dfCodesheet.ICD9 <- data.table(dfCodesheet.ICD9)
setkey(dfCodesheet.ICD9,"ICD9")

###### OPCS4
dfCodesheet.OPCS4.coding240 <- data.frame(fread("data/OPCS4.coding240.tsv"))[,c(1,2)]
names(dfCodesheet.OPCS4.coding240) <- c("OPCS4","text")
dfCodesheet.OPCS4 <- dfCodesheet.OPCS4.coding240
dfCodesheet.OPCS4$ICD10 <- NA
dfCodesheet.OPCS4$ICD9 <- NA
dfCodesheet.OPCS4$READ <- NA
dfCodesheet.OPCS4$CTV3 <- NA
dfCodesheet.OPCS4$n_20003 <- NA
dfCodesheet.OPCS4$source="OPCS4"
dfCodesheet.OPCS4$text <- str_replace_all(dfCodesheet.OPCS4$text,"[^/[:^punct:]]", "") 
dfCodesheet.OPCS4 <- data.table(dfCodesheet.OPCS4)
setkey(dfCodesheet.OPCS4,"OPCS4")

#####  FINAL OBJECT. 
LstdfCodesheets <- list(
  READ = dfCodesheet.READ[,c("READ", "CTV3", "ICD10", "ICD9", "OPCS4", "n_20003", "text","source")],
  CTV3 = dfCodesheet.CTV3[,c("READ", "CTV3", "ICD10", "ICD9", "OPCS4","n_20003", "text","source")],
  ICD10 = dfCodesheet.ICD10[,c("READ", "CTV3", "ICD10", "ICD9", "OPCS4","n_20003", "text","source")],
  ICD9 = dfCodesheet.ICD9[,c("READ", "CTV3", "ICD10", "ICD9", "OPCS4","n_20003", "text","source")],
  OPCS4 = dfCodesheet.OPCS4[,c("READ", "CTV3", "ICD10", "ICD9", "OPCS4","n_20003", "text","source")],
  ALL = rbind(dfCodesheet.READ,dfCodesheet.CTV3,dfCodesheet.ICD10,dfCodesheet.ICD9,dfCodesheet.OPCS4)
)
########################################
##### FUNCTIONS to convert different codings. 
########################################

convert.coding<- function(Str,
                          from.code="READ.CODE",
                          to.code="UKB.Coding",
                          lookuptable=dfCodesheetREAD_SR.Coding,ignore.case=FALSE){
  # Str<-"f3...,f36z.,f31"
  Str<-as.character(Str)
  Str=gsub(pattern = ".",replacement = "\\.",Str,fixed=T) ### INMPORTANT FOR READ CODEES!!!!!
  if(is.na(Str)) { return(NA)}
  VctStr<-unlist(strsplit(Str,","))
  #VctRXstrings<-strsplit(df[!is.na(df$n_20003_),]$n_20003_,",")[[1]]
  c <- paste(unique(unlist(lapply(VctStr,  function(x)
    lookuptable[,get(to.code)] [ grep(paste("^", x,sep=""),lookuptable[,get(from.code)] ,ignore.case=ignore.case )]
  ))),collapse=",")
  
  return(c)
}


add.description.to.codes <- function(Str,code.id="UKB.Coding",description.id="Meaning",description.lookuptable=dfCodesheetREAD_SR.Coding,ignore.case=FALSE,firstcodeonly=TRUE) {
  if(is.na(Str) | Str =="NA"){return(Str)}
  Str<-as.character(Str)
  if(is.na(Str)) { return(NA)}
  
  VctStr<-unlist(strsplit(Str,","))
  
  
  
  c<- sapply(VctStr,  function(x){
    x.d <- description.lookuptable[,get(description.id)] [ grep(paste("^", x,sep=""),description.lookuptable[,get(code.id)] ,ignore.case=ignore.case )]
    if(length(x.d)==0){return(paste0(x," (NA)"))}
    x.d <- str_replace_all(x.d,  "[^/[:^punct:]]", "") # replace all symbols to not mess up downstream things.
    if(firstcodeonly==TRUE){
      x.d <- x.d[1]
    }
    x.d <- paste0(x.d,collapse=" /")
    x.d <- paste0(x," (",x.d,")")
    x.d
  },USE.NAMES = F)
  
  c <- paste(unique(unname(c)),collapse=",")
  return(c)
}

add.description.to.vectorofcodes <- function(vctcodes=c("G551.","G55.."),code.id="CTV3",description.id="text",description.lookuptable=LstdfCodesheets$CTV3,ignore.case=FALSE,firstcodeonly=TRUE) {
  #c <- unlist(sapply(vctcodes,na.omit,USE.NAMES = F))
  if(length(vctcodes)==0){return(c)}
  
  c<-as.character(vctcodes)
  df_c <- description.lookuptable[.(c)] #[[code.id]]
  df_c <- df_c %>% select(c=eval(code.id),text) %>% unique %>% filter(!is.na(c)) %>% dplyr::group_by(c) %>% mutate(text = paste0(text, collapse = "/")) %>% unique %>% mutate(c_text=paste0(c," (",text,")"))
  df_c<- data.frame(df_c)
  #c <- paste(unique(unname(c)),collapse=",")
  return(df_c)
}


expand_clean_codes <- function(col=df$ICD10CODES,
                               from.code="ALT_CODE",
                               description.id='DESCRIPTION',
                               lookuptable = dfCodesheet.icd10_lkp,
                               add_description=T){
  # col=df$ICD10CODES
  # lookuptable=LstdfCodesheets$ICD10
  # from.code="ICD10"
  # description.id='text'
  
  
  col <- PreProcessDfDefinitions( df=data.frame(col=col,tmp=rep("NA",length(col))),VctAllColumns = c("col","tmp"),VctColstoupper="col")[,1]
  
  to.code="self"
  lookuptable$self <- lookuptable[,get(from.code)]
  #c <- unlist(lapply(col, convert.coding,from.code=from.code,to.code=to.code,lookuptable=lookuptable))
  
  c1 <- paste(col, unlist(lapply(col, convert.coding,from.code=from.code,to.code=to.code,lookuptable=lookuptable)),sep=",")
  c2 <- unlist(lapply(c1,function(x) {  x = unique(strsplit(x,"," )[[1]]); if(length(x)==1 & x[1] =="NA"){ return("NA")} else{ return( paste(x[x != "NA"],collapse=",") )} }))
  
  if(add_description==T){
    c2 <- sapply( c2, add.description.to.codes,
                  code.id=from.code,
                  description.id=description.id,
                  description.lookuptable=lookuptable,USE.NAMES = F)
  }
  return(c2)
}

lookup_list_in_df <- function(lst,df.lookup){
  # lookup list in df, as long as names are corresponding to df headers..
  df.all<-data.frame()
  for (i in 1:length(lst)){
    lookup=na.omit(lst[i][[1]])
    d <- df.lookup[ df.lookup[,get(names(lst[i]) )] %in% lookup ,] %>% unique()
    if(nrow(d)>1) d$source = names(lst[i])
    df.all<-rbind(df.all,d)
    
  }
  df.all <- df.all %>% unique()
  return(df.all)
}

annotate_codes <- function(c,LstdfCodesheets=LstdfCodesheets){
  lst_lookup_annotated<-list()
  for (i in 1:length(c)){
    c_anno <- add.description.to.vectorofcodes(vctcodes = c[i][[1]],code.id = names(c[i]),description.id = "text",description.lookuptable = LstdfCodesheets[[names(c[i])]])
    lst_lookup_annotated[[names(c[i])]] = c_anno
    #add.description.to.vectorofcodes(lst_lookup$ICD10,code.id = "ICD10",description.id = "text",description.lookuptable = LstdfCodesheets$ICD10)
  }
  return(lst_lookup_annotated)
}

# suggest new codes; 
lookup_codes <- function(codes=row, LstdfCodesheets=LstdfCodesheets,expand_input=T ,fromcodes = c("ICD10","ICD9","READ","CTV3","OPCS4"),cols=c( "ICD10CODES", "ICD9CODES", "READCODES", "CTV3CODES","OPCS4CODES")){
  # codes=row_exp$ICD10CODES
  # df.lookup=LstdfCodesheets$ALL #dfCodesheet.ICD10
  # cols=c( "ICD10CODES", "ICD9CODES", "READCODES", "CTV3CODES","OPCS4CODES")
  # fromcodes <- c("ICD10","ICD9","READ","CTV3","OPCS4")
  input=codes
  if(expand_input){
    for(code in fromcodes){
      colname <- glue::glue("{code}CODES")
      codes[[colname]] <- expand_clean_codes(col = codes[[colname]] , from.code=code,description.id='text',lookuptable = LstdfCodesheets[code][[1]],add_description=F)
    }
    # row$ICD10CODES <- expand_clean_codes(col =row$ICD10CODES, from.code="ICD10",description.id='text',lookuptable = LstdfCodesheets["ICD10"][[1]],add_description=F)
    # row$ICD9CODES <- expand_clean_codes(col =row$ICD9CODES, from.code="ICD9",description.id='text',lookuptable = LstdfCodesheets["ICD9"][[1]],add_description=F)
    # row$READCODES <- expand_clean_codes(col =row$READCODES, from.code="READ",description.id='text',lookuptable = LstdfCodesheets["READ"][[1]],add_description=F)
    # row$CTV3CODES <- expand_clean_codes(col =row$CTV3CODES, from.code="CTV3",description.id='text',lookuptable =  LstdfCodesheets["CTV3"][[1]],add_description=F)
    # row$OPCS4CODES <- expand_clean_codes(col =row$OPCS4CODES, from.code="OPCS4",description.id='text',lookuptable =  LstdfCodesheets["OPCS4"][[1]],add_description=F)
  }
  
  r <- PreProcessDfDefinitions( df=codes,VctAllColumns = cols ,VctColstoupper=c())
  c <- sapply(r[cols],function(x) unique(unlist(strsplit(x,","))))
  names(c) <- sub(pattern = "CODES","",names(c))
  df_lookup <- lookup_list_in_df(lst = c,df.lookup = LstdfCodesheets$ALL)
  lst_lookup <- apply(df_lookup[,fromcodes,with=FALSE],2,function(x) c(na.omit(unique(x)) ))
  #lst_lookup annotated:
  lst_lookup_anno <- annotate_codes(c=lst_lookup,LstdfCodesheets = LstdfCodesheets)
  
  df_lookup_anno <- df_lookup
  for(col in fromcodes){
    df_lookup_anno[[col]] <-   lst_lookup_anno[[col]]$c_text[match(df_lookup[,get(col)] ,lst_lookup_anno[[col]]$c)]
  }
  
  return(list(input=input,
              input_c=c,
              lst_lookup=lst_lookup,
              df_lookup=df_lookup,
              lst_lookup_anno=lst_lookup_anno,
              df_lookup_anno=df_lookup_anno))
}

########################################
##### STAND ALONE EXAMPLE FOR 1 ROW. 
########################################
# ICD10, icd9, read2,ct3,
dfDefinitions_file="/Users/niek/repos/ukpheno/data/dfDefinitions.tsv"
df = data.frame(fread(dfDefinitions_file))
df <- ProcessDfDefinitions(df=df,fill_dependencies = F)
df <- df[,c("TRAIT","DESCRIPTION", "ICD10CODES","ICD9CODES","READCODES","CTV3CODES","OPCS4CODES")]

# DO IIT FOR 1 ROW suggest codes for 1 selected row. 
irow=13#12
row <- df[irow,]
#row$OPCS4CODES <- "K02"


codes.lookup <- lookup_codes(codes = row,LstdfCodesheets=LstdfCodesheets,expand_input=T)


#### SHINY: 


library(shiny)
library(DT)
# Define UI for app that draws a histogram ----
ui <- fluidPage(
  # App title ----
  titlePanel("Hello Shiny!"),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Slider for the number of bins ----
      textAreaInput(inputId="iICD10", label="ICD10", value = "I421,I422", width = NULL, placeholder = NULL),
      textAreaInput(inputId="iICD9", label="ICD9", value = "4251", width = NULL, placeholder = NULL),
      textAreaInput(inputId="iREAD", label="READ", value = "", width = NULL, placeholder = NULL),
      textAreaInput(inputId="iCTV3", label="CTV3", value = "", width = NULL, placeholder = NULL),
      textAreaInput(inputId="iOPCS4", label="OPCS4", value = "", width = NULL, placeholder = NULL),
      checkboxInput(inputId = "iExpandcodes", "Expand codes, e.g. I50 -> I501,I502, etc.  ", FALSE),
      actionButton("goButton", "Go!"),
      HTML("<br><br>note; this is for exploration, translations are not reliable. will add some other options. medication lookups are not supported (yet).")
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      #verbatimTextOutput("oICD10")
      #DT::dataTableOutput("table_input")
      
      tabsetPanel(
        tabPanel("ICD10 ", DT::dataTableOutput("table_oICD10")),
        tabPanel("ICD9",  DT::dataTableOutput("table_oICD9")),
        tabPanel("READ",  DT::dataTableOutput("table_oREAD")),
        tabPanel("CTV3",  DT::dataTableOutput("table_oCTV3")),
        tabPanel("OPCS4",  DT::dataTableOutput("table_oOPCS4")),
        tabPanel("raw_data",  DT::dataTableOutput("table_oraw"))
        
        
        
      )
    )
  )
)


server <- function(input, output) {

  
    observeEvent(input$goButton, {
    showModal(modalDialog( "please wait" ,easyClose = FALSE,footer=NULL))
    row <- data.frame(ICD10CODES=input$iICD10,
               ICD9CODES=input$iICD9,
               READCODES=input$iREAD,
               CTV3CODES=input$iCTV3,
               OPCS4CODES=input$iOPCS4)
    
  
    codes.lookup <- lookup_codes(codes = row,LstdfCodesheets=LstdfCodesheets,expand_input=input$iExpandcodes)
    
    
    convert_lookup_to_df <- function(codes,input_c=row_exp$ICD10CODES){
      data.frame(codes=codes$c,codes$text, new=codes$c %in% input_c)
      # codes_ <- sapply(strsplit(unname(codes)," "),function(x) x[1]) 
      # df <- data.frame(codes=names(codes), expanded_codes=unname(codes),new=!codes_ %in% strsplit(input,split = ",")[[1]])
      # df
    }
    
    output$table_oICD10 = DT::renderDataTable({
      convert_lookup_to_df(codes=codes.lookup$lst_lookup_anno$ICD10,input_c=codes.lookup$input_c$ICD10)
    })
    output$table_oICD9 = DT::renderDataTable({
      convert_lookup_to_df(codes=codes.lookup$lst_lookup_anno$ICD9,input_c=codes.lookup$input_c$ICD9)
    })
    output$table_oREAD = DT::renderDataTable({
      convert_lookup_to_df(codes=codes.lookup$lst_lookup_anno$READ,input_c=codes.lookup$input_c$READ)
    })
    output$table_oCTV3 = DT::renderDataTable({
      convert_lookup_to_df(codes=codes.lookup$lst_lookup_anno$CTV3,input_c=codes.lookup$input_c$CTV3)
    })
    output$table_oOPCS4 = DT::renderDataTable({
      convert_lookup_to_df(codes=codes.lookup$lst_lookup_anno$OPCS4,input_c=codes.lookup$input_c$OPCS4)
    })
   # data.frame(x="123",t="asd")
    removeModal()
  })

  output$table_oraw = DT::renderDataTable({
    LstdfCodesheets$ALL
  })
  #output$oICD10 <- renderText({ })
}
shinyApp(ui, server)



