# ukbpheno
[![DOI](https://zenodo.org/badge/241869442.svg)](https://zenodo.org/badge/latestdoi/241869442)


![ukbpheno_concept](https://user-images.githubusercontent.com/9621370/170734618-092b9fb3-3a3d-41c3-bcf9-8ccf83a4e860.jpg)


## Description
ukbpheno is an R package for efficiently munging the files provided by UK Biobank to generate data tables of with unified format for further analysis such as making dichotomous phenotypes for UKbio and a composite time-to-event variable combining record level data (HESIN/GP/cancer registry) and main dataset (self reports i.e. nurse interview / touchscreen). Aim of the package is to define binary phenotype data for different types of longitudinal data analysis (e.g. GWAS analysis, cox regressions, baseline tables) in a standardized and reproducible manner. The package can also be used for data exploration with efficient subsetting of the main dataset and visualization functions.

Please check out the [wiki](https://github.com/niekverw/ukbpheno/wiki) for short tutorials on downloading the data as well as usage of the package.

## Installation

```R 
devtools::install_github("niekverw/ukbpheno")
```


## Basic Usage

```R
library(data.table)
library(ukbpheno)

# the directory with datafiles
pheno_dir <-"mydata/ukb99999/"

# main dataset 
fukbtab <- paste(pheno_dir,"ukb99999.tab",sep="")

# meta data file
fhtml <- paste(pheno_dir,"ukb99999.html",sep="")

# hospital inpatient data
fhesin <- paste(pheno_dir,"hesin.txt",sep="")
fhesin_diag <- paste(pheno_dir,"hesin_diag.txt",sep="")
fhesin_oper <- paste(pheno_dir,"hesin_oper.txt",sep="")

# GP data
fgp_clinical <- paste(pheno_dir,"gp_clinical.txt",sep="")
fgp_scripts <- paste(pheno_dir,"gp_scripts.txt",sep="")

# harmonize the data without any definition
lst.harmonized.data<-harmonize_ukb_data(f.ukbtab = fukbtab,f.html = fhtml,f.gp_clinical = fgp_clinical,f.gp_scripts = fgp_scripts,f.hesin = fhesin,f.hesin_diag = fhesin_diag,f.hesin_oper=fhesin_oper,allow_missing_fields = TRUE)
```

## Ascertainment of health outcomes    
Health outcomes are ascertained using data from linkage with national registries (e.g. primary /secondary care) or self report. Full definitions are described in https://bit.ly/3KrMsYD

- Coronary artery disease (doi: 10.1161/CIRCRESAHA.117.312086) 
  - Ischemic heart diseases diagnosis codes 
  - Myocardial infarction diagnosis codes 
  - Coronary Artery Bypass Graft operation codes 
  - Percutaneous Coronary Intervention operation codes 
- Heart failure due to ischemia vs no heart failure after ischemia
    - heart failure among participants with coronary artery disease
    - exclude any participant with cardiomyopathy diagnosis from controls
  

```R
# definition table included in the package 
fdefinitions <- system.file("extdata", "definitions_cardiometabolic_traits.tsv", package="ukbpheno")
# data setting file included in the package
fdata_setting <- system.file("extdata", "data.settings.tsv", package="ukbpheno")
dfData.settings <-fread(fdata_setting)
# process the definition table based on data setting
dfDefinitions_processed_expanded<-read_defnition_table(fdefinitions,fdata_setting,dir.code.map=system.file("extdata", package="ukbpheno"))
# harmonize data
lst.harmonized.data<-harmonize_ukb_data(f.ukbtab = fukbtab,f.html = fhtml,dfDefinitions=dfDefinitions_processed_expanded,f.gp_clinical = fgp_clinical,f.gp_scripts = fgp_scripts,f.hesin = fhesin,f.hesin_diag = fhesin_diag,f.hesin_oper=fhesin_oper,allow_missing_fields = TRUE)

# to identify cases/controls status for CAD  
trait<-"Cad"
df_reference_dt_v0<-lst.harmonized.data$dfukb[,c("identifier","f.53.0.0")]
# read withdrawal list, individuals to be removed from the analysis
f_particip_withdraw<-paste(pheno_dir,"w12345_20210809.csv",sep="")
df_withdrawal<-fread(f_particip_withdraw)
df_reference_dt_v0<-df_reference_dt_v0[! identifier  %in% df_withdrawal$V1]
lst.Cad.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT==trait), lst.harmonized.data$lst.data,dfData.settings, df_reference_date=df_reference_dt_v0)
# summary of diagnosis per participant including case/control status before/after the reference date (baseline visit) as well as the corresponding time-to-event information
View(lst.Cad.case_control$df.casecontrol)

# HF in CAD
trait<-"HfInCad"

# the reference date is the date of diagnosis of CAD
lst.HfInCad.case_control <- get_cases_controls(definitions=dfDefinitions_processed_expanded %>% filter(TRAIT==trait), lst.harmonized.data$lst.data,dfData.settings, vct.identifiers=df_reference_dt_v0$identifier)
```


![dotplot4readme](https://user-images.githubusercontent.com/9621370/151220378-1ade1fa5-8e38-469e-9b9d-aa74138b8be0.png)

Figure above: Relative contribution of different data sources to selected cardiovascular diseases


## Code lookup with shiny app 
Required:
- the code maps (Excel workbook) provided by UK Biobank Showcase Resource 592. 
- R library "optparse" 

```shell
cd ../ukbpheno/inst/util
# show input options 
Rscript shiny.lookup_codes.R --help
# to start the app
Rscript shiny.lookup_codes.R --fcoding_xls path_to_download/all_lkps_maps_v3.xlsx
```



## Citation
If you use ukbpheno, please cite [Yeung, M. W., van der Harst, P., & Verweij, N. (2022). ukbpheno v1.0: An R package for phenotyping health-related outcomes in the UK Biobank. STAR Protocols, 3(3), 101471.](https://star-protocols.cell.com/protocols/1733#article-info)
