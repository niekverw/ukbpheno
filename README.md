# ukbpheno

## Description
ukbpheno is an R package for efficiently munging the files provided by UK Biobank to generate data tables of with unified format for further analysis such as making dichotomous phenotypes for UKbio and a composite time-to-event variable combining record level data (HESIN/GP/cancer registry) and main dataset (self reports i.e. nurse interview / touchscreen). The package can also be used to efficiently extract required columns from a huge main dataset for exploration which is aided by visualization functionalities included in the package. 


## Installation

```R 
devtools::install("niekverw/ukbpheno")
```


## Usage

```R
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

# harmonize the data
lst.harmonized.data<-harmonize_ukb_data(f.ukbtab = fukbtab,f.html = fhtml,f.gp_clinical = fgp_clinical,f.gp_scripts = fgp_scripts,f.hesin = fhesin,f.hesin_diag = fhesin_diag,f.hesin_oper=fhesin_oper,allow_missing_fields = TRUE)

```
![dotplot4readme](https://user-images.githubusercontent.com/9621370/151220378-1ade1fa5-8e38-469e-9b9d-aa74138b8be0.png)

## Code lookup with shiny app 
Required: the code maps (Excel workbook) provided by UK Biobank Showcase Resource 592. 

```shell
cd ../ukbpheno/inst/util
Rscript shiny.lookup_codes.R --fcoding_xls path_to_download/all_lkps_maps_v3.xlsx
```
