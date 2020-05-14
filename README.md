# ukbpheno

`run.r` is the base script that I use to test things; functions in different scripts, can easily convert it into R package. 

Currently it only reads in the data and converts data into dataframes with 4-5 colums: Identifier, event-date, event-code, event (yes/no), episode-duration (optional) 


## data types
### Time to event data
- HESIN: ICD9, ICD10, OPCS3, OPCS4. decide how to optionize with primary and secondary diagnoses. 
- Primary care: clinical (ReadV2, CTV3) scripts(readv2, BNF, DMD) 
- Death records 
- Self report: non-cancer (20002 + 20009), cancer and operation. 
- medication nurse interview (only first event available)

### Other data, e.g. yes/no, categorical without dates. 
- touchscreen variables
- abnormal biomarker cutoffs? 



## 
Given a visit-date or other date (e.g. date of  diagnosis) and some codes, calculate:
- days to first event in history
- days to first event in future
- if participant is positive +/- 5? days from the visit (e.g. relevant to see if individual is on medication)
- days to death (primary cause)
- days to death (secondary cause)
- some number on severity? 
- some number about if data looks strange? 


--------

-  age of first diagnosis. 
-  censoring dates (e.g. cancer have different censoring dates). can we make further variables that summarize full censoring columns that can be used in cox/kaplan for age-of-diagnosis or (new-onset) event from baseline ? 
-  how to merge data where you have no days to event. 
-  more complicated queries, e.g. need at least X codes for primary care to be considered? 
-  more complex combinations e.g. specific combination of codes like operation + diagnosis.  
-  more complex phenotype schemes, e.g. exclude asthma from COPD cases. 
-  scoring cases/intersections e.g. to define high confidence cases vs low confidence vs controls. 
- functions that can help with defining a disease, e.g. suggest codes based on others. 

