# ukbpheno

## data types
### Time to event data
- HESIN: ICD9, ICD10, OPCS3, OPCS4. decide how to optionize with primary and secondary diagnoses. 
- Death records, decide how to optionize primary and secondary diagnoses. 
- Self report: non-cancer (20002 + 20009), cancer and operation. Decide how to deal with different answers across visits ( e.g. take the mean  and set event_dt to NA if  >10 years apart.)

### Other data, e.g. yes/no, categorical without dates. 
- medication nurse interview
- touchscreen variables
- abnormal biomarker cutoffs.  



## functions
- Function that reads tsv with definitions 

- Function that converts hes, primary care and self reported tables into a list of dataframes with |n_eid| date| event_dt | event_dur , where event_dur is optional (HESIN).

- Function that pulls out the earliest date that someone is diagnosed before and after the reference date (e.g. visit of interest) for each data-time-to-event data source
  - note that for self report we should only pull out ealiest historical date
  - Input: Disease definition (e.g. code lists of ICD10,ICD9,etc) 
  - Output: list of dataframes per source of diagnosis with 
    1) HXn: 1/0 if disease in history, or number of events for HESIN only
    2) HXd days until event in history 
    3) HXt: time of episode duration in history (for HESIN only, NA otherwise)
    4) FUn: 1/0 if disease in future, or number of events for HESIN only
    5) FUd: days until event in future 
    6) FUt: time of episode duration in history (for HESIN only, NA otherwise)

- Function that pulls out non time-to-event data, e.g. answers that are made during the visit or biomarker cutoffs. 

- Function that merges different data sources into single dataframe, where we can decide that HESIN may cause recurrent events. 
  - output: single dataframe
  
- Function that exports dataframes to STATA format. 

- Functions to wrap above functions. 

- Functions to summarize the overlap between data sources. 





--------

- think about age of first diagnosis. 
- think about censoring dates (e.g. cancer have different censoring dates). can we make further variables that summarize full censoring columns that can be used in cox/kaplan for age-of-diagnosis or (new-onset) event from baseline ? 
- think about how to merge data where you have no days to event. 
- think about more complicated merges, e.g. need at least X codes in hesin or primary care? 



