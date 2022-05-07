rm(list=ls())
rm(list=ls())
library(RODBC)
library(data.table)
library("dplyr")
dir.create("dbmart_extract/")
dir.output="dbmart_extract/"
odbc <- odbcConnect("XXX", "XXX", .rs.api.askForPassword("password"))

##### Select for PCR tested and U07.1 annotated patients to define cohort 
covidpcr = sqlQuery(odbc, paste("SELECT * FROM COVID19_Mart.Analytics.COVIDCohort"))
length(unique(covidpcr$Patient_Num))
covidU07.1 = sqlQuery(odbc, paste("SELECT * FROM COVID19_Mart.Analytics.Diagnoses
                                  WHERE DiagnosisCD IN ('ICD10:U07.1')"))
covidU09.9 = sqlQuery(odbc, paste("SELECT * FROM COVID19_Mart.Analytics.Diagnoses
                                  WHERE DiagnosisCD IN ('ICD10:U09.9')"))
length(unique(covidU07.1$Patient_Num))
length(unique(c(covidpcr$Patient_Num,covidU07.1$Patient_Num)))
save(covidpcr,covidU07.1,covid.patient.num,file=paste0("covid_cohort_",length(covid.patient.num),".RData"))


for(ii in c(1:70,72:115)){
  print(ii)
  if(ii==1){start=1;end=6000
  } else if(ii>1 & ii<115){start=((ii-1)*6000+1);end=(ii*6000)
  } else if(ii==115){start=((ii-1)*6000+1);end=length(covid.patient.num)}
  
  
    dbmart = sqlQuery(odbc,paste0("SELECT f.patient_num, f.encounter_num, f.start_date, c.concept_cd, nval_num, tval_char, 
                                 v.location_cd,  c.concept_path, c.name_char  
                                 FROM COVID19_Mart.RPDR.OBSERVATION_FACT_RPDR f 
                                 INNER JOIN COVID19_Mart.RPDR.VISIT_DIMENSION v ON f.patient_num = v.patient_num AND f.encounter_num = v.encounter_num
                                 INNER JOIN COVID19_Mart.RPDR.CONCEPT_DIMENSION c ON f.concept_cd = c.concept_cd
                                 WHERE f.patient_num IN ('",paste(covid.patient.num[start:end],collapse = "','"),"')"))
    
    dbmart$patient_num <- as.character(dbmart$patient_num)
    dbmart$encounter_num <- as.character(dbmart$encounter_num)
    dbmart$start_date <- as.POSIXct(dbmart$start_date, "%Y-%m-%d")
    dbmart$concept_cd <- as.character(dbmart$concept_cd)
    dbmart$tval_char <- as.character(dbmart$tval_char)
    dbmart$location_cd <- as.character(dbmart$location_cd)
    dbmart$concept_path <- as.character(dbmart$concept_path)
    dbmart$name_char <- as.character(dbmart$name_char)
    dbmart$nval_num <- as.numeric(dbmart$nval_num)
    
    dbmart$concept <- "OTHER"
    dbmart$concept <- ifelse(dbmart$concept_path %like% "\\Diagnosis","Diagnosis",dbmart$concept)
    dbmart$concept <- ifelse(dbmart$concept_path %like% "\\Medication","Medication",dbmart$concept)
    dbmart$concept <- ifelse(dbmart$concept_path %like% "\\Procedures","Procedure",dbmart$concept)
    dbmart$concept <- ifelse(dbmart$concept_path %like% "\\LabTests","LabTests",dbmart$concept)
    dbmart$concept <- ifelse(dbmart$concept_path %like% "\\HealthHistory","HealthHistory",dbmart$concept)
    
    dbmart = dbmart[!duplicated(paste0(dbmart$patient_num,dbmart$encounter_num,dbmart$concept_cd,dbmart$start_date)), ]
    
    saveRDS(dbmart,file=paste0(dir.output,"dbmart_",ii,".RDS"))
    ## for some reason 62 did not work
}





