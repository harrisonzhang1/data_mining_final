# Data Mining Final: Data Mining to Find Patterns of Long COVID in Electronic Health Records

This repository contains the analysis code and final PDF report of the Data Mining Final Project.

## Project Description
There is a general consensus that much work needs to be done to learn more about long COVID as public health scientists and healthcare providers seek to develop treatment policies for the disease. Since electronic health records (EHRs)—which have information recorded by physicians during routine interactions with patients—contain a wealth of information covering a wide range of patient populations and would not require significant data collection efforts, there is interest in examining whether existing EHR data can be leveraged to find disease patterns in post-COVID patients to better understand the wide spectrum of conditions that emerge after an initial infection. In this analysis, I linked MGB healthcare system EHR data with COVID-19 testing data to first identify COVID-19 patients and then cluster these patients using their post-infection diagnosis codes to identify patterns of potential long COVID. Using K=5 clusters, I was able to identify 5 subpopulations of patients with different combinations of new onset diseases that developed after initial infection.

## File Descriptions
Final_Project_Report.pdf: this PDF document contains the full report of the project. 

Final_Project_Report.docx: this Word document contains the full report of the project. 

analysis_main.R: this contains the R code used to link, clean, explore, analyze, and visualize the datasets. 

step1_data_extraction.R: this contains the R code used to extract data from SQL databases and save them as .csv formats for downstream analysis.

phe_str.Rdata: this contains the mapping file used to add names and body system groups to diagnosis codes used in the project. 

## Data Availability Statement
Data used in this study are only available to those with an approved institutional review board protocol. 

Project instructions: https://leewtai.github.io/courses/data_mining/homeworks/proj3.html
