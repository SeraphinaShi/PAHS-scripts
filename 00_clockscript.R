#===========================
# This R script load in the required data sets for the clocks calculation.
# Author: Seraphina Shi, Lars Van Der Laan
#===========================

library(data.table) #install.packages("data.table")
library(dplyr)
setwd("/Users/seraphinashi/Desktop/Projects/PAHS/Xuanwei Study coal and DNA metylation study")
#---------------------------
# cpg data
#---------------------------

# Store the cpg matrix as a data.table with the first column being the cpg names/ids (usually the rownames but maybe not always)
# This should be the beta/cpg matrix as stored in the RData file.
dat_cpg <- data.table::fread("CGR data/methBeta/RothmanAsianMethStudy7.normBeta.txt")
dim(dat_cpg)
colnames(dat_cpg)[1] <- "CpGName"


# Download the necessary horvath cpgs from the calculator site
horvathCPG <- data.table(CpGName = read.csv("https://dnamage.genetics.ucla.edu/sites/all/files/tutorials/datMiniAnnotation3.csv")$Name)

# This "trick" will subset dat_cpg to only rows/cpgs that appear in horvathCPG and it will add "NA" values for missing cpgs as required for the horvath clock calculator
horvath_data <- dat_cpg[horvathCPG, on = "CpGName"]
 


#---------------------------
# phenotype data
#---------------------------
# pheno is a phenotype file with ids and sheet is a file that maps the phenotype ids to the cpg data ids. 
# The old sheet file may already be merged with the pheno file. In which case, this can be simplified.

pheno_old <- read.csv("Covariates file/LEX_methylation_phenotype_cov_final.csv")
# pheno <- merge(pheno, sheet, by.x = "id", by.y = "id")  # May not be necessary
pheno_old$TargetID <- gsub("X","",pheno_old$TargetID)
pheno_old$grp <- "LEX"

#----------------
# new phenotype files 
pheno_new <- read.csv("Seraphina's folder/data/new_exposure_files/LEX_subjects.csv")
pheno_new$SID <- paste0("CQ ", pheno_new$Visit, pheno_new$SbjctD)
pheno_new$grp <- "LEX"

library("readxl")
pheno_control <- read_excel("Seraphina's folder/data/new_exposure_files/LCW_CleanControls_Methylation_subjectID.xlsx")
pheno_control$SID <- paste("CR", pheno_control$ids_LCW)

childExp <- read.csv("Seraphina's folder/data/new_exposure_files/LEX_fuelCHLD.csv")
childExp$SID <- paste0("CQ ", childExp$Visit, childExp$SbjctD)

cumExp <- read.csv("Seraphina's folder/data/new_exposure_files/LEX_fuelCUM.csv")
cumExp$SID <- paste0("CQ ", cumExp$Visit, cumExp$SbjctD)

fuelStove <- read.csv("Seraphina's folder/data/new_exposure_files/LEX_fuelstoveCUR.csv")
fuelStove$SID <- paste0("CQ ", fuelStove$Visit, fuelStove$SbjctD)

pheno <- merge(pheno_old, pheno_new, by = "SID",  all = T)
pheno <- merge(pheno, pheno_control, by = "SID", all = T)
pheno <- merge(pheno, childExp %>% select(brthFuel, SID), by = "SID", all = T)
pheno <- merge(pheno, cumExp %>% select(cumFuel, SID), by = "SID", all = T)
pheno <- merge(pheno, fuelStove %>% select(curFuel, curStove, SID), by = "SID", all = T)


pheno <- pheno %>% select(TargetID, Sample_ID, Sample_Plate, SID, SbjctD, SbjctD_LEX, Visit, VISIT,
                          ids_LCW, grp.x, grp.y, grp, CGR_Internal_CONTROL, Replicate, Failed_CGR_QC,
                          age.x, age.y, Gender, BMI, bmi, ses, edu, education_cat, 
                          county, loc_N, loc_E, dist_min, inpatient_sms, diagnosis_code_sms, enroll_age_sms, diagnosis_eng_sms, 
                          econ_any_items, curr_commune,
                          brthFuel, cumFuel, curFuel, fuel, fueltype, curStove, stovetype, laibin_coal, coal_area,
                          everything())

# clean up group
pheno$Group <- pheno$grp
pheno$Group[is.na(pheno$Group)] <- pheno$grp.x[is.na(pheno$Group)]
pheno$Group[is.na(pheno$Group)] <- pheno$grp.y[is.na(pheno$Group)]
pheno <- pheno %>% select(-c(grp, grp.x, grp.y))

# clean up age
# The pheno file should have a column called exactly "Age" with the chronological ages
identical(pheno$age.x[(! is.na(pheno$age.x)) & (! is.na(pheno$age.y))], pheno$age.y[(! is.na(pheno$age.x)) & (! is.na(pheno$age.y))])
pheno$age <- pheno$age.x
pheno$age[is.na(pheno$age)] <- pheno$age.y[is.na(pheno$age)]
pheno$age[is.na(pheno$age)] <- pheno$enroll_age_sms[is.na(pheno$age)]
pheno <- pheno %>% select(-c(age.x, age.y, enroll_age_sms))

# clean up BMI
pheno$BMI[is.na(pheno$BMI)] <- pheno$bmi[is.na(pheno$BMI)]
pheno <- pheno %>% select(-bmi)

 
# remove empty columns
empty_columns <- colSums(is.na(pheno) | pheno == "") == nrow(pheno)
pheno <- pheno[, !empty_columns]

# clean up fuel types and stove types
pheno$brthFuel[pheno$brthFuel == 1] = "Wood"
pheno$brthFuel[pheno$brthFuel == 2] = "Smokeles"
pheno$brthFuel[pheno$brthFuel == 3] = "Mix"
pheno$brthFuel[pheno$brthFuel == 4] = "Smoky"
pheno$brthFuel[pheno$Group == "LCW Clean Ctrl"] = "Gas/Elect"

pheno$cumFuel[pheno$cumFuel == 1] = "Mix"
pheno$cumFuel[pheno$cumFuel == 2] = "Smoky"
pheno$cumFuel[pheno$Group == "LCW Clean Ctrl"] = "Gas/Elect"

pheno$curFuel[pheno$curFuel == 1] = "Wood_and_or_Plant"
pheno$curFuel[pheno$curFuel == 2] = "Smokeles"
pheno$curFuel[pheno$curFuel == 3] = "Smoky"
pheno$curFuel[pheno$Group == "LCW Clean Ctrl"] = "Gas/Elect"

pheno$curStove[pheno$curStove == 1] = "Firepit_and_unventilated"
pheno$curStove[pheno$curStove == 2] = "Ventilated"
pheno$curStove[pheno$curStove == 3] = "Portable_stove"
pheno$curStove[pheno$curStove == 4] = "Mix"
pheno$curStove[pheno$Group == "LCW Clean Ctrl"] = "Gas/Elect"


pheno <- pheno %>% select(TargetID, Sample_ID, Sample_Plate, SID, SbjctD, Visit, 
                          ids_LCW, Group, CGR_Internal_CONTROL, Replicate, Failed_CGR_QC,
                          age, Gender, BMI, ses, edu, education_cat, 
                          county, loc_N, loc_E, dist_min, inpatient_sms, diagnosis_code_sms, diagnosis_eng_sms, 
                          econ_any_items, curr_commune,
                          brthFuel, cumFuel, curFuel, fuel, fueltype, curStove, stovetype, laibin_coal, coal_area,
                          everything())

# The pheno file should have a column called exactly "Female" that is 1 for "female" and 0 for "male".
pheno$Female <- as.numeric(pheno$Gender == "F") 

# The pheno file should have a column called exactly "Age" with the chronological ages
pheno$Age <- pheno$age

# remove controls that are not clean
# pheno <- pheno[!(grepl("CR" ,pheno$SID) & pheno$Group == "LEX"), ]

# remove males
# dat <- dat[dat$Gender != "M", ]




keep_cols <- intersect(colnames(horvath_data), pheno$TargetID) # Find common observations between pheno and beta mat
horvath_data <- horvath_data[,c("CpGName", keep_cols), with = F] # If there are extra observations not found in both then remove them.
pheno <- pheno[na.omit(match(keep_cols, pheno$TargetID)), ] # Make sure row order of phenofile matches the column order of beta/cpg matrix (check ids)

# missing values should be fine, careful 

keep <- which(! is.na(pheno$Age))
pheno <- pheno[keep, ]
horvath_data <- horvath_data[, c(1, keep+1), with=F]

#---------------------------
# output data
#---------------------------
write.csv(horvath_data, "Seraphina's folder/data/0.calculate_clocks/horvathcalculator_beta_input.csv", row.names = F) # To be fed into online calculator
write.csv(pheno, "Seraphina's folder/data/0.calculate_clocks/horvathcalculator_pheno_input.csv", row.names = F)# To be fed into online calculator
 
