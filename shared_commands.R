#========================
# This file loads libraries and defines common variables 
#========================

# load library
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gtsummary)
library(grid)
library(gridExtra)
library(stringr)

# 6 clock names
clock_names <- c("Horvath", "Hannum", 
                 "PhenoAge","Skin&Blood", 
                 "GrimAge", "DNAmTL",
                 "IEAA", "EEAA")
clock_names_extra <- c("DNAmAge Horvath (years)", "DNAmAge Hannum (years)", 
                       "DNAm PhenoAge (years)","DNAMaGE Skin&Blood (years)",
                       "DNAm GrimAge (years)", "DNAm TL (kb)")

# define 6 clocks 
clocks_vars <- c("DNAmAge", "DNAmAgeHannum", 
                 "DNAmPhenoAge", "DNAmAgeSkinBloodClock", 
                 "DNAmGrimAge", "DNAmTL")

# define 8 epigenetic age accelerations 
EAAs_names <- c("Horvath EAA", "Hannum EAA", "PhenoAge EAA", 
                "Skin&Blood EAA", "GrimAge EAA", "DNAmTL",
                "IEAA", "EEAA")

EAAs_vars <- c("AgeAccelerationResidual", "AgeAccelerationResidualHannum", "AgeAccelPheno",
          "DNAmAgeSkinBloodClockAdjAge", "AgeAccelGrim", "DNAmTLAdjAge", "IEAA", "EEAA")



# define colors for each clock
colors_clock <- brewer.pal(n = 8, name = 'Dark2')
names(colors_clock) <- EAAs_names


# exposure variables names
ambient_exp_names <- c("bap", "pm25", "ANY", "BPE", "BaA", "BbF", "BkF", 
                       "CHR", "DBA", "FLT", "FLU", "IPY", "NAP", "PHE", "PYR")
ambient_exp_vars <- paste0(ambient_exp_names, "_air")

ambient_exp_label <- c("Indoor benzo[a]pyrene (BaP) concentration",
                        "Indoor PM2.5 concentration",
                        "Indoor acenaphthylene (ANY) concentration",
                        "Indoor benzo[ghi]-perylene (BPE) concentration",
                        "Indoor benzo[a]anthracene (BaA) concentration",
                        "Indoor benzo[b]fluoranthene (BbF) concentration",
                        "Indoor benzo[k]flouranthene (BkF) concentration",
                        "Indoor chrysene(CHR) concentration",
                        "Indoor dibenz(ah)anthracene(DBA) concentration",
                        "Indoor fluoranthene (FLT) concentration",
                        "Indoor fluorine (FLU) concentration",
                        "Indoor indeno(1,2,3-cd)pyrene (IPY) concentration",
                        "Indoor naphthalene (NAP) concentration",
                        "Indoor phenanthrene (PHE) concentration",
                        "Indoor pyrene(PYR) concentration")

urinary_exp_names <- c("Benzanthracene_Chrysene", "Naphthalene", "2.Methylnaphthalene", 
                       "1.Methylnaphthalene", "Acenaphthene", "Phenanthrene_Anthracene", 
                       "Fluoranthene", "Pyrene")
urinary_exp_vars <- paste0(c("Benzanthracene_Chrysene", "Naphthalene", "Methylnaphthalene_2", 
                             "Methylnaphthalene_1", "Acenaphthene", "Phenanthrene_Anthracene", 
                             "Fluoranthene", "Pyrene"), "_urine")

urinary_exp_label <- c("Benzo(a)anthracene and Chrysene concentration",
                       "Naphthalene concentration",
                       "2-Methylnaphthalene concentration",
                       "1-Methylnaphthalene concentration",
                       "Acenaphthene concentration",
                       "Phenanthrene and Anthracene concentration",
                       "Fluoranthene concentration",
                       "Pyrene concentration")




var_cat <- c("predictedGender", 
             "curFuel", "brthFuel", "cumFuel", "fueltype", "fuel",
             "curStove",  "stovetype", 
             "education_cat", "county")

var_cont <- c("age","BMI")



var_name_list <- list("ID" = c("SID", "SbjctD", "Visit"),
                      "other_info" = c("Group", "Gender"), 
                      "clocks" = clocks_vars,
                      "EEAs" = EAAs_vars,
                      "confounders" = c("Age", "county", "BMI", "ses", "edu"),
                      "fuel_exp" = c("curFuel", "brthFuel", "cumFuel"),
                      "stove_exp" = "curStove", 
                      "ambient_exp" = ambient_exp_vars,
                      "urinary_exp" = urinary_exp_vars,
                      "cluster_exp" = list("clusCUR6" = c("CUR6_BC_PAH6", "CUR6_PAH31", "CUR6_NkF", "CUR6_PM_RET",  "CUR6_NO2", "CUR6_SO2"),
                                           "clusCHLD5" = c("CHLD5_X7",  "CHLD5_X33", "CHLD5_NkF", "CHLD5_NO2", "CHLD5_SO2"),
                                           "clusCUM6" = c("CUM6_BC_NO2_PM", "CUM6_PAH36", "CUM6_DlP",  "CUM6_NkF" ,"CUM6_RET", "CUM6_SO2"),
                                           "clusMEAS6" = c("MEAS6_BC_PM_RET", "MEAS6_X31","MEAS6_X5","MEAS6_DlP","MEAS6_NkF","MEAS6_NO2_SO2"),
                                           "clusURI5" = c("URI5_NAP_1M_2M", "URI5_ACE","URI5_FLU_PHE", "URI5_PYR","URI5_CHR")))

