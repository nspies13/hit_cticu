##### Set Up Pipeline Options #####
library(targets)
library(tarchetypes)
library(tidyverse)

# Set target options
tar_option_set(
  packages = c("tidyverse", "qs", "here"),
  format = "qs",
  deployment = "main"
)

# Set Figure Defaults
theme_ns <- theme(text = element_text(family = "Helvetica"),
                  title = element_text(face = "bold", size = 14),
                  plot.title = element_text(hjust = 0.5),
                  axis.title = element_text(size = 10, face = "italic"),
                  legend.title = element_text(face = "bold", size = 12),
                  axis.line = element_line(),
                  axis.ticks = element_blank(),
                  axis.text = element_text(size = 8, face = "italic"),
                  panel.grid = element_blank(),
                  panel.background = element_blank(),
                  strip.text = element_text(size = 8, face = "italic"),
                  strip.background = element_blank())
theme_set(theme_ns)

color_map <- c("Negative" = "gray60", "Positive" = "#8B0000", "See Interpretation" = "#8B0000")
alpha_map <- c("Negative" = 1, "Positive" = 1, "See Interpretation" = 0.8)

# Set up parallel processing
library(future.callr)
plan(callr, workers = 32)

# Load Helper Functions
tar_source()

# Define Pipeline
##### Load Data #####
load <- list(
  
  tar_target(admissions, read_delim("../Data/29851 ZaydmanSpies -Heparin-induced Thrombocytopenia 20220525 Admissions.txt")),
  tar_target(allergies, read_delim("../Data/29851 ZaydmanSpies -Heparin-induced Thrombocytopenia 20220525 Allergens.txt")),
  tar_target(demographics, read_delim("../Data/29851 ZaydmanSpies -Heparin-induced Thrombocytopenia 20220525 Demographics.txt")),
  tar_target(lda, read_delim("../Data/29851 ZaydmanSpies -Heparin-induced Thrombocytopenia 20220525 LinesDrainsAirwaysWounds.txt")),
  tar_target(measurements, read_delim("../Data/29851 ZaydmanSpies -Heparin-induced Thrombocytopenia 20220525 Measurements.txt")),
  tar_target(pmh, read_delim("../Data/29851 ZaydmanSpies -Heparin-induced Thrombocytopenia 20220525 Medical History.txt")),
  tar_target(diagnoses, read_delim("../Data/29851 ZaydmanSpies -Heparin-induced Thrombocytopenia 20220525 Prior 6 months Diagnoses.txt")),
  tar_target(labs, read_delim("../Data/29851 ZaydmanSpies -Heparin-induced Thrombocytopenia 20220525 Prior 6 months Labs.txt")),
  tar_target(mar, read_delim("../Data/29851 ZaydmanSpies -Heparin-induced Thrombocytopenia 20220525 Prior 6 months Medication Administration.txt")),
  tar_target(med_orders, read_delim("../Data/29851 ZaydmanSpies -Heparin-induced Thrombocytopenia 20220525 Prior 6 months Medication Orders.txt")),
  tar_target(procedures, read_delim("../Data/29851 ZaydmanSpies -Heparin-induced Thrombocytopenia 20220525 Prior 6 months Procedures.txt")),
  tar_target(order_times, read_delim("../Data/HIT_order_times.txt") %>% rename(Epic_MRN = MRN))
    
)

##### Exploratory Data Analysis #####
eda <- list(
  
  tar_target(demo_input, convertDemographics(demographics))
  
)

##### Find HIT Work-ups ##### 
tested <- list(
  
  tar_target(hit_labs, labs %>% filter(MED_PT_TYPE == "Hospital Encounter" & TEST_NAME %in% c("HIT ANTIBODIES W/REFLEX TO SEROTONIN RELEASE ASSAY (SRA)", "SEROTONIN RELEASE ASSAY", "HEPARIN PLATELET FACTOR 4 (HIT) ANTIBODY"))),
  tar_target(hit_patients, hit_labs %>% select(Epic_MRN) %>% distinct() %>% pluck(1)), 
  tar_target(hit_mar, mar %>% filter(Epic_MRN %in% hit_patients & grepl("HEPARIN|BIVAL|FONDA|LOVEN", MED_NAME))),
  tar_target(hit_med_orders, med_orders %>% filter(Epic_MRN %in% hit_patients & (grepl("HEPARIN|BIVAL|FONDA|LOVEN", MED_NAME) | (THERAPEUTIC_CLASS_NAME == "ANTICOAGULANTS")))),
  tar_target(hit_patient_diagnoses, diagnoses %>% filter(Epic_MRN %in% hit_patients))
  
)

##### Capture Decision Making #####
decisions <- list(
  
  tar_target(lia_results, getLIAresults(hit_labs)),
  tar_target(sra_results, getSRAresults(hit_labs)),
  tar_target(combined_results, combineResults(lia_results, sra_results)),
  tar_target(order_decisions, getMARchanges(combined_results, hit_mar, hit_labs))

)

##### Make Figures #####
figures <- list(
  
  #tar_target(table_1, makeDemographicsTable(demographics)),
  tar_target(lia_positivity_rates, makePositivityRatesFigures(combined_results, order_times)),
  tar_target(lia_repeat, makeRepeatLIAFigures(combined_results, lia_results)),
  tar_target(sra_repeat, makeRepeatSRAFigures(combined_results, sra_results)),
  tar_target(sra_by_lia, makeSRAbyLIAplot(combined_results))
  
)

##### Run Pipeline #####
list(load, eda, tested, figures, decisions)
