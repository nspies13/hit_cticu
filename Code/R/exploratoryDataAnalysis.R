convertDemographics <- function(demographics = demographics){
  
  demo_input <- 
    demographics %>% 
    mutate(AGE = as.numeric(difftime(as_datetime("2022-03-01"), BIRTH_DATE, units = "days") / 365.25)) %>%
    mutate(RACE = case_when(RACE == "White" ~ "White", RACE == "Black or African American" ~ "Black", T ~ "Other/Missing")) %>%
    mutate(across(c(ETHNICITY, SEX, RACE), ~as_factor(.)))
  
  demo_input
  
}

makeDemographicsTable <- function(demo_input = demo_input){
  
  library(gtsummary)
  
  demo_table <- 
    demo_input %>% select(AGE, SEX, RACE) %>% tbl_summary(missing = "no") %>%
    bold_labels() %>% italicize_labels() %>% modify_header(label = "**Demographic**") %>% 
    modify_footnote(all_stat_cols() ~ "Median (IQR) for Age, otherwise N (%)")
  
  demo_table
  
}

makeLabResultsTable <- function(hit_labs = hit_labs, demo_input = demo_input){
  
  lab_table_input <- 
    hit_labs %>% 
      group_by(Epic_MRN, COMPONENT_NAME) %>% 
      arrange(COLLECTION_DTTM) %>% 
      slice_head(n = 1) %>%
      filter(COMPONENT_NAME %in% c("HIT AB", "UNFRACTIONATED HEPARIN RESULT")) %>%
      select(Epic_MRN, COMPONENT_NAME, LAB_RSLT_VALUE) %>% pivot_wider(c(Epic_MRN), names_from = COMPONENT_NAME, values_from = LAB_RSLT_VALUE) %>%
      rename(LIA_RESULT = `HIT AB`, SRA_RESULT = `UNFRACTIONATED HEPARIN RESULT`) %>% 
      left_join(demo_input, .) %>% 
      mutate(LIA_TESTED = ifelse(is.na(LIA_RESULT), F, T), SRA_TESTED = ifelse(is.na(SRA_RESULT), F, T))
    
  gt_tested <- 
    lab_table_input %>%
      select(AGE, SEX, RACE, LIA_TESTED) %>% tbl_summary(missing = "no", by = LIA_TESTED) %>% add_p() %>%
      bold_labels() %>% italicize_labels() %>% italicize_levels() %>% modify_header(label = "", stat_1 = "*Tested*", stat_2 = "*Not Tested*", p.value = "*p-value*") %>% 
      modify_spanning_header(all_stat_cols() ~ "***LIA Ordered***")
    
  gt_results <- 
    lab_table_input %>% 
      select(AGE, SEX, RACE, LIA_RESULT) %>% tbl_summary(missing = "no", by = LIA_RESULT) %>% add_p() %>% bold_p() %>%
      bold_labels() %>% italicize_labels() %>%  italicize_levels() %>% modify_header(label = "", stat_1 = "*Negative*", stat_2 = "*Positive*", p.value = "*p-value*") %>% 
      modify_spanning_header(all_stat_cols() ~ "***LIA Result***") 
  
  gt_sra <- 
    lab_table_input %>% 
      filter(SRA_RESULT != "See Interpretation") %>%
      select(AGE, SEX, RACE, SRA_RESULT) %>% tbl_summary(missing = "no", by = SRA_RESULT) %>% add_p() %>%
      bold_labels() %>% italicize_labels() %>% italicize_levels() %>%  modify_header(label = "", stat_1 = "*Negative*", stat_2 = "*Positive*", p.value = "*p-value*") %>% 
      modify_spanning_header(all_stat_cols() ~ "***SRA Result***")
  
  table_1 <- tbl_merge(list(gt_tested, gt_results, gt_sra), tab_spanner = c("***LIA Ordered***", "***LIA Results***", "***SRA Results***"))
    
  table_1
  
}
