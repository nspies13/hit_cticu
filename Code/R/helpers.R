plotTAT <- function(hit_labs){
  
  gg_input <- 
    hit_labs %>% filter(COMPONENT_NAME == "HIT AB LIA UNITS") %>% 
      transmute(TAT = as.numeric(difftime(LAB_RSLT_DTTM, COLLECTION_DTTM, units = "mins")), RESULT = ifelse(LAB_RSLT_VALUE < 1, "Negatve", "Positive")) %>%
      mutate(TAT = ifelse(TAT > 480, 480, TAT))
 
  quantile(gg_input$TAT, probs = c(.50, .95, .995))
  
  ggplot(gg_input, aes(TAT)) + 
    stat_ecdf(geom = "step") + 
      annotate(geom = "rect", color = "red", fill = NA, ymin = -1000, ymax = 0.5, xmin = -1000, xmax = quantile(gg_input$TAT, probs = c(.50))[[1]]) +
      annotate(geom = "rect", color = "darkred", fill = NA, ymin = -1000, ymax = 0.95, xmin = -1000, xmax = quantile(gg_input$TAT, probs = c(.95))[[1]]) +
      coord_cartesian(xlim = c(0,240), ylim = c(0,1))
    
   
}

getMARchanges <- function(combined_results, hit_mar, hit_labs){
  
  lia <- 
    hit_labs %>% 
      filter(COMPONENT_NAME == "HIT AB LIA UNITS") %>%
      arrange(COLLECTION_DTTM) %>% 
      group_by(Epic_MRN) %>% 
      slice_head(n = 1) %>% 
      select(Epic_MRN, COLLECTION_DTTM, LAB_RSLT_DTTM, LAB_RSLT_VALUE) %>%
      mutate(BINARY = ifelse(LAB_RSLT_VALUE < 1, "Negative", "Positive"))
  
  mar <- left_join(hit_mar, lia, by = "Epic_MRN")
  
  data <- 
    mar %>% 
      arrange(MED_ADMIN_DOSE_START_DATE, MED_ADMIN_DOSE_STOP_DATE) %>% 
      transmute(
        Patient = as_factor(Epic_MRN),
        Medication = MED_NAME,
        Order = ORDER_ID,
        Start = as.numeric(difftime(MED_ADMIN_DOSE_START_DATE, COLLECTION_DTTM, units = "hours")), 
        Stop = as.numeric(difftime(MED_ADMIN_DOSE_STOP_DATE, COLLECTION_DTTM, units = "hours"))) %>% 
      filter(!grepl("SYRINGE|FLUSH", Medication)) %>% 
      mutate(Medication = case_when(
        grepl("HEPARIN", Medication) ~ "Heparin",
        grepl("BIVAL", Medication) ~ "Bivalirudin",
        grepl("FONDA", Medication) ~ "Fondaparinux"))
  
  rows <- seq(-24, 24)
  
  gg_input <- 
    expand_grid(Patient = unique(data$Patient), Time = rows) %>%
      left_join(data) 
  
  gg_in <- 
    gg_input %>% 
      mutate(Med = ifelse(Start < Time & Time < Stop, Medication, NA))
  
  gg_input <- 
    gg_in %>% 
      distinct(Patient, Time, Med) %>%
      group_by(Patient, Time) %>%
      na.omit()
  
  gg_order <- gg_input %>% pivot_wider(names_from = Time, values_from = Med, values_fn = last, values_fill = "None")

  hep_counts <- apply(gg_order, 1, function(x) sum(grepl("epari", x)))
  bival_counts <- apply(gg_order, 1, function(x) sum(grepl("ival", x)))
  
  gg_order_in <-
    bind_cols(gg_order, hep_counts = hep_counts, bival_counts = bival_counts) %>%
      arrange(bival_counts, desc(hep_counts)) %>%
      mutate(Patient = factor(Patient), Patient = fct_reorder2(Patient, bival_counts, -hep_counts))
  
  anyHepPrior <-
    gg_order_in %>% 
      transmute(across(matches("-"), ~grepl("Heparin", .x))) %>% 
      apply(1, function(x) any(x))
  
  anyHepPost <-
    gg_order_in %>% 
      transmute(across(c('0','1','2','3','4','5','6','7','8'), ~grepl("Heparin", .x))) %>% 
      apply(1, function(x) any(x))
  
  anyBivalPrior <-
    gg_order_in %>% 
      transmute(across(matches("-"), ~grepl("Bival", .x))) %>% 
      apply(1, function(x) any(x))
  
  anyBivalPost <-
    gg_order_in %>% 
      transmute(across(c('0','1','2','3','4','5','6','7','8'), ~grepl("Bival", .x))) %>% 
      apply(1, function(x) any(x))
  
  nonePost <- 
    gg_order_in %>% 
      transmute(across(c('0','1','2','3','4','5','6','7','8'), ~grepl("None", .x))) %>%
      apply(1, function(x) all(x[2:10]))
    
  priorHepCount <- sum(anyHepPrior, na.rm = T)
  postHepCount <- sum(anyHepPost, na.rm = T)
  
  priorBivalCount <- sum(anyBivalPrior, na.rm = T)
  postBivalCount <- sum(anyBivalPost, na.rm = T)
  
  nonePostCount <- sum(nonePost, na.rm = T)
  
  sum(anyHepPrior & (nonePost | anyBivalPost), na.rm = T)
  
  lia_timings <- 
    lia %>%
      mutate(Patient = factor(Epic_MRN), LIA_time = as.numeric(difftime(LAB_RSLT_DTTM, COLLECTION_DTTM, units = "hours")))
  
  gg_long <- 
    gg_order_in %>%
      pivot_longer(cols = !matches("Patient|count|Order"), names_to = "Time", values_to = "Med") %>%
      left_join(lia_timings %>% select(Patient, LIA_time, BINARY)) %>%
      mutate(Time = as.numeric(Time))
    
  order <- reorder(gg_long$Patient, gg_long$bival_counts)
  
  gg_long <- 
    gg_long %>% 
      mutate(Patient = factor(Patient, levels = unique(order)))
  
  gg_long %>%
    ggplot(aes(Time, fct_rev(Patient), color = Med)) + 
      geom_line(linewidth = 0.9) +
      geom_vline(xintercept = 0) + 
      scale_color_manual(values = c("None" = "white", "Heparin" = "grey60", "Bivalirudin" = "#da0000")) +
      coord_cartesian(xlim = c(-24, 12)) + 
      scale_x_continuous(breaks = seq(-24, 12, by = 12)) + 
      xlab("Hours From LIA Order") + ylab("Patient") + 
      annotate("text", x = 10, y = 8, label = "Bivalirudin", fontface = "bold.italic", size = 6) +
      annotate("text", x = 10, y = 140, label = "Heparin", fontface = "bold.italic", size = 6) +
      annotate("text", x = 10, y = 55, label = "None", fontface = "bold.italic", size = 6) +
      theme(axis.title = element_text(size = 12), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = "none", axis.line.y = element_blank())
  ggsave("../Figures/order_timeline.pdf", height = 5, width = 8)
  
  sra_pos <- combined_results %>% filter(TEST_NUM == 1 & TEST_NAME == "SRA" & BINARY == "Positive") %>% select(Epic_MRN) %>% pluck(1) 
  
  gg_long_results <- 
    gg_long %>% 
      mutate(Patient = fct_reorder(Patient, hep_counts))
  
  gg_pos <- 
    gg_long_results %>%
      filter(BINARY == "Positive") %>%
      ggplot(aes(Time, fct_rev(Patient), color = Med)) + 
      geom_line(linewidth = 1.05) +
      geom_point(aes(LIA_time, fct_rev(Patient)), color = "black", size = 0.5, shape = 20) +
      geom_vline(xintercept = 0) + 
      scale_color_manual(values = c("None" = "white", "Heparin" = "grey60", "Bivalirudin" = "#da0000")) +
      coord_cartesian(xlim = c(-24, 12)) + 
      scale_x_continuous(breaks = seq(-24, 12, by = 12)) + 
      annotate("text", x = 9, y = 5, label = "Bivalirudin", fontface = "bold.italic", size = 4) +
      ggtitle("Positive") + xlab("Hours From LIA Order") + ylab("Patient") + 
      theme(title = element_text(size = 12, face = "italic"), axis.title.y = element_blank(), axis.text.y = element_blank(), strip.text = element_text(size = 8), axis.ticks = element_blank(), legend.position = "none", axis.line.y = element_blank())
  
  gg_neg <- 
    gg_long_results %>%
      filter(BINARY == "Negative") %>%
      ggplot(aes(Time, fct_rev(Patient), color = Med)) + 
      geom_line(linewidth = 1) +
      geom_point(aes(LIA_time, fct_rev(Patient)), color = "black", size = 0.5, shape = 20) +
      geom_vline(xintercept = 0) + 
      scale_color_manual(values = c("None" = "white", "Heparin" = "grey60", "Bivalirudin" = "#da0000")) +
      coord_cartesian(xlim = c(-24, 12)) + 
      scale_x_continuous(breaks = seq(-24, 12, by = 12)) + 
      annotate("text", x = 9, y = 110, label = "Heparin", fontface = "bold.italic", size = 4) +
      annotate("text", x = 9, y = 30, label = "None", fontface = "bold.italic", size = 4) +
      ggtitle("Negative") + xlab("Hours From LIA Order") + ylab("Patient") + 
      theme(title = element_text(size = 12, face = "italic"), axis.title = element_blank(), axis.text = element_blank(), strip.text = element_text(size = 8), axis.ticks = element_blank(), legend.position = "none", axis.line.y = element_blank())
    
  ggpubr::ggarrange(gg_neg, gg_pos, nrow = 2, ncol = 1, heights = c(0.75, 0.25), labels = "AUTO")
  ggsave("../Figures/order_timeline_by_lia.pdf", height = 6, width = 4)
  
  
}

plotMARchanges <- function(mar_nested = mar_nested){

  MRN_order <- mar_nested %>% arrange(heparin_prior_12, heparin_prior_24, heparin_prior_48, bival_prior_12, bival_prior_24, bival_prior_48, heparin_next_12, heparin_next_24, heparin_next_48, bival_next_12, bival_next_24) %>% select(Epic_MRN)
  
  mar_gg_input_heparin <- pivot_longer(mar_nested %>% select(!matches("bival"), -data), cols = matches("heparin"), names_to = "TIME", values_to = "HEPARIN") %>% mutate(TIME = gsub("heparin_", "", TIME))
  mar_gg_input_bival <- pivot_longer(mar_nested %>% select(!matches("heparin"), -data), cols = matches("bival"), names_to = "TIME", values_to = "BIVAL") %>% mutate(TIME = gsub("bival_", "", TIME))
  
  mar_gg_input <- 
    full_join(mar_gg_input_heparin, mar_gg_input_bival) %>% 
    mutate(Epic_MRN = factor(Epic_MRN, levels = MRN_order[[1]]), TIME = factor(TIME, levels = c("prior_48", "prior_24", "prior_12", "next_12", "next_24", "next_48"))) %>%
    mutate(HEPARIN = replace_na(HEPARIN, FALSE), BIVAL = replace_na(BIVAL, FALSE)) %>% 
    mutate(ANTICOAGULATION = case_when(
      HEPARIN & BIVAL ~ "BOTH",
      HEPARIN & !BIVAL ~ "HEPARIN",
      !HEPARIN & BIVAL ~ "BIVAL",
      !HEPARIN & !BIVAL ~ "NEITHER"
    )) 
  
  ggplot(mar_gg_input, aes(x = TIME, y = Epic_MRN, fill = ANTICOAGULATION)) + 
    geom_tile() + 
    scale_fill_manual(values = c("BOTH" = "gray50", "HEPARIN" = viridis(9, option = "E")[[3]], "BIVAL" = viridis(11, option = "E")[[8]], "NEITHER" = "gray70", `NA` = "gray50")) + 
    scale_x_discrete(limits = c("prior_12", "next_12", "next_24", "next_48")) + 
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
  
  
}

getOrderChanges <- function(trigger = hit_collect, orders = hit_med_orders){
  
  orders <- left_join(orders, trigger, by = "Epic_MRN") %>% mutate(time_until_HIT_order = difftime(COLLECTION_DTTM, MED_ORDER_DATE, units = "hours"))
  
  orders_within_24 <- orders %>% filter(abs(time_until_HIT_order) < 24)
  
  orders_prior_48h <- 
    orders %>% 
      distinct %>% 
        filter(difftime(MED_ORDER_DATE, COLLECTION_DTTM, units = "hours") < 48 & MED_ORDER_DATE < COLLECTION_DTTM) 
  
  orders_prior_24h <- 
    orders %>% 
      distinct %>% 
        filter(difftime(MED_ORDER_DATE, COLLECTION_DTTM, units = "hours") < 24 & MED_ORDER_DATE < COLLECTION_DTTM) 
  
  orders_prior_12h <- 
    orders %>% 
      distinct %>% 
       filter(difftime(MED_ORDER_DATE, COLLECTION_DTTM, units = "hours") < 12 & MED_ORDER_DATE < COLLECTION_DTTM) 
  
  orders_next_48h <- 
    orders %>% 
     distinct %>% 
      filter(difftime(MED_ORDER_DATE, COLLECTION_DTTM, units = "hours") < 48 & difftime(MED_ORDER_DATE, COLLECTION_DTTM, units = "hours") > 24 & MED_ORDER_DATE > COLLECTION_DTTM) 
    
  orders_next_24h <- 
    orders %>% 
      distinct %>% 
        filter(difftime(MED_ORDER_DATE, COLLECTION_DTTM, units = "hours") < 24 & difftime(MED_ORDER_DATE, COLLECTION_DTTM, units = "hours") > 12 & MED_ORDER_DATE > COLLECTION_DTTM) 
  
  orders_next_12h <- 
    orders %>% 
      distinct %>% 
        filter(difftime(MED_ORDER_DATE, COLLECTION_DTTM, units = "hours") < 12 & MED_ORDER_DATE > COLLECTION_DTTM) 
  
  orders_nested_48 <- orders_prior_48h %>% select(Epic_MRN, MED_NAME, MED_ORDER_DATE, COLLECTION_DTTM) %>% nest_by(Epic_MRN)
  orders_nested_24 <- orders_prior_24h %>% select(Epic_MRN, MED_NAME, MED_ORDER_DATE, COLLECTION_DTTM) %>% nest_by(Epic_MRN)
  orders_nested_12 <- orders_prior_12h %>% select(Epic_MRN, MED_NAME, MED_ORDER_DATE, COLLECTION_DTTM) %>% nest_by(Epic_MRN)
  
  orders_nested_next_48 <- orders_next_48h %>% select(Epic_MRN, MED_NAME, MED_ORDER_DATE, COLLECTION_DTTM) %>% nest_by(Epic_MRN)
  orders_nested_next_24 <- orders_next_24h %>% select(Epic_MRN, MED_NAME, MED_ORDER_DATE, COLLECTION_DTTM) %>% nest_by(Epic_MRN)
  orders_nested_next_12 <- orders_next_12h %>% select(Epic_MRN, MED_NAME, MED_ORDER_DATE, COLLECTION_DTTM) %>% nest_by(Epic_MRN)
  
  orders_nested <- orders_prior_48h %>% select(Epic_MRN, MED_NAME, MED_ORDER_DATE, COLLECTION_DTTM) %>% nest_by(Epic_MRN)
  
  orders_nested$heparin_prior_48 <- orders_nested_48$data %>% map(1) %>% map(~ifelse(grepl("HEPARIN", .), T, F)) %>% map(~any(.)) %>% unlist
  orders_nested$heparin_prior_24 <- orders_nested_24$data %>% map(1) %>% map(~ifelse(grepl("HEPARIN", .), T, F)) %>% map(~any(.)) %>% unlist
  orders_nested$heparin_prior_12 <- orders_nested_12$data %>% map(1) %>% map(~ifelse(grepl("HEPARIN", .), T, F)) %>% map(~any(.)) %>% unlist
  
  orders_nested$bival_prior_48 <- orders_nested_48$data %>% map(1) %>% map(~ifelse(grepl("BIVAL", .), T, F)) %>% map(~any(.)) %>% unlist
  orders_nested$bival_prior_24 <- orders_nested_24$data %>% map(1) %>% map(~ifelse(grepl("BIVAL", .), T, F)) %>% map(~any(.)) %>% unlist
  orders_nested$bival_prior_12 <- orders_nested_12$data %>% map(1) %>% map(~ifelse(grepl("BIVAL", .), T, F)) %>% map(~any(.)) %>% unlist

  orders_nested <- orders_nested_next_48$data %>% map(1) %>% map(~ifelse(grepl("HEPARIN", .), T, F)) %>% map(~any(.)) %>% unlist %>% bind_cols(Epic_MRN = orders_nested_next_48$Epic_MRN, heparin_next_48 = .) %>% left_join(orders_nested, .)
  orders_nested <- orders_nested_next_24$data %>% map(1) %>% map(~ifelse(grepl("HEPARIN", .), T, F)) %>% map(~any(.)) %>% unlist %>% bind_cols(Epic_MRN = orders_nested_next_24$Epic_MRN, heparin_next_24 = .) %>% left_join(orders_nested, .)
  orders_nested <- orders_nested_next_12$data %>% map(1) %>% map(~ifelse(grepl("HEPARIN", .), T, F)) %>% map(~any(.)) %>% unlist %>% bind_cols(Epic_MRN = orders_nested_next_12$Epic_MRN, heparin_next_12 = .) %>% left_join(orders_nested, .)
  
  orders_nested <- orders_nested_next_48$data %>% map(1) %>% map(~ifelse(grepl("BIVAL", .), T, F)) %>% map(~any(.)) %>% unlist %>% bind_cols(Epic_MRN = orders_nested_next_48$Epic_MRN, bival_next_48 = .) %>% left_join(orders_nested, .)
  orders_nested <- orders_nested_next_24$data %>% map(1) %>% map(~ifelse(grepl("BIVAL", .), T, F)) %>% map(~any(.)) %>% unlist %>% bind_cols(Epic_MRN = orders_nested_next_24$Epic_MRN, bival_next_24 = .) %>% left_join(orders_nested, .)
  orders_nested <- orders_nested_next_12$data %>% map(1) %>% map(~ifelse(grepl("BIVAL", .), T, F)) %>% map(~any(.)) %>% unlist %>% bind_cols(Epic_MRN = orders_nested_next_12$Epic_MRN, bival_next_12 = .) %>% left_join(orders_nested, .)
    
  orders_nested 
  
}

plotOrderChanges <- function(orders_nested = orders_nested){
  
  MRN_order <- orders_nested %>% arrange(heparin_prior_12, heparin_prior_24, heparin_prior_48, bival_prior_12, bival_prior_24, bival_prior_48, heparin_next_12, heparin_next_24, heparin_next_48, bival_next_12, bival_next_24) %>% select(Epic_MRN)
  
  mar_gg_input_heparin <- pivot_longer(mar_nested %>% select(!matches("bival"), -data), cols = matches("heparin"), names_to = "TIME", values_to = "HEPARIN") %>% mutate(TIME = gsub("heparin_", "", TIME))
  mar_gg_input_bival <- pivot_longer(mar_nested %>% select(!matches("heparin"), -data), cols = matches("bival"), names_to = "TIME", values_to = "BIVAL") %>% mutate(TIME = gsub("bival_", "", TIME))
  
  mar_gg_input <- 
    full_join(mar_gg_input_heparin, mar_gg_input_bival) %>% 
    mutate(Epic_MRN = factor(Epic_MRN, levels = MRN_order[[1]]), TIME = factor(TIME, levels = c("prior_48", "prior_24", "prior_12", "next_12", "next_24", "next_48"))) %>%
    mutate(HEPARIN = replace_na(HEPARIN, FALSE), BIVAL = replace_na(BIVAL, FALSE)) %>% 
    mutate(ANTICOAGULATION = case_when(
      HEPARIN & BIVAL ~ "BOTH",
      HEPARIN & !BIVAL ~ "HEPARIN",
      !HEPARIN & BIVAL ~ "BIVAL",
      !HEPARIN & !BIVAL ~ "NEITHER"
    )) 
  
  ggplot(mar_gg_input, aes(x = TIME, y = Epic_MRN, fill = ANTICOAGULATION)) + 
    geom_tile() + 
    scale_fill_manual(values = c("BOTH" = "black", "HEPARIN" = viridis(9, option = "D")[[3]], "BIVAL" = viridis(9, option = "D")[[6]], "NEITHER" = "gray50", `NA` = "gray50")) + 
    scale_x_discrete(limits = c("prior_12", "next_12", "next_24", "next_48")) + 
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
  
}

getOrderDecisions <- function(order_changes = order_changes){
  
  order_changes %>% mutate(heparin_stopped = ifelse(heparin_prior_24 & !heparin_next_24, T, F))
  
}

getLIAresults <- function(hit_labs){
  
  LIA <- 
    hit_labs %>% 
      filter(COMPONENT_NAME == "HIT AB LIA UNITS") %>% 
      group_by(Epic_MRN) %>% 
      select(Epic_MRN, COLLECTION_DTTM, LAB_RSLT_VALUE) %>% 
      distinct() %>% 
      mutate(
        BINARY = ifelse(LAB_RSLT_VALUE > 1 | grepl(">", LAB_RSLT_VALUE), "Positive", "Negative"),
        TEST_NUM = row_number()
        )
  
  LIA
  
}

getSRAresults <- function(hit_labs){
  
  SRA <- 
    hit_labs %>% 
      filter(COMPONENT_NAME == "UNFRACTIONATED HEPARIN RESULT") %>% 
      group_by(Epic_MRN) %>% 
      select(Epic_MRN, COLLECTION_DTTM, LAB_RSLT_VALUE) %>% 
      distinct() %>% 
      mutate(
        TEST_NUM = row_number(),
        BINARY = LAB_RSLT_VALUE
        )
    
  SRA
  
}
    
combineResults <- function(LIA = lia_results, SRA = sra_results){
  
  results <- 
    bind_rows(
      LIA %>% mutate(TEST_NAME = "LIA"),
      SRA %>% mutate(TEST_NAME = "SRA")
      ) %>% 
    group_by(Epic_MRN) %>% 
    arrange(Epic_MRN, COLLECTION_DTTM)
  
  results 
  
}

makePositivityRatesFigures <- function(results = combined_results, order_times = order_times){
  
  pos_rate <- 
    results %>% 
      filter(TEST_NUM == 1) %>%
      group_by(TEST_NAME) %>%
      count(BINARY) %>% mutate(prop = round(n / sum(n), digits = 2), BINARY = factor(BINARY, levels = c("See Interpretation", "Positive", "Negative"))) %>%
        ggplot(aes(x = TEST_NAME, y = prop, fill = BINARY, alpha = BINARY)) +
          geom_col(width = .50) + 
          geom_text(data = . %>% filter(BINARY != "See Interpretation"), aes(label = prop), fontface = "bold.italic", color = "white", size = 8, position = position_stack(vjust = 0.5)) +
          scale_fill_manual(values = color_map) +
          scale_alpha_manual(values = alpha_map) +
          ylab("Proportion") + ggtitle("Positivity Rates of First Test") +
          theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 16, face = "bold.italic", color = "black"), axis.text.y = element_blank(), legend.position = "none")
  ggsave("../Figures/PositivityRates.svg", pos_rate, width = 8, height = 8)
  
  pos_count <- 
    results %>% 
      filter(TEST_NUM == 1) %>%
      group_by(TEST_NAME) %>%
      count(BINARY) %>% mutate(prop = round(n / sum(n), digits = 2), BINARY = factor(BINARY, levels = c("See Interpretation", "Positive", "Negative"))) %>%
        ggplot(aes(x = TEST_NAME, y = n, fill = BINARY, alpha = BINARY)) +
          geom_col(width = .50) + 
          geom_text(data = . %>% filter(BINARY == "Negative" | TEST_NAME == "LIA"), aes(label = n), color = "white", size = 8, position = position_stack(vjust = 0.5)) +
          annotate("text", x = "SRA", y = 40, label = "7*", color = "black", size = 8) +
          scale_fill_manual(values = color_map) +
          scale_alpha_manual(values = alpha_map) +
          ylab("Count") + ggtitle("Result Counts of First Test") +
          theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 16, face = "bold.italic", color = "black"), axis.text.y = element_blank(), legend.position = c(0.75, 0.88), legend.title = element_blank())
  ggsave("../Figures/PositivityCounts.svg", pos_count, width = 8, height = 8)
  
  ggpubr::ggarrange(pos_count, pos_rate, nrow = 1) 
  ggsave(here("../Figures/PositivityFigure.svg"), width = 8, height = 4)
  
}

makeRepeatLIAFigures <- function(results = combined_results, LIA = lia_results){
  
  LIA_wide <- 
    LIA %>% 
      pivot_wider(id_cols = Epic_MRN, names_from = TEST_NUM, values_from = BINARY, names_prefix = "LIA_")
  
  LIA_first_neg <- 
    LIA_wide %>% filter(LIA_1 == "Negative")
  
  first_neg_count <- LIA_first_neg %>% group_by(LIA_1, LIA_2) %>% count()
  LIA_neg_then_pos <- LIA_first_neg %>% filter(LIA_2 == "Positive")
  
  ggplot() + 
    geom_bar(data = LIA %>% filter(TEST_NUM == 1), aes(x = 1, fill = BINARY), position = "stack") + 
    geom_bar(data = LIA_first_neg, aes(x = 2, y = after_stat(prop), fill = LIA_2), position = "stack")
  
  library(ggalluvial)
  
  ggplot() + 
    scale_fill_manual(values = color_map) +
    geom_flow(data = LIA %>% filter(TEST_NUM < 3),
              aes(x = TEST_NUM, stratum = BINARY, alluvium = Epic_MRN,
                  fill = BINARY, label = BINARY), stat = "alluvium", lode.guidance = "backward", aes.flow = "backward") +
    geom_stratum(data = LIA %>% filter(TEST_NUM < 3),
                 aes(x = TEST_NUM, stratum = BINARY, alluvium = Epic_MRN,
                     fill = BINARY, label = BINARY), na.rm = T, color = NA) +
    geom_text(data = LIA %>% filter(TEST_NUM < 2),
              aes(x = TEST_NUM, stratum = BINARY, alluvium = Epic_MRN,
                  fill = BINARY, label = BINARY), stat = "stratum", color = "white", size = 4, fontface = "bold.italic", na.rm = T) +
    geom_flow(data = LIA %>% filter(TEST_NUM > 1 & TEST_NUM < 4),
              aes(x = TEST_NUM, stratum = BINARY, alluvium = Epic_MRN,
                  fill = BINARY, label = BINARY), stat = "alluvium", lode.guidance = "forward", aes.flow = "backward", na.rm = T) +
    geom_stratum(data = LIA %>% filter(TEST_NUM > 1 & TEST_NUM < 4),
                 aes(x = TEST_NUM, stratum = BINARY, alluvium = Epic_MRN,
                     fill = BINARY, label = BINARY), na.rm = T, color = NA) +
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("First", "Second", "Third")) + 
    ylab("Count") + 
    theme(legend.position = "none", axis.title.x = element_blank())
    ggsave(here("../Figures/RepeatAlluvialPlot.svg"), width = 8, height = 8)
  
}

makeRepeatSRAFigures <- function(results = combined_results, SRA = sra_results){
  
  SRA_wide <- 
    SRA %>% 
    pivot_wider(id_cols = Epic_MRN, names_from = TEST_NUM, values_from = BINARY, names_prefix = "SRA_")
  
  SRA_first_neg <- 
    SRA_wide %>% filter(SRA_1 == "Negative")
  
  first_neg_count <- SRA_first_neg %>% group_by(SRA_1, SRA_2) %>% count()
  SRA_neg_then_pos <- SRA_first_neg %>% filter(SRA_2 == "Positive")
  
  library(ggalluvial)
  
  ggplot() + 
    geom_bar(data = SRA %>% filter(TEST_NUM == 1), aes(x = 1, fill = BINARY), position = "stack") + 
    geom_bar(data = SRA_first_neg, aes(x = 2, y = after_stat(prop), fill = SRA_2), position = "stack")
  
  ggplot() + 
    geom_flow(data = SRA,
              aes(x = TEST_NUM, stratum = BINARY, alluvium = Epic_MRN,
                  fill = BINARY, label = BINARY), stat = "alluvium", lode.guidance = "backward", aes.flow = "backward") +
    geom_stratum(data = SRA,
                 aes(x = TEST_NUM, stratum = BINARY, alluvium = Epic_MRN,
                     fill = BINARY, label = BINARY, alpha = BINARY), na.rm = T, color = NA) +
    geom_text(data = SRA,
              aes(x = TEST_NUM, stratum = BINARY, alluvium = Epic_MRN,
                  fill = BINARY, label = BINARY), stat = "stratum", color = "white", size = 4, fontface = "bold.italic", na.rm = T) +
    scale_x_continuous(breaks = c(1, 2), labels = c("First", "Second")) + 
    scale_fill_manual(values = color_map) +
    scale_alpha_manual(values = alpha_map) +
    ylab("Count") + 
    theme(legend.position = "none", axis.title.x = element_blank())
  ggsave(here("../Figures/RepeatAlluvialPlotSRA.svg"), width = 8, height = 8)
  
}

makeSRAbyLIAplot <- function(results = combined_results){
  
  
  f <- function(x, count, group) {
    data.frame(x, count, group) %>%
      add_count(x, wt = count) %>%
      group_by(x, group) %>%
      mutate(prop = count / n) %>%
      pull(prop)
  }
  
  input <- 
    results %>% 
      pivot_wider(id_cols = c(Epic_MRN, COLLECTION_DTTM), names_from = TEST_NAME, values_from = LAB_RSLT_VALUE) %>%
      mutate(LIA = gsub(">16.00", "9", LIA))
  
  breaks = c(0, 1, 2, 4, 6, 8, 10)
  
  ggplot(input %>% filter(SRA != "NA" & SRA !="See Interpretation"), aes(x = as.numeric(as.character(LIA)))) +
    geom_histogram(aes(group = factor(SRA), fill = factor(SRA), 
                       y = f(..x.., ..count.., ..group..)), breaks = breaks, color = "gray90") + 
    stat_bin(aes(y = f(..x.., ..count.., ..group..),
                 label = ifelse((..count.. < 10 & f(..x.., ..count.., ..group..) < .05),"", paste(..count..,"\n","(",scales::percent(f(..x.., ..count.., ..group..), accuracy = 1), ")", sep = "")), group = SRA), 
             breaks = breaks, geom = "text", position = position_stack(vjust = 0.5), size = 6, color = "white", fontface = "bold.italic") +
    stat_bin(aes(y = f(..x.., ..count.., ..group..), label = ..count..), 
             breaks = breaks, geom = "text", vjust = -.5, size = 5, fontface = "italic") + guides(fill=guide_legend("SRA Result")) +
    scale_fill_manual(values=c("gray60", "#8B0000", "#8B0000")) + geom_segment(aes(x = 1, xend = 1, y = 0, yend = 1)) + 
    coord_cartesian(xlim = c(0, 10), ylim = c(0,1.05), clip = "off") + scale_x_continuous(name = "LIA units", breaks = breaks, labels = c("0", "1", "2", "4", "6", ">8", "")) +
    ggtitle("SRA Positivity by LIA Result") + ylab("Count") +
    theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold.italic"), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  ggsave(here("../Figures/SRAbyLIA.svg"), width = 8, height = 8)
  
}

getRecentPlatelets <- function(labs = labs, trigger = hit_collect){
  
  prior_results <- full_join(labs, trigger %>% rename(LIA_COLLECTION_DTTM = COLLECTION_DTTM)) %>% filter(COLLECTION_DTTM < LIA_COLLECTION_DTTM)
  
  prior_platelets <- 
    prior_results %>% 
      filter(COMPONENT_NAME == "PLATELET") %>% 
      arrange(Epic_MRN, desc(COLLECTION_DTTM)) %>% 
      mutate(hours_until_LIA = -1 * difftime(COLLECTION_DTTM, LIA_COLLECTION_DTTM, units = "hours")) %>%
      filter(hours_until_LIA < 240) %>% mutate(LAB_RSLT_VALUE = as.numeric(LAB_RSLT_VALUE))
  
  prior_platelets <- 
    prior_platelets %>% 
      group_by(Epic_MRN) %>% 
      mutate(max_platelet_10d = max(LAB_RSLT_VALUE, na.rm = T), percent_of_max = LAB_RSLT_VALUE/max_platelet_10d)

  
  
}

makeTimelineTibble <- function(labs = hit_labs, mar = hit_mar, orders = hit_med_orders, order_times = order_times){
  
  first_lia_orders <- order_times %>% rename(Epic_MRN = MRN) %>% arrange(Epic_MRN, ORIG_ORDER_DT_TM) %>% group_by(Epic_MRN) %>% slice_head(n = 1) %>% select(Epic_MRN, ORIG_ORDER_DT_TM)
  first_lia <- left_join(first_lia_orders, hit_labs) %>% filter(COMPONENT_NAME == "HIT AB") %>% arrange(ORIG_ORDER_DT_TM) %>% group_by(Epic_MRN) %>% slice_head(n = 1)
  first_lia <- first_lia %>% mutate(diff = as.numeric(difftime(ORIG_ORDER_DT_TM, COLLECTION_DTTM, units = "hours"))) %>% filter()
  
  first_sra <- labs %>% filter(COMPONENT_NAME == "UNFRACTIONATED HEPARIN RESULT") %>% group_by(Epic_MRN) %>% arrange(COLLECTION_DTTM) %>% slice_head(n = 1)
  first_sra_time <- first_sra %>% select(Epic_MRN, SRA_TIME = COLLECTION_DTTM)
  
  first_times <- left_join(first_lia, first_sra)
  
  orders_nested <- getOrderChanges(first_lia, orders)
  
}

tmp <- function(){
  
  LIA %>% 
    ggplot(aes(TEST_NUM, fill = BINARY)) + 
    geom_histogram(bins = 10) + 
    scale_x_continuous(breaks = seq(1:10))
  
  counts <- 
    LIA %>% 
    pivot_wider(id_cols = Epic_MRN, names_from = TEST_NUM, values_from = BINARY, names_prefix = "Test_") %>% 
    group_by(Test_1, Test_2, Test_3, Test_4, Test_5) %>% 
    count() %>% ungroup() %>% mutate(prop = n / sum(n))
  
  ggplot(LIA, aes(x = TEST_NUM, fill = BINARY)) +
    geom_bar(position = "stack")
  
  
  first_pos <- counts %>% filter(Test_1) %>% group_by(Test_2) %>% tally(n)
  first_two_pos <- counts %>% filter(Test_1 & Test_2) %>% group_by(Test_3) %>% tally(n)
  
  first_neg <- counts %>% filter(!Test_1) %>% group_by(Test_2) %>% tally(n)
  first_two_neg <- counts %>% filter(!Test_1 & !Test_2) %>% group_by(Test_3) %>% tally(n)
  
  ggplot() + 
    geom_bar(data = LIA, aes(x = TEST_NUM), position = "stack")
  
  
}
  