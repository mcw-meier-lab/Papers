###############################################################################################################
#                       Example Code for original article submitted to Brain Communications                   #
#     Profiling central nervous system-injury proteins in multiple mild traumatic brain injury cohorts        #
#                             Corresponding Author: Timothy B. Meier (tmeier@mcw.edu)                         #
#                           Code Author: Mitchell J. Andersson (mandersson@mcw.edu)                           #
#                                          DATE: August 10th, 2026                                            #
###############################################################################################################

# TABLE OF CONTENTS
# 1. SPORTS-RELATED CONCUSSION COHORT
# 1.1. A priori Between-Subjects Comparison
# 1.2. Discovery Between-Subjects Comparison
# 1.3. A priori T-Tests - Within Subjects vs Baseline
# 1.4. Discovery T-Tests - Within Subjects vs Baseline

# 2. PEDIATRIC COHORT
# 2.1. A priori Between-Subjects Comparison (at 1 week & 4 months)
# 2.2. Discovery Between-Subjects Comparison (at 1 week & 4 months)
# 2.3. Sensitivity Analysis: Exclude CT/MRI+
# 2.4. Sensitivity analysis: Assessing biomarker concentration among pediatric patients < 2 days

# 3. OLDER ADULT COHORT
# 3.1. A priori Between-Subjects Comparison (at 1 week)
# 3.2. Discovery Between-Subjects Comparison (at 1 week)
# 3.3. Comparing mTBI patients at 4 months versus controls - a priori
# 3.4. Comparing mTBI patients at 4 months versus controls - discovery
# 3.5. Sensitivity Analysis: Exclude CT/MRI+

# LOAD PACKAGES
library(plyr)
library(tidyverse)   # Data manipulation (dplyr, tidyr, purrr)
library(gtsummary)   # Demographic summary tables
theme_gtsummary_compact() # Compact styling theme for tables
library(DT)          # Interactive datatables
library(broom)

#####################################################################################################################
######################################## 1. SPORTS-RELATED CONCUSSION COHORT ########################################
#####################################################################################################################

# POST-IMPORT AND DATA CLEANING #

# src_df = SPORTS-RELATED CONCUSSION DATAFRAME CONSISTING OF ATHLETES WHO SUSTAINED AN mTBI and MATCHED CONTACT CONTROLS. 
#          type distinguishes groups, subject_id is participant ID
# ped_df = PEDIATRIC mTBI DATAFRAME CONSISTING OF PEDIATRIC mTBI PATIENTS AND MATCHED CONTROLS. 
#          dx_fac distinguishes groups, URSI is participant ID
# adults_matched_df = OLDER ADULT mTBI DATAFRAME CONSISTING OF OLDER ADULT mTBI PATIENTS AND MATCHED CONTROLS. 
#          dx_fac distinguishes groups, URSI is participant ID


# DEFINE A PRIORI AND DISCOVERY BIOMARKERS

target_biomarkers <- c("GFAP", "IL6", "Tau", "NfL", "pTau181", "S100B", "UCHL1")

discovery_vars <- strsplit("ACHE,AGRN,ANXA5,APOE,APOE4,ARSA,A38,AB40,AB42,
                           BACE1,BASP1,BDNF,
                           CALB2,CCL11,CCL13,CCL17,CCL2,CCL22,CCL26,CCL3,CCL4,
                           CD40LG,CD63,CHI3L1,CHIT1,CNTN2,CRH,CRP,CSF2,CST3,CX3CL1,
                           CXCL1,CXCL10,CXCL8,DDC,ENO2,FABP3,FCN2,FGF2,FLT1,FOLR1,
                           GDF15,GDI1,GDNF,GOT1,HBA1,HTT,
                           ICAM1,IFNG,IGF1R,IGFBP7,IL10,IL12p70,IL13,IL15,IL16,
                           IL17A,IL18,IL1B,IL2,IL33,IL4,IL5,IL6R,IL7,IL9,
                           KDR,KLK6,MDH1,MME,MSLN,NEFH,NGF,NPTX1,NPTX2,NPTXR,NPY,NRGN,
                           OligoSNCA,PARK7,PDGFRB,PDLIM5,PGF,PGK1,POSTN,PRDX6,PSEN1,
                           pSNCA129,pTau217,pTau231,pTDP43409,PTN,REST,RUVBL2,
                           S100A12,SAA1,SFRP1,SFTPD,SLIT2,SMOC1,SNAP25,SNCA,SNCB,SOD1,SQSTM1,
                           TAFA5,TARDBP,TEK,TIMP3,TNF,TREM1,TREM2,UBB,
                           VCAM1,VEGFA,VEGFD,VGF,VSNL1,YWHAG,YWHAZ", split = ",")[[1]]

###############################################
## 1.1. A Priori Between-Subjects Comparison ##
###############################################

# WELCH T-TEST FUNCTION
t_test_fun<-function(x, y) {
  ttest=(t.test(x ~ y))}

# FILTER DATAFRAME TO INCLUDE ONLY mTBI & CONTACT CONTROLS AT BASELINE, 6-HOURS, 48-HOURS & 15 DAYS POST-INJURY
dat_apriori <- src_df %>%  
  filter(((type == "mTBI" | type ==  "Contact Control") & (visit_fac == "Baseline" | visit_fac == "6-Hour" | visit_fac == "48-Hour" | visit_fac == "Day 15"))) %>% 
  select(subject_id,visit_fac, type, GFAP,IL6,NfL,Tau,pTau181,S100B,UCHL1) %>% 
  pivot_longer(c(-subject_id,-visit_fac, -type))

# EXAMPLE MODEL
dat_apriori %>% filter(visit_fac == "Baseline" & name == "GFAP") %>% 
  t.test(value ~ type, data = .)

# EXAMPLE CHECK VALID VALUES PER GROUP
dat_apriori %>% 
  filter(visit_fac == "Baseline" & name == "IL6") %>% view()

dat_apriori %>% 
  filter(visit_fac == "Baseline" & name == "GFAP") %>% 
  select(type, value) %>% tbl_summary(by = type) %>% add_overall(last = T)

# RUN MODELS ACROSS ALL BIOMARKERS AT EACH type# RUN MODELS ACROSS ALL BIOMARKERS AT EACH VISIT
ttest_models_apriori <- ddply(dat_apriori, .(visit_fac, name), summarise,z=t_test_fun(value,type)$statistic,
                              pval=t_test_fun(value,type)$p.value,
                              cc_mean=t_test_fun(value,type)$estimate[1],
                              mtbi_mean=t_test_fun(value,type)$estimate[2],
                              n_total = sum(!is.na(value)))

# CALCULATE FOLD CHANGE
ttest_models_apriori$log2fc <- ttest_models_apriori$mtbi_mean - ttest_models_apriori$cc_mean

DT::datatable(ttest_models_apriori, 
              extensions = c('FixedColumns',"FixedHeader"),
              options = list(scrollX = TRUE, 
                             paging=FALSE,
                             fixedHeader=TRUE))

################################################
## 1.2. Discovery Between-Subjects Comparison ##
################################################

# FILTER DATAFRAME TO INCLUDE ONLY mTBI & CONTACT CONTROLS AT BASELINE, 6-HOURS, 48-HOURS & 15 DAYS POST-INJURY
dat_discovery <- src_df %>%  
  filter(((type == "mTBI" | type ==  "Contact Control") & (visit_fac == "Baseline" | visit_fac == "6-Hour" | visit_fac == "48-Hour" | visit_fac == "Day 15"))) %>% 
  select(subject_id, visit_fac, type, all_of(discovery_vars)) %>% 
  pivot_longer(c(-subject_id, -visit_fac, -type))

# EXAMPLE MODEL
dat_discovery %>% filter(visit_fac == "Baseline" & name == "ACHE") %>% 
  t.test(value ~ type, data = .)

# RUN MODELS ACROSS ALL BIOMARKERS AT EACH VISIT
ttest_models_discovery <- ddply(dat_discovery, .(visit_fac, name), summarise,z=t_test_fun(value,type)$statistic,
                                pval=t_test_fun(value,type)$p.value,
                                cc_mean=t_test_fun(value,type)$estimate[1],
                                mtbi_mean=t_test_fun(value,type)$estimate[2],
                                n_total = sum(!is.na(value)))

# CALCULATE FOLD CHANGE
ttest_models_discovery$log2fc <- ttest_models_discovery$mtbi_mean - ttest_models_discovery$cc_mean

# FDR-ADJUST TESTS BY VISIT
ttest_bl <- ttest_models_discovery %>% filter(visit_fac=='Baseline') %>% 
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))
ttest_6h <- ttest_models_discovery %>% filter(visit_fac=='6-Hour') %>% 
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))
ttest_48h <- ttest_models_discovery %>% filter(visit_fac=='48-Hour') %>% 
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))
ttest_d15 <- ttest_models_discovery %>% filter(visit_fac=='Day 15') %>% 
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))

m <- rbind(ttest_bl,ttest_6h,ttest_48h,ttest_d15)

DT::datatable(ttest_models_discovery, 
              extensions = c('FixedColumns',"FixedHeader"),
              options = list(scrollX = TRUE, 
                             paging=FALSE,
                             fixedHeader=TRUE))

#########################################################
## 1.3. A Priori T-Tests - Within Subjects vs Baseline ##
#########################################################

# PAIRED T-TEST FUJNCTION
paired_t_test_fun <- function(x, y) {
  t.test(x, y, paired = TRUE)}

# SELECT BIOMARKERS AND RESHAPE DATAFRAME
pair_dat_apriori <- src_df %>%
  select(subject_id, visit_fac, type, GFAP, IL6, NfL, Tau, pTau181, S100B, UCHL1) %>%
  filter(type == "mTBI" & visit_fac %in% c("Baseline", "6-Hour", "48-Hour", "Day 15")) %>%
  pivot_longer(cols = -c(subject_id, visit_fac, type)) %>%
  filter(!is.na(value)) %>%
  pivot_wider(names_from = visit_fac, values_from = value, id_cols = c(subject_id, name, type)) %>%
  select(-type) %>%
  pivot_longer(cols = c("6-Hour", "48-Hour", "Day 15"), names_to = "visit_fac") %>%
  filter(!is.na(value))

# EXAMPLE MODEL
pair_dat_apriori %>% 
  filter(visit_fac == "Day 15" & name == "GFAP") %>% 
  t.test(value-Baseline ~ 1, data = .)
# OR
pair_dat_apriori %>% 
  filter(visit_fac == "Day 15" & name == "GFAP") %>% 
  with(t.test(value, Baseline, paired = TRUE))

# EXAMPLE CHECK VALID VALUES PER GROUP
src_df %>% 
  filter((visit_fac == "6-Hour" | visit_fac == "Baseline") & type == "mTBI") %>% 
  select(subject_id, type, visit_fac, GFAP) %>% view()

pair_dat_apriori %>% 
  filter(visit_fac == "6-Hour" & name == "GFAP") %>% 
  mutate(m_diff = value-Baseline) %>%  
  select(m_diff) %>% tbl_summary()

# RUN PAIRED T-TESTS FOR EACH BIOMARKER AT EACH TIMEPOINT VS BASELINE
paired_ttest_models_apriori <- ddply(
  pair_dat_apriori,
  .(visit_fac, name),
  summarise,
  t = paired_t_test_fun(value, Baseline)$statistic,
  pval = paired_t_test_fun(value, Baseline)$p.value,
  log2fc = paired_t_test_fun(value, Baseline)$estimate[1],
  n_total = sum(!is.na(value)))

DT::datatable(
  paired_ttest_models_apriori,
  extensions = c('FixedColumns', 'FixedHeader'),
  options = list(scrollX = TRUE, paging = FALSE, fixedHeader = TRUE))

############################################################
## # 1.4. Discovery T-Tests - Within Subjects vs Baseline ##
############################################################

# SELECT BIOMARKERS AND RESHAPE DATAFRAME
pair_dat_discovery <- src_df %>%
  select(subject_id, visit_fac, type, all_of(discovery_vars)) %>%
  filter(type == "mTBI" & visit_fac %in% c("Baseline", "6-Hour", "48-Hour", "Day 15")) %>%
  pivot_longer(cols = -c(subject_id, visit_fac, type)) %>%
  filter(!is.na(value)) %>%
  pivot_wider(names_from = visit_fac, values_from = value, id_cols = c(subject_id, name, type)) %>%
  select(-type) %>%
  pivot_longer(cols = c("6-Hour", "48-Hour", "Day 15"), names_to = "visit_fac") %>%
  filter(!is.na(value))

# EXAMPLE MODEL
pair_dat_discovery %>% 
  filter(visit_fac == "6-Hour" & name == "ACHE") %>% 
  t.test(value-Baseline ~ 1, data = .)


# RUN PAIRED T-TESTS FOR EACH BIOMARKER AT EACH TIMEPOINT VS BASELINE
paired_ttest_models_discovery <- ddply(
  pair_dat_discovery,
  .(visit_fac, name),
  summarise,
  t = paired_t_test_fun(value, Baseline)$statistic,
  pval = paired_t_test_fun(value, Baseline)$p.value,
  log2fc = paired_t_test_fun(value, Baseline)$estimate[1],
  n_total = sum(!is.na(value)))


# FDR-ADJUST TESTS BY VISIT
paired_ttest_discovery_6h <- paired_ttest_models_discovery %>% filter(visit_fac=='6-Hour') %>% 
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))

paired_ttest_discovery_48h <- paired_ttest_models_discovery %>% filter(visit_fac=='48-Hour') %>% 
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))

paired_ttest_discovery_d15 <- paired_ttest_models_discovery %>% filter(visit_fac=='Day 15') %>% 
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))

paired_ttest_models_discovery <- rbind(paired_ttest_discovery_6h,
                                       paired_ttest_discovery_48h,
                                       paired_ttest_discovery_d15)


DT::datatable(
  paired_ttest_models_discovery,
  extensions = c('FixedColumns', 'FixedHeader'),
  options = list(scrollX = TRUE, paging = FALSE, fixedHeader = TRUE))

#####################################################################################################################
############################################# 2. PEDIATRIC mTBI COHORT ##############################################
#####################################################################################################################

###############################################
## 2.1. A Priori Between-Subjects Comparison ##
###############################################

# SELECT BIOMARKERS AND RESHAPE DATAFRAME
dat_apriori_ped <- ped_df %>%  
  filter(((dx_fac == "mTBI" | dx_fac ==  "Healthy Controls") & (visit_fac == "1 week" | visit_fac == "4 months"))) %>% 
  select(URSI, visit_fac, dx_fac, GFAP,IL6,NfL,Tau,pTau181,S100B,UCHL1) %>% 
  pivot_longer(c(-URSI, -visit_fac, -dx_fac))

# EXAMPLE MODEL
dat_apriori_ped %>% filter(visit_fac == "1 week" & name == "GFAP") %>% 
  t.test(value ~ dx_fac, data = .)

# EXAMPLE CHECK VALID VALUES PER GROUP
dat_apriori_ped %>% 
  filter(visit_fac == "4 months" & name == "GFAP") %>% view()

dat_apriori_ped %>% 
  filter(visit_fac == "4 months" & name == "GFAP") %>% 
  select(dx_fac, value) %>% tbl_summary(by = dx_fac) %>% add_overall(last = T)

# RUN T-TESTS FOR EACH BIOMARKER AT EACH TIMEPOINT
ttest_models_apriori_ped <- ddply(dat_apriori_ped, .(visit_fac, name), summarise,z=t_test_fun(value,dx_fac)$statistic,
                                  pval=t_test_fun(value,dx_fac)$p.value,
                                  cc_mean=t_test_fun(value,dx_fac)$estimate[1],
                                  mtbi_mean=t_test_fun(value,dx_fac)$estimate[2],
                                  n_total = sum(!is.na(value)))

# CALCULATE FOLD CHANGE
ttest_models_apriori_ped$log2fc <- ttest_models_apriori_ped$mtbi_mean - ttest_models_apriori_ped$cc_mean

DT::datatable(ttest_models_apriori_ped, 
              extensions = c('FixedColumns',"FixedHeader"),
              options = list(scrollX = TRUE, 
                             paging=FALSE,
                             fixedHeader=TRUE))

################################################
## 2.2. Discovery Between-Subjects Comparison ##
################################################

# SELECT BIOMARKERS AND RESHAPE DATAFRAME
dat_discovery_ped <- ped_df %>%  
  filter(((dx_fac == "mTBI" | dx_fac ==  "Healthy Controls") & (visit_fac == "1 week" | visit_fac == "4 months"))) %>% 
  select(visit_fac, dx_fac, all_of(discovery_vars)) %>% 
  pivot_longer(c(-visit_fac, -dx_fac))

# EXAMPLE MODEL
dat_discovery_ped %>% filter(visit_fac == "1 week" & name == "ACHE") %>% 
  t.test(value ~ dx_fac, data = .)

# RUN T-TESTS FOR EACH BIOMARKER AT EACH TIMEPOINT
ttest_models_discovery_ped <- ddply(dat_discovery_ped, .(visit_fac, name), summarise,z=t_test_fun(value,dx_fac)$statistic,
                                    pval=t_test_fun(value,dx_fac)$p.value,
                                    cc_mean=t_test_fun(value,dx_fac)$estimate[1],
                                    mtbi_mean=t_test_fun(value,dx_fac)$estimate[2],
                                    n_total = sum(!is.na(value)))

# CALCULATE FOLD CHANGE
ttest_models_discovery_ped$log2fc <- ttest_models_discovery_ped$mtbi_mean - ttest_models_discovery_ped$cc_mean

# APPLY FDR-CORRECTION
ttest_models_discovery_ped_1week <- ttest_models_discovery_ped %>% filter(visit_fac == "1 week") %>% 
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))

ttest_models_discovery_ped_4months <- ttest_models_discovery_ped %>% filter(visit_fac == "4 months") %>% 
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))

ttest_models_discovery_ped <- rbind(ttest_models_discovery_ped_1week,ttest_models_discovery_ped_4months)

DT::datatable(ttest_models_discovery_ped, 
              extensions = c('FixedColumns',"FixedHeader"),
              options = list(scrollX = TRUE, 
                             paging=FALSE,
                             fixedHeader=TRUE))

#################################################
##  2.3. Sensitivity Analysis: Exclude CT/MRI+ ##
#################################################

# SUBSET DATAFRAME TO THOSE WITHOUT POSITIVE IMAGING
ped_df_ct_mri_neg <- ped_df %>%
  filter(!URSI %in% (filter(., ct_mri_pos == 1) %>% pull(URSI)))

# SELECT BIOMARKERS AND RESHAPE DATAFRAME
dat_apriori_ped_sens <- ped_df_ct_mri_neg %>%  
  filter(((dx_fac == "mTBI" | dx_fac ==  "Healthy Controls") & (visit_fac == "1 week" | visit_fac == "4 months"))) %>% 
  select(URSI, visit_fac, dx_fac, GFAP,IL6,NfL,Tau,pTau181,S100B,UCHL1) %>% 
  pivot_longer(c(-URSI,-visit_fac, -dx_fac))

# EXAMPLE MODEL
dat_apriori_ped_sens %>% filter(visit_fac == "1 week" & name == "Tau") %>% 
  t.test(value ~ dx_fac, data = .)

# RUN T-TESTS FOR EACH BIOMARKER AT EACH TIMEPOINT
ttest_models_apriori_ped_sens <- ddply(dat_apriori_ped_sens, .(visit_fac, name), summarise,z=t_test_fun(value,dx_fac)$statistic,
                                       pval=t_test_fun(value,dx_fac)$p.value,
                                       cc_mean=t_test_fun(value,dx_fac)$estimate[1],
                                       mtbi_mean=t_test_fun(value,dx_fac)$estimate[2],
                                       n_total = sum(!is.na(value)))

# CALCULATE FOLD CHANGE
ttest_models_apriori_ped_sens$log2fc <- ttest_models_apriori_ped_sens$mtbi_mean - ttest_models_apriori_ped_sens$cc_mean

DT::datatable(ttest_models_apriori_ped_sens, 
              extensions = c('FixedColumns',"FixedHeader"),
              options = list(scrollX = TRUE, 
                             paging=FALSE,
                             fixedHeader=TRUE))

# COMPARE TO PRIMARY ANALYSES
ttest_models_apriori_ped %>%
  left_join(
    ttest_models_apriori_ped_sens %>% 
      select(visit_fac, name, log2fc_sens = log2fc, pval_sens = pval),
    by = c("visit_fac", "name")) %>% 
  select(visit_fac, name, log2fc, pval, log2fc_sens, pval_sens) %>%
  mutate(log2fc_abs_diff = abs(log2fc_sens - log2fc)) %>% 
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x)))

# SELECT BIOMARKERS AND RESHAPE DATAFRAME
dat_discovery_ped_sens <- ped_df_ct_mri_neg %>%  
  filter(((dx_fac == "mTBI" | dx_fac ==  "Healthy Controls") & (visit_fac == "1 week" | visit_fac == "4 months"))) %>% 
  select(URSI,visit_fac, dx_fac, all_of(discovery_vars)) %>% 
  pivot_longer(c(-URSI,-visit_fac, -dx_fac))

# RUN T-TESTS FOR EACH BIOMARKER AT EACH TIMEPOINT
ttest_models_discovery_ped_sens <- ddply(dat_discovery_ped_sens, .(visit_fac, name), summarise,z=t_test_fun(value,dx_fac)$statistic,
                                         pval=t_test_fun(value,dx_fac)$p.value,
                                         cc_mean=t_test_fun(value,dx_fac)$estimate[1],
                                         mtbi_mean=t_test_fun(value,dx_fac)$estimate[2],
                                         n_total = sum(!is.na(value)))

# CALCULATE FOLD CHANGE
ttest_models_discovery_ped_sens$log2fc <- ttest_models_discovery_ped_sens$mtbi_mean - ttest_models_discovery_ped_sens$cc_mean

# APPLY FDR-CORRECTION
ttest_models_discovery_ped_sens_1week <- ttest_models_discovery_ped_sens %>% filter(visit_fac == "1 week") %>% 
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))

ttest_models_discovery_ped_sens_4months <- ttest_models_discovery_ped_sens %>% filter(visit_fac == "4 months") %>% 
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))

ttest_models_discovery_ped_sens <- rbind(ttest_models_discovery_ped_sens_1week, ttest_models_discovery_ped_sens_4months)


DT::datatable(ttest_models_discovery_ped_sens, 
              extensions = c('FixedColumns',"FixedHeader"),
              options = list(scrollX = TRUE, 
                             paging=FALSE,
                             fixedHeader=TRUE))

# COMPARE TO PRIMARY ANALYSES
ttest_models_discovery_ped %>%
  left_join(
    ttest_models_discovery_ped_sens %>% 
      select(visit_fac, name, log2fc_sens = log2fc, pval_fdr_sens = pval_fdr),
    by = c("visit_fac", "name")) %>% 
  select(visit_fac, name, log2fc, pval_fdr, log2fc_sens, pval_fdr_sens) %>%
  mutate(
    log2fc_abs_diff = abs(log2fc_sens - log2fc)) %>% 
  arrange(desc(log2fc_abs_diff)) %>% 
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x)))

##################################################################################################
##  2.4. Sensitivity analysis: Assessing biomarker expression among pediatric patients < 2 days ##
##################################################################################################

# The objective of these analyses is to enable better comparisons with the SRC 6 hour group. 
# Comparing the days post injury variable between mTBI pediatric patients who arrived within 
# 2 days of injury with SRC athletes at the 6-hour time point, the time is nearly identical 
# and does not differ statistically.

# MEAN & SD OF TIME TO POST-INJURY BLOOD DRAW IN SRC GROUP
src_df %>%
  filter(visit_fac == "6-Hour", type == "mTBI") %>%
  mutate(DPI = as.numeric(DPI)) %>%
  summarise(mean_DPI = mean(DPI, na.rm = TRUE),
            sd_DPI = sd(DPI, na.rm = TRUE),
            n_src = sum(!is.na(DPI)))

# SUBSET PEDIATRIC COHORT TO THOSE SEEN WITHIN 3 DAYS POST-INJURY
ped_df_sens_superacute <- ped_df %>% 
  filter(visit_fac == "1 week") %>% 
  filter((dx_fac == "mTBI" & DPI < 3) | (dx_fac == "Healthy Controls"))

ped_df_sens_superacute %>%
  mutate(DPI = as.numeric(DPI)) %>%
  summarise(mean_DPI = mean(DPI, na.rm = TRUE),
            sd_DPI = sd(DPI, na.rm = TRUE),
            n_ped = sum(!is.na(DPI)))

# COMPARE PEDIATRIC PATIENTS SEEN EARLIER VERSUS LATER ON DEMOGRAPIC VARIABLES
ped_mtbi_sens_comparison_df <- ped_df %>%
  filter(dx_fac == "mTBI" & ((DPI < 3) | (DPI > 3 & DPI < 15))) %>%
  mutate(
    group = case_when(
      DPI < 3 ~ "mTBI_superacute",
      DPI > 3 & DPI < 15 ~ "mTBI_acute"))

ped_mtbi_sens_comparison_df %>% 
  select(group,age_num,female_fac,race_fac,ethnicity_fac,
         NumPrevInj_binary, WRAT, locpta) %>% 
  tbl_summary(by = group,
              statistic = list(all_continuous() ~ "{mean} ({sd})")) %>% 
  add_p(pvalue_fun = label_style_pvalue(digits = 3))

# DOES THE TIME TO BLOOD DRAW STATISTICALLY DIFFER BETWEEN SPORT 6-HOUR TIMEPOINT AND EARLY PEDIATRIC GROUP
t.test(
  ped_raw %>% filter(DX == 1, DPI < 3) %>% pull(DPI),
  src_df %>% filter(visit_fac == "6-Hour", type == "mTBI") %>% pull(DPI))

#### DEFINE BIOMARKER LISTS ####
superacute_apriori_biomarkers <- c(
  "Tau", "NfL",
  "GFAP", "pTau181",
  "IL6", "S100B", "UCHL1")

# DISCOVERY MARKERS DERIVED FROM PRIMARY RESULTS #
superacute_discovery_biomarkers <- c(
  "ANXA5", "PGK1", "PRDX6", "PSEN1", "TARDBP",
  "MDH1", "S100A12", "YWHAG",
  "CSF2", "FABP3", "GOT1", "HBA1", "IL16",
  "IL18", "pTau231", "REST")

# EXAMPLE A PRIORI MODEL
ped_df_sens_superacute %>% 
  t.test(GFAP ~ dx_fac, data = .)

# EXAMPLE DISCOVERY MODEL
ped_df_sens_superacute %>% 
  t.test(HBA1 ~ dx_fac, data = .)


#### RUN T-TESTS FOR A PRIORI BIOMARKERS ####
t_results_apriori_superacute <- map_dfr(superacute_apriori_biomarkers, function(name) {
  formula <- as.formula(paste(name, "~ dx_fac"))
  
  t.test(formula, data = ped_df_sens_superacute) %>%
    tidy() %>%
    mutate(name = name)}) %>%
  select(name, estimate1, estimate2, statistic, p.value) %>%
  mutate(
    log2fc = estimate2 - estimate1,
    visit_fac = "0-2 Days^")

#### RUN T-TESTS FOR DISCOVERY BIOMARKERS ####
t_results_discovery_superacute <- map_dfr(superacute_discovery_biomarkers, function(name) {
  formula <- as.formula(paste(name, "~ dx_fac"))
  
  t.test(formula, data = ped_df_sens_superacute) %>%
    tidy() %>%
    mutate(name = name)}) %>%
  select(name, estimate1, estimate2, statistic, p.value) %>%
  mutate(
    log2fc = estimate2 - estimate1,
    visit_fac = "0-2 Days^",
    pval_fdr = p.adjust(p.value, method = "fdr"))

# View results
t_results_apriori_superacute %>%
  mutate(diffexpressed = case_when(
    log2fc > 0 & p.value < 0.05 ~ "UP", # CHANGED TO pval
    log2fc < 0 & p.value < 0.05 ~ "DOWN",
    p.value >= 0.05             ~ "NO")) %>% 
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x)))

t_results_discovery_superacute  %>%
  mutate(diffexpressed = case_when(
    log2fc > 0 & pval_fdr < 0.05 ~ "UP",
    log2fc < 0 & pval_fdr < 0.05 ~ "DOWN",
    pval_fdr >= 0.05             ~ "NO")) %>% 
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x)))



#####################################################################################################################
############################################## 3. OLDER ADULTS COHORT ###############################################
#####################################################################################################################

###############################################
## 3.1. A Priori Between-Subjects Comparison ##
###############################################

# SELECT BIOMARKERS AND RESHAPE DATAFRAME
dat_apriori_adults <- adults_df_matched %>%  
  filter(((dx_fac == "mTBI" | dx_fac ==  "Healthy Controls") & (visit_fac == "1 week"))) %>% 
  select(visit_fac, dx_fac, GFAP,IL6,NfL,Tau,pTau181,S100B,UCHL1) %>% 
  pivot_longer(c(-visit_fac, -dx_fac))

# EXAMPLE MODEL
dat_apriori_adults %>% filter(visit_fac == "1 week" & name == "GFAP") %>% 
  t.test(value ~ dx_fac, data = .)

# EXAMPLE CHECK VALID VALUES PER GROUP
dat_apriori_adults %>% 
  filter(visit_fac == "1 week" & name == "UCHL1") %>% view()

dat_apriori_adults %>% 
  filter(visit_fac == "1 week" & name == "UCHL1") %>% 
  select(dx_fac, value) %>% tbl_summary(by = dx_fac) %>% add_overall(last = T)

# RUN T-TESTS FOR EACH BIOMARKER AT 1 WEEK
ttest_models_apriori_adults <- ddply(dat_apriori_adults, .(visit_fac, name), summarise,z=t_test_fun(value,dx_fac)$statistic,
                                     pval=t_test_fun(value,dx_fac)$p.value,
                                     cc_mean=t_test_fun(value,dx_fac)$estimate[1],
                                     mtbi_mean=t_test_fun(value,dx_fac)$estimate[2],
                                     n_total = sum(!is.na(value)))

# CALCULATE FOLD CHANGE
ttest_models_apriori_adults$log2fc <- ttest_models_apriori_adults$mtbi_mean - ttest_models_apriori_adults$cc_mean

DT::datatable(ttest_models_apriori_adults, 
              extensions = c('FixedColumns',"FixedHeader"),
              options = list(scrollX = TRUE, 
                             paging=FALSE,
                             fixedHeader=TRUE))

################################################
## 3.2. Discovery Between-Subjects Comparison ##
################################################

# SELECT BIOMARKERS AND RESHAPE DATAFRAME
dat_discovery_adults <- adults_df_matched %>%  
  filter(((dx_fac == "mTBI" | dx_fac ==  "Healthy Controls") & (visit_fac == "1 week"))) %>% 
  select(visit_fac, dx_fac, all_of(discovery_vars)) %>% 
  pivot_longer(c(-visit_fac, -dx_fac))

# EXAMPLE MODEL
dat_discovery_adults %>% filter(visit_fac == "1 week" & name == "ACHE") %>% 
  t.test(value ~ dx_fac, data = .)

# RUN T-TESTS FOR EACH BIOMARKER AT 1 WEEK
ttest_models_discovery_adults <- ddply(dat_discovery_adults, .(visit_fac, name), summarise,z=t_test_fun(value,dx_fac)$statistic,
                                       pval=t_test_fun(value,dx_fac)$p.value,
                                       cc_mean=t_test_fun(value,dx_fac)$estimate[1],
                                       mtbi_mean=t_test_fun(value,dx_fac)$estimate[2],
                                       n_total = sum(!is.na(value)))

# CALCULATE FOLD CHANGE
ttest_models_discovery_adults$log2fc <- ttest_models_discovery_adults$mtbi_mean - ttest_models_discovery_adults$cc_mean

# APPLY FDR-CORRECTION
ttest_models_discovery_adults <- ttest_models_discovery_adults %>%
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))

DT::datatable(ttest_models_discovery_adults, 
              extensions = c('FixedColumns',"FixedHeader"),
              options = list(scrollX = TRUE, 
                             paging=FALSE,
                             fixedHeader=TRUE))


#########################################################################
## 3.3. Comparing mTBI patients at 4 months versus controls - a priori ##
#########################################################################

# SELECT BIOMARKERS AND RESHAPE DATAFRAME
# FOR THIS ANALYSIS, WE COMPARE CONTROLS AT 1 WEEK TO mTBI PATIENTS AT 4 MONTHS
dat_apriori_adults_followup <- adults_df_matched %>%
  filter(
    (dx_fac == "mTBI" & visit_fac == "4 months") |
      (dx_fac == "Healthy Controls" & visit_fac == "1 week")) %>%
  select(visit_fac, dx_fac, GFAP, IL6, NfL, Tau, pTau181, S100B, UCHL1) %>%
  pivot_longer(cols = c(GFAP, IL6, NfL, Tau, pTau181, S100B, UCHL1))

# EXAMPLE MODEL
dat_apriori_adults_followup %>% filter(name == "Tau") %>% 
  t.test(value ~ dx_fac, data = .)

# EXAMPLE CHECK VALID VALUES PER GROUP
dat_apriori_adults_followup %>% 
  filter(name == "UCHL1") %>% view()

dat_apriori_adults_followup %>% 
  filter(name == "UCHL1") %>% 
  select(dx_fac, value) %>% tbl_summary(by = dx_fac) %>% add_overall(last = T)


# RUN T-TESTS FOR EACH BIOMARKER
ttest_models_apriori_adults_followup <- ddply(dat_apriori_adults_followup, .(name), summarise,
                                              z = t_test_fun(value, dx_fac)$statistic,
                                              pval = t_test_fun(value, dx_fac)$p.value,
                                              cc_mean = t_test_fun(value, dx_fac)$estimate[1],
                                              mtbi_mean = t_test_fun(value, dx_fac)$estimate[2],
                                              n_total = sum(!is.na(value)))

# CALCULATE FOLD CHANGE
ttest_models_apriori_adults_followup$log2fc <- 
  ttest_models_apriori_adults_followup$mtbi_mean - ttest_models_apriori_adults_followup$cc_mean


DT::datatable(ttest_models_apriori_adults_followup, 
              extensions = c('FixedColumns', "FixedHeader"),
              options = list(scrollX = TRUE, paging = FALSE, fixedHeader = TRUE))

##########################################################################
## 3.4. Comparing mTBI patients at 4 months versus controls - discovery ##
##########################################################################

# SELECT BIOMARKERS AND RESHAPE DATAFRAME
# FOR THIS ANALYSIS, WE COMPARE CONTROLS AT 1 WEEK TO mTBI PATIENTS AT 4 MONTHS
dat_discovery_adults_followup <- adults_df_matched %>%
  filter(
    (dx_fac == "mTBI" & visit_fac == "4 months") |
      (dx_fac == "Healthy Controls" & visit_fac == "1 week")) %>%
  select(visit_fac, dx_fac, all_of(discovery_vars)) %>%
  pivot_longer(cols = c(-visit_fac, -dx_fac))

# EXAMPLE MODEL
dat_discovery_adults_followup %>% filter(name == "ACHE") %>% 
  t.test(value ~ dx_fac, data = .)

# RUN T-TESTS FOR EACH BIOMARKER 
ttest_models_discovery_adults_followup <- ddply(dat_discovery_adults_followup, .(name), summarise,
                                                z = t_test_fun(value, dx_fac)$statistic,
                                                pval = t_test_fun(value, dx_fac)$p.value,
                                                cc_mean = t_test_fun(value, dx_fac)$estimate[1],
                                                mtbi_mean = t_test_fun(value, dx_fac)$estimate[2],
                                                n_total = sum(!is.na(value)))

# CALCULATE FOLD CHANGE
ttest_models_discovery_adults_followup$log2fc <- 
  ttest_models_discovery_adults_followup$mtbi_mean - ttest_models_discovery_adults_followup$cc_mean

# APPLY FDR CORRECTION
ttest_models_discovery_adults_followup <- ttest_models_discovery_adults_followup %>%
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))

DT::datatable(ttest_models_discovery_adults_followup, 
              extensions = c('FixedColumns', "FixedHeader"),
              options = list(scrollX = TRUE, paging = FALSE, fixedHeader = TRUE))

#################################################
##  3.5. Sensitivity Analysis: Exclude CT/MRI+ ##
#################################################

adults_df_matched_ct_mri_neg <- adults_df_matched %>%
  filter(!URSI %in% (filter(., ct_mri_pos == 1) %>% pull(URSI)))

dat_apriori_adults_sens <- adults_df_matched_ct_mri_neg %>%  filter(((dx_fac == "mTBI" | dx_fac ==  "Healthy Controls") & (visit_fac == "1 week"))) %>% select(visit_fac, dx_fac, GFAP,IL6,NfL,Tau,pTau181,S100B,UCHL1) %>% pivot_longer(c(-visit_fac, -dx_fac))

ttest_models_apriori_adults_sens <- ddply(dat_apriori_adults_sens, .(visit_fac, name), summarise,z=t_test_fun(value,dx_fac)$statistic,
                                          pval=t_test_fun(value,dx_fac)$p.value,
                                          cc_mean=t_test_fun(value,dx_fac)$estimate[1],
                                          mtbi_mean=t_test_fun(value,dx_fac)$estimate[2],
                                          n_total = sum(!is.na(value)))

ttest_models_apriori_adults_sens$log2fc <- ttest_models_apriori_adults_sens$mtbi_mean - ttest_models_apriori_adults_sens$cc_mean

DT::datatable(ttest_models_apriori_adults_sens, 
              extensions = c('FixedColumns',"FixedHeader"),
              options = list(scrollX = TRUE, 
                             paging=FALSE,
                             fixedHeader=TRUE))

dat_discovery_adults_sens <- adults_df_matched_ct_mri_neg %>%  
  filter(((dx_fac == "mTBI" | dx_fac ==  "Healthy Controls") & (visit_fac == "1 week"))) %>% 
  select(visit_fac, dx_fac, all_of(discovery_vars)) %>% 
  pivot_longer(c(-visit_fac, -dx_fac))

ttest_models_discovery_adults_sens <- ddply(dat_discovery_adults_sens, .(visit_fac, name), summarise,z=t_test_fun(value,dx_fac)$statistic,
                                            pval=t_test_fun(value,dx_fac)$p.value,
                                            cc_mean=t_test_fun(value,dx_fac)$estimate[1],
                                            mtbi_mean=t_test_fun(value,dx_fac)$estimate[2],
                                            n_total = sum(!is.na(value)))

ttest_models_discovery_adults_sens$log2fc <- ttest_models_discovery_adults_sens$mtbi_mean - ttest_models_discovery_adults_sens$cc_mean


ttest_models_discovery_adults_sens <- ttest_models_discovery_adults_sens %>%
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))

DT::datatable(ttest_models_discovery_adults_sens, 
              extensions = c('FixedColumns',"FixedHeader"),
              options = list(scrollX = TRUE, 
                             paging=FALSE,
                             fixedHeader=TRUE))

ttest_models_discovery_adults %>%
  left_join(
    ttest_models_discovery_adults_sens %>% 
      select(visit_fac, name, log2fc_sens = log2fc, pval_fdr_sens = pval_fdr),
    by = c("visit_fac", "name")) %>% 
  select(visit_fac, name, log2fc, pval_fdr, log2fc_sens, pval_fdr_sens) %>%
  mutate(
    log2fc_abs_diff = abs(log2fc_sens - log2fc)) %>% 
  arrange(desc(log2fc_abs_diff)) %>% slice_head(n = 50)

dat_apriori_adults_followup_sens <- adults_df_matched_ct_mri_neg %>%
  filter(
    (dx_fac == "mTBI" & visit_fac == "4 months") |
      (dx_fac == "Healthy Controls" & visit_fac == "1 week")) %>%
  select(visit_fac, dx_fac, GFAP, IL6, NfL, Tau, pTau181, S100B, UCHL1) %>%
  pivot_longer(cols = c(GFAP, IL6, NfL, Tau, pTau181, S100B, UCHL1))

# Run t-tests
ttest_models_apriori_adults_followup_sens <- ddply(dat_apriori_adults_followup_sens, .(name), summarise,
                                                   z = t_test_fun(value, dx_fac)$statistic,
                                                   pval = t_test_fun(value, dx_fac)$p.value,
                                                   cc_mean = t_test_fun(value, dx_fac)$estimate[1],
                                                   mtbi_mean = t_test_fun(value, dx_fac)$estimate[2],
                                                   n_total = sum(!is.na(value)))

ttest_models_apriori_adults_followup_sens$log2fc <- 
  ttest_models_apriori_adults_followup_sens$mtbi_mean - ttest_models_apriori_adults_followup_sens$cc_mean

DT::datatable(ttest_models_apriori_adults_followup_sens, 
              extensions = c('FixedColumns', "FixedHeader"),
              options = list(scrollX = TRUE, paging = FALSE, fixedHeader = TRUE))

ttest_models_apriori_adults_followup %>% 
  left_join(ttest_models_apriori_adults_followup_sens %>% 
              select(name, log2fc_sens = log2fc, pval_sens = pval), by = "name") %>% 
  select(name, log2fc, pval, log2fc_sens, pval_sens) %>% 
  mutate(log2fc_abs_diff = abs(log2fc_sens - log2fc)) %>% 
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x)))

dat_discovery_adults_followup_sens <- adults_df_matched_ct_mri_neg %>%
  filter(
    (dx_fac == "mTBI" & visit_fac == "4 months") |
      (dx_fac == "Healthy Controls" & visit_fac == "1 week")) %>%
  select(visit_fac, dx_fac, all_of(discovery_vars)) %>%
  pivot_longer(cols = c(-visit_fac, -dx_fac))

# Run t-tests
ttest_models_discovery_adults_followup_sens <- ddply(dat_discovery_adults_followup_sens, .(name), summarise,
                                                     z = t_test_fun(value, dx_fac)$statistic,
                                                     pval = t_test_fun(value, dx_fac)$p.value,
                                                     cc_mean = t_test_fun(value, dx_fac)$estimate[1],
                                                     mtbi_mean = t_test_fun(value, dx_fac)$estimate[2],
                                                     n_total = sum(!is.na(value)))

ttest_models_discovery_adults_followup_sens$log2fc <- 
  ttest_models_discovery_adults_followup_sens$mtbi_mean - ttest_models_discovery_adults_followup_sens$cc_mean

# FDR correction
ttest_models_discovery_adults_followup_sens <- ttest_models_discovery_adults_followup_sens %>%
  mutate(pval_fdr = p.adjust(pval, method = "fdr"))

DT::datatable(ttest_models_discovery_adults_followup_sens, 
              extensions = c('FixedColumns', "FixedHeader"),
              options = list(scrollX = TRUE, paging = FALSE, fixedHeader = TRUE))

ttest_models_discovery_adults_followup %>% 
  left_join(ttest_models_discovery_adults_followup_sens %>% 
              select(name, log2fc_sens = log2fc, pval_sens_orig = pval, pval_fdr_sens = pval_fdr), by = "name") %>% 
  select(name, log2fc, pval_fdr, log2fc_sens, pval_sens_orig, pval_fdr_sens) %>% 
  mutate(log2fc_abs_diff = abs(log2fc_sens - log2fc)) %>% 
  arrange(desc(log2fc)) %>% 
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x)))
