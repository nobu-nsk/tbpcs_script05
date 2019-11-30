library(haven)
library(labelled)
library(tidyverse)
library(magrittr)
library(survey)
library(binom)
library(finalfit)
library(pivottabler)
library(knitr)
library(kableExtra)
library(expss)

## functions ####
numstyle <- function(x){
  sapply(x, function(y) prettyNum(y,big.mark =" ",digits=2,nsmall=ifelse(y<9.95,1,0)))
}
summary_survey <- function(.data, dependent, explanatory, cont = "mean", cont_range = FALSE, 
                           p = FALSE, column = FALSE, digits = c(1, 1, 3, 1), fit_id = FALSE){
  # Define variable type
  explanatory_type = .data$variables %>% 
    select(explanatory) %>% 
    purrr::map(is.numeric)
  
  # Hypothesis test
  if(p){
    p_tests = explanatory %>% 
      purrr::map2(explanatory_type,
                  ~if(!.y){
                    survey::svychisq(as.formula(paste0("~", ., "+", dependent)), .data)$p.value %>% 
                      p_tidy(digits[3])
                  } else if (.y & cont == "mean"){
                    survey::svyttest(as.formula(paste0(.x, "~", dependent)), .data)$p.value %>% 
                      p_tidy(digits[3])
                  } else if (.y & cont == "median"){
                    survey::svyranktest(as.formula(paste0(.x, "~", dependent)), .data)$p.value %>% 
                      p_tidy(digits[3])
                  }
      )
  }
  
  # Output table
  explanatory %>% 
    purrr::map2(explanatory_type,
                ~ if(!.y){
                  survey::svytable(as.formula(paste0("~", .x, "+", dependent)), .data) %>% 
                    as.data.frame(stringsAsFactors = FALSE) %>% 
                    { if(column) {
                      dplyr::group_by(., !! sym(dependent)) %>% 
                        dplyr::mutate(
                          total = sum(Freq),
                          prop = 100 * Freq / total
                        )
                    } else { 
                      dplyr::group_by_at(., 1) %>% 
                        dplyr::mutate(
                          total = sum(Freq),
                          prop = 100 * Freq / total
                        )
                    }
                    } %>% 
                    dplyr::mutate(
                      value = paste0(Freq %>% round_tidy(digits[4]), " (", 
                                     prop %>% round_tidy(digits[4]), ")")
                    ) %>%
                    dplyr::select(-total, -prop, -Freq) %>% 
                    tidyr::pivot_wider(names_from = !! dependent, values_from = value) %>% 
                    dplyr::mutate(
                      label = names(.)[1]
                    ) %>% 
                    dplyr::rename(levels = 1)
                } else {
                  { if(cont == "mean") {
                    survey::svyby(as.formula(paste0("~", .x)), as.formula(paste0("~", dependent)), .data, svymean, na.rm = TRUE) %>%
                      dplyr::mutate(
                        value = paste0(!! sym(.x) %>% round_tidy(digits[1]), " (", 
                                       se %>%  round_tidy(digits[2]), ")") 
                      ) %>% 
                      dplyr::select(-c(.x, se)) %>% 
                      tidyr::pivot_wider(names_from = !! dependent, values_from = value) %>% 
                      dplyr::mutate(
                        label = .x,
                        levels = "Mean (SE)"
                      )
                  } else if(cont == "median"){
                    survey::svyby(as.formula(paste0("~", .x)), as.formula(paste0("~", dependent)), .data, 
                                  svyquantile, quantiles = 1:3/4, ci = TRUE, na.rm = TRUE) %>% 
                      dplyr::rename(Q1 = 2,
                                    Q2 = 3,
                                    Q3 = 4) %>% 
                      dplyr::mutate(
                        IQR = Q3 - Q1) %>% 
                        { if(cont_range){
                          dplyr::mutate(., 
                                        value = paste0(Q2 %>% round_tidy(digits[1]), " (", 
                                                       Q1 %>%  round_tidy(digits[2]), " to ",
                                                       Q3 %>% round_tidy(digits[2]), ")")
                          )                        
                        } else {
                          dplyr::mutate(.,
                                        value = paste0(Q2 %>% round_tidy(digits[1]), " (", 
                                                       IQR %>%  round_tidy(digits[2]), ")")
                          )
                        }} %>% 
                      dplyr::select(-c(2:8)) %>%                      
                      tidyr::pivot_wider(names_from = !! dependent, values_from = value) %>% 
                      dplyr::mutate(
                        label = .x,
                        levels = "Median (IQR)"
                      )
                  }
                  }
                }) %>% 
    
    # Add hypothesis test
                { if(p){
                  purrr::map2_df(., p_tests,
                                 ~ mutate(.x,
                                          p = .y)
                  )} else {
                    dplyr::bind_rows(.)
                  }} %>%
    dplyr::select(label, levels, dplyr::everything()) %>% 
    as.data.frame() %>%
    { if(fit_id){
      levels_id = .$levels
      drop = levels_id %in% c("Mean (SE)", "Median (IQR)")
      levels_id[drop] = ""
      mutate(., 
             fit_id = paste0(label, levels_id),
             index = 1:dim(.)[1])
    } else {
      .
    }} %>% 
    rm_duplicate_labels()
}
summary_cont <- function(.data, dependent = NULL, explanatory, ...){
  sum_mean <- summary_factorlist(.data, dependent, explanatory, cont="mean", cont_cut=2, ...) %>% 
    rowid_to_column("order")
  sum_meadian <- 
    summary_factorlist(.data, dependent, explanatory, cont="median", cont_range=TRUE, cont_cut=2, ...) %>% 
    mutate(label="") %>% 
    rowid_to_column("order") %>% mutate(order=order+0.5)
  rbind(sum_mean,sum_meadian) %>% arrange(order) %>% select(-order)
}
meansd <- function(x,digits=1){paste0(round(mean(x),digits)," (",round(sd(x),digits),")")}
mediqr <- function(x){paste0(round(median(x),1)," (",round(quantile(x,0.25),1),"-",round(quantile(x,0.75),1),")")}
n_perc <- function(x,digits=1){paste0(sum(x)," (",round(mean(x)*100,digits),")")}

#### load dataset ####
d <- read_dta("dta/TBPCS_iso3_2019_imputed.dta")
#d <- read_dta("dta/TZA_TBPCS_2019_imputed - usd 27.11.19.dta")

#### variable manupulations ####
d <- d %>% 
  replace_na(list(visit_before=0, n_dot_visits=0, n_pickup_total=0, s_fu=0)) %>% 
  mutate(total_visits=rowSums(dplyr::select(., visit_before, n_dot_visits, n_pickup_total, s_fu), na.rm = TRUE)) %>% 
#  select(visit_before,n_dot_visits,n_pickup_total,s_fu,total_visits) %>% 
  mutate(age_c = cut(age,breaks=c(0,15,25,35,45,55,65,Inf),right=FALSE,
                     labels = c("0-14","15-24","25-34","35-44","45-54","55-64",">65"))) %>% 
  mutate(mdr=factor(mdr,labels=c("DS-TB","DR-TB"))) %>% 
  mutate(cc1_f=factor(cc1, labels=c('No','Yes'))) %>% 
  mutate(educ_level=factor(educ_level)) %>% 
#  mutate(below_poverty=NA) # placeholder for poverty indicator 
  mutate(below_poverty=factor(below_poverty)) %>% 
  mutate(hh_quintile=factor(hh_quintile, 
                            labels=c("Poorest","Second","Third","Fourth","Wealthiest"))) %>% 
  mutate(hh_quintile_r=fct_rev(hh_quintile)) %>% 
  mutate(prev_hosp=factor(prev_hosp, labels=c('No','Yes'))) 
#  mutate(___=factor(___, levels=c('__','__'))) %>% 
# mutate(cc1=factor(cc1, labels=c('No','Yes'))) %>% 
# mutate(hh_quintile=factor(hh_quintile,labels=c("Poorest","Second","Third","Fourth","Wealthiest"))) %>% 
# mutate(hh_quintile_r=rev(hh_quintile)) %>% 
# mutate(prev_hosp=factor(prev_hosp)) %>%

d = apply_labels(d,
                 age_c = "Age group",
                 prev_hosp = "Hospitalized during current phase",
                 sex = "Sex",
                 mdr = "Drug-suceptible status",
                 hiv = "Recorded HIV status",
                 hh_quintile_r = "Income quintile",
                 self_admin = "Treatment modality"
)


# arsenal: setting for tabulation parameters for arsenal package  ####
my_controls <- tableby.control(
  test = T,
  total = T,
  numeric.test = "kwt", cat.test = "chisq",
  numeric.stats = c("meansd", "medianq1q3"), 
  cat.stats = c("countpct"),
  stats.labels = list(
    meansd = "Mean (SD)",
    medianq1q3 = "Median (Q1, Q3)",
    #    range = "Min - Max",
    Nmiss2 = "Missing"
  )
)

# descriptive stat table stored in tab1 object
tab1 <- tableby( ~ bw + lbw + gender + mage + weeks + premie + visits + marital + gained_kg + smoking,
                 data = d,
                 control = my_controls
)

#
tab1_sum <- summary(tab1, title = "Table 1. Characteristics", digits=1)
write2html(tab1_sum,"Table1.html")
write2word(tab1_sum,"Table1.docx")



d %>% group_by(f)

table(d$facility)
table(d$district)

table(d$cc1)

mean(d$cc1)


## svydesing settings #######
## adjustment 
## 
## survey design setting 
## 
## srs / cluster_design / cluster_design_wt

# simple binomial proportion 
t <- table(d$cc1) %>% addmargins(margin = 1)
binom.confint(t[2],t[3], method='asymptotic')

# assuming simple random sampling (theoritically results should be same as above)
srs <- svydesign(id=~1, data=d, probs = 1) ## assuming simple random sampling 
svytotal(~cc1,srs)
svytable(~cc1,srs)
svyciprop(~I(cc1==1), design=srs, method=c('asin'))
cbind((t <- svymean(~total_cost_hc,srs)),confint(t))

# setting cluster sampling design (district as cluster)
cluster_design <- svydesign(id=~district, data=d, probs = 1) # set district as cluster sampling unit 
svytotal(~cc1,cluster_design)
svytable(~cc1,cluster_design)
svyciprop(~I(cc1==1), design=cluster_design, method=c('asin'))
cbind((t <- svymean(~total_cost_hc,cluster_design)),confint(t))

svyby(~cc1,~mdr,cluster_design,svyciprop,vartype="ci")

# setting cluster sampling design (district as cluster) and weighted for non-response 
#All TB cases notified in 2018: 75828
#DR-TB Patients started on treatment in 2018: 409
#DS-TB in 2018: 75828-409=75419
# DR-TB: 0.539%  DS-TB: 99.461%
tball <- 75828
csize <- 26
d <- d %>% 
  group_by(district) %>% summarize(actualsize=n()) %>% 
  mutate(nc = round(actualsize/csize,0)) %>% 
  mutate(sampprov = actualsize/(csize*nc)) %>% 
  mutate(wt = 1/sampprov) %>% select(district, sampprov,wt) %>% # create data.frame with district and weight
  left_join(d,.) # join with original d 
d %>% ggplot(aes(x=sampprov)) + geom_histogram() #check distribution of sampling probabilities 

# setting cluster sampling design (district as cluster) with wt for non-response adjustment (under- or over-enrolment)
cluster_design_wt <- svydesign(id=~district, data=d, weights=~wt) 

svyciprop(~I(cc1==1), design=cluster_design_wt, method=c('asin'))
cbind((t <- svymean(~total_cost_hc,cluster_design_wt)),confint(t))
svyquantile(~total_cost_hc,cluster_design_wt,quantiles=c(0.5,0.25,0.75))


## Table 2 (svy) ####
## 
dependent = "mdr"
explanatory = c("sex","tb_type","phase","hiv","income_hh_pre_annual")
explanatory = c("sex","tb_type","phase","hiv","treat_group_nmdr","treat_group_mdr","self_admin",
                "income_hh_pre_annual_reported","income_hh_pre_annual","expend_hh_annual",
                "below_poverty")
Table2_svy <- 
  d %>% 
  mutate(sex = ff_label(sex, "Sex")) %>% 
  #   hiv = ff_lable(hiv,"Recorded HIV status")
  #  ) %>% 
  svydesign(ids=~1,data=.) %>% 
  summary_survey(dependent, c(expl_1,expl_2,expl_3), cont = "mean", cont_range = FALSE, p = FALSE, 
                 column = TRUE, digits = c(1, 1, 3, 0)) %>% 
  rename("Variable"=1," "=2,"DS-TB"=3,"DR-TB"=4)

kable(Table2_svy, "html", caption = "Table 2. Descriptive statistics of survey sample by MDR status",
      align = 'llcc', na_include = TRUE) %>%
  kable_styling(bootstrap_options = c("striped"), full_width=F,position='left') %T>%
  cat(., file = "Table 2. Descriptive statistics.html")

## Table 2 ####
## Tabulation by using 
## finalfit::summary_factorlist 
##
expl_1 = c("sex","age","age_c","tb_type","phase","hiv","self_admin")#,"treat_group_nmdr","treat_group_mdr")
expl_2 = c("income_hh_pre_annual_reported","income_hh_pre_annual","expend_hh_annual")
expl_3 = c("below_poverty","educ_level")

Table2_1 <- 
  d %>% 
  mutate(
    age = ff_label(age, "Age (years)"),
    age_c = ff_label(age_c, "Age group"),
    sex = ff_label(sex, "Sex"),
         tb_type = ff_label(tb_type, "Type of TB"),
         phase = ff_label(phase, "Treatment phase"),
         hiv = ff_label(hiv, "Recorded HIV status"),
         self_admin = ff_label(self_admin, "DOT")) %>% 
  summary_factorlist(dependent,expl_1,cont = "mean",cont_range=TRUE,p=FALSE,total_col = TRUE,
                     na_include=TRUE,na_to_missing=TRUE,digits = c(1, 1, 3, 1))
Table2_2_mean <- 
  d %>% 
  summary_factorlist(dependent,expl_2,cont = "mean",p=FALSE,total_col = TRUE,digits = c(0, 0, 3, 0))
Table2_2_mean

Table2_2_median <- 
  d %>% 
  summary_factorlist(dependent,expl_2,cont = "median",cont_range=TRUE,p=FALSE,total_col = TRUE,
                     digits = c(0, 0, 3, 0)) %>% 
  mutate(label="")
Table2_2_median
table_2_2 <- rbind(Table2_2_mean[1,],Table2_2_median[1,],
                 Table2_2_mean[2,],Table2_2_median[2,],
                 Table2_2_mean[3,],Table2_2_median[3,])
Table2_3 <- 
  d %>% 
  mutate(below_poverty = ff_label(below_poverty, "Below poverty line")) %>% 
  summary_factorlist(dependent,expl_3,p=FALSE,total_col = TRUE,na_include=TRUE,na_to_missing=TRUE,
                     digits = c(1, 1, 3, 1))

Table2 <- rbind(Table2_1,table_2_2,Table2_3) %>% remove_rownames()
  
kable(Table2, "html", caption = "Table 2. Descriptive statistics of survey sample by drug-resistance status",
      align = 'llccc', na_include = TRUE) %>%
  kable_styling(bootstrap_options = c("striped"), full_width=F,position='left') %T>%
  cat(., file = "Table 2. Descriptive statistics.html")


## #### Table 2
## Tabulation by using 
## finalfit::summary_survey
##

expl_1 = c("sex","age","age_c","tb_type","phase","hiv","self_admin")#,"treat_group_nmdr","treat_group_mdr")
expl_2 = c("income_hh_pre_annual_reported","income_hh_pre_annual","expend_hh_annual")
expl_3 = c("below_poverty","educ_level")

Table2a_1 <- 
  cluster_design_wt %>% 
  summary_survey(dependent,expl_1,cont = "mean",digits = c(1, 1, 3, 1))
#                 total_col = TRUE,na_include=TRUE,na_to_missing=TRUE)
Table2a_2_mean <- 
  cluster_design_wt %>% 
  summary_survey(dependent,expl_2,cont = "mean",digits = c(0, 0, 3, 0))

Table2a_2_median <- 
  cluster_design_wt %>% 
  summary_survey(dependent,expl_2,cont = "median",cont_range=TRUE,digits = c(0, 0, 3, 0)) %>% 
  mutate(label="")
table_2a_2 <- rbind(Table2a_2_mean[1,],Table2a_2_median[1,],
                   Table2a_2_mean[2,],Table2a_2_median[2,],
                   Table2a_2_mean[3,],Table2a_2_median[3,])
Table2a_3 <- 
  cluster_design_wt %>% 
  summary_survey(dependent,expl_3,p=FALSE,digits = c(1, 1, 3, 1))

Table2a <- rbind(Table2a_1,table_2a_2,Table2a_3) %>% remove_rownames()

kable(Table2a, "html", caption = "Table 2. Descriptive statistics of survey sample by drug-resistance status",
      align = 'llccc', na_include = TRUE) %>%
  kable_styling(bootstrap_options = c("striped"), full_width=F,position='left') %T>%
  cat(., file = "Table 2a. Descriptive statistics (weighted).html")

## Table 3 ####
## Tabulation by using 
## finafit::summary_factorlist 
##

d <- d %>% mutate(
  current_hosp   = ff_label(current_hosp,  "Hospitalized at time of interview"),
  prev_hosp      = ff_label(prev_hosp,     "Hospitalized during current phase"),
  t_stay_current = ff_label(t_stay_current,"Days hospitalized during current phase"),
  total_visits = ff_label(total_visits,"Number of visits per episode: Total"),
  n_dot_visits = ff_label(n_dot_visits,"Number of visits: DOT"),
  s_fu = ff_label(s_fu,"Number of visits: follow-up"),
  n_pickup_total = ff_label(n_pickup_total,"Number of visits: Drug pick-up"),
  visit_before = ff_label(visit_before,"Number of visits: Pre-diagnosis"),
  #  ___ = ff_label(_____,"Number of visits pre-diagnosis (non-public facility), mean"),
  duration_int = ff_label(duration_int,"Intensive phase (month)"),
  duration_cont = ff_label(duration_cont,"Continuation phase (month)"),
  weeks_before_tx = ff_label(weeks_before_tx,"Treatment delay (week)")
)

dependent = "mdr"
expl_1 = c("current_hosp","prev_hosp","t_stay_current")
expl_2 = c("total_visits","n_dot_visits","s_fu","n_pickup_total","visit_before")
expl_3 = c("duration_int","duration_cont")
expl_4 = c("weeks_before_tx")

Table3_1 <-   
  d %>% 
  summary_factorlist(dependent,expl_1,cont = "mean",cont_range=TRUE,p=FALSE, total_col = TRUE,
                     na_include=TRUE,na_to_missing=TRUE,digits = c(1, 1, 3, 1)) %>% 
  ff_remove_ref() %>% mutate(levels=ifelse(levels=="1","n (%)",levels))

Table3_2 <- d %>% 
  summary_cont(dependent,expl_2,p=FALSE,total_col = TRUE,digits = c(1, 1, 3, 0))

Table3_3 <- 
  d %>% summary_factorlist(dependent,expl_3,cont_cut=3,cont="mean",p=FALSE,total_col = TRUE,digits = c(1, 1, 3, 0))

Table3_4 <- 
  summary_cont(d, dependent,expl_4,p=FALSE,total_col = TRUE,digits = c(1, 1, 3, 0))

Table3 <- rbind(Table3_1,Table3_2,Table3_3,Table3_4) %>% remove_rownames() %>% 
  rename(" "=1," "=2)
kable(Table3, "html", caption = "Table 3. Model of care for survey sample",
      align = 'llccc', na_include = TRUE) %>%
  kable_styling(bootstrap_options = c("striped"), full_width=F,position='left') %>% 
  pack_rows("Hospitalisation",1,3) %>% 
  pack_rows("Ambulatory care",4,13) %>% 
  pack_rows("Treatment duration",14,15) %>% 
  pack_rows("Treatment delay",16,17) %T>%
  cat(., file = "Table 3. Model of care for survey sample.html")

## Table 4. Hours lost by seeking and accessing care ####
## Tabulation by using 
## pivottabler:: PivotTable
##
d <- d %>% mutate(t_before=replace_na(t_before,0),      # replace na in summary time variables 
                  t_current=replace_na(t_current,0)) #%>% select(t_before,t_current)

pt <- PivotTable$new()
pt$addData(d)
pt$addRowDataGroups("phase", totalCaption="Row total")
pt$addColumnDataGroups("mdr", totalCaption="Total")
pt$defineCalculation(calculationName="Number", caption="n",summariseExpression="n()")
pt$defineCalculation(calculationName="Mean", caption="Mean (SD)",summariseExpression="meansd(t_before)")
pt$defineCalculation(calculationName="Median", caption="Median (IQR)",summariseExpression="mediqr(t_before)")
pt$renderPivot()
table1 <- pt$asMatrix()

pt <- PivotTable$new()
pt$addData(d)
pt$addRowDataGroups("phase", totalCaption="Row total")
pt$addColumnDataGroups("mdr", totalCaption="Total")
pt$defineCalculation(calculationName="Number", caption="n", summariseExpression="n()")
pt$defineCalculation(calculationName="Mean", caption="Mean(SD)",summariseExpression="meansd(t_current)")
pt$defineCalculation(calculationName="Median", caption="Median(IQR)",summariseExpression="mediqr(t_current)")
pt$renderPivot()
table2 <- pt$asMatrix()

combined <- rbind(table1,table2)[c(3,8,9),] %>% as.data.frame.matrix(stringsAsFactors=FALSE)
colnames(combined) <- table1[2,]
combined[,1] <- c("Prediagnosis","Intensive phase","Continuation phase")

kable(combined, "html", caption = "Table 4. Hours lost by patient while seeking and accessing care",
      align = 'rccccccccc', na_include = TRUE) %>%
  kable_styling(bootstrap_options = c("striped"), full_width=F,position='left') %>% 
  #  pack_rows("Hours lost by patient",1,3) %>% 
  add_header_above(c(" " = 1, "DS-TB" = 3, "DR-TB" = 3, "Total" = 3)) %T>%
  cat(., file = "Table 4. Hours lost by patient while seeking and accessing care.html")

## alternative method ######
t1 <- d %>% filter(mdr==0,phase=='phase1') %>% select(t_before) %>% 
  summarise_all(list(n=length,mean=meansd,median=mediqr))
t2 <- d %>% filter(mdr==1,phase=='phase1') %>% select(t_before) %>% 
  summarise_all(list(n=length,mean=meansd,median=mediqr))
t3 <- d %>% filter(phase=='phase1') %>% select(t_before) %>% 
  summarise_all(list(n=length,mean=meansd,median=mediqr))

t4 <- d %>% filter(mdr==0) %>% group_by(phase) %>% select(t_current) %>% 
  summarise_all(list(n=length,mean=meansd,median=mediqr))
t5 <- d %>% filter(mdr==1) %>% group_by(phase) %>% select(t_current) %>% 
  summarise_all(list(n=length,mean=meansd,median=mediqr))
t6 <- d %>% filter() %>% group_by(phase) %>% select(t_current) %>% 
  summarise_all(list(n=length,mean=meansd,median=mediqr))

Table4_alt <- rbind(cbind(t1,t2,t3),cbind(t4[,-1],t5[,-1],t6[,-1]))
kable(Table4, "html", caption = "Table 4. Hours lost by seeking and accessing care (alt)",
      align = 'llccc', na_include = TRUE) %>%
  kable_styling(bootstrap_options = c("striped"), full_width=F,position='left')


## Table 5. Estimated costs borne by households affected by TB ####
## Tabulation by using 
##   dplyr
##

mean_sd_med_iqr <- function(x,digits=1){
  paste0(
    round(mean(x),digits)," (",round(sd(x),digits),")__",
    round(median(x),1)," (",round(quantile(x,0.25),1),"-",round(quantile(x,0.75),1),")"
  )
}

cost_cat <- d %>% 
  mutate(cat_medical=cat_before_med+cat_current_med,
         cat_nmedical=cat_before_nmed+cat_current_accomodation+cat_current_travel+
           cat_current_food+cat_current_nutri) %>% 
  select(mdr,cat_before_med,cat_before_nmed,cat_before_indirect,
         cat_current_med,cat_current_travel,cat_current_accomodation,
         cat_current_food,cat_current_nutri,
         cat_medical,
         cat_nmedical,
         cat_current_indirect,
         total_cost_hc)
cost_cat %>% 
  mutate(mdr="Total") %>% 
  group_by(mdr) %>% 
  summarize_all(mean_sd_med_iqr) %>% t() -> cost_total
cost_cat %>% 
  group_by(mdr) %>% 
  summarize_all(mean_sd_med_iqr) %>% t() -> cost_ds_dr
Table5 <- cbind(cost_ds_dr,cost_total)[-1,] %>% as.data.frame.matrix(stringsAsFactors = FALSE) %>%
  separate(V1,into=c("C1","C2"), sep='__') %>% 
  separate(V2,into=c("C3","C4"), sep='__') %>% 
  separate(V3,into=c("C5","C6"), sep='__') %>% 
  rownames_to_column() %>% 
  rename("Category"=1,"Mean (SD)"=2,"Median (IQR)"=3,"Mean (SD)"=4,"Median (IQR)"=5,"Mean (SD)"=6,"Median (IQR)"=7)

kable(Table5, "html", caption = "Table 5. Estimated costs borne by households affected by TB",
      align = 'rcccccc', na_include = TRUE) %>%
  kable_styling(bootstrap_options = c("striped"), full_width=F,position='left') %>% 
  #  pack_rows("Hours lost by patient",1,3) %>% 
  add_header_above(c(" " = 1, "DS-TB" = 2, "DR-TB" = 2, "Total" = 2)) %>% 
  pack_rows("Pre-diagnosis",1,3) %>% 
  pack_rows("Post-diagnosis",4,8) %>% 
  pack_rows("Subtotal",9,11) %T>%
  pack_rows("Total",12,12) %T>%
  cat(., file = "Table 5. Estimated costs borne by households affected by TB.html")

## Table 6. Reported coping mechanisms and social consequences ####
## Tabulation by using 
##  dplyr
##
coping_quint <- d %>% 
  mutate(impact1=as.numeric(cope_impact=='cope_impact1'),impact2=as.numeric(cope_impact=='cope_impact2'),
         impact3=as.numeric(cope_impact=='cope_impact3'),impact4=as.numeric(cope_impact=='cope_impact4')) %>% 
  select(hh_quintile,borrow,asset_sale,coping,
         social_food_insec,social_divorce,social_lossofjob,social_dropout,social_exclusion,
         days_lost,impact1,impact2,impact3,impact4) 
coping_quint %>% 
  group_by(hh_quintile) %>% 
  summarise_all(n_perc) %>% t() -> coping_brk
coping_quint %>% 
  group_by(hh_quintile="Total") %>% 
  summarise_all(n_perc) %>% t() -> coping_tot
Table6 <- cbind(coping_brk,coping_tot)[-1,] %>% as.data.frame.matrix() %>% rownames_to_column()
quint_lab <- paste0(c("Lowest","Second","Third","Fourth","Highest","Total"),
                    " (N=",c(table(d$hh_quintile),length(d$hh_quintile)),")")
names(Table6) <- c(" ",quint_lab)

kable(Table6, "html", caption = "Table 6. Reported dissaving mechanisms and social consequences",
      align = 'rcccccc', na_include = TRUE) %>%
  kable_styling(bootstrap_options = c("striped"), full_width=F,position='left') %>% 
  #  pack_rows("Hours lost by patient",1,3) %>% 
  add_header_above(c(" " = 1, "Income Quintile" = 5, " "=1)) %>% 
  pack_rows("Coping strategies",1,3) %>% 
  pack_rows("Social consequences",4,9) %>% 
  pack_rows("Perceived impact",10,13) %T>%
#  pack_rows("Total",12,12) %T>%
  cat(., file = "Table 6. Reported dissaving mechanisms and social consequences.html")





str(Table6)




table(d$cope_impact)

d$cope_impact==


%>% 

  
  
    summary

table(d$cope_impact)

-> coping_brk



coping_quint



cope_impact2


=="cope_impact2"))
         
         ,(impact==2))
         
         ,impact3,impact4)

table(d$cope_impact)

d$borrow
         asset_sale1
         borrow1
         coping1
         
         
         
         
         social_food_insec1
         


## ####
## Figure 
##        
##

# var=d$pct1
# from=0
# to=2
# graph=TRUE

threshold_graph <- function(var, from=0, to=2, threshold=c(1:6/10), graph=TRUE){
  x=seq(from=from,to=to,by=1/100)
  y=NULL
  n=length(var)
  for(i in x){
    y=c(y,sum(var>i)/n)
  }
  thcurve <- tibble(x=x,y=y)
  cc_by_threshold <- thcurve[x %in% threshold,] %>% rename(threshold=1,cc=2) %>% 
  mutate(thr_lab = factor(paste0(threshold*100,"%")))
    
  pathdata <- 
    cc_by_threshold %>% mutate(x1=0,y1=cc,x2=threshold,y2=cc,x3=threshold,y3=0) %>% 
    pivot_longer(x1:y3,names_to=c(".value","order"),names_pattern = "(.)(.)") %>% 
    mutate(size=ifelse(threshold==0.2,1.2,0.6))
  
  p <- 
    ggplot(thcurve,aes(x=x,y=y)) + geom_line(size=0.6,col='darkblue') +
    scale_color_brewer(palette = "Spectral",direction = 1) +
    scale_x_continuous(breaks=0:4/2,labels=c("0%","50%","100%","150%","200%"),limits=c(0,2))+
    scale_y_continuous(breaks=0:10/10,labels=paste0(0:10*10,"%"))+
    labs(title="Proportion of households faced catastrophic costs by different threshold",
         x="Threadshold as percentage of annual household income/expenditure/consumption (%)",
         y="Proportion of households faced catastrophic costs (%)") + 
    scale_size_identity() +
    geom_path(data=pathdata,aes(x=x,y=y,col=thr_lab,size=size), 
              arrow = arrow(ends = "first", type = "open",angle=15,length=unit(0.1, "inches")))
  if(graph==TRUE) print(p)
  print(cc_by_threshold[,c('thr_lab','cc')])
  invisible(list(plot=p, thcurve=thcurve, cc_by_threshold=cc_by_threshold))
}

#gh <- read_dta(file = "../tbpcs_ghana/ghana_pcs_clean.dta")
mmr <- read_dta(file = "../tbpcs_myanmar/mmr_pcs_clean.dta")
pct1_mmr <- mmr$pct1[!is.na(mmr$pct1)]
threshold_graph(gh$pct1, from=0, to=2, graph=TRUE)
ggsave("Ghana_threshold_pct1.pdf")

gh %>% select(pct1) %>% 
  arrange(pct1) %>% summary() 
  mutate(rn=row_number())
  ggplot(aes(x=pct1)) + stat_ecdf(geom = "step")


d %>% select(cat_before_med,cat_current_med,
             cat_before_nmed,
             cat_current_travel, cat_current_accomodation, cat_current_food,cat_current_nutri, 
             cat_before_indirect,cat_current_indirect,cat_indirect) %>% View()


d$c_nmed_before_travel
d$c_nmed_before_food
d$c_nmed_before_accommodation
d$c_nmed_before_dis

denom <- d$expend_hh_annual
tc    <- d$total_cost_hc
pct   <- tc/denom

d %>% select(pct1) %>% arrange(pct1) %>% mutate(rn=row_number()) %>% 
  ggplot() + geom_step(aes(x=rn,y=pct1))

d %>% select(pct1) %>% arrange(pct1) %>% mutate(rn=row_number()) %>% 
ggplot(d, aes(pct1)) + stat_ecdf(geom = "step")

summary(d$pct1)
summary(pct)

summary(denom)
summary(tc)

d %>% filter(income_hh_pre_annual < 1800000) %>% 
  ggplot(aes(x=pct1)) + stat_ecdf(geom = "step")
ggplot(d, aes(d$income_hh_pre_annual)) + geom_histogram()

d %>% filter(income_hh_pre_annual < 1800000) %>% 


boxplot(d$income_hh_pre_annual_reported)


library(RColorBrewer)
display.brewer.pal(6, 'YlOrRd')
brewer.pal(6, 'YlOrRd')
display.brewer.all()

d$pct1

rm(list=ls())





##
##
##
##
dependent = c("phase")
expl_1 = c("t_before","t_current")
expl_2 = c("total_visits","n_dot_visits","s_fu","n_pickup_total","visit_before")
expl_3 = c("duration_int","duration_cont")
expl_4 = c("weeks_before_tx")

Table3_1 <-   
  d %>% 
  summary_cont(dependent,expl_1,total_col = TRUE,
               na_include=TRUE,na_to_missing=TRUE,digits = c(1, 1, 3, 1)) 
Table3_2 <- d %>% 
  summary_cont(dependent,expl_2,p=FALSE,total_col = TRUE,digits = c(1, 1, 3, 0))

Table3_3 <- 
  d %>% summary_factorlist(dependent,expl_3,cont_cut=3,cont="mean",p=FALSE,total_col = TRUE,digits = c(1, 1, 3, 0))

Table3_4 <- 
  summary_cont(d, dependent,expl_4,p=FALSE,total_col = TRUE,digits = c(1, 1, 3, 0))

Table3 <- rbind(Table3_1,Table3_2,Table3_3,Table3_4) %>% remove_rownames() %>% 
  rename(" "=1," "=2)
kable(Table3, "html", caption = "Table 3. Model of care for survey sample",
      align = 'llccc', na_include = TRUE) %>%
  kable_styling(bootstrap_options = c("striped"), full_width=F,position='left') %>% 
  pack_rows("Hospitalisation",1,3) %>% 
  pack_rows("Ambulatory care",4,13) %>% 
  pack_rows("Treatment duration",14,15) %>% 
  pack_rows("Treatment delay",16,17) %T>%
  cat(., file = "Table 3. Model of care for survey sample.html")








## Table 8. Risk factors of experiencing catastrophic costs ####
## 
##

fullmodel <- glm(cc1_f ~ age_c + sex + mdr + hiv + hh_quintile + prev_hosp + self_admin +district +facility +weeks_before_tx,
                 data=d, family=binomial())
k_p1=qchisq(0.1,1,lower.tail=FALSE) # define k value corresponding to p=0.05 in Chi-square distribution (one side)
sw_backward <- step(fullmodel,k=k_p1) # backward selection from the full model
sw_backward

dependent = 'cc1_f'
explanatory = c('age_c','sex','mdr','hiv','hh_quintile_r','prev_hosp','self_admin')
explanatory_final = c('age_c','sex','hh_quintile_r','prev_hosp','self_admin')
t_or <- d %>% 
  finalfit(dependent, explanatory, explanatory_final, digits=c(2,2,3)) %>% 
  ff_remove_p() #%>% fit2df()

t_bivar <- d %>%
  summary_factorlist(dependent, explanatory,p=TRUE, column=FALSE,total_col = TRUE, fit_id=FALSE)

#t_ffmerged <- ff_merge(t_bivar,t_or,last_merge = TRUE) does not work
Table8 <- cbind(t_bivar[,c(1,2,5,4,6)],t_or[,c(5,6)]) %>% remove_rownames()
kable(Table8, "html", caption = "Table 8. Risk factors for catastrophic costs") %>%
  kable_styling(bootstrap_options = c("striped"), full_width=F) %T>%
  cat(., file = "Table 8. Risk factors for catastrophic costs.html")





t_bivar_or <- cbind(t_bivar[,c(1,2,5,4,6)],t_or[,c(5,6)])
t_bivar_or
kable(t_bivar_or, "html", caption = "Table 8. Risk factors of experiencing catastrophic costs") %>%
  kable_styling(bootstrap_options = c("striped"), full_width=F) %T>%
  cat(., file = "Table 8. Risk factors of experiencing catastrophic costs.html")


hh_quintile_r=fct_rev(hh_quintile),


str(d$current_hosp)
str(d$prev_hosp)


#sw_backward <- step(fullmodel,k=k_p025) # backward selection from the full model
#summary(sw_backward)
#sw_forward <- step(priorimodel,direction="forward", k=k_p025, # foward selection 
sw_forward <- step(priorimodel,direction="forward", k=k_p1, # foward selection 
                   scope=list(lower=formula(priorimodel),     # starting with priori model 
                              upper=formula(fullmodel)))      # a possibility to extend to full model 





'facility_type'
'prev_hosp')

table(d$facility_type)
table(d$prev_hosp)
table(d$current_hosp)


d

  d$self_admin
  
d %>% select(self_admin,self_admin_int)
d$facility_type
#出力
explanatory       = c("age_study_c", "sex", "area3", "maincaregiver_relationship", "mother_education", "father_education", 
                      "mother_occupation",  "father_occupation", "ses_pc1", "n_watersource_drink", 
                      "water_source_drink_treat", "water_source_bath_treat", "n_watersource_wash", "water_source_wash_treat",
                      "water_contacts", "water_contact2_playing", "stunting", "underweight", "breastfeeding_b",
                      "combined_anemia", "malaria_rdt", "hiv_mother_3c", "hiv_child_3c")
explanatory_multi = c("age_study_c", "sex", "area3")


t_or <- d %>% 
  finalfit(dependent, explanatory, explanatory_multi,digits = c(2,2,3),column = FALSE)
t_bivar <- d %>%
  summary_factorlist(dependent, explanatory, p = TRUE, column = FALSE,total_col = TRUE)

t_bivar_or <- cbind(t_bivar[,c(1,2,5,4,6)],t_or[,c(5,6)])
kable(t_bivar_or, "html", caption = "Table #. Prevalence of Kato-Katz positivity by risk factor") %>%
  kable_styling(bootstrap_options = c("striped", full_width=F)) %T>%
  cat(., file = "Table #. Kato-Katz_tabulation(Finafit).html")
write_csv(t_bivar_or, "Table #. KAto-Katz_tabulation(Finafit).csv")


