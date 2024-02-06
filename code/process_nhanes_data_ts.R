#----------------------------------------------------------------------------
# Process the NHANES Physical Activity Data

# This script includes *some* modifications of code from Andrew Leroux:
# rnhanesdata: NHANES Accelerometry Data Pipeline. https://github.com/andrew-leroux/rnhanesdata.
# Reference: Leroux et al (2019, Stat in Biosciences). https://link.springer.com/article/10.1007/s12561-018-09229-9
# Reference: Smirnova et al. (2019, J Gerontol A Biol Sci Med Sci). https://www.ncbi.nlm.nih.gov/pubmed/31504213
#----------------------------------------------------------------------------
# install rnhanesdata package from github
devtools::install_github("andrew-leroux/rnhanesdata")

#rm(list = ls())
# Aggregate activity counts into intervals of 'k_min' minutes
k_min = 5 
#----------------------------------------------------------------------------
# Import the data:
sapply(
  c("rnhanesdata", "devtools","magrittr","dplyr", "forcats", "reshape2", "ggplot2", "MetBrewer", "gridExtra"                
  ), function(x) require(x, character.only=TRUE)
)

source("code/helper_functions.R")

# Create a (local) temporary directory 
# where lab measurement (cholesterol, blood presure) data will be downloaded from the CDC website 
# and then loaded into R. These files need to be downloaded separately as 
# the raw files associated with these lab measurements are not included in the rnhanesdata package.
dir_tmp = tempfile()
dir.create(dir_tmp)

if (!dir.exists(dir_tmp)){
  dir.create(dir_tmp, showWarnings = FALSE)
}
dl_file = function(url) {
  bn = basename(url) 
  destfile = file.path(dir_tmp, bn)
  if (!file.exists(destfile)) {
    out = download.file(url, destfile = destfile, mode="wb")
  }
  stopifnot(file.exists(destfile))
}

## download the lab measurement data for the cohort 2005-2006
# Total Cholesterol: LBXTC
dl_file("https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/TCHOL_D.XPT")
# HDL Cholesterol: LBDHDD
dl_file("https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/HDL_D.XPT")
# Blood Pressure, up to 4 measurements per person: BPXSY1 , BPXSY2, BPXSY3 and BPXSY4
dl_file("https://wwwn.cdc.gov/Nchs/Nhanes/2005-2006/BPX_D.XPT")


varnames <- c("LBXTC","LBXHDD","LBDHDD",           ## 1. cholesterol. Note LBXHDD and LBDHDD are the same variable, 
              ##    but different names for 2003-2004 and 2005-2006 cohorts
              "BPXSY1","BPXSY2","BPXSY3", "BPXSY4" ## 2. blood pressure measurements
)

## load and merge the lab data
lab_data <- process_covar(varnames=varnames,localpath=dir_tmp)

## combine waves
CVMarkers <- bind_rows(lab_data$Covariate_C, lab_data$Covariate_D)

rm(list=c("lab_data","dir_tmp","varnames"))


## load the data
data("PAXINTEN_D") # '05-06 minute level PAXINTEN_D data (1 = Sunday)
data("Flags_D")
data("Covariate_D")

## re-code activity counts which are considered "non-wear" to be 0
## this doesn't impact many data points, most estimated non-wear times correspond to 0 counts
PAXINTEN_D[,paste0("MIN",1:1440)] <-PAXINTEN_D[,paste0("MIN",1:1440)]*Flags_D[,paste0("MIN",1:1440)]

## Merge covariate and accelerometry data
## note that both PAXINTEN_* and Covariate_* have a column
## called "SDDSRVYR" indicating which NHANES wave the data is associated with.
## To avoid duplicating this column in the merged data, we add this variable to the "by"
## argument in left_join()
AllAct <- left_join(PAXINTEN_D, Covariate_D, by=c("SEQN", "SDDSRVYR"))
AllFlags <- left_join(Flags_D, Covariate_D, by=c("SEQN", "SDDSRVYR"))

#merge with cardiovascular markers 
AllAct <- left_join(AllAct, CVMarkers, by = "SEQN")
AllFlags <- left_join(AllFlags, CVMarkers, by = "SEQN")
rm(list=c("CVMarkers"))

## Create Age in years using the age at examination (i.e. when participants wore the device)
AllAct$Age <- AllFlags$Age <- AllAct$RIDAGEEX/12

## Re-level comorbidities to assign refused/don't know as not having the condition
## Note that in practice this does not affect many individuals, but it is an assumption we're making.
levels(AllAct$CHD)    <- levels(AllFlags$CHD)    <- list("No" = c("No","Refused","Don't know"), "Yes" = c("Yes"))
levels(AllAct$CHF)    <- levels(AllFlags$CHF)    <- list("No" = c("No","Refused","Don't know"), "Yes" = c("Yes"))
levels(AllAct$Stroke) <- levels(AllFlags$Stroke) <- list("No" = c("No","Refused","Don't know"), "Yes" = c("Yes"))
levels(AllAct$Cancer) <- levels(AllFlags$Cancer) <- list("No" = c("No","Refused","Don't know"), "Yes" = c("Yes"))
levels(AllAct$Diabetes) <- levels(AllFlags$Diabetes) <- list("No" = c("No","Borderline", "Refused","Don't know"), "Yes" = c("Yes"))

## Re-level education to have 3 levels and categorize don't know/refused to be missing
levels(AllAct$EducationAdult) <- levels(AllFlags$EducationAdult) <- 
  list("Less than high school" = c("Less than 9th grade", "9-11th grade"),
       "High school" = c("High school grad/GED or equivalent"),
       "More than high school" = c("Some College or AA degree", "College graduate or above"))

## Re-level alcohol consumption to include a level for "missing"
levels(AllAct$DrinkStatus) <- levels(AllFlags$DrinkStatus) <- c(levels(AllAct$DrinkStatus), "Missing alcohol")
AllAct$DrinkStatus[is.na(AllAct$DrinkStatus)] <- AllFlags$DrinkStatus[is.na(AllAct$DrinkStatus)] <- "Missing alcohol"

# systolic blood pressure calculation 
AllAct$SYS <- AllFlags$SYS <-  round(apply(AllAct[,c("BPXSY1","BPXSY2","BPXSY3", "BPXSY4")],
                                           1,mean, na.rm= TRUE))

## Re-order columns so that activity and wear/non-wear flags are the last 1440 columns of our two
## data matrices. This is a personal preference and is absolutely not necessary.
act_cols <- which(colnames(AllAct) %in% paste0("MIN",1:1440))
oth_cols <- which(!colnames(AllAct) %in% paste0("MIN",1:1440))
AllAct   <- AllAct[,c(oth_cols,act_cols)]
AllFlags <- AllFlags[,c(oth_cols,act_cols)]
rm(list=c("act_cols","oth_cols"))

# And remove these:
rm(PAXINTEN_D, Covariate_D)

## make dataframe with one row per individual to create table 1.
## Remove columns associated with activity to avoid any confusion.
table_dat <- AllAct[!duplicated(AllAct$SEQN), 
                    -which(colnames(AllAct) %in% 
                             c(paste0("MIN",1:1440),"SBout","ABout","SATP","ASTP"))]

## subset based on our age inclusion/exclusion criteria
## note that individuals age 85 and over are coded as NA
#number of individuals excluded due to subset selection
# table_dat <- subset(table_dat, !(!(Age %in% seq(35,85, by= .5)) | is.na(Age)))
table_dat <- subset(table_dat, !(Age < 35 | is.na(Age)))

## get the SEQN (id variable) associated with individuals with fewer than 3 days accelerometer wear time
## with at least 10 hours OR had their data quality/device calibration flagged by NHANES
keep_inx       <- exclude_accel(AllAct, AllFlags)
Act_Analysis   <- AllAct[keep_inx,]
Flags_Analysis <- AllFlags[keep_inx,]
nms_rm         <- unique(c(Act_Analysis$SEQN[-which(Act_Analysis$SEQN %in% names(table(Act_Analysis$SEQN))[table(Act_Analysis$SEQN)>=3])],
                           setdiff(AllAct$SEQN, Act_Analysis$SEQN))
)
rm(list=c("keep_inx"))


## Additional inclusion/exclusion criteria.
## Aside from mortality or accelerometer weartime, the only missingness is in
## BMI, Education, SYS, total cholesterol, LBXTC and HDL cholesterol, LBDHDD.
criteria_vec <- c("(is.na(table_dat$BMI_cat))",         # missing BMI
                  "(is.na(table_dat$EducationAdult))",  # missing education
                  "(table_dat$SEQN %in% nms_rm)",       # too few "good" days of accelerometery data
                  "(is.na(table_dat$SYS) | (is.na(table_dat$LBXTC)) | (is.na(table_dat$LBDHDD)) )" #missing lab measures
)

## add in column indicating exclusion:
##   Exclude = 1 indicates an individual does not meet our inclusion criteria
##   Exclude = 0 indicates an individual does meet our inclusion criteria
eval(parse(text=paste0("table_dat$Exclude <- as.integer(", paste0(criteria_vec,collapse="|"), ")")))

## Create our dataset for analysis with one row per subject
## containing only those subjects who meet the inclusion criteria.
data_analysis  <- subset(table_dat, Exclude == 0)

# Remove the extraneous variables from data analysis 
extra_var = c("PAXCAL","PAXSTAT", "WEEKDAY", "SDDSRVYR", "WTMEC2YR", 
              "SDMVPSU", "SDMVSTRA", "WTINT2YR", "RIDAGEMN", "RIDAGEEX", "RIDAGEYR",
              "Exclude", "BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4")
data_analysis = data_analysis %>% dplyr::select(-extra_var)

# Rename some variables:
data_analysis = data_analysis %>% rename(
  "Total Cholesterol" = "LBXTC",
  "HDL Cholesterol" = "LBDHDD",
  "Systolic Blood Pressure" = "SYS",
)

# Combine some categories of these variables and remove those with mobility problems:
data_analysis = data_analysis %>% mutate(
  Race = case_when(
    Race == "White" ~ "White",
    Race == "Black" ~ "Black",
    Race == "Other" ~ "Other",
    (Race == "Mexican American") | (Race == "Other Hispanic") ~ "Hispanic"
  ) %>% fct_relevel("White", "Black", "Hispanic", "Other"),
  
  Education = case_when(
    EducationAdult == "Less than high school" ~ "< HS",
    EducationAdult == "High school" ~ "= HS",
    EducationAdult == "More than high school" ~ "> HS"
  ) %>% fct_relevel("< HS", "= HS", "> HS"),
  
  MobilityProblem = case_when(
    MobilityProblem == "No Difficulty" ~ "No",
    MobilityProblem == "Any Difficulty" ~ "Yes"
  ) %>% fct_relevel("No", "Yes"),
) %>% filter(
  MobilityProblem == "No"
)

## Get activity/flag data for only those included participants AND days.
## Since we've already removed the "bad" days from Act_Analysis and Act_Flags,
## we need only subset based on subject ID now
Act_Analysis   <- subset(Act_Analysis, SEQN %in% data_analysis$SEQN)
Flags_Analysis <- subset(Flags_Analysis, SEQN %in% data_analysis$SEQN)

## Sort the SEQN in the ascending order in both Act_Analysis and Flags_Analysis
Act_Analysis <- Act_Analysis[sort.int(Act_Analysis$SEQN, index.return = T)$ix,]
Flags_Analysis <- Flags_Analysis[sort.int(Flags_Analysis$SEQN, index.return = T)$ix,]

rm(list=c("AllAct","AllFlags","nms_rm"))
#----------------------------------------------------------------------------
# Construct the subject-specific data
#----------------------------------------------------------------------------
# Number of individuals:
n = length(data_analysis$SEQN)

# Minutes in the day:
times = 1:1440

# Aggregate into intervals:
all_times = rep(seq(k_min, length(times), by = k_min), 
                each = k_min)

# Unique observation points:
tau = sort(unique(all_times)); 

# Compute total activity count by time-of-day:
Act = matrix(NA, nrow = n, ncol = length(tau));

# And number of curve observations per weekday/weekend 
n_wkday = n_wkend = numeric(n)

# Loop over individuals:
for(i in 1:n){
  # Indices of the ith individual:
  inds = which(data_analysis$SEQN[i] == Act_Analysis$SEQN)
  
  # Total activity for individual i in each minute (possibly across multiple days):
  Act_i = colSums(Act_Analysis[inds,paste0("MIN",times)], na.rm=TRUE)
  
  # Aggregate activitiy in k_min-minute intervals:
  Act[i,] = sapply(tau, function(tau_j) sum(Act_i[which(all_times == tau_j)]))
  
  # Count number of weekdays and weekends for individual i:
  # (1 = Sun, 2 = Mon, 3 = Tues, 4 = Weds, 5 = Thurs, 6 = Fri, 7 = Sat)
  wkday_i = Act_Analysis[inds,]$WEEKDAY
  n_wkday[i] = sum((wkday_i>=2) & (wkday_i <= 6))
  n_wkend[i] = sum((wkday_i==1) | (wkday_i == 7))
}

# Add weekday/weekend variables to dataset:
data_analysis = data_analysis %>% 
  mutate(n_wkday = n_wkday,
         n_wkend = n_wkend)

# Subset to complete cases:
sub_subj = complete.cases(Act) & 
  complete.cases(data_analysis)

# Predictor data:
data_analysis = data_analysis[sub_subj,] %>% dplyr::select(- c("SEQN"))

levels(data_analysis$Education) <- list("lessHS" = "< HS", "equalHS" = "= HS", "greaterHS" = "> HS")
names(data_analysis) <- gsub(" ", "_", names(data_analysis))


Yminute <- data.frame(t(Act_Analysis[,grep("MIN",names(Act_Analysis))]))
minutegroup <- 10 #aggregate by x minute intervals
grp <- (1:nrow(Yminute) - 1)%/%minutegroup

Y <- apply(Yminute, 2, function(w) aggregate(w, list(grp), function(z) sum(z, na.rm=FALSE))$x)

# remove subject with NA observations (just 1?)
Y <- Y[,rep(sub_subj, table(Act_Analysis$SEQN))]
nomiss <- apply(Y, 2, function(x) sum(is.na(x)))==0
Y <- Y[,nomiss]

Y <- sqrt(Y)

# start at 4am
Y <- Y[c((nrow(Y)/6 + 1):nrow(Y)),]

W <- Act_Analysis$SEQN[rep(sub_subj, table(Act_Analysis$SEQN))]
W <- W[nomiss]


data_analysis2 <- rowrep(data_analysis,  table(Act_Analysis$SEQN[rep(sub_subj, table(Act_Analysis$SEQN))]))[nomiss,]
data_analysis2$Weekday <- Act_Analysis$WEEKDAY[rep(sub_subj, table(Act_Analysis$SEQN))][nomiss]
data_analysis2$Weekend <- data_analysis2$Weekday %in% c(6,7)


# Construct the matrix of predictors:
X = model.matrix( ~ 
                    BMI + 
                    Race + Gender + Age +
                    Education +
                    DrinksPerWeek +
                    SmokeCigs +
                    Diabetes + CHF + CHD + Cancer + Stroke +
                    HDL_Cholesterol +
                    Total_Cholesterol +
                    Systolic_Blood_Pressure + Weekend, data = data_analysis2)
# 
# X[,-1] <- scale(X[,-1])

