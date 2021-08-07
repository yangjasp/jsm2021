######
####Example Multi-wave Survey under Two-Phase Design in R
######

######
### Load Packages - Only uses optimall, survey, and survival  ------------
#######

# install.packages("optimall")
# install.packages("survey")
# install.packages("survival")

library(optimall)
library(survey)
library(survival)

######
### Load Data from survival package, pre-process to match hypothetical example -------
#####

nwt <- subset(survival::nwtco, select =  - c(in.subcohort, study))
nwt$stage <- ifelse(nwt$stage %in% c(1,2), 0, 1) # dichotomize stage for simplicity
names(nwt)[names(nwt) == "seqno"] <- "id"
names(nwt)[names(nwt) == "instit"] <- "local"
names(nwt)[names(nwt) == "histol"] <- "central"
names(nwt)[names(nwt) == "edrel"] <- "trelapse"
names(nwt)[names(nwt) == "rel"] <- "relapse"

######
### Set Up Phase 1 Data
######
phase1 <- subset(nwt, select = -central) # Suppose "central" is initially unavailable

# Stratify: 

phase1$strata <- phase1$stage # Initialize strata column

phase1 <- split_strata(data = phase1, strata = "strata", split = NULL,
                       split_var = "age", split_at = 0.5,  
                       type = "global quantile") # Divide all strata at global median

table(phase1$new_strata)

# Stratify further by splitting the largest stratum (not shown in slides):
phase1 <- split_strata(data = phase1, strata = "new_strata", split = "0.age_[0,37]",
                       split_var = "age", split_at = 0.5, trunc = 3,  
                       type = "local quantile") # Split specific stratum at its local median.

# Rename strata that were split twice
phase1$new_strata[phase1$new_strata == "0.age_[0,37].age_(18,37]"] <- "0.age_(18,37]"
phase1$new_strata[phase1$new_strata == "0.age_[0,37].age_[0,18]"] <- "0.age_[0,18]"

table(phase1$new_strata) # Final strata

############
#### Build Naive Model for Phase 2, 
#### Wave 1 allocation (using local histology instead of central lab): -----------------
#############

fit1 <- survival::coxph(Surv(trelapse, relapse) ~ 
                          local + age + 
                          stage, 
                        data = phase1)

# Get influence functions and add them as column to Phase 1 data.
phase1$infl<-resid(fit1,type="dfbeta",
                   weighted=FALSE)[,1]

########
#### Initialize the multi-wave object and set up metadata
########

# Initialize. Give it a title and specify phase 2 will have 2 waves.
NwtSurvey <- new_multiwave(
  phases = 2,  waves = c(1, 2),
  phase1 = phase1, 
  metadata = list(title = "NWT Survey"))

#### Metadata

# Phase 1 Metadata
get_data(NwtSurvey, phase =  1, slot = "metadata") <- 
  list(description = "4028 observations of 10 variables. 
  Includes influence functions for naive model 
       with error-prone local histology.")

# View diagram
multiwave_diagram(NwtSurvey)

# Phase 2 metadata. Specify strata and id for any functions we later use on phase 2
# in apply_multiwave()
get_data(NwtSurvey, phase =  2, slot = "metadata") <- 
  list(description = "Two waves of 300 central histology samples each.",
       allocation =  "Influence functions for",
       strata = "new_strata",
       id = "id")

###########
###### Phase 2, Wave 1 (allocate based on naive model): ---------------
###########

#### General Steps for Any Wave:

# 1. Specify metadata.
# 2. Allocate with optimum_allocation() or allocate_wave(), place result in "design" slot
#    of the wave.
# 3. Select samples according to the design with apply_multiwave(fun = "sample_strata") for SRS
#    within strata according to design.
# 4. Gather Data.
# 5. Merge collected data into data from the previous wave (or phase 1, if current wave is wave 1) 
#    with merge_samples(). The resulting merged data will be in the "data" slot of the current wave
#    in the Multiwave object.
# 6. Re-estimate influence functions with most recent data if necessary.

# So, for Wave 1:

# 1. Specify metadata (not shown in slides)
get_data(NwtSurvey, phase =  2, wave =  1, slot = "metadata") <- 
  list(description = "Wave 1 - 300 Samples")

# 2a. Use optimum_allocation to allocate the wave. In our case, we 
#    allocate based on influence functions from naive model we fit earlier.

#    Note: we could have instead used 'apply_multiwave(fun = "optimum_allocation") here.
#    If we had, the result would automatically be placed in the "design slot

design <- optimum_allocation(data = phase1, strata = "new_strata",
                             y = "infl", nsample = 300,
                             method = "WrightII") # Recall that "infl" is our variable of interest.
design

# 2b. Place the design in the design slot of phase 2, wave 1 of the multiwave object.
get_data(NwtSurvey, phase = 2, wave = 1, slot = "design") <- design 

# 3. Select samples based on this design.

set.seed(500)
NwtSurvey <- apply_multiwave(NwtSurvey, phase = 2, wave = 1, 
                             fun = "sample_strata", 
                             design_strata = "strata",
                             n_allocated = "stratum_size")

get_data(NwtSurvey, phase = 2, wave = 1, slot = "samples") # View result. 300 unique ids.

# 4. Gather data. In real-life, this is actual survey data collection. 
#    But for our hypothetical example, we just have to pull 'central' out of the 
#    nwt data set for the ids we selected for wave 1 sampling.

get_data(NwtSurvey, phase = 2, wave = 1, slot = "sampled_data") <- subset(nwt, select = c(id, central))[as.character(nwt$id) %in% get_data(NwtSurvey, phase = 2, wave = 1, slot = "samples"),]

# 5. Merge collected data into full data with "merge_samples":

# Merge. Indicator for being sampled in phase 2 is in new "sampled_phase2" column
NwtSurvey <- apply_multiwave(NwtSurvey, phase = 2, wave = 1,
                             fun = "merge_samples", 
                             sampled_ind = "sampled_phase2") 

table(is.na(get_data(NwtSurvey, phase = 2, wave = 1, slot  =  "data")$central)) # 300 samples.

# 6. Re-estimate influence function using a weighted Cox model with wave 1 data.
#    Since we collected 300 samples of central, we can use this, out actual, 
#    variable of interest, in the model instead of the error-prone "local" that we
#    used before.

# Use the survey package to run a Cox model with our survey design.
svydesign <- twophase(id = list(~1,~1), strata = list(NULL, ~ new_strata),
                      subset = ~ sampled_phase2,
                      data = get_data(NwtSurvey, phase = 2, wave = 1, slot  =  "data"),
                      method="simple")

fit2 <- svycoxph(Surv(trelapse, relapse) ~ central + age + stage, design = svydesign)	
summary(fit2)

# Add the updated influence functions to the data of Wave 1
infl_wave1 <- resid(fit2,type="dfbeta",weighted=FALSE)[,1]
infl_wave1_df <- data.frame(id = names(infl_wave1),
                            infl_wave1 = infl_wave1)


get_data(NwtSurvey, phase = 2, wave = 1, slot = "data") <- 
  merge(get_data(NwtSurvey, phase = 2, wave = 1, slot = "data"), infl_wave1_df, 
        by = "id", all.x = TRUE)

###########
#### Phase 2, Wave 2
##########

# 1. Specify Metadata
get_data(NwtSurvey, phase =  2, wave =  2, slot = "metadata") <- 
  list(description = "Wave 2 - 300 Samples")

# 2. Allocate wave with apply_multiwave(fun = "allocate_wave") based on wave 1 influence functions

# Allocate wave takes into account that strata have already been sampled in previous waves.
NwtSurvey <- apply_multiwave(NwtSurvey, phase =  2, wave =  2,
                             strata = "new_strata",
                             fun =  "allocate_wave",
                             y = "infl_wave1", nsample = 300,
                             already_sampled = "sampled_phase2")

# Result automatically placed in design slot because we used allocate_wave
get_data(NwtSurvey,  phase = 2, wave = 2, slot = "design")

# 3. Randomly select samples based on allocation

set.seed(200)
NwtSurvey <- apply_multiwave(NwtSurvey, phase = 2, wave = 2, 
                             fun = "sample_strata", 
                             design_strata = "strata",
                             already_sampled = "sampled_phase2")

get_data(NwtSurvey, phase = 2, wave = 2, slot = "samples")

# 4. "Collect" data

get_data(NwtSurvey, phase = 2, wave = 2, slot = "sampled_data") <- subset(nwt, select = c(id, central))[as.character(nwt$id) %in% get_data(NwtSurvey, phase = 2, wave = 2, slot = "samples"),]

# 5. Merge collected data back into full data with merge_samples().

NwtSurvey <- apply_multiwave(NwtSurvey, phase = 2, wave = 2,
                             fun = "merge_samples", 
                             sampled_ind = "sampled_phase2")

table(is.na(get_data(NwtSurvey, phase = 2, wave = 2, slot  =  "data")$central)) # 600 total samples

multiwave_diagram(NwtSurvey)

######
#### Make model with our 600 samples
######
svydesign <- twophase(id = list(~1,~1), strata = list(NULL, ~ new_strata),
                      subset = ~ sampled_phase2,
                      data = get_data(NwtSurvey, phase = 2, wave = 2, slot  =  "data"),
                      method="simple")

fit_final <- svycoxph(Surv(trelapse, relapse) ~ central + age + stage,
                      design = svydesign)	
summary(fit_final)

# Compare to SRS
test <- get_data(NwtSurvey, phase = 2, wave = 2, slot = "data")

set.seed(350)
ids <- sample(nwt$id, size = 600, replace = FALSE)
temp <- nwt %>%
  dplyr::filter( id %in% ids) %>%
  dplyr::mutate(test_ind = 1) %>%
  dplyr::select(id, central_test = central, test_ind)

test <- merge(test, temp, by = "id", all.x = TRUE)
test$test_ind <- ifelse(is.na(test$test_ind), 0, 1)

svydesign_test <- twophase(id = list(~1,~1), strata = list(NULL, ~ new_strata),
                      subset = ~ test_ind,
                      data = test,
                      method="simple")

fit3 <- svycoxph(Surv(trelapse, relapse) ~ central_test + age + stage, 
                 design = svydesign_test)	
summary(fit3)
