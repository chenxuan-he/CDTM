# import the necessary packages
library(stringr)
library(NLP)
library(tm)
library(readxl)
library(Matrix)
library(Rcpp)
library(numDeriv)
library(Rcgmin)
source("functions.R")

# Step 1: read the data in
# Load cut abstract, see _stemBeforeWDM.R for details
data <- read.csv("data_cut.csv")

# Step 2: Pre-process the data
source("data_preprocess.R")
lower.thresh <- 200
covariates_name <- c("Journal", "Num.author", 
                     "USA", "UK", "PRC", "France", "Canada", "Germany")
time_slices <- "Time_slices"
# change the data into digital numbers, refer to function data_preprocess() for details
out <- data_preprocess(data,
                       covariates_name = covariates_name,
                       text_name = "Abstract", 
                       time_slices = time_slices,
                       split_by = "\\|", 
                       word.lower.thresh=lower.thresh)

# Step 3: Initialize all the settings
# Set the topic numbers
source("settings_init.R")
K <- 8
# obs_variance refers to pre-set variance $\varsigma^2$
obs_variance <- .005
# chain_variance refers to pre-set variance $\sigma^2$, 
# where $\tau_t \mid \tau_{t-1} \sim N(\tau_{t-1}, \sigma^2 I)
chain_variance <- .005
verbose <- 1
# covariates model
prevalence =~ Journal+Num.author+USA+UK+PRC+France+Canada+Germany
settings <- settings_init(out, K, prevalence = prevalence, verbose = verbose,
                          obs_variance = obs_variance, chain_variance = chain_variance)

# Step 4: Initialize all the parameters
source("cdtm_init.R")
model <- cdtm_init(out$documents, settings)


# Step 5: Training
source("cdtm_control.R")
# To source estep, cpp environment should be enabled
source("estep.R")
source("CDTMmu.R")
source("CDTMsigma.R")
source("CDTMbeta_var.R")
source("CDTMoptbeta.R")
source("CDTMconvergence.R")

# update mean/var fwd_mean/var in the global state
mean.global <- model$mean
var.global <- model$var
zeta.global <- model$zeta
obs.global_byT <- NULL
obs.global_byK <- NULL
model_cdtm <- cdtm_control(documents=out$documents,
                                 vocab=out$vocab, 
                                 time_slices = out$meta$Time_slices, 
                                 settings=settings, 
                                 model=model)
