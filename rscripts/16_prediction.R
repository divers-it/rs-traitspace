rm(list=ls())

####
# Read in data ----
####

# load data set
full_df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# load imputed data set
full_df <-read.csv("outputs/imputed_with_phylo.csv", row.names=1, stringsAsFactors = TRUE)

# load untransformed data set
full_df <- readRDS(file = here::here("outputs/5_df_filt.rds"))

####
# Predicting flower sex ----
####

####
## Decision trees ----
#### 
# https://bradleyboehmke.github.io/HOML/DT.html

# remove flower sex traits
df <- subset(full_df, select=-c(SexualSystem,Mating))

# remove NA
df <- df[!is.na(df$FlowerSex),]

# make a decision tree
# to understand control options
# https://cran.r-project.org/web/packages/rpart/vignettes/longintro.pdf
flowersex_dt1 <- rpart::rpart(
  formula = FlowerSex ~ .,
  data    = df,
  control = rpart::rpart.control(xval = 100,
                          minsplit = 20, # default = 20
                          minbucket = 7, # default = minsplit / 3
                          maxcompete = 4, # default = 4
                          xva = 100,
                          maxsurrogate = 5, # default = 5
                          usesurrogate = 2#, # default = 2
                          # cp = # threshold complexity
                          )
)

# look at decision tree output
flowersex_dt1

# plot decision tree
rpart.plot::rpart.plot(flowersex_dt1)

# save plot with rules
pdf("figures/16_flowersex_dt1.pdf",width=10,height=10)
rpart.plot::rpart.plot(flowersex_dt1, type = 3, clip.right.labs = FALSE, branch = .3, under = TRUE)
dev.off()

# Cross-validation error
# lower x-axis is cost complexity value, upper x-axis is number of terminal nodes
# it’s common to use the smallest tree within 1 standard error (SE) of the 
# minimum CV error (represented by the dashed line).
rpart::plotcp(flowersex_dt1)

# rpart cross validation results
flowersex_dt1$cptable

# caret cross validation results
# Note: only works without missing data
if(sum(is.na(df))<1){
  
  flowersex_dt3 <- caret::train(
    FlowerSex ~ .,
    data = df,
    method = "rpart",
    trControl = caret::trainControl(method = "cv", number = 100),
    tuneLength = 60
  )
  
  # look at model results
  flowersex_dt3
  
  # plot cross-validation results
  # Lower α values (deeper trees) help to minimize errors.
  ggplot(flowersex_dt3)
  
  # To measure feature importance, the reduction in the loss function (e.g., SSE)
  # attributed to each variable at each split is tabulated.
  # 
  # When using caret, feature importance values are standardized so that the most important 
  # feature has a value of 100 and the remaining features are scored based on 
  # their relative reduction in the loss function.
  vip::vip(flowersex_dt3, num_features = 10, bar = FALSE)

}

####
## Random forests ----
####
# https://bradleyboehmke.github.io/HOML/random-forest.html

# load packages
library(dplyr)
library(h2o)

# number of features
n_features <- length(setdiff(names(df), "FlowerSex"))

# Stratified sampling with the rsample package
set.seed(123)
split <- rsample::initial_split(df, prop = 0.7, 
                                strata = "FlowerSex")
df_train  <- rsample::training(split)
df_test   <- rsample::testing(split)

# train a default random forest model
flowersex_rf1 <- ranger::ranger(
  FlowerSex ~ ., 
  data = df,
  mtry = floor(n_features / 3),
  respect.unordered.factors = "order",
  na.action = "na.learn",
  seed = 123
)

# look at output
flowersex_rf1

# get OOB accuracy
default_acc <- flowersex_rf1$prediction.error

# create hyperparameter grid
hyper_grid <- expand.grid(
  mtry = floor(n_features * c(.05, .15, .25, .333, .4)),
  min.node.size = c(1, 3, 5, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.5, .63, .8),                       
  acc = NA                                               
)

# execute full cartesian grid search
for(i in seq_len(nrow(hyper_grid))) {
  # fit model for ith hyperparameter combination
  fit <- ranger::ranger(
    formula         = FlowerSex ~ ., 
    data            = df_train, 
    num.trees       = n_features * 10,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$min.node.size[i],
    replace         = hyper_grid$replace[i],
    sample.fraction = hyper_grid$sample.fraction[i],
    verbose         = FALSE,
    seed            = 123,
    respect.unordered.factors = 'order',
  )
  # export OOB error 
  hyper_grid$acc[i] <- fit$prediction.error
}

# assess top 10 models
hyper_grid %>%
  arrange(acc) %>%
  mutate(perc_gain = (default_acc - acc) / default_acc * 100) %>%
  head(10)

# once you’ve identified the optimal parameter values from the grid search, 
# you will want to re-run your model with these hyperparameter values. 
# You can also crank up the number of trees, which will help create more stables
# values of variable importance.

# re-run model with impurity-based variable importance
#
# NOTE: unsure what this is, but many sources say it is biased compared to permutation-based
rf_impurity <- ranger::ranger(
  formula = FlowerSex ~ ., 
  data = df_train, 
  num.trees = 2000,
  mtry = 7,
  min.node.size = 1,
  sample.fraction = .63,
  replace = FALSE,
  importance = "impurity",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed  = 123
)

# re-run model with permutation-based variable importance
# 
# In the permutation-based approach, for each tree, the OOB sample is passed down
# the tree and the prediction accuracy is recorded. Then the values for each variable
# (one at a time) are randomly permuted and the accuracy is again computed. 
# The decrease in accuracy as a result of this randomly shuffling of feature 
# values is averaged over all the trees for each predictor. The variables with
# the largest average decrease in accuracy are considered most important.
rf_permutation <- ranger::ranger(
  formula = FlowerSex ~ ., 
  data = df_train, 
  num.trees = 2000,
  mtry = 7,
  min.node.size = 1,
  sample.fraction = .63,
  replace = FALSE,
  importance = "permutation",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed  = 123
)

# make individual plots
p1 <- vip::vip(rf_impurity, num_features = 15, bar = FALSE)
p2 <- vip::vip(rf_permutation, num_features = 15, bar = FALSE)

# Typically, you will not see the same variable importance order between the two options;
# however, you will often see similar variables at the top of the plots (and also the bottom).
gridExtra::grid.arrange(p1, p2, nrow = 1)

####
### Using h2o to fit rf model and tune hyperparameters ----
####

# initiate h2o session
h2o.no_progress()
h2o.init(max_mem_size = "2g")

# convert training data to h2o object
train_h2o <- as.h2o(df_train)

# set the response column
response <- "FlowerSex"

# set the predictor names
predictors <- setdiff(colnames(df_train), response)

h2o_rf1 <- h2o.randomForest(
  x = predictors, 
  y = response,
  training_frame = train_h2o, 
  ntrees = n_features * 10,
  seed = 123
)

# check output and accuracy
h2o_rf1

# hyperparameter grid
hyper_grid <- list(
  mtries = floor(n_features * c(.05, .15, .25, .333, .4)),
  min_rows = c(1, 3, 5, 10),
  max_depth = c(10, 20, 30),
  sample_rate = c(.55, .632, .70, .80)
)

# random grid search strategy
search_criteria <- list(
  strategy = "RandomDiscrete",
  stopping_metric = "mse",
  stopping_tolerance = 0.001,   # stop if improvement is < 0.1%
  stopping_rounds = 10,         # over the last 10 models
  max_runtime_secs = 60*5      # or stop search after 5 min.
)

# perform grid search 
random_grid <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_random_grid",
  x = predictors, 
  y = response, 
  training_frame = train_h2o,
  hyper_params = hyper_grid,
  ntrees = n_features * 10,
  seed = 123,
  stopping_metric = "RMSE",   
  stopping_rounds = 10,           # stop if last 10 trees added 
  stopping_tolerance = 0.005,     # don't improve RMSE by 0.5%
  search_criteria = search_criteria
)

# collect the results and sort by our model performance metric of choice
random_grid_perf <- h2o.getGrid(
  grid_id = "rf_random_grid", 
  sort_by = "accuracy", 
  decreasing = FALSE
)

# summarize results
random_grid_perf@summary_table

# shut down session
h2o.shutdown(prompt=FALSE)

####
# Predicting mating system ----
####

####
## Decision trees ----
#### 
# https://bradleyboehmke.github.io/HOML/DT.html

# remove flower sex traits
df <- subset(full_df, select=-c(SexualSystem,FlowerSex))

# remove NA
df <- df[!is.na(df$Mating),]

# make a decision tree
mating_dt1 <- rpart::rpart(
  formula = Mating ~ .,
  data    = df
)

# look at decision tree output
mating_dt1

# plot decision tree
rpart.plot::rpart.plot(mating_dt1)

# save plot with rules
pdf("figures/16_mating_dt1.pdf",width=10,height=10)
rpart.plot::rpart.plot(mating_dt1, type = 3, clip.right.labs = FALSE, branch = .3, under = TRUE)
dev.off()

# Cross-validation error
# lower x-axis is cost complexity value, upper x-axis is number of terminal nodes
# it’s common to use the smallest tree within 1 standard error (SE) of the 
# minimum CV error (represented by the dashed line).
rpart::plotcp(mating_dt1)

# rpart cross validation results
mating_dt1$cptable

# caret cross validation results
# Note: only works without missing data
if(sum(is.na(df))<1){
  
  mating_dt3 <- caret::train(
    Mating ~ .,
    data = df,
    method = "rpart",
    trControl = caret::trainControl(method = "cv", number = 100),
    tuneLength = 60
  )
  
  # look at model results
  mating_dt3
  
  # plot cross-validation results
  # Lower α values (deeper trees) help to minimize errors.
  ggplot(mating_dt3)
  
  # To measure feature importance, the reduction in the loss function (e.g., SSE)
  # attributed to each variable at each split is tabulated.
  # 
  # When using caret, feature importance values are standardized so that the most important 
  # feature has a value of 100 and the remaining features are scored based on 
  # their relative reduction in the loss function.
  vip::vip(mating_dt3, num_features = 10, bar = FALSE)
  
}

####
## Random forests ----
####
# https://bradleyboehmke.github.io/HOML/random-forest.html

# load packages
library(dplyr)
library(h2o)

# number of features
n_features <- length(setdiff(names(df), "Mating"))

# Stratified sampling with the rsample package
set.seed(123)
split <- rsample::initial_split(df, prop = 0.7, 
                                strata = "Mating")
df_train  <- rsample::training(split)
df_test   <- rsample::testing(split)

# train a default random forest model
Mating_rf1 <- ranger::ranger(
  Mating ~ ., 
  data = df,
  mtry = floor(n_features / 3),
  respect.unordered.factors = "order",
  seed = 123
)

# look at output
Mating_rf1

# get OOB accuracy
default_acc <- Mating_rf1$prediction.error

# create hyperparameter grid
hyper_grid <- expand.grid(
  mtry = floor(n_features * c(.05, .15, .25, .333, .4)),
  min.node.size = c(1, 3, 5, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.5, .63, .8),                       
  acc = NA                                               
)

# execute full cartesian grid search
for(i in seq_len(nrow(hyper_grid))) {
  # fit model for ith hyperparameter combination
  fit <- ranger::ranger(
    formula         = Mating ~ ., 
    data            = df_train, 
    num.trees       = n_features * 10,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$min.node.size[i],
    replace         = hyper_grid$replace[i],
    sample.fraction = hyper_grid$sample.fraction[i],
    verbose         = FALSE,
    seed            = 123,
    respect.unordered.factors = 'order',
  )
  # export OOB error 
  hyper_grid$acc[i] <- fit$prediction.error
}

# assess top 10 models
hyper_grid %>%
  arrange(acc) %>%
  mutate(perc_gain = (default_acc - acc) / default_acc * 100) %>%
  head(10)

# once you’ve identified the optimal parameter values from the grid search, 
# you will want to re-run your model with these hyperparameter values. 
# You can also crank up the number of trees, which will help create more stables
# values of variable importance.

# re-run model with impurity-based variable importance
#
# NOTE: unsure what this is, but many sources say it is biased compared to permutation-based
rf_impurity <- ranger::ranger(
  formula = Mating ~ ., 
  data = df_train, 
  num.trees = 2000,
  mtry = 7,
  min.node.size = 1,
  sample.fraction = .63,
  replace = FALSE,
  importance = "impurity",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed  = 123
)

# re-run model with permutation-based variable importance
# 
# In the permutation-based approach, for each tree, the OOB sample is passed down
# the tree and the prediction accuracy is recorded. Then the values for each variable
# (one at a time) are randomly permuted and the accuracy is again computed. 
# The decrease in accuracy as a result of this randomly shuffling of feature 
# values is averaged over all the trees for each predictor. The variables with
# the largest average decrease in accuracy are considered most important.
rf_permutation <- ranger::ranger(
  formula = Mating ~ ., 
  data = df_train, 
  num.trees = 2000,
  mtry = 7,
  min.node.size = 1,
  sample.fraction = .63,
  replace = FALSE,
  importance = "permutation",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed  = 123
)

# make individual plots
p1 <- vip::vip(rf_impurity, num_features = 15, bar = FALSE)
p2 <- vip::vip(rf_permutation, num_features = 15, bar = FALSE)

# Typically, you will not see the same variable importance order between the two options;
# however, you will often see similar variables at the top of the plots (and also the bottom).
gridExtra::grid.arrange(p1, p2, nrow = 1)

####
### Using h2o to fit rf model and tune hyperparameters ----
####

# initiate h2o session
h2o.no_progress()
h2o.init(max_mem_size = "2g")

# convert training data to h2o object
train_h2o <- as.h2o(df_train)

# set the response column
response <- "Mating"

# set the predictor names
predictors <- setdiff(colnames(df_train), response)

h2o_rf1 <- h2o.randomForest(
  x = predictors, 
  y = response,
  training_frame = train_h2o, 
  ntrees = n_features * 10,
  seed = 123
)

# check output and accuracy
h2o_rf1

# hyperparameter grid
hyper_grid <- list(
  mtries = floor(n_features * c(.05, .15, .25, .333, .4)),
  min_rows = c(1, 3, 5, 10),
  max_depth = c(10, 20, 30),
  sample_rate = c(.55, .632, .70, .80)
)

# random grid search strategy
search_criteria <- list(
  strategy = "RandomDiscrete",
  stopping_metric = "mse",
  stopping_tolerance = 0.001,   # stop if improvement is < 0.1%
  stopping_rounds = 10,         # over the last 10 models
  max_runtime_secs = 60*5      # or stop search after 5 min.
)

# perform grid search 
random_grid <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_random_grid",
  x = predictors, 
  y = response, 
  training_frame = train_h2o,
  hyper_params = hyper_grid,
  ntrees = n_features * 10,
  seed = 123,
  stopping_metric = "RMSE",   
  stopping_rounds = 10,           # stop if last 10 trees added 
  stopping_tolerance = 0.005,     # don't improve RMSE by 0.5%
  search_criteria = search_criteria
)

# collect the results and sort by our model performance metric of choice
random_grid_perf <- h2o.getGrid(
  grid_id = "rf_random_grid", 
  sort_by = "accuracy", 
  decreasing = FALSE
)

# summarize results
random_grid_perf@summary_table

# shut down session
h2o.shutdown(prompt=FALSE)


####
# Predicting reproductive system ----
####

# read in original sexual systems
ors <- read.csv("outputs/original_reproductive_systems.csv")

#check order
ors$species == rownames(full_df)

####
## Decision trees ----
#### 
# https://bradleyboehmke.github.io/HOML/DT.html

# remove flower sex traits
df <- subset(full_df, select=-c(SexualSystem,Mating,FlowerSex))

# add reproductive system
df$RS <- as.factor(ors$RS)

# remove unknown reproductive system
df <- df[df$RS!="unknown",]
df$RS <- droplevels(df$RS)

# remove NA
df <- df[!is.na(df$RS),]

# make a decision tree
# to understand control options
# https://cran.r-project.org/web/packages/rpart/vignettes/longintro.pdf
RS_dt1 <- rpart::rpart(
  formula = RS ~ .,
  data    = df,
  control = rpart::rpart.control(xval = 100,
                                 minsplit = 20, # default = 20
                                 minbucket = 7, # default = minsplit / 3
                                 maxcompete = 4, # default = 4
                                 xva = 100,
                                 maxsurrogate = 5, # default = 5
                                 usesurrogate = 2#, # default = 2
                                 # cp = # threshold complexity
  )
)

# look at decision tree output
RS_dt1

# plot decision tree
rpart.plot::rpart.plot(RS_dt1)

# save plot with rules
pdf("figures/16_RS_dt1.pdf",width=10,height=10)
rpart.plot::rpart.plot(RS_dt1, type = 3, clip.right.labs = FALSE, branch = .3, under = TRUE)
dev.off()

# Cross-validation error
# lower x-axis is cost complexity value, upper x-axis is number of terminal nodes
# it’s common to use the smallest tree within 1 standard error (SE) of the 
# minimum CV error (represented by the dashed line).
rpart::plotcp(RS_dt1)

# rpart cross validation results
RS_dt1$cptable

# caret cross validation results
# Note: only works without missing data
if(sum(is.na(df))<1){
  
  RS_dt3 <- caret::train(
    RS ~ .,
    data = df,
    method = "rpart",
    trControl = caret::trainControl(method = "cv", number = 100),
    tuneLength = 60
  )
  
  # look at model results
  RS_dt3
  
  # plot cross-validation results
  # Lower α values (deeper trees) help to minimize errors.
  ggplot(RS_dt3)
  
  # To measure feature importance, the reduction in the loss function (e.g., SSE)
  # attributed to each variable at each split is tabulated.
  # 
  # When using caret, feature importance values are standardized so that the most important 
  # feature has a value of 100 and the remaining features are scored based on 
  # their relative reduction in the loss function.
  vip::vip(RS_dt3, num_features = 10, bar = FALSE)
  
}

####
## Random forests ----
####
# https://bradleyboehmke.github.io/HOML/random-forest.html

# load packages
library(dplyr)
library(h2o)

# number of features
n_features <- length(setdiff(names(df), "RS"))

# Stratified sampling with the rsample package
set.seed(123)
split <- rsample::initial_split(df, prop = 0.7, 
                                strata = "RS")
df_train  <- rsample::training(split)
df_test   <- rsample::testing(split)

# train a default random forest model
RS_rf1 <- ranger::ranger(
  RS ~ ., 
  data = df_train,
  mtry = floor(n_features / 3),
  respect.unordered.factors = "order",
  seed = 123
)

# look at output
RS_rf1

# get OOB accuracy
default_acc <- RS_rf1$prediction.error

# create hyperparameter grid
hyper_grid <- expand.grid(
  mtry = floor(n_features * c(.05, .15, .25, .333, .4)),
  min.node.size = c(1, 3, 5, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.5, .63, .8),                       
  acc = NA                                               
)

# execute full cartesian grid search
for(i in seq_len(nrow(hyper_grid))) {
  # fit model for ith hyperparameter combination
  fit <- ranger::ranger(
    formula         = RS ~ ., 
    data            = df_train, 
    num.trees       = n_features * 10,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$min.node.size[i],
    replace         = hyper_grid$replace[i],
    sample.fraction = hyper_grid$sample.fraction[i],
    verbose         = FALSE,
    seed            = 123,
    respect.unordered.factors = 'order',
  )
  # export OOB error 
  hyper_grid$acc[i] <- fit$prediction.error
}

# assess top 10 models
hyper_grid %>%
  arrange(acc) %>%
  mutate(perc_gain = (default_acc - acc) / default_acc * 100) %>%
  head(10)

# once you’ve identified the optimal parameter values from the grid search, 
# you will want to re-run your model with these hyperparameter values. 
# You can also crank up the number of trees, which will help create more stables
# values of variable importance.

# re-run model with impurity-based variable importance
#
# NOTE: unsure what this is, but many sources say it is biased compared to permutation-based
rf_impurity <- ranger::ranger(
  formula = RS ~ ., 
  data = df_train, 
  num.trees = 2000,
  mtry = 7,
  min.node.size = 1,
  sample.fraction = .63,
  replace = FALSE,
  importance = "impurity",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed  = 123
)

# re-run model with permutation-based variable importance
# 
# In the permutation-based approach, for each tree, the OOB sample is passed down
# the tree and the prediction accuracy is recorded. Then the values for each variable
# (one at a time) are randomly permuted and the accuracy is again computed. 
# The decrease in accuracy as a result of this randomly shuffling of feature 
# values is averaged over all the trees for each predictor. The variables with
# the largest average decrease in accuracy are considered most important.
rf_permutation <- ranger::ranger(
  formula = RS ~ ., 
  data = df_train, 
  num.trees = 2000,
  mtry = 7,
  min.node.size = 1,
  sample.fraction = .63,
  replace = FALSE,
  importance = "permutation",
  respect.unordered.factors = "order",
  verbose = FALSE,
  seed  = 123
)

# make individual plots
p1 <- vip::vip(rf_impurity, num_features = 15, bar = FALSE)
p2 <- vip::vip(rf_permutation, num_features = 15, bar = FALSE)

# Typically, you will not see the same variable importance order between the two options;
# however, you will often see similar variables at the top of the plots (and also the bottom).
gridExtra::grid.arrange(p1, p2, nrow = 1)

####
### Using h2o to fit rf model and tune hyperparameters ----
####

# initiate h2o session
h2o.no_progress()
h2o.init(max_mem_size = "2g")

# convert training data to h2o object
train_h2o <- as.h2o(df_train)

# set the response column
response <- "RS"

# set the predictor names
predictors <- setdiff(colnames(df_train), response)

h2o_rf1 <- h2o.randomForest(
  x = predictors, 
  y = response,
  training_frame = train_h2o, 
  ntrees = n_features * 10,
  seed = 123
)

# check output and accuracy
h2o_rf1

# hyperparameter grid
hyper_grid <- list(
  mtries = floor(n_features * c(.05, .15, .25, .333, .4)),
  min_rows = c(1, 3, 5, 10),
  max_depth = c(10, 20, 30),
  sample_rate = c(.55, .632, .70, .80)
)

# random grid search strategy
search_criteria <- list(
  strategy = "RandomDiscrete",
  stopping_metric = "mse",
  stopping_tolerance = 0.001,   # stop if improvement is < 0.1%
  stopping_rounds = 10,         # over the last 10 models
  max_runtime_secs = 60*5      # or stop search after 5 min.
)

# perform grid search 
random_grid <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_random_grid",
  x = predictors, 
  y = response, 
  training_frame = train_h2o,
  hyper_params = hyper_grid,
  ntrees = n_features * 10,
  seed = 123,
  stopping_metric = "RMSE",   
  stopping_rounds = 10,           # stop if last 10 trees added 
  stopping_tolerance = 0.005,     # don't improve RMSE by 0.5%
  search_criteria = search_criteria
)

# collect the results and sort by our model performance metric of choice
random_grid_perf <- h2o.getGrid(
  grid_id = "rf_random_grid", 
  sort_by = "accuracy", 
  decreasing = FALSE
)

# summarize results
random_grid_perf@summary_table

# shut down session
h2o.shutdown(prompt=FALSE)