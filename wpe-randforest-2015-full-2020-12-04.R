# 2015 Predicted Woody Cover and Biophysical Variable Random Forest
# Austin Rutherford
# arutherford@email.arizona.edu
# 2020-12-09

# Load packages
library(tidyverse)
library(tidymodels)
library(tidypredict)
library(ranger)
library(vip)
library(doParallel)

# Read in data
pwc_data <- read_csv('woody_landsat_2015/2015_allpts_NDVI_samples_full.csv')

# Check out the data
summary(pwc_data)
glimpse(pwc_data)

# Replace -9999 (missing data) value with NA, 
# Set negative 2015 NDVI Landsat predicted woody cover (pwc) values to 0,
# Divide pwc by 100 to convert percent to decimal,
pwc_clean <- pwc_data %>%
  mutate(clay = na_if(clay, -9999),
         depclay = na_if(depclay, -9999),
         fallmean = na_if(fallmean, -9999),
         fallsum = na_if(fallsum, -9999),
         ppt = na_if(ppt, -9999),
         slope = na_if(slope, -9999),
         springmean = na_if(springmean, -9999),
         springsum = na_if(springsum, -9999),
         summermean = na_if(summermean, -9999),
         summersum = na_if(summersum, -9999),
         tmax = na_if(tmax, -9999),
         washdist = na_if(washdist, -9999),
         wintermean = na_if(wintermean, -9999),
         wintersum = na_if(wintersum, -9999),
         summertmax = na_if(summertmax, -9999),
         summertmin = na_if(summertmin, -9999),
         springtmax = na_if(springtmax, -9999),
         springtmin = na_if(springtmin, -9999),
         wintertmax = na_if(wintertmax, -9999),
         wintertmin = na_if(wintertmin, -9999),
         falltmax = na_if(falltmax, -9999),
         falltmin = na_if(falltmin, -9999),
         pwc = replace(pwc, which(pwc<0),NA),
         pwc = pwc/100,
         pwcclass = as.factor(pwcclass),
         ecosite = as.factor(ecosite))

# Look at data again
glimpse(pwc_clean)

# Drop NAs and create data subset (half of the data)
pwc_samps <- pwc_clean %>% drop_na()

# Set up training and testing split
set.seed(4242)

# Train and test split
data_split <- initial_split(pwc_samps, prop = 0.80, strata = 'pwcclass')
pwc_train <- training(data_split)
pwc_test <- testing(data_split)

# Preprocess data
pwc_rec <- recipe(pwc ~ ., data = pwc_train) %>%
  step_rm(pwcclass, point, easting, northing) %>% # remove from model
  update_role(point, new_role = "id variable") %>% # make point an id for back tracking
  step_dummy(ecosite) %>% # use dummy variables for factors in model
  step_center(all_predictors()) %>% #mean center data, m = 0
  step_scale(all_predictors()) %>% #scale so sd = 1, finish normalizing
  step_nzv(all_predictors()) #remove zero variance (no info) cols

# Apply recipe preprocessing to training data
pwc_prepped <- prep(pwc_rec, training = pwc_train) # preps data, applies recipe

# Run (bake) prepped preprocessng to training data to see the number of final dummy variables
pwc_train_bake <- bake(pwc_prepped, new_data = pwc_train)

# Setup our model (using Random Forest, ranger package/engine)
rf_mod <- rand_forest(mtry = tune(), min_n = tune(),
                      trees = 1000) %>%
  set_mode("regression") %>%
  set_engine("ranger", importance = "impurity")

# Build workflow to pair model and cross validation and tuning with data preprocessing
pwc_wflow <- workflow() %>% 
  add_model(rf_mod) %>% 
  add_recipe(pwc_rec)

# Set up cross validation
folds <- vfold_cv(pwc_train, v = 5, repeats = 5) #5-fold (v) cross validation, only do it once (no repeats)
folds

# Set up the initial tuning grid for finding best mtry and min_n (hyperparameters)
# Match with number of predictor variables (31) following preprocessing/dummy variable creation from ecosites divided by three
rf_param <-
  pwc_wflow %>%
  parameters() %>%
  update(mtry = mtry(range = c(1L, 10L)))

rf_tune_grid <- grid_regular(rf_param, levels = 10)

rf_tune_grid

# Training RF model/tuning
# First, initiate parallel processing
num_cores <- detectCores() -1
cl <- makeCluster(num_cores, outfile = "", type = "FORK")
clusterEvalQ(cl, library(tidymodels))
registerDoParallel(cl)

# Run random forest models to tune hyperparameters
rf_fit_results <- tune_grid(pwc_wflow, 
                            resamples = folds,
                            grid = rf_tune_grid,
                            metrics = metric_set(rmse, mae, rsq))

# Save an RDS file for future use if needed, save that memory
saveRDS(rf_fit_results, "./randforest_rds/first_model_RF_full.rds")

# Graphs rmse for all min_n and mtry
rf_param_plot <- rf_fit_results %>%
  collect_metrics() %>%
  dplyr::filter(.metric == "rmse") %>%
  dplyr::select(mean, min_n, mtry) %>%
  tidyr::pivot_longer(min_n:mtry,
                      values_to = "value",
                      names_to = "parameter") %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "RMSE")

ggsave('Graphs/rf_param_plot.png', plot = rf_param_plot)

# Look at best random forest models based on the Mean Absolute Error (and/or RMSE)
show_best(rf_fit_results, metric = "mae")
show_best(rf_fit_results, metric = "rmse")

# Update tuning grid based on best range of min_n and mtry
rf_grid_update <- grid_regular(
  mtry(range = c(8, 10)),
  min_n(range = c(2, 6)),
  levels = 5)

# Tune random forest with targets min_n and mtry ranges
rf_fit_update <- tune_grid(pwc_wflow, 
                           resamples = folds,
                           grid = rf_grid_update,
                           metrics = metric_set(rmse, mae, rsq))

# Graphs rmse for targeted mtry and min_n
rf_param_plot_tune <- rf_fit_update  %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, 100*mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "RMSE")

ggsave('Graphs/rf_param_plot_tune.png', plot = rf_param_plot_tune)

# Pick the best set of hyperparameters
show_best(rf_fit_update, metric = "mae")

# Pull out parameters of the best model based on RMSE (for prediction), mtry = 8, min_n = 2
pwc_rf_best <- rf_fit_update %>% 
  select_best(metric = "rmse")

# Finalize the model with the parameters of the best selected model
final_rf <- finalize_model(rf_mod, pwc_rf_best)

# Rerun final model and output the variable importance based on the GINI index, point is ID only
pwc_vip <- final_rf %>%
  fit(pwc ~ .,
      data = pwc_train_bake) %>%
  vip(geom = "point", num_features = 31) # variable importance graph

ggsave('Graphs/pwc_vip_31.png', plot = pwc_vip)

# Finalize Workflow
final_wf <- workflow() %>%
  add_recipe(pwc_rec) %>%
  add_model(final_rf)

# Function for unregistrering cluster for parallel processing, problems with 'ranger'?  
#(https://stackoverflow.com/questions/25097729/un-register-a-doparallel-cluster)
# unregister <- function() {
#   env <- foreach:::.foreachGlobals
#   rm(list=ls(name=env), pos=env)
# }

# Apply last workflow to training and testing data
final_res <- final_wf %>%
  last_fit(data_split)

# Look at final model performance metrics (mae, rmse, rsq)
final_res %>%
  collect_metrics()

# Save an RDS file for future use if needed, save that memory
saveRDS(final_res, "./randforest_rds/final_model_RF_full.rds")

# Compare predicted pwc values to test pwc
final_res %>% 
  collect_predictions(summarize = FALSE) %>% arrange(desc(pwc))

# Graph of predicted values to testing pwc
rf_plot <- final_res %>%
  collect_predictions() %>%
  ggplot(aes(pwc, .pred)) +
  geom_point(size = 0.5, alpha = 0.5) +
  labs(y = 'Random Forest Predicted Woody Cover (%)', x = 'Landsat NDVI Predicted Woody Cover (%)') +
  ylim(c(0, 0.60))+
  xlim(c(0, 0.80))+
  scale_color_manual(values = c("gray80", "darkred"))+
  geom_abline(intercept = 0, slope = 1, color = "blue", size = 1)+
  stat_smooth(method = "lm", formula = y ~ x, color = "red")+
  ggpmisc::stat_poly_eq(formula = y ~ x, 
                        aes(label =  paste(stat(eq.label),
                                           stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
                        parse = TRUE)+
  theme_bw()


# Save comparison plot of full random forest model
ggsave('Graphs/RF_Model_Full.png', plot = rf_plot)
