#Load libraries
library(tidyverse)

###########################################

# Define folder path
path = "~/Research/missing_surrogate/output/vary_prop/mar_1"
# Get a list of all CSV files in the folder
file_list <- list.files(path = path, pattern = "*.csv", full.names = TRUE)

# Load and combine all CSV files into one data frame
all_sim_data <- lapply(file_list, function(file) {
  # Read each file
  data <- read.csv(file)
  # Extract the file name
  data$filename <- basename(file)
  return(data)
}) %>%
  bind_rows()

# Set sim setting columns as factor
all_sim_data  <- all_sim_data  %>%
  mutate(across(c(prop_obs_s0, prop_obs_s1), as.factor))

# Find all columns relevant to parameter
models <- c("gold","complete", "ipw.ipw", "ipw.complete", "smle.ipw", "smle.complete")
parameter_model <- as.vector(outer("R.", models, paste0))

# Reshape the data from wide to long format
long_all_sim_data <- all_sim_data %>%
  select(prop_obs_s0, prop_obs_s1, starts_with("R")) %>%
  pivot_longer(cols = all_of(parameter_model),  # Adjust for relevant columns
               names_to = "Model_parameter",
               values_to = "Estimation")

# Make name for each sim settings
long_all_sim_data <-  long_all_sim_data |>
  mutate(group = paste(prop_obs_s0, prop_obs_s1, sep="-")) 
# Filter for group of interest before running this code
# Draw boxplot of simulation results for each methods across settings
long_all_sim_data |>
  ggplot(aes(x = Model_parameter, y = Estimation)) + 
  geom_boxplot() +
  labs(title = "Estimated R for difference proportions of data observed (control - treatment)",
       x = "Models",
       y = "R estimation") +
  geom_hline(yintercept=0.5, linetype='dotted', col = 'red') +
  facet_wrap(~group, ncol=3) 

