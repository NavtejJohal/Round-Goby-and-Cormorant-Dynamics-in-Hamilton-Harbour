library(tidyverse)
library(fixr)
library(glmmTMB)
set.seed(134)

###############################
## Initial Data Orginization ##
###############################
gobyData <- read_csv("Data/CPUE_BaitedOnly.csv")
cormorantData <- read_csv("Data/CormorantHWRaw.csv") |> 
  fix_column_names()
WQData <- read_csv("Data/WQ_2004-2020.csv") #water quality data

summed_Goby_Data <- gobyData |> 
  group_by(Year, Site) |> 
  summarise(
    across((cpue_Total),
           list(
             mean = ~ mean(.x, na.rm = TRUE)
           ),
           .names = "{.col}_{.fn}"
    ),
    .groups = "drop")

combData <- summed_Goby_Data %>%
  left_join(cormorantData, by = "Year")

summed_WQData <- WQData |> 
  group_by(Year, Site) |> 
  summarise(
    across(c(Secchi_cm, Oxygen_mgmL, Temp_C, pH),
           list(
             mean = ~ mean(.x, na.rm = TRUE)
           ),
           .names = "{.col}_{.fn}"
    ),
    .groups = "drop")


combData <- combData %>%
  left_join(summed_WQData, by = c("Year", "Site"))
#################
#################
### 2002-2022 ###
#################
#################

# Z transforming
nests <- combData |>
  distinct(Year, Active_Nests) |>
  arrange(Year)

nest_stats <- nests |>
  summarise(
    nest_mean = mean(Active_Nests, na.rm = TRUE),
    nest_sd = sd(Active_Nests, na.rm = TRUE)
  )

combData <- combData |>
  mutate(
    Active_Nests_z = (Active_Nests - nest_stats$nest_mean) /
      nest_stats$nest_sd
  )


combData$cpue_z <- as.numeric(scale(combData$cpue_Total_mean))
combData$Secchi_z <- as.numeric(scale(combData$Secchi_cm_mean))
combData$Oxygen_z <- as.numeric(scale(combData$Oxygen_mgmL_mean))
combData$pH_z <- as.numeric(scale(combData$pH_mean))
combData$Temp_z <- as.numeric(scale(combData$Temp_C_mean))
combData$Year_z <- as.numeric(scale(combData$Year))
combData$Year <- as.factor(combData$Year) #making sure year is a factor
combData$Site <- as.factor(combData$Site)
##########################
## Parameter estimation ##
##########################
# Estimating site + year variation 
site_sd <- 0.5
year_sd <- 0.4

residualvar <- sd(combData$cpue_z)^2 - year_sd^2 - site_sd^2
residualsd <- sqrt(residualvar)
################
## Simulation ##
################
#z scaled simulation
sim1 <- simulate_new(~ 
                       Active_Nests_z + 
                       Year_z + 
                       (1|Year) + 
                       (1|Site),
                     newdata = combData,
                     newparams = list(
                       beta = c(0,-0.1,-0.4), #intercept,active nests,year
                       theta = c(log(year_sd), log(site_sd)),#1|year 0.5,1|Site, 0.5
                       betadisp = 2*log(residualsd)#residual 0.77
                     )
)

# Plots
combData <- combData |> 
  mutate(
    cpue_zsim1 = sim1[[1]]
  )

sim1_goby_year_plot <- combData |>
  
  ggplot(aes(x=Year_z)) +
  
  geom_point(
    aes(y=cpue_z),
    colour = "red",
    alpha = 0.5) +
  
  geom_point(
    aes(y=cpue_zsim1),
    colour = "blue",
    alpha = 0.5) +
  
  theme_classic() +
  
  geom_smooth(
    aes(y=cpue_z),
    colour = "red",
    method = lm,
    alpha = 0.5) +
  
  #facet_wrap(~ Site) +
  
  geom_smooth(
    aes(y=cpue_zsim1),
    colour = "blue",
    method = lm,
    alpha = 0.5) +
  
  labs(
    subtitle = "Red = observed; Blue = simulated",
    x = "z transformed year",
    y = "abundance (z transformed cpue)"
  )

print(sim1_goby_year_plot)

sim1_goby_DCCO_plot <- combData |>
  
  ggplot(aes(x = Active_Nests_z)) +
  
  geom_point(
    aes(y = cpue_z),
    colour = "red",
    alpha = 0.5
  ) +
  
  geom_point(
    aes(y = cpue_zsim1),
    colour = "blue",
    alpha = 0.5
  ) +
  
  geom_smooth(
    aes(y = cpue_z),
    method = "lm",
    se = FALSE,
    colour = "red"
  ) +
  
  geom_smooth(
    aes(y = cpue_zsim1),
    method = "lm",
    se = FALSE,
    colour = "blue"
  ) +
  
  theme_bw() +
  
  labs(
    x = "Standardized Active Cormorant Nests",
    y = "Standardized Goby Abundance",
    title = "Observed and simulated standardized goby abundance vs standardized cormorant nest counts",
    subtitle = "Red = observed; Blue = simulated"
  )


print(sim1_goby_DCCO_plot)

####################
## Power_Analysis ##
####################
library(broom.mixed)
library(purrr)
library(patchwork)
library(performance)
# Extract one observation per year for predator abundance
nests <- nests |>
  mutate(
    Active_Nests_z = (Active_Nests - nest_stats$nest_mean) /
      nest_stats$nest_sd
  )

# nest statistics
nest_mean <- mean(nests$Active_Nests, na.rm = TRUE)
nest_sd <- sd(nests$Active_Nests, na.rm = TRUE)

# Year statistics
#year_mean <- mean(nests$Year, na.rm = TRUE)
year_mean <- 4800
year_sd <- sd(nests$Year, na.rm = TRUE)

# make year a dbl
nests <- nests |>
  mutate(
    Year = as.numeric(Year)
  )

# Estimate trend in nest count
nest_trend <- lm(
  Active_Nests ~ Year,
  data = nests
)

# Estimate variation around nest count trend
stable_sd <- sd(
  residuals(nest_trend),
  na.rm = TRUE
)

# Generate one 40-year annual cormorant dataset
make_annual_data <- function(n_years) {
  
  observed_years <- nrow(nests)
  
  additional_years <- n_years - observed_years
  
  if (additional_years < 0) {
    stop("Requested years cannot be fewer than observed years.")
  }
  
  if (additional_years > 0) {
    
    future_data <- tibble(
      Year = max(nests$Year) +
        seq_len(additional_years),
      
      Active_Nests = round(
        pmax(
          0,
          rnorm(
            additional_years,
            mean = nest_mean,
            sd = stable_sd
          )
        )
      )
    )
    
    annual_extended <- bind_rows(
      nests |>
        select(Year, Active_Nests),
      
      future_data
    )
    
  } else {
    
    annual_extended <- nests|>
      select(Year, Active_Nests)
  }
  
  annual_extended|>
    mutate(
      Active_Nests_z =
        (Active_Nests - nest_mean) / nest_sd,
      
      Year_z =
        (Year - year_mean) / year_sd,
      
      Year_f = factor(Year)
    )
}

# Generate the 40-year dataset once
annual_40_full <- make_annual_data(40)
#plot it
ggplot(
  annual_40_full,
  aes(x = Year, y = Active_Nests)
) +
  geom_line() +
  geom_point() +
  geom_vline(
    xintercept = max(nests$Year),
    linetype = "dashed"
  ) +
  labs(
    title = "Observed and simulated predator abundance",
    x = "Year",
    y = "Active predator nests"
  ) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5000))

# Extract the six sites
sites <- levels(combData$Site)

# Function to create site-year datasets
make_design <- function(n_years) {
  
  annual_subset <- annual_40_full %>%
    slice_head(n = n_years)
  
  crossing(
    Year = annual_subset$Year,
    Site = sites
  ) %>%
    left_join(
      annual_subset,
      by = "Year"
    ) %>%
    mutate(
      Site = factor(Site, levels = sites),
      Year_f = factor(Year),
      Fish_Abundance_z = 0
    )
}

# Create all four sample sizes
designs <- list(
  "20 years" = make_design(20),
  "25 years" = make_design(25),
  "30 years" = make_design(30),
  "40 years" = make_design(40)
)

# Check sample sizes
map_int(designs, nrow)

# model and simulation formula
model_formula <- Fish_Abundance_z ~
  Active_Nests_z + Year_z +
  (1 | Year_f) + (1 | Site)

sim_formula <- ~
  Active_Nests_z + Year_z +
  (1 | Year_f) + (1 | Site)

# parameters
## active nest fixed effect
beta_scenarios <- tibble(
  beta_true = c(0, -0.10, -0.20, -0.30)
)
## Fixed effect for year
beta_year <- -0.4

# #Gaussian residual SD
residual_sd <- 0.70

# Random effect SDs
year_sd <- 0.40
site_sd <- 0.50

run_simulation <- function(
    design,
    beta_true,
    nsim = 10,
    seed = 123
) {
  
  # Define parameters
  params <- list(
    
    beta = c(
      0,
      beta_true,
      beta_year
    ),
    
    betadisp = 2 * log(residual_sd),
    
    theta = c(
      log(year_sd),
      log(site_sd)
    )
  )
  
  # Simulate datasets
  simulated_values <- glmmTMB::simulate_new(
    object = sim_formula,
    newdata = design,
    family = gaussian,
    newparams = params,
    nsim = nsim,
    seed = seed
  )
  
  # Fit the model to each simulated dataset
  results <- map_dfr(seq_len(nsim), function(i) {
    
    sim_data <- design
    
    sim_data$Fish_Abundance_z <-
      simulated_values[[i]]
    
    fit <- tryCatch(
      suppressWarnings(
        glmmTMB(
          model_formula,
          family = gaussian,
          data = sim_data
        )
      ),
      error = function(e) NULL
    )
    
    # Handle failed model fits
    if (is.null(fit)) {
      
      return(tibble(
        replicate = i,
        estimate = NA_real_,
        std.error = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        converged = FALSE
      ))
    }
    
    converged <- (
      fit$fit$convergence == 0 &&
        isTRUE(fit$sdr$pdHess)
    )
    
    # Extract fixed effects and Wald confidence intervals
    tidy_fit <- tryCatch(
      broom.mixed::tidy(
        fit,
        effects = "fixed",
        component = "cond",
        conf.int = TRUE,
        conf.level = 0.95,
        conf.method = "wald"
      ),
      error = function(e) NULL
    )
    
    if (is.null(tidy_fit)) {
      
      return(tibble(
        replicate = i,
        estimate = NA_real_,
        std.error = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        converged = FALSE
      ))
    }
    
    predator_result <- tidy_fit %>%
      filter(term == "Active_Nests_z")
    
    if (nrow(predator_result) == 0) {
      
      return(tibble(
        replicate = i,
        estimate = NA_real_,
        std.error = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        converged = FALSE
      ))
    }
    
    tibble(
      replicate = i,
      estimate = predator_result$estimate[1],
      std.error = predator_result$std.error[1],
      conf.low = predator_result$conf.low[1],
      conf.high = predator_result$conf.high[1],
      converged = converged
    )
  })
  
  # Keep successful model fits
  successful <- results %>%
    filter(
      converged,
      is.finite(estimate),
      is.finite(conf.low),
      is.finite(conf.high)
    )
  
  # Calculate performance measures
  summary_results <- successful %>%
    summarise(
      
      n_successful = n(),
      
      convergence_rate = n_successful / nsim,
      
      mean_estimate = mean(estimate),
      
      bias = mean(estimate - beta_true),
      
      RMSE = sqrt(
        mean((estimate - beta_true)^2)
      ),
      
      mean_CI_width = mean(
        conf.high - conf.low
      ),
      
      coverage = mean(
        conf.low <= beta_true &
          conf.high >= beta_true
      ),
      
      power = if (beta_true < 0) {
        mean(conf.high < 0)
      } else if (beta_true > 0) {
        mean(conf.low > 0)
      } else {
        NA_real_
      }
    )
  
  list(
    replicates = results,
    summary = summary_results
  )
}

all_replicates <- list()
all_summaries <- list()

counter <- 1

for (design_name in names(designs)) {
  
  current_design <- designs[[design_name]]
  
  for (j in seq_len(nrow(beta_scenarios))) {
    
    current_beta <- beta_scenarios$beta_true[j]
    
    # Run simulations
    output <- run_simulation(
      design = current_design,
      beta_true = current_beta,
      nsim = 1000, # number of simmulations
      seed = 1000 + counter
    )
    
    # Save individual simulation results
    all_replicates[[counter]] <-
      output$replicates %>%
      mutate(
        sample_size = design_name,
        n_years = n_distinct(current_design$Year),
        beta_true = current_beta
      )
    
    # Save summary statistics
    all_summaries[[counter]] <-
      output$summary %>%
      mutate(
        sample_size = design_name,
        n_years = n_distinct(current_design$Year),
        beta_true = current_beta
      )
    
    counter <- counter + 1
  }
}

# Combine results
replicate_results <- bind_rows(all_replicates)

power_summary <- bind_rows(all_summaries)

#check summary
power_summary %>%
  select(
    sample_size,
    beta_true,
    n_successful,
    convergence_rate,
    mean_estimate,
    bias,
    RMSE,
    mean_CI_width,
    coverage,
    power
  )

plot_data <- power_summary %>%
  mutate(
    beta_group = factor(
      beta_true,
      levels = c(0, -0.10, -0.20, -0.30),
      labels = c(
        "Beta = 0.00",
        "Beta = -0.10",
        "Beta = -0.20",
        "Beta = -0.30"
      )
    )
  )
p_beta <- ggplot(
  plot_data,
  aes(
    x = n_years,
    y = mean_estimate,
    colour = beta_group,
    group = beta_group
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_line(
    aes(y = beta_true),
    linetype = "dashed"
  ) +
  labs(
    title = "Estimated beta vs. sample size",
    x = "Number of years",
    y = "Mean estimated beta",
    colour = "True beta"
  ) +
  theme_classic()

p_beta

p_width <- ggplot(
  plot_data,
  aes(
    x = n_years,
    y = mean_CI_width,
    colour = beta_group,
    group = beta_group
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(
    title = "Confidence interval width vs. sample size",
    x = "Years",
    y = "Mean 95% CI width",
    colour = "True beta"
  ) +
  theme_classic()

p_width

p_coverage <- ggplot(
  plot_data,
  aes(
    x = n_years,
    y = coverage,
    colour = beta_group,
    group = beta_group
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    title = "Coverage vs. sample size",
    x = "Years",
    y = "Coverage probability",
    colour = "True beta"
  ) +
  theme_classic()

p_coverage

p_power <- ggplot(
  plot_data,
  aes(
    x = n_years,
    y = power,
    colour = beta_group,
    group = beta_group
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(
    yintercept = 0.80,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    limits = c(0, 1)
  ) +
  labs(
    title = "Power vs. sample size",
    x = "Years",
    y = "Statistical power",
    colour = "True beta"
  ) +
  theme_classic()

p_power

p_bias <- ggplot(
  power_summary,
  aes(
    x = n_years,
    y = bias,
    colour = factor(beta_true),
    group = beta_true
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Bias of cormorant effect estimates",
    x = "Years",
    y = "Bias",
    colour = "True Beta"
  ) +
  theme_classic()

p_bias

p_rmse <- ggplot(
  power_summary,
  aes(
    x = n_years,
    y = RMSE,
    colour = factor(beta_true),
    group = beta_true
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(
    title = "RMSE of cormorant effect estimates",
    x = "Years",
    y = "Root Mean Square Error (RMSE)",
    colour = "True Beta"
  ) +
  theme_classic()

p_rmse
# Save individual simulation results
write.csv(
  replicate_results,
  "replicate_results.csv",
  row.names = FALSE
)

# Save summary statistics
write.csv(
  power_summary,
  "power_summary.csv",
  row.names = FALSE
)

# Save R objects
saveRDS(
  list(
    replicate_results = replicate_results,
    power_summary = power_summary
  ),
  "power_analysis_results.rds"
)
#save plots
ggsave("beta_plot.png",p_beta)
ggsave("CI_width_plot.png",p_width)
ggsave("coverage_plot.png",p_coverage)
ggsave("power_plot.png",p_power)
ggsave("bias_plot.png",p_bias)
ggsave("rmse_plot.png",p_rmse)
