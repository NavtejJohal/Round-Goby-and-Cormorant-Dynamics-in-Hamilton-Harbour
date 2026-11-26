# Functions

get_goby_data <- function(goby_file){
  read_csv(goby_file)
}

get_cormorant_data <- function(cormorant_file){
  cormorantData <- read_csv(cormorant_file) |> 
    fix_column_names()
}

combine_data <- function(goby_data, cormorant_data) {
  summed_Goby_Data <- goby_data |> 
    group_by(Year, Site) |> 
    summarise(
      across((cpue_Total),
             list(
               mean = ~ mean(.x, na.rm = TRUE)
             ),
             .names = "{.col}_{.fn}"
      ),
      .groups = "drop")
  
  summed_Goby_Data |>
    left_join(cormorant_data, by = "Year")
}

get_sites <- function(comb_data){
  comb_data$Site <- factor(comb_data$Site)
  sites <- levels(comb_data$Site)
  return(sites)
}

fix_nest_data <- function(comb_data) {
  comb_data |>
    distinct(Year, Active_Nests) |>
    arrange(Year)
}

model_cormorant_trend <- function(nest_data){
# Estimate trend in nest count
  nest_trend <- lm(
    Active_Nests ~ Year,
    data = nest_data
  )
  
  # Estimate variation around nest count trend
  sd(
    residuals(nest_trend),
    na.rm = TRUE
  )
}

make_annual_cormorant_data <- function(n_years, nest_mean, sd, nest_data) { #(number of total years, constant mean, sd around that mean)
  
  observed_years <- nrow(nest_data) # of already observed years
  
  additional_years <- n_years - observed_years # anumber mount of extra years to generate
  
  # If statement to avoid requests of less years
  if (additional_years < 0) { 
    stop("Requested years cannot be fewer than observed years.") # Stops the function
  }
  
  if (additional_years > 0) {
    #create dataset to store additional years
    future_data <- tibble( 
      Year = max(nest_data$Year) + seq_len(additional_years),
      Active_Nests = round( # round so there are not decimal nest counts
        pmax(
          0,
          rnorm(
            additional_years,
            mean = nest_mean,
            sd = sd
          )
        )
      )
    )
    
    annual_extended <- bind_rows(nest_data, future_data)
    
  } else {
    
    annual_extended <- nest_data
  }
  
  annual_extended |>
    mutate(
      Active_Nests_z =
        (Active_Nests - mean(Active_Nests)) / sd(Active_Nests),
      Year_z =
        (Year - mean(Year)) / sd(Year),
      Year_f = factor(Year)
    )
}

make_output_design <- function(n_years, cormorantData, sites) {
  
  annual_subset <- cormorantData |> # adjusts cormorant data to the number of years of goby data needed
    slice_head(n = n_years)
  
  crossing(
    Year = annual_subset$Year,
    Site = sites
  ) |>
    left_join(
      annual_subset,
      by = "Year"
    ) |>
    mutate(
      Site = factor(Site, levels = sites),
      Year_f = factor(Year),
      Fish_Abundance_z = 0
    )
}

# Create datasets for all 4 year amounts
create_four_designs <- function(proposedAvgNests = 4000, cormorant_residual_sd, sites, nest_data){ #4000 avg nests
  list(
    "21 years" = make_output_design(21, make_annual_cormorant_data(21,proposedAvgNests,cormorant_residual_sd, nest_data), sites),
    "25 years" = make_output_design(25, make_annual_cormorant_data(25,proposedAvgNests,cormorant_residual_sd, nest_data), sites),
    "30 years" = make_output_design(30, make_annual_cormorant_data(30,proposedAvgNests,cormorant_residual_sd, nest_data), sites),
    "40 years" = make_output_design(40, make_annual_cormorant_data(40,proposedAvgNests,cormorant_residual_sd, nest_data), sites)
  )
}

run_simulation <- function(
    design,
    sim_formula = ~
      Active_Nests_z + Year_z +
      (1 | Year_f) + (1 | Site),
    model_formula = Fish_Abundance_z ~
      Active_Nests_z + Year_z +
      (1 | Year_f) + (1 | Site),
    beta_true,
    beta_year = -0.4,
    residual_sd = 0.7,
    ryear_sd = 0.4,
    rsite_sd = 0.5,
    nsim = 10,
    seed = 123
) {
  
  # Define parameters
  params <- list(
    beta = c(
      0, #intercept
      beta_true,
      beta_year
    ),
    betadisp = 2 * log(residual_sd),
    theta = c(
      log(ryear_sd),
      log(rsite_sd)
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
  
  ## Fit the model to each simulated dataset
  # Save the results of each simulation
  results <- map_dfr(seq_len(nsim), function(i) {
    
    sim_data <- design
    
    sim_data$Fish_Abundance_z <-
      simulated_values[[i]]
    
    fit <-
      glmmTMB(
        model_formula,
        family = gaussian,
        data = sim_data
      )
    
    tidy_fit <-
      broom.mixed::tidy(
        fit,
        effects = "fixed",
        component = "cond",
        conf.int = TRUE,
        conf.level = 0.95,
        conf.method = "wald"
      )
    # Save the effect of cormorants on nests
    cormorant_result <- tidy_fit |>
      filter(term == "Active_Nests_z")
    
    tibble(
      replicate = i,
      estimate = cormorant_result$estimate[1],
      std.error = cormorant_result$std.error[1],
      conf.low = cormorant_result$conf.low[1],
      conf.high = cormorant_result$conf.high[1],
    )
  })
  
  # Calculate analysis measures: estimate mean, bias, RMSE, CI width, power
  summary_results <- results |>
    summarise(
      
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
      
      power = mean(conf.high < 0)
    )
  
  list(
    replicates = results,
    summary = summary_results
  )
}

run_power_analysis <- function(
    designs,
    beta_cormorant = c(0,-0.1,-0.2,-0.3), #list of beta values
    beta_year,
    residual_sd,
    ryear_sd,
    rsite_sd,
    nsim = 10
){
  ## prepping running loop of simulations
  # Creating result lists
  all_replicates <- list()
  all_summaries <- list()
  
  #counter for changing seed (gives same result each power analysis can remove)
  counter <- 1
  
  for (design_name in names(designs)) {
    
    current_design <- designs[[design_name]]
    
    for (j in 1:length(beta_cormorant)) {
      
      current_beta <- beta_cormorant[j]
      
      # Run simulations and model fits
      output <- run_simulation(
        design = current_design,
        beta_true = current_beta,
        nsim = 10, # number of simulations
        seed = 1000 + counter
      )
      
      # Save individual simulation results
      all_replicates[[counter]] <-
        output$replicates |>
        mutate(
          sample_size = design_name,
          n_years = length(unique(current_design$Year)),
          beta_true = current_beta
        )
      
      # Save summary statistics
      all_summaries[[counter]] <-
        output$summary |>
        mutate(
          sample_size = design_name,
          n_years = length(unique(current_design$Year)),
          beta_true = current_beta
        )
      
      counter <- counter + 1
    }
  }
  return(list(bind_rows(all_replicates), bind_rows(all_summaries))) # output all replicate data and summary data
}

plot_power_analysis <- function(power_analysis_summary = power_analysis[[2]]){
  plot_data <- power_analysis_summary |>
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
  
  print(p_coverage)
  
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
  
  print(p_power)
  
  p_bias <- ggplot(
    plot_data,
    aes(
      x = n_years,
      y = bias,
      colour = beta_group,
      group = beta_group
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
  
  print(p_bias)
  
  p_rmse <- ggplot(
    plot_data,
    aes(
      x = n_years,
      y = RMSE,
      colour = beta_group,
      group = beta_group
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
  
  print(p_rmse)
  
  return((p_power + p_coverage)/(p_bias + p_rmse))
}
