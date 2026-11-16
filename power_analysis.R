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
WQData <- read_csv("Data/WQ_2004-2020.csv")

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
  
  theme_bw() +
  
  geom_smooth(
    aes(y=cpue_z),
    colour = "red",
    method = lm,
    alpha = 0.5) +
  
  #  facet_wrap(~ Site) +
  
  geom_smooth(
    aes(y=cpue_zsim1),
    colour = "blue",
    method = lm,
    alpha = 0.5) +
  
  labs(
    subtitle = "Red = observed; Blue = simulated"
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