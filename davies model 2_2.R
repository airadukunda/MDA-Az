
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#Modeling the Population-Level Impact of Azithromycin Mass Drug Administration on the Emergence, Persistence, and Spread of Antimicrobial-Resistant Echerchia Coli    #>> 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 
# This model use the demographic and social contacts patterns
# from 3.1.EpiSignalDetetion_1_OneYearContact .
# In the above mentioned script, you need to change
# the WHO EMRO countries by WHO Afro countries to get Tanzania
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#Work envronment
getwd()
setwd("C:/Disk F/4.Oxford Modelling for Global Health/Afox_Ubuntu/Afox Placement with Ben Cooper")
getwd()
#Packages
pacman::p_load(deSolve, viridis, ggplot2, tidyr, dplyr, readr) 
# Control parameters
save.fig <- FALSE
# 1. Demographic parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.a. Age structure (101 groups: 0,1, 100+)
age_groups <- c(as.character(0:99), "100+")
# Number of age groups
n_age <- length(age_groups)
# Aging matrix
A <- n_age  # 101 age groups (0,...100+)
dd <- rep(1, A)  # 
ageing <- t(diff(diag(dd), lag = 1) / (1 * 365.25)) # 1 year age groups
ageing <- cbind(ageing, rep(0, A))  # No ageing from last compartment
# Population structure
setwd("C:/Disk F/4.Oxford Modelling for Global Health/Placement project disk")

#Population structure 2000-2023 

Population_emro_2023 <- read_csv("Population_emro_2023_1yearage.csv")
table(Population_emro_2023$Country)
head(Population_emro_2023)
Tanzania_pop <- as.data.frame(Population_emro_2023 %>%
                                filter(Year == "2023") %>%
                                filter(Country == "United Republic of Tanzania"))                     
head(Tanzania_pop)
# Population structure in thousands 
Tanzania_pop_in_thousands <- as.data.frame(Tanzania_pop %>%
                                             mutate(Population_age = Population_age * 1000,
                                                    Annual_population = Annual_population * 1000))

# 1.c.Visualization of Tanzania's population structure
pacman::p_load(ggplot2, dplyr, plotly)

pyramyd <- ggplot(Tanzania_pop_in_thousands, aes(x = Age_Category, y = Population_age)) +
  geom_col(fill = "steelblue") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +  
  labs(title = "Tanzania's population structure (2023)",
       x = "Age", y = "Population( in million)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(pyramyd)
#  Shorten the name
popstruc <- Tanzania_pop_in_thousands
#  Number of age groups 
A <- n_age <- nrow(Tanzania_pop_in_thousands)  # Fixed: use nrow instead of length
#  Births by age of mother

popbirth <- read.csv("3.U.1.Birth_1year_Afro.csv", header = TRUE)

table(popbirth$Year)
print(popbirth)
table(popbirth$Country)
str(popbirth)
#  Head(popbirth)
popbirth <- popbirth %>%
  filter(Country == "United Republic of Tanzania") %>%
  filter(Year == 2023)
#
names(popbirth)
pyramyd_birth <- ggplot(popbirth, aes(x = Age, y = Births)) +
  geom_col(fill = "blue") +
  labs(title = "Birth in Tanzania (2023)",
       x = "Age", y = "Births") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(pyramyd_birth)

#Convert from 1000s per 1 year period to per person per day
group_durations_years <- rep(1, n_age)
group_durations_days <- group_durations_years * 365
print(popbirth)
popbirth[, 5] <- 1000 * popbirth[, 5] / (1 * popstruc[, 5] * 365.25)
# 1.g. Natural mortality per person per year
mortality <- read.csv("3.U.1.EMRO_mortality_by_age_group_1yearage.csv", header = TRUE)

table(mortality$Country)

popmort <- as.data.frame(mortality %>% 
                           filter(Country == "United Republic of Tanzania") %>%
                           filter(Year == 2023))

pyramyd_mort <- ggplot(popmort, aes(x = Age, y = Percentage)) +
  geom_col(fill = "red") +
  labs(title = "Mortality in Tanzania (2023)",
       x = "Age", y = "Deaths per 1000 pop") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(pyramyd_mort)

# Convert mortality from 1000s per 1 year period to per person per day
group_durations_years <- rep(1, n_age)
group_durations_days <- group_durations_years * 365

popmort[, 5] <- 1000 * popmort[, 5] / (1 * popstruc[, 5] * 365.25)
mort <- popmort[, 5]

#1.h.Contact matrix  and visualization
(m_contact_1y <- as.matrix(read.csv("3.U.1.contact_Tanzania_1y.csv")))
dim(m_contact_1y)

for(i in 1:n_age){
  for(j in 1:n_age){
    m_contact_1y[i,j]<-m_contact_1y[i,j]/25
  }
}
colnames(m_contact_1y ) <-c(as.character(0:99), "100+")
rownames(m_contact_1y ) <- c(as.character(0:99), "100+")
m_contact_1y 
#Visualization of my contact matrix
pacman::p_load(ggplot2,reshape2)
#data frame
#df <- melt(m_contact_1y )
df <- reshape2::melt(m_contact_1y)
colnames(df) <- c("Contactee", "Contactor", "Contacts")
#plot contact matrix
(p<-ggplot(df, aes(x = Contactor, y = Contactee, fill = Contacts)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    labs(title = "Contact Matrix Heatmap for Tanzania")
)
# Back into work environment
setwd("C:/Disk F/4.Oxford Modelling for Global Health/Afox_Ubuntu/Afox Placement with Ben Cooper")

# Indices for compartments 
Xindex <- 1:n_age                        # Uninfected, untreated
Sindex <- (1*n_age+1):(2*n_age)          # drug-sensitive, untreated
Rindex <- (2*n_age+1):(3*n_age)          # drug-resistant, untreated
Srindex <- (3*n_age+1):(4*n_age)         # drug-sensitive, treated
Rsindex <- (4*n_age+1):(5*n_age)         # drug-resistant, treated
#...............................................................................
# Considering how to add MDA as an intervention
 # mda_start_times <- c(365, 79, 4380, 4745)
 (mda_start_times<-(2:10)*365.25)
 (mda_duration <- 30)
 
 mda_active <- function(time, mda_starts, duration) {
   any(sapply(mda_starts, function(start) {
     time >= start && time < (start + duration)
   }))
 }
#..............................................................................

 bacteria.odes <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Extract compartments 
    X <- state[Xindex]   # Uninfected, untreated
    S <- state[Sindex]   # drug-sensitive, untreated
    R <- state[Rindex]   # drug-resistant, untreated
    Sr <- state[Srindex] # drug-sensitive, treated
    Rs <- state[Rsindex] # drug-resistant, treated
    # Total population
    N <- X + S + R + Sr + Rs 
    # Births
    b1 <- sum(popbirth[, 5] * N)
    births <- rep(0, n_age)
    births[1] <- b1   
    #
    S.tot <- S + Sr  # Susceptible co-colonised total
    R.tot <- R + Rs  # Resistance co-colonised total
    
    
    lamda.S <- beta.S * (m_contact_1y %*% S.tot)   # vectorized force of infection
    lamda.R <- beta.R * (m_contact_1y %*% R.tot)
    
    #lamda.S <- beta.S * (m_contact_1y %*% (S.tot/N))
    #lamda.R <- beta.R * (m_contact_1y %*% (R.tot/N))
    
    #Intervention : MDA implementation
    
    #...........................................................................
    #    conidering how to add the   MDA intervention into the odes
    mda_targeted_ages <- 0:5         # Targeting 
    azt <- rep(0, n_age)             # Initialize azt vector
    azt[mda_targeted_ages] <- 1      # Ages targeted 
    b <- ifelse(mda_active(t, mda_start_times, mda_duration), (1*12 / 365.25), 0)
    a <- b * azt                     # Apply a only on targeted ages
    #............................................................................
    
    # ODEs system
    
    dX <- births + u.S * S + u.R * R + u.c * (Sr + Rs) - lamda.S * X - lamda.R * X +  ageing %*% X - mort * X  # + MDA * S # MDA * (1-kappa) * S  # Successfully treated return to X 
    
    dS <- lamda.S * X - S * u.S - S * a - k * lamda.R * S +  ageing %*% S - mort * S                           # - MDA  * S                         # Remove all in S
    
    dR <- lamda.R * X - u.R * R - k * lamda.S * R + a * (Sr + Rs) +  ageing %*% R - mort * R                            # MDA * kappa * S # Treated sensitive may become resistant 
    
    dSr <- k * lamda.R * S - Sr * u.c - Sr * a +  ageing %*% Sr - mort * Sr  
    
    dRs <- k * lamda.S * R - Rs * u.c - Rs * a + ageing %*% Rs - mort * Rs  #a=epsilon
    
    list(c(dX, dS, dR, dSr, dRs))
  })
}
bacteria.solve <- function(t, state, parameters) {
  parameters[["beta.R"]] <- parameters[["beta.S"]] * (1 - parameters[["c"]])
  out_Tanzania <- as.data.frame(ode(state, t, bacteria.odes, parameters))
  return(out_Tanzania)
}


# Parameters
parms.orig <- list(
  beta.S = 5,              # transmission of sensitive         : (β = 5 month−1)   
  u.S = 1,                 # Clearance sensitive (natural)     : (u = 1 month−1)
  u.R = 1,                 # Clearance resistant (natural)     : (u = 1 month−1) :lower than susceptible?
  u.c = 1,                 # Clearance co-colonised (natural)
  a = 0.16,                # Clearance sensitive (drug-induced)
  k=  0.05,                # The efficiency of co-colonisation : (k= 0.25, 0.5, and 1.0)
  m_contact = m_contact_1y, # Social contacts per day
  mda_cycle= 365,
  mda_duration = 90,
  mda_cov =  0.8,
  kappa =    0    # 0.05  # Proportion that develop/select resistance (Assumed)
)
# Convert daily
parms.orig[1:5] <- lapply(parms.orig[1:5], function(x) x*12/365.25)# convert to daily

# MDA parameters
#parms.orig["mda_cov"]= 0
#parms.orig["mda_cycle"]mda_cycle = 365        # Period between 2 mda
#parms.orig["mda_duration"] = 2*30             # mda campaign duration in days

parms.orig[["c"]] <- 0.1  # cost of resistance on transmission
parms.orig[["k"]] <- 0.5  # efficiency of co-colonisation
# Initial conditions
state.orig <- c(
  X = rep(0.95, n_age),  # Uninfected, untreated
  S = rep(0.025, n_age), # drug-sensitive, untreated
  R = rep(0.025, n_age), # drug-resistant, untreated
  Sr = rep(0, n_age),    # drug-sensitive, treated
  Rs = rep(0, n_age)     # drug-resistant, treated
)
tvec_1 <- seq(0, 365.25 *2,1)       #Pre-MDA
tvec_5 <- seq(0, 365.25 *6,1)       #5 years of MDA
tvec_10_a <- seq(0, 365.25 *11,1)   #10 years of MDA
tvec_10_b <- seq(0, 365.25 *11,1)   #10 years No MDA, i will need a = 10
# Run model  
parms <- parms.orig
#parms.orig_0 <- parms.orig
#parms.orig_0["mda_cov"]= 0
#parms_0<-parms.orig_0
state <- state.orig
start <- Sys.time()
out_1_Tanzania <- bacteria.solve(tvec_1, state, parms)
out_5_Tanzania <- bacteria.solve(tvec_5, state, parms)
out_10_a_Tanzania <- bacteria.solve(tvec_10_a, state, parms)
# For no MDA scenario,set a =0 or start times very high
(mda_start_times<-12*365.25)
out_10_b_Tanzania <- bacteria.solve(tvec_10_b, state, parms)
end<-Sys.time()
end-start
#
results_1_Tanzania <- as.data.frame(out_1_Tanzania)
results_5_Tanzania <- as.data.frame(out_5_Tanzania)
results_10_a_Tanzania <- as.data.frame(out_10_a_Tanzania)
results_10_b_Tanzania <- as.data.frame(out_10_b_Tanzania)

# Column names 
compartment_names <- c("X", "S", "R", "Sr", "Rs")
col_names <- c("time")
for(comp in compartment_names) {
  for(age in age_groups) {
    col_names <- c(col_names, paste0(comp, "_", age))
  }
}
colnames(results_1_Tanzania) <- col_names
colnames(results_5_Tanzania) <- col_names
colnames(results_10_a_Tanzania) <- col_names
colnames(results_10_b_Tanzania) <- col_names
names(results_1_Tanzania)
names(results_5_Tanzania)
names(results_10_a_Tanzania)
names(results_10_b_Tanzania)
# Plot model
cols <- viridis(5)

par(mfrow = c(1, 1))
#Total across age groups: Here i will use Index
X_total <- rowSums(out_1_Tanzania[,  Xindex + 1])  # Here i added +1 because first column is time
S_total <- rowSums(out_1_Tanzania[,  Sindex + 1])
R_total <- rowSums(out_1_Tanzania[,  Rindex + 1]) 
Rs_total <- rowSums(out_1_Tanzania[, Rsindex + 1])
Sr_total <- rowSums(out_1_Tanzania[, Srindex + 1])
#
R_total_10_a <- rowSums(out_10_a_Tanzania[,  Rindex + 1]) #+rowSums(out_10_a_Tanzania[,  Rsindex + 1])  
R_total_10_b <- rowSums(out_10_b_Tanzania[,  Rindex + 1]) #+rowSums(out_10_b_Tanzania[,  Rsindex + 1]) 
Rs_total_10_a <- rowSums(out_10_a_Tanzania[,  Rsindex + 1]) 
Rs_total_10_b <- rowSums(out_10_b_Tanzania[,  Rsindex + 1]) 


#Visualization
plot(
  out_1_Tanzania$time, X_total,
  type = "l",
  col = cols[1],
  las = 1,
  xaxs = "i", yaxs = "i",
  ylim = c(0, 100),#max(c(X_total, S_total, R_total,Rs_total, Sr_total))),
  bty = "n",
  lwd = 3.5,
  xlab = "Day", ylab = "Proportion",
  main = "Population level Bacterial Dynamics over time"
)
lines(out_1_Tanzania$time, S_total, col = cols[2], lwd = 3.5)
lines(out_1_Tanzania$time, R_total, col = cols[3], lwd = 3.5)
lines(out_1_Tanzania$time, Rs_total, col = cols[4], lwd = 3.5)
lines(out_1_Tanzania$time, Sr_total, col = cols[5], lwd = 3.5)
legend(
  "topright",
  bty = "n",
  col = c(cols),
  legend = c("X_total", "S_total", "R_total","Rs_total", "Sr_total"),
  lty = 1,
  lwd = 3.5,
  ncol = 2
)
#Visualization for the ten years
#a.Visualization for the 10 years R
plot(
  out_10_a_Tanzania$time, R_total_10_a,
  type = "l",
  col = cols[1],
  las = 1,
  xaxs = "i", yaxs = "i",
  ylim = c(0,max(c(R_total_10_a,R_total_10_b))),
  bty = "n",
  lwd = 3.5,
  xlab = "Day", ylab = "Proportion (R)",
  main = "Population level Bacterial Dynamics over time"
)
abline(v = 365*2, col = "red", lwd = 2, lty = 2)
# Add text above or next to the vertical line
text(
  x = 1000,
  y = max(c(R_total_10_a, R_total_10_b)) * 0.90,   # vertical position (95% of ymax)
  labels = "MDA start",
  pos = 3,                  # 3 = above the point (can use 1,2,4)
  cex = 1,
  col = "red"
)
#lines(out_10_a_Tanzania$time, Rs_total_10_a, col = cols[2], lwd = 3.5)
lines(out_10_b_Tanzania$time, R_total_10_b, col = cols[4], lwd = 3.5)
#lines(out_10_b_Tanzania$time, Rs_total_10_b, col = cols[4], lwd = 3.5)
legend(
  "topright",
  bty = "n",
  col = c(cols[c(1,4)]),
  #col = c("red","blue"),
  legend = c("MDA", "No MDA"),
  lty = 1,
  lwd = 3.5,
  ncol = 2
)

#b.Visualization for the 10 years Rs
plot(
  out_10_b_Tanzania$time, Rs_total_10_a,
  type = "l",
  col = cols[1],
  las = 1,
  xaxs = "i", yaxs = "i",
  ylim = c(0,100),#max(c(Rs_total_10_a, Rs_total_10_b))),
  bty = "n",
  lwd = 3.5,
  xlab = "Day", ylab = "Proportion(Rs)",
  main = "Population level Bacterial Dynamics over time"
)
#lines(out_10_a_Tanzania$time, Rs_total_10_a, col = cols[2], lwd = 3.5)
lines(out_10_b_Tanzania$time, R_total_10_b, col = cols[3], lwd = 3.5)
#lines(out_10_b_Tanzania$time, Rs_total_10_b, col = cols[4], lwd = 3.5)
legend(
  "topright",
  bty = "n",
  col = c(cols[c(1,3)]),
  legend = c("MDA","No MDA"),
  lty = 1,
  lwd = 3.5,
  ncol = 2
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                    Additional-visualization                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
pacman::p_load(data.table)   # This package will allow us to reshape data set faster
results_1_Tanzania<- as.data.table(results_1_Tanzania)
#Long format : Here i will be using melt to be faster
results_1_Tanzania_long <- melt(
  results_1_Tanzania,
  id.vars = "time",
  variable.name = "variable",
  value.name = "value"
)
# data.table 
results_1_Tanzania_long <- as.data.table(results_1_Tanzania_long)
#
results_1_Tanzania_long[, c("compartment", "age_group") := tstrsplit(variable, "_", fixed = TRUE)]
results_1_Tanzania_long[, `:=`(
  age_group = factor(age_group, levels = age_groups),
  time_years = time / 365,
  time_days = time
)]
results_1_Tanzania_long[, variable := NULL]
#Structure of the long format
names(results_1_Tanzania_long)
head(results_1_Tanzania_long)
table(results_1_Tanzania_long$compartment)
#Here i will need to speed the visualization 
pacman::p_load(data.table,ggplot2,scales)
# Data.table 
setDT(results_1_Tanzania_long)
#Transformations in one data.table
# Proportions (in one step)
results_1_Tanzania_long[as.integer(age_group) <= 101, `:=`(
  total_by_age = sum(value), 
  proportion = value / sum(value) * 100
), by = age_group]
# Pre-aggregate data
Final_1_Tanzania_summary <- results_1_Tanzania_long[
  as.integer(age_group) <= 101,# Exclude 100 + age
  .(total_value = sum(value)), 
  by = .(age_group, compartment)
][, proportion := total_value / sum(total_value) * 100, by = age_group]
# Plot-aggregated data
names(Final_1_Tanzania_summary)
dim(Final_1_Tanzania_summary)
table(Final_1_Tanzania_summary$compartment)
plot_a <- ggplot(Final_1_Tanzania_summary, 
                 aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "stack",col= NA) +#or dodge
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = " Population level Bacterial Dynamics over time in Tanzania", x = "Age", y = "Percent ") +
  scale_fill_manual(values = c("X" = "#00c4aa", "S" = "#e573f3", 
                               "R" = "#00b3f4", "Sr" = "#9b9602", "Rs" = "#fc726c"))
print(plot_a)
plot_b <- ggplot(Final_1_Tanzania_summary, 
                 aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "dodge",col= NA) +#or dodge
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = " Population level Bacterial Dynamics over time in Tanzania", x = "Age", y = "Percent ") +
  scale_fill_manual(values = c("X" = "#00c4aa", "S" = "#e573f3", 
                               "R" = "#00b3f4", "Sr" = "#9b9602", "Rs" = "#fc726c"))
print(plot_b)
table(Final_1_Tanzania_summary$compartment)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#  MDA visualization                                                           #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ##############################################################################
# 4  Scenarios: 1, 5, 10 Years MDA and 10 No MDA
################################################################################
# 
# Age-Stratified Comparison of Three Time Scenarios: 1, 5, and 10 Years
# 
pacman::p_load(data.table, ggplot2, scales, dplyr, patchwork)
# 
# Data for all three scenarios
# 
# Results for each scenario
process_scenario <- function(results_data_Tanzania, scenario_name) {
  # Column names 
  compartment_names <- c("X", "S", "R", "Sr", "Rs")
  col_names <- c("time")
  for(comp in compartment_names) {
    for(age in age_groups) {
      col_names <- c(col_names, paste0(comp, "_", age))
    }
  }
  colnames(results_data_Tanzania) <- col_names
  
  # Conversion in  data.table
  results_dt <- as.data.table(results_data_Tanzania)
  
  # Long format
  results_long <- melt(
    results_dt,
    id.vars = "time",
    variable.name = "variable",
    value.name = "value"
  )
  
  # I Will need Split compartment and age
  results_long[, c("compartment", "age_group") := tstrsplit(variable, "_", fixed = TRUE)]
  results_long[, `:=`(
    age_group = factor(age_group, levels = age_groups),
    time_years = time / 365,
    time_days = time,
    scenario = scenario_name
  )]
  results_long[, variable := NULL]
  
  # Calculate proportions by age
  results_long[as.integer(age_group) <= 101, `:=`(
    total_by_age = sum(value), 
    proportion = value / sum(value) * 100
  ), by = .(age_group, scenario)]
  
  # Aggregated data
  final_summary <- results_long[
    as.integer(age_group) <= 101,
    .(total_value = sum(value)), 
    by = .(age_group, compartment, scenario)
  ][, proportion := total_value / sum(total_value) * 100, by = .(age_group, scenario)]
  
  return(final_summary)
}

# Process all three scenarios
scenario_1yr_Tanzania <- process_scenario(results_1_Tanzania, "Pre-MDA")
scenario_5yr_Tanzania <- process_scenario(results_5_Tanzania, "5 Years MDA")
scenario_10yr_Tanzania <- process_scenario(results_10_a_Tanzania, "10 Years MDA")
scenario_10yr_no_MDA_Tanzania <- process_scenario(results_10_b_Tanzania, "10 Years No MDA")

# All scenarios
all_scenarios_Tanzania <- rbind(scenario_1yr_Tanzania, scenario_5yr_Tanzania, scenario_10yr_Tanzania,scenario_10yr_no_MDA_Tanzania)

# Factor levels for scenarios
all_scenarios_Tanzania$scenario <- factor(all_scenarios_Tanzania$scenario, 
                                          levels = c("Pre-MDA", "5 Years MDA", "10 Years MDA","10 Years No MDA"))

table(all_scenarios_Tanzania$scenario)

all_scenarios_Tanzania <-all_scenarios_Tanzania|>
    filter(scenario ==c("10 Years MDA","10 Years No MDA"))
  
#I want to remove one scenario
#all_scenarios_Tanzania <-all_scenarios_Tanzania|>
#  filter(scenario !="10 Years No MDA")
# 
table(all_scenarios_Tanzania$compartment)
all_scenarios_Tanzania <-all_scenarios_Tanzania|>
    filter(compartment==c("R","Rs"))
  
# Plot 1: Stacked bar plots comparing all three scenarios
# 
plot_stacked <- ggplot(all_scenarios_Tanzania, 
                       aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "stack", col = NA) +
  facet_wrap(~scenario, ncol = 2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  ) +
  labs(
    title = "Population level Bacterial Dynamics over time in Tanzania",
    subtitle = "",
    x = "Age", 
    y = "Percent",
    fill = "Compartment"
  ) +
  scale_fill_manual(values = c(
    "X" = "#00c4aa", 
    "S" = "#e573f3", 
    "R" = "#00b3f4", 
    "Sr" = "#9b9602", 
    "Rs" = "#fc726c"
  ))

print(plot_stacked)

#
#  Plot 2: Dodged bar plots comparing all three scenarios
# 
plot_dodged <- ggplot(all_scenarios_Tanzania, 
                      aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "dodge", col = NA) +
  facet_wrap(~scenario, ncol = 2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  ) +
  labs(
    title = "Population level Bacterial Dynamics over time in Tanzania",
    subtitle = "",
    x = "Age", 
    y = "Percent",
    fill = "Compartment"
  ) +
  scale_fill_manual(values = c(
    "X" = "#00c4aa", 
    "S" = "#e573f3", 
    "R" = "#00b3f4", 
    "Sr" = "#9b9602", 
    "Rs" = "#fc726c"
  ))

print(plot_dodged)
# 
# Plot 3: Side-by-side comparison (stacked) - all scenarios in one row
#
plot_horizontal <- ggplot(all_scenarios_Tanzania, 
                          aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "stack", col = NA) +
  facet_wrap(~scenario, ncol = 2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  ) +
  labs(
    title = "Population level Bacterial Dynamics over time in Tanzania",
    subtitle = "",
    x = "Age", 
    y = "Percent",
    fill = "Compartment"
  ) +
  scale_fill_manual(values = c(
    "X" = "#00c4aa", 
    "S" = "#e573f3", 
    "R" = "#00b3f4", 
    "Sr" = "#9b9602", 
    "Rs" = "#fc726c"
  ))

print(plot_horizontal)

# 
# Plot 4: Focus on Resistance (R) only across scenarios
# 
resistance_only <- all_scenarios_Tanzania[compartment == "R"]

plot_resistance <- ggplot(resistance_only, 
                          aes(x = age_group, y = proportion, fill = scenario)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  ) +
  labs(
    title = "Resistant Bacteria (R) across age groups in Tanzania",
    subtitle = "",
    x = "Age", 
    y = "Percent Resistant",
    fill = "Scenario"
  ) +
  scale_fill_viridis_d(option = "plasma", end = 0.8)

print(plot_resistance)

#
# Plot 5: Total Resistance (R + Rs) comparison
# 
total_resistance <- all_scenarios_Tanzania[compartment %in% c("R", "Rs"), 
                                           .(proportion = sum(proportion)), 
                                           by = .(age_group, scenario)]

plot_total_resistance <- ggplot(total_resistance, 
                                aes(x = age_group, y = proportion, 
                                    color = scenario, group = scenario)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 0.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  ) +
  labs(
    title = "Total resistance (R + Rs) across age groups in  Tanzania",
    subtitle = "",
    x = "Age", 
    y = "Percent with Resistance",
    color = "Scenario"
  ) +
  scale_color_viridis_d(option = "plasma", end = 0.8)

print(plot_total_resistance)

# 
# Plot 6: Individual compartments across scenarios (small multiples)
# 
plot_by_compartment <- ggplot(all_scenarios_Tanzania, 
                              aes(x = age_group, y = proportion, fill = scenario)) +
  geom_col(position = "dodge") +
  facet_wrap(~compartment, scales = "free_y", ncol = 2) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  ) +
  labs(
    title = "Compartment-specific Dynamics across scenarios",
    subtitle = "",
    x = "Age", 
    y = "Percent",
    fill = "Scenario"
  ) +
  scale_fill_viridis_d(option = "plasma", end = 0.8)

print(plot_by_compartment)

# Save plots (optional)
#
if(save.fig) {
  ggsave("comparison_stacked.png", plot_stacked, width = 12, height = 10, dpi = 300)
  ggsave("comparison_dodged.png", plot_dodged, width = 12, height = 10, dpi = 300)
  ggsave("comparison_horizontal.png", plot_horizontal, width = 16, height = 6, dpi = 300)
  ggsave("resistance_comparison.png", plot_resistance, width = 12, height = 6, dpi = 300)
  ggsave("total_resistance_line.png", plot_total_resistance, width = 12, height = 6, dpi = 300)
  ggsave("compartment_multiples.png", plot_by_compartment, width = 12, height = 10, dpi = 300)
  cat("\nPlots saved successfully!\n")
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Resistance consumption curve at k = 0.5
avec <- seq(0, 2, 0.1)
rvec <- rep(NA, length(avec))

for (ii in seq_along(avec)) {
  parms[["a"]] <- avec[ii] / 365.25
  out_Tanzania <- bacteria.solve(tvec_1, state, parms)
  
  # Resistance proportion 
  R_final <- rowSums(out_Tanzania[nrow(out_Tanzania), Rindex + 1]) + rowSums(out_Tanzania[nrow(out_Tanzania), Rsindex + 1])
  X_final <- rowSums(out_Tanzania[nrow(out_Tanzania), Xindex + 1])
  
  rvec[ii] <- R_final / (1 - X_final)
}

# Plot resistance vs consumption
plot(
  avec, rvec,
  type = "l",
  col = cols[1],
  las = 1,
  xaxs = "i", yaxs = "i",
  ylim = c(0, 1),
  bty = "n",
  lwd = 2.5,
  xlab = "Consumption", ylab = "Resistance Proportion",
  main = "Resistance and consumption"
)

