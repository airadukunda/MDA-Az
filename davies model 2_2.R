
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
ageing <- t(diff(diag(dd), lag = 1) / (1 * 365.25))
ageing <- cbind(ageing, rep(0, A))  # No ageing from last compartment
# Population structure
setwd("C:/Disk F/4.Oxford Modelling for Global Health/Placement project disk")

#Population structure 2000-2023                             
Population_emro_2023 <- read_csv("Population_emro_2023_1yearage.csv")
table(Population_emro_2023$Country)
head(Population_emro_2023)
Afghanistan_pop <- as.data.frame(Population_emro_2023 %>%
                                   filter(Year == "2023") %>%
                                   filter(Country == "Afghanistan"))                     
head(Afghanistan_pop)
# Population structure in thousands 
Afghanistan_pop_in_thousands <- as.data.frame(Afghanistan_pop %>%
                                                mutate(Population_age = Population_age * 1000,
                                                       Annual_population = Annual_population * 1000))

# 1.c. Visualization of Afghanistan's population structure
pacman::p_load(ggplot2, dplyr, plotly)

pyramyd <- ggplot(Afghanistan_pop_in_thousands, aes(x = Age_Category, y = Population_age)) +
  geom_col(fill = "steelblue") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +  
  labs(title = "Afghanistan's population structure (2023)",
       x = "Age", y = "Population( in million)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(pyramyd)
#  Shorten the name
popstruc <- Afghanistan_pop_in_thousands
#  Number of age groups 
A <- n_age <- nrow(Afghanistan_pop_in_thousands)  # Fixed: use nrow instead of length
#  Births by age of mother
popbirth <- read.csv("3.U.1.Birth_1year_emro.csv", header = TRUE)
table(popbirth$Country)
#  Head(popbirth)
popbirth <- popbirth %>%
  filter(Country == "Afghanistan") %>%
  filter(Year == 2023)
#
pyramyd_birth <- ggplot(popbirth, aes(x = Age, y = Birth)) +
  geom_col(fill = "blue") +
  labs(title = "Birth in Afghanistan (2023)",
       x = "Age", y = "Births") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(pyramyd_birth)

# Convert from 1000s per 1 year period to per person per day
group_durations_years <- rep(1, n_age)
group_durations_days <- group_durations_years * 365

popbirth[, 6] <- 1000 * popbirth[, 6] / (1 * popstruc[, 5] * 365.25)
# 1.g. Natural mortality per person per year
mortality <- read.csv("3.U.1.EMRO_mortality_by_age_group_1yearage.csv", header = TRUE)

popmort <- as.data.frame(mortality %>% 
                           filter(Country == "Afghanistan") %>%
                           filter(Year == 2023))

pyramyd_mort <- ggplot(popmort, aes(x = Age, y = Percentage)) +
  geom_col(fill = "red") +
  labs(title = "Mortality in Afghanistan (2023)",
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
(m_contact_1y <- as.matrix(read.csv("3.U.1.contact_Afghanistan_1y.csv")))
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
    labs(title = "Contact Matrix Heatmap for Afghanistan")
)

# Back into work environment
setwd("C:/Disk F/4.Oxford Modelling for Global Health/Afox_Ubuntu/Afox Placement with Ben Cooper")

# Indices for compartments 
Xindex <- 1:n_age                        # Uninfected, untreated
Sindex <- (1*n_age+1):(2*n_age)          # drug-sensitive, untreated
Rindex <- (2*n_age+1):(3*n_age)          # drug-resistant, untreated
Srindex <- (3*n_age+1):(4*n_age)         # drug-sensitive, treated
Rsindex <- (4*n_age+1):(5*n_age)         # drug-resistant, treated

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
    b1 <- sum(popbirth[, 6] * N)
    births <- rep(0, n_age)
    births[1] <- b1
    
    S.tot <- S + Sr  # Susceptible co-colonised total
    R.tot <- R + Rs  # Resistance co-colonised total
    
    #lamda.S <- beta.S *(m_contact_1y %*% sum(S.tot)) # 
    #lamda.R <- beta.R *(m_contact_1y %*%  sum(R.tot)) # 
    
    lamda.S <- beta.S * (m_contact_1y %*% S.tot)   # vectorized force of infection
    lamda.R <- beta.R * (m_contact_1y %*% R.tot)
    
    
  
    #Intervention : MDA implementation
    
    #mda_targeted_ages <- 0:5        # Targeting 
    #azt <- rep(0, n_age)            # Initialize azt vector
    #azt[mda_targeted_ages] <- 1 
    
    #Frequency of MDA (once per year according to the WHO guidlnes)
    #time_in_cycle <- t %% mda_cycle
    # mda_active <- (time_in_cycle < mda_duration) * 1
    # mda starting after 2025 (25 years from start year 2000)
    # mda_rate <- ifelse(((t > 25*365) & mda_active), mda_cov, 0)
    
    # ODEs
    dX <- births + u.S * S + u.R * R + u.c * (Sr + Rs) - lamda.S * X - lamda.R * X + 
      ageing %*% X - mort * X  
    
    dS <- lamda.S * X - S * u.S - S * a - k * lamda.R * S + 
      ageing %*% S - mort * S  
    
    dR <- lamda.R * X - u.R * R - k * lamda.S * R + a * (Sr + Rs) + 
      ageing %*% R - mort * R  
    
    dSr <- k * lamda.R * S - Sr * u.c - Sr * a + 
      ageing %*% Sr - mort * Sr  
    
    dRs <- k * lamda.S * R - Rs * u.c - Rs * a + 
      ageing %*% Rs - mort * Rs  
    
    list(c(dX, dS, dR, dSr, dRs))
  })
}
bacteria.solve <- function(t, state, parameters) {
  parameters[["beta.R"]] <- parameters[["beta.S"]] * (1 - parameters[["c"]])
  out <- as.data.frame(ode(state, t, bacteria.odes, parameters))
  return(out)
}


# Parameters
parms.orig <- list(
  beta.S = 5,              # transmission of sensitive
  u.S = 1,                 # Clearance sensitive (natural)
  u.R = 1,                 # Clearance resistant (natural) - lower than susc
  u.c = 1,                 # Clearance co-colonised (natural)
  a = 0.16,                # Clearance sensitive (drug-induced)
  m_contact = m_contact_1y #contacts parameters (per day)
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
tvec <- seq(0, 365.25 *1, 1) # 20 years
# Run model
parms <- parms.orig
state <- state.orig
start <- Sys.time()
out <- bacteria.solve(tvec, state, parms)
end<-Sys.time()
end-start
#
results_Afghanistan <- as.data.frame(out)
# Column names 
compartment_names <- c("X", "S", "R", "Sr", "Rs")

col_names <- c("time")
for(comp in compartment_names) {
  for(age in age_groups) {
    col_names <- c(col_names, paste0(comp, "_", age))
  }
}
colnames(results_Afghanistan) <- col_names
names(results_Afghanistan)

# Plot model
cols <- viridis(5)

par(mfrow = c(1, 1))
#Total across age groups: Here i will use Index
X_total <- rowSums(out[,  Xindex + 1])  # Here i added +1 because first column is time
S_total <- rowSums(out[,  Sindex + 1])
R_total <- rowSums(out[,  Rindex + 1]) 
Rs_total <- rowSums(out[, Rsindex + 1])
Sr_total <- rowSums(out[, Srindex + 1])
#Visualization
plot(
  out$time, X_total,
  type = "l",
  col = cols[1],
  las = 1,
  xaxs = "i", yaxs = "i",
  ylim = c(0, max(c(X_total, S_total, R_total,Rs_total, Sr_total))),
  bty = "n",
  lwd = 2.5,
  xlab = "Day", ylab = "Proportion",
  main = "Population level Bacterial Dynamics over time"
)
lines(out$time, S_total, col = cols[2], lwd = 2.5)
lines(out$time, R_total, col = cols[3], lwd = 2.5)
lines(out$time, Rs_total, col = cols[4], lwd = 2.5)
lines(out$time, Sr_total, col = cols[5], lwd = 2.5)

legend(
  "topright",
  bty = "n",
  col = c(cols),
  legend = c("X_total", "S_total", "R_total","Rs_total", "Sr_total"),
  lty = 1,
  lwd = 1.2,
  ncol = 2
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Additional-visualization 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(data.table)   # This package will allow us to reshape dataset faster
results_Afghanistan_0<- as.data.table(results_Afghanistan)
#Long format : Here i will be using melt to be faster
results_Afghanistan_long <- melt(
  results_Afghanistan_0,
  id.vars = "time",
  variable.name = "variable",
  value.name = "value"
)
# data.table 
results_Afghanistan_long <- as.data.table(results_Afghanistan_long)
#
results_Afghanistan_long[, c("compartment", "age_group") := tstrsplit(variable, "_", fixed = TRUE)]
results_Afghanistan_long[, `:=`(
  age_group = factor(age_group, levels = age_groups),
  time_years = time / 365,
  time_days = time
)]
results_Afghanistan_long[, variable := NULL]
#Structure of the long format
names(results_Afghanistan_long)
head(results_Afghanistan_long)
table(results_Afghanistan_long$compartment)
#Here i will need to speed the visualization 
pacman::p_load(data.table,ggplot2,scales)
# Data.table 
setDT(results_Afghanistan_long)
#Transformations in one data.table
# Proportions (in one step)
results_Afghanistan_long[as.integer(age_group) <= 101, `:=`(
  total_by_age = sum(value), 
  proportion = value / sum(value) * 100
), by = age_group]
# Pre-aggregate data
Final_Afghanistan_summary <- results_Afghanistan_long[
  as.integer(age_group) <= 100,# Exclude 100 + age
  .(total_value = sum(value)), 
  by = .(age_group, compartment)
][, proportion := total_value / sum(total_value) * 100, by = age_group]
# Plot-aggregated data
names(Final_Afghanistan_summary)
dim(Final_Afghanistan_summary)
table(Final_Afghanistan_summary$compartment)
plot <- ggplot(Final_Afghanistan_summary, 
                    aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "stack",col= NA) +#or dodge
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = " Population level Bacterial Dynamics over time in Afghanistan", x = "Age", y = "Percent ") +
  scale_fill_manual(values = c("X" = "#00c4aa", "S" = "#e573f3", 
                               "R" = "#00b3f4", "Sr" = "#9b9602", "Rs" = "#fc726c"))
print(plot)
plot1 <- ggplot(Final_Afghanistan_summary, 
               aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "dodge",col= NA) +#or dodge
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = " Population level Bacterial Dynamics over time in Afghanistan", x = "Age", y = "Percent ") +
  scale_fill_manual(values = c("X" = "#00c4aa", "S" = "#e573f3", 
                               "R" = "#00b3f4", "Sr" = "#9b9602", "Rs" = "#fc726c"))
print(plot1)
table(Final_Afghanistan_summary$compartment)

# Resistance consumption curve at k = 0.5
avec <- seq(0, 2, 0.1)
rvec <- rep(NA, length(avec))

for (ii in seq_along(avec)) {
  parms[["a"]] <- avec[ii] / 365.25
  out <- bacteria.solve(tvec, state, parms)
  
  # Resistance proportion 
  R_final <- rowSums(out[nrow(out), Rindex + 1]) + rowSums(out[nrow(out), Rsindex + 1])
  X_final <- rowSums(out[nrow(out), Xindex + 1])
  
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

