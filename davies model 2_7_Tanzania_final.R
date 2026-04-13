#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# Modeling the Population-Level Impact of Azithromycin Mass Drug Administration on the Emergence, Persistence, and Spread of Antimicrobial-Resistant Echerchia Coli    #>>
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# This model use the demographic and social contacts patterns
# from 3.1.EpiSignalDetetion_1_OneYearContact .
# In the above mentioned script,  i will need to change
# the WHO EMRO countries by WHO Afro countries to get Tanzania
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Work environment
getwd()
setwd("C:/Disk F/4.Oxford Modelling for Global Health/Afox_Ubuntu/Afox Placement with Ben Cooper")
getwd()
# Packages
pacman::p_load(deSolve, viridis, ggplot2, tidyr, dplyr, readr)
# Section I. Routine data collected at global and state level over time
# Source: WHO:https://data.who.int/indicators/i/918081E/745F475

# Global distribution and incidence of multidrug resistant and ESBL producing Escherichia coli: an observational study of the ATLAS dataset
# Source : https://www.sciencedirect.com/science/article/pii/S2213398425002398

# Pfizer – Antimicrobial Testing Leadership and Surveillance (ATLAS)
# Source : https://www.amrindustryalliance.org/case-study/antimicrobial-testing-leadership-and-surveillance-atlas/
pacman::p_load(readr)
data <- read_csv("E.coli_resistance_C3_ALL_LATEST.csv")
print(data)
dim(data)
names(data)
names(data)[5] <- "Period"
names(data)[11] <- "Country"
names(data)[12] <- "Percentage"
names(data)
# Filter data
pacman::p_load(dplyr)
table(data$Country)
data$Country[data$Country == "United Republic of Tanzania"] <- "Tanzania"
data1 <- data %>%
  filter(Country %in% c("World", "Nigeria", "Namibia ", "Rwanda", "Sudan", "Tanzania", "Zambia", "Uganda", "Malawi")) %>%
  filter(Country %in% c("World", "Niger", "Tanzania", "Malawi"))
# Visualisation
pacman::p_load(ggplot2)
ggplot(data1, aes(x = Period, y = Percentage, color = Country, group = Country)) +
  # geom_line(linewidth = 1) +
  geom_point(size = 3) +
  # geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) +
  scale_x_continuous(breaks = seq(min(data1$Period),
    max(data1$Period),
    by = 1
  )) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20)
  ) +
  labs(
    title = "Proportion of bloodstream infection due to Escherichia coli resistant to  C3 (%)",
    subtitle = "Annual estimates by country",
    # title = "Proportion of bloodstream infection due to Escherichia coli resistant to third-generation cephalosporins (%)",
    x = "Period",
    y = "Percentage(%)",
    # y = "Rate per 100,000",
    color = "Country"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.text = element_text(size = 10)
  )
#
g_0 <- ggplot(data1, aes(x = Period, y = Percentage, color = Country, group = Country)) +
  geom_ribbon(
    data = subset(data1, Country == "World"),
    aes(
      x = Period,
      ymin = 0,
      ymax = Percentage
    ),
    inherit.aes = FALSE,
    fill = "#C2A5CF",
    # fill = "grey80",
    alpha = 0.4
  ) +
  geom_point(size = 3) +
  # geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(min(data1$Period), max(data1$Period), by = 1)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20)
  ) +
  labs(
    title = "Proportion of bloodstream infection due to Escherichia coli resistant to C3 (%)",
    subtitle = "WHO data, Official estimate updated 7 May 2025",
    x = "Period",
    y = "Percentage (%)",
    color = "Country"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
g_0
ggsave("Figure 1.Propportion of bloodstream infection due to E.Coli resistant to c3.png",
  plot = last_plot(),
  bg = "white",
  width = 10,
  height = 8,
  dpi = 300
)
# Control parameters
save.fig <- FALSE
# 1. Demographic parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.a. Age structure (101 groups: 0,1, 100+)
age_groups <- c(as.character(0:99), "100+")
# Number of age groups
n_age <- length(age_groups)
# Aging matrix
A <- n_age # 101 age groups (0,...100+)
dd <- rep(1, A) #
ageing <- t(diff(diag(dd), lag = 1) / (1 * 365.25)) # 1 year age groups
ageing <- cbind(ageing, rep(0, A)) # No ageing from last compartment
# Population structure
setwd("C:/Disk F/4.Oxford Modelling for Global Health/Placement project disk")
# Population structure 2000-2023
Population_emro_2023 <- read_csv("Population_emro_2023_1yearage.csv")
table(Population_emro_2023$Country)
head(Population_emro_2023)
Tanzania_pop <- as.data.frame(Population_emro_2023 %>%
  filter(Year == "2023") %>%
  filter(Country == "United Republic of Tanzania"))
head(Tanzania_pop)
# Population structure in thousands
Tanzania_pop_in_thousands <- as.data.frame(Tanzania_pop %>%
  mutate(
    Population_age = Population_age * 1000,
    Annual_population = Annual_population * 1000
  ))
# 1.c.Visualization of Tanzania's population structure
pacman::p_load(ggplot2, dplyr, plotly)
#
pyramyd <- ggplot(Tanzania_pop_in_thousands, aes(x = Age_Category, y = Population_age)) +
  geom_col(fill = "steelblue") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +
  labs(
    title = "Tanzania's population structure (2023)",
    x = "Age", y = "Population"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(pyramyd)
# Shorten the name
popstruc <- Tanzania_pop_in_thousands
#  Number of age groups
A <- n_age <- nrow(Tanzania_pop_in_thousands) # Fixed: use nrow instead of length
#  Births by age of mother
popbirth <- read.csv("3.U.1.Birth_1year_Afro.csv", header = TRUE)
#
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
  labs(
    title = "Birth in Tanzania (2023)",
    x = "Age", y = "Births"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
print(pyramyd_birth)
# Convert from 1000s per 1 year period to per person per day
group_durations_years <- rep(1, n_age)
#
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
  # scale_x_continuous(
  # limits = c(0, 100),
  # breaks = seq(0, 100, by = 10))+
  labs(
    title = "Mortality in Tanzania (2023)",
    x = "Age", y = "Deaths per 1000 pop"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
print(pyramyd_mort)
# Convert mortality from 1000s per 1 year period to per person per day
group_durations_years <- rep(1, n_age)
group_durations_days <- group_durations_years * 365
#
popmort[, 5] <- 1000 * popmort[, 5] / (1 * popstruc[, 5] * 365.25)
mort <- popmort[, 5] # for dynamic population
# 1.h.Contact matrix  and visualization
(m_contact_1y_Tanzania <- as.matrix(read.csv("3.U.1.contact_Tanzania_1y.csv")))
dim(m_contact_1y_Tanzania)
colSums(ageing)
#
# for (i in 1:n_age) {
#  for (j in 1:n_age) {
#    m_contact_1y_Tanzania[i, j] <- m_contact_1y_Tanzania[i, j] / 25
#  }
# }
m_contact_1y_Tanzania <- m_contact_1y_Tanzania / 5

colnames(m_contact_1y_Tanzania) <- c(as.character(0:99), "100+")
rownames(m_contact_1y_Tanzania) <- c(as.character(0:99), "100+")
m_contact_1y_Tanzania
# Visualization of my contact matrix
pacman::p_load(ggplot2, reshape2)
# data frame
# df <- melt(m_contact_1y_Tanzania )
df <- reshape2::melt(m_contact_1y_Tanzania)
colnames(df) <- c("Contactee", "Contactor", "Contacts")
# Plot contact matrix
(p <- ggplot(df, aes(x = Contactor, y = Contactee, fill = Contacts)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "Contact Matrix Heatmap for Tanzania")
)
# Back into work environment
setwd("C:/Disk F/4.Oxford Modelling for Global Health/Afox_Ubuntu/Afox Placement with Ben Cooper")
# Indices for compartments
# Xindex <- (0*n_age+1):(1*n_age)       # Uninfected, untreated
Xindex <- 1:(1 * n_age) # Uninfected, untreated
Sindex <- (1 * n_age + 1):(2 * n_age) # drug-sensitive, untreated
Rindex <- (2 * n_age + 1):(3 * n_age) # drug-resistant, untreated
Srindex <- (3 * n_age + 1):(4 * n_age) # drug-sensitive, treated
Rsindex <- (4 * n_age + 1):(5 * n_age) # drug-resistant, treated
Dindex <- (5 * n_age + 1):(6 * n_age) # Cummulative deaths
CumIncRindex <- (6 * n_age + 1):(7 * n_age) # Cummulative resistance
AMRDindex <- (7 * n_age + 1):(8 * n_age) # Cummulative resistance
# MDA  intervention..............................................................
# mda_start_times <- c(365, 79, 4380, 4745)
(mda_start_times <- (0:50) * 365.25)
(mda_duration <- 30)

# mda_active <- function(time, mda_starts, duration) {
#  any(sapply(mda_starts, function(start) {
#    time >= start && time < (start + duration)
#  }))
# }
mda_active <- function(time, mda_starts, duration) {
  any(time >= mda_starts & time < (mda_starts + duration))
}
# ODE system.....................................................................
bacteria.odes <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Extract compartments
    X <- state[Xindex] # Uninfected, untreated
    S <- state[Sindex] # drug-sensitive, untreated
    R <- state[Rindex] # drug-resistant, untreated
    Sr <- state[Srindex] # drug-sensitive, treated
    Rs <- state[Rsindex] # drug-resistant, treated
    D <- state[Dindex] # cumulative deaths
    CumIncR <- state[CumIncRindex] # cumulative resistances
    AMRD <- state[AMRDindex] # AMR related deaths
    # Total population
    N <- X + S + R + Sr + Rs
    #
    S.tot <- S + Sr # Susceptible co-colonised total
    R.tot <- R + Rs # Resistance co-colonised total
    #
    # lamda.S <- beta.S * (m_contact_1y_Tanzania %*% (S.tot/N)) #Between host transmission
    # lamda.R <- beta.R * (m_contact_1y_Tanzania %*% (R.tot/N))
    lamda.S <- beta.S * (m_contact %*% (S.tot / N))
    lamda.R <- beta.R * (m_contact %*% (R.tot / N))

    # Intervention : MDA implementation..........................................
    # is_mda <- use_mda&&mda_active(t, mda_start_times, mda_duration)  # called ONCE
    # mda_starts <- parameters[["mda_start_times"]]
    is_mda <- mda_active(t, mda_start_times, mda_duration) #

    # Option B
    # b    <- ifelse(is_mda, a + tau, tau)
    # a_t  <- b * azt
    # bc      <- ifelse(is_mda, a.C + tau, tau)
    # a.C_t   <- bc * azt
    # Option C
    b <- ifelse(is_mda, a, 0)
    a_t <- b * azt + tau
    bc <- ifelse(is_mda, a.C, 0)
    a.C_t <- bc * azt + tau
    # Mortality
    mort_eff <- mort
    # if (is_mda) mort_eff[0:4+1] <- mort[0:4+1] * (1 - theta)
    # if (use_mda=="TRUE") mort_eff[0:4+1] <- mort[0:4+1] * (1 - theta)
    if (use_mda) {
      mort_eff[1:5] <- mort[1:5] * (1 - theta)
    }
    # ...........................................................................
    # /tau<--log(1-parms.orig$a.use)      #a.use_1:daily rate of antibiotics use           :option A
    # //p_treated = parms.orig$a.use_p/1000 # parms.orig$a.use_p=ddd/1000/d                 :option B
    # //tau<-p_treated/parms.orig$d         # daily antibiotic use rate using ddd/1000/day  :option B

    # Caluclation of a_t
    # //b <- ifelse(mda_active(t, mda_start_times, mda_duration),a + (a.use.eff*tau), (a.use.eff*tau)) #:option A
    # //b <- ifelse(mda_active(t, mda_start_times, mda_duration),a + tau, tau)                          #:option B
    # //b <- ifelse(mda_active(t, mda_start_times, mda_duration), a, 0)                                #Option C: sean
    # //a_t <- b * azt                                                                                  #:option A,b
    # //a_t <- b * azt + tau                                                                           #option C: Sean

    # Caluclation of a.C_t
    # //bc <- ifelse(mda_active(t, mda_start_times, mda_duration),a.C + (a.use.eff*tau), (a.use.eff*tau))#:option A
    # //bc <- ifelse(mda_active(t, mda_start_times, mda_duration),a.C + tau, tau)                         #:option B
    # //bc <- ifelse(mda_active(t, mda_start_times, mda_duration), a.C, 0)                               #:option c
    # //a.C_t <- bc * azt                                                                                 #:option A,B
    # //a.C_t <- bc * azt + tau :                                                                        #Option C

    # Mortality reduction during MDA
    # //mort_eff <- mort  # start from baseline
    # Pulse for MDA period
    # //if (mda_active(t, mda_start_times, mda_duration)) {
    # // mort_eff[mda_targeted_ages] <- mort[mda_targeted_ages] * (1 - theta)
    # //}
    # ...........................................................................
    # Births
    births <- rep(0, n_age)
    births[1] <- sum(popbirth[, 5] * N) # for dynamic population
    # total_deaths <- sum(mort_eff * N  )  #for static population
    total_deaths <- sum(mort_eff * N + amrd_rate * (R + Rs)) # for static population
    births[1] <- total_deaths # for static population
    # Browser()
    # ...........................................................................
    # ODEs system # Here i added a_t
    dX <- births + (u.S + a_t) * S + u.R * R + u.C * (Sr + Rs) - (lamda.S + lamda.R) * X + ageing %*% X - mort_eff * X #
    #
    dS <- lamda.S * X - u.S * S - k * lamda.R * S - a_t * S + ageing %*% S - mort_eff * S #
    #
    dR <- lamda.R * X - u.R * R - k * lamda.S * R + a.C_t * (Sr + Rs) + ageing %*% R - (mort_eff + amrd_rate) * R #
    #
    dSr <- k * lamda.R * S - Sr * u.C - a.C_t * Sr + ageing %*% Sr - mort_eff * Sr #
    #
    dRs <- k * lamda.S * R - Rs * u.C - a.C_t * Rs + ageing %*% Rs - (mort_eff + amrd_rate) * Rs #

    # Counting
    dD <- mort_eff * X + mort_eff * S + mort_eff * R + mort_eff * Sr + mort_eff * Rs + amrd_rate * (R + Rs) # Cummulative                    #
    dCumIncR <- lamda.R * X + k * lamda.R * S # Cummulative incidence of New resistant infections :a.Uninfected person get resistant strains #
    dAMRD <- amrd_rate * (R + Rs) # Cummulative : Mortality attributable to AMR
    # dCumRs <- k*lamba.S*R                      # Cummulative
    list(c(dX, dS, dR, dSr, dRs, dD, dCumIncR, dAMRD))
  })
}

bacteria.solve <- function(t, state, parameters) {
  parameters[["beta.R"]] <- parameters[["beta.S"]] * (1 - parameters[["c"]])
  out_Tanzania <- as.data.frame(ode(state, t, bacteria.odes, parameters))
  return(out_Tanzania)
}
# Parameters---------------------------------------------------------------------
parms.orig <- list(
  # Pathogen parameters..........................................................#Recommended,Range
  beta.S = 0.03, # 5,       # Transmission of sensitive   : (β = 5 month−1)    :[0.04-0.08][0.03-0.10
  u.S = 1, # Clearance sensitive (natural):(u = 1 month−1)          :[0.01-0.05]
  u.R = 1, # Clearance resistant (natural)     : (u = 1 month−1)    :[0.008-0.02] lower than susceptible?
  u.C = 1, # Clearance co-colonised (natural)  : (u = 1 month−1)    :[0.01-0.04]
  k = 0.5, # The efficiency of co-colonisation : (k = 0.25,0.5,1.0)
  c = 0.20, # The fitness cost : (c = 0-10%)
  # MDA Azithromycin parameters..................................................
  a = 0.16, # Clearance sensitive (drug-induced)                     :[0.05-0.10]
  a.C = 0.16, # Clearance co-colonised (drug-induced)                  :[0.03-0.08]
  mda_cycle = 365, # MDA frequency
  mda_duration = 30, # MDA duration
  mda_cov = 0.6, # MDA coverage
  # theta   =  0,           # (for static population) under-five mortality reduction
  theta = 0.13, # (for dynamic population) under-five mortality reduction due to MDA
  # Baseline antibiotic use parameters------------------------------------------
  a.use = 0.06, # Antibiotic use in % (~0.01–0.05:routine use,0.05 – 0.20: High use communitie,0.20 – 0.80 (short period) :MDA  )
  a.use.eff = 0.05, # Antibiotic effect: 0.005−0.05 (assumed)
  #---DDD/1000/day---------#
  a.use_p = 23.1, # ddd/1000/day  in general population in 2018 CrI: 25[23.1-26.9]
  a.use_c = 36.9, # Antibiotic use in % in under five in 2018 CrI:36.9[31.9-42.4]
  d = 7, # Duration of antibiotics treatment in days  7[5-10]
  # Country specif social contact patterns---------------------------------------
  m_contact = m_contact_1y_Tanzania, # Social contacts per day
  # Others parameters------------------------------------------------------------
  kappa = 0, # 0.05                #Proportion that develop/select resistance (Assumed),
  # amrd_rate =  0           # 27.3/(100000/365), #AMR related mortality per person per day,
  amrd_rate = (27.3 / 100000 / 365) / (0.15 * 0.9) # Sean:0.9(colonization),: 0.15(Prevalence)
)
# MDA rate calculation: Exponential decay
(parms.orig$r_mda <- -log(1 - parms.orig$mda_cov) / parms.orig$mda_duration)
parms.orig
parms.orig[1:4]
# Convert daily
# [1:5]
parms.orig[2:4] <- lapply(parms.orig[2:4], function(x) x * 12 / 365.25) # Daily
parms.orig[1:4]
# Adjusted clearance rates: Esther et al
# parms.orig["u.S"]= 0.0098
# parms.orig["u.R"]= 0.0098
# parms.orig["u.C"]= 0.0098
# βR=−ln(1−p)/D where p is the daily probability of acquisition
(prob <- 3.92 / 100) # 3.92% from Rebecca et al
(ar <- log(1 / (1 - prob)))
# [1] 0.03998901
# MDA intervention into the ODE: Pre-computation
# a.Exposure to antibiotic
# //tau<-parms.orig$a.use*parms.orig$a.use.eff     #Option
# //tau<--log(1-parms.orig$a.use)       #a.use_1:daily rate of antibiotics use       :option A
p_treated <- parms.orig$a.use_p / 1000 # parms.orig$a.use_p=ddd/1000/d                 :option B
tau <- p_treated / parms.orig$d # daily antibiotic use rate using ddd/1000/day  :option B
# b.Age and coverage
mda_targeted_ages <- 1:5 #     Index of Targeted ages
azt <- rep(0, n_age) #     Initialize azt vector
azt[mda_targeted_ages] <- mda_cov # 0.8 #MDA coverage
# c.Derrived  parameters
parms.orig$tau <- tau
parms.orig$azt <- azt
parms.orig$use_mda <- TRUE
parms.orig$mda_start_times <- mda_start_times #
parms <- parms.orig
parms
# ~~~~~~~~~No-annual MDA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
parms_noMDA <- parms # No annual MDA  parameters
parms_noMDA$r_mda <- 0
parms_noMDA$a <- 0
parms_noMDA$a.C <- 0
parms_noMDA$theta <- 0
parms_noMDA$rho <- 0
parms_noMDA$use_mda <- FALSE
parms_noMDA
# parms_noMDA<-parms.orig

# Inside bacteria.odes:
# is_mda <- use_mda && mda_active(t, mda_start_times, mda_duration)
# Initial conditions------------------------------------------------------------
names(Tanzania_pop_in_thousands)
head(Tanzania_pop_in_thousands)
dim(Tanzania_pop_in_thousands)
initP <- Tanzania_pop_in_thousands[, 5]
sum(initP)
# initS <- Tanzania_pop_in_thousands[,5]
initX <- 0.95 * Tanzania_pop_in_thousands[, 5]
initS <- 0.025 * Tanzania_pop_in_thousands[, 5]
initR <- 0.025 * Tanzania_pop_in_thousands[, 5]
initSr <- 0 * Tanzania_pop_in_thousands[, 5]
initRs <- 0 * Tanzania_pop_in_thousands[, 5]
initD <- 0 * Tanzania_pop_in_thousands[, 5]
initCumIncR <- 0 * Tanzania_pop_in_thousands[, 5]
initAMRD <- 0 * Tanzania_pop_in_thousands[, 5]
# Combine initial states
state.orig <- c(initX, initS, initR, initSr, initRs, initD, initCumIncR, initAMRD)
# States
tvec_1_a <- seq(0, 1 * 365.25, 1) # 1 Year-MDA
tvec_1_b <- seq(0, 1 * 365.25, 1) # Pre-MDA
tvec_5_a <- seq(0, 5 * 365.25, 1) # 5 years of MDA
tvec_5_b <- seq(0, 5 * 365.25, 1) # 5 years of MDA
tvec_10_a <- seq(0, 10 * 365.25, 1) # 10 years of MDA
tvec_10_b <- seq(0, 10 * 365.25, 1) # 10 years No MDA, i will need a = 10
# Run model

state <- state.orig
start <- Sys.time()
# Baseline: Equilibrium at 70 years
tvec_0_b <- seq(0, 110 * 365.25, 1)
parms_noMDA$mda_start_times <- numeric(0) # (mda_start_times<-(0:70)*365.25)
start <- Sys.time()
out_0_b_Tanzania <- bacteria.solve(tvec_0_b, state, parms_noMDA)
end <- Sys.time()
print(end - start)
last10 <- tail(out_0_b_Tanzania, 3650)
# Total population stability
total_pop <- rowSums(last10[, c(Xindex, Sindex, Rindex, Srindex, Rsindex)])
# plot(total_pop, type="l")
# Max difference
max(abs(diff(total_pop)))
state <- as.numeric(out_0_b_Tanzania[nrow(out_0_b_Tanzania), -1])
# print(state )
# Modelling different scenarios
parms
parms_noMDA
# Simulations#--------------------------------------------------------------------
start <- Sys.time()
# ~~~~~~~~~No MDA and annual MDA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
(parms$mda_start_times <- (0:1) * 365.25) # 1 year MDA
out_1_a_Tanzania <- bacteria.solve(tvec_1_a, state, parms) # Annual
parms_noMDA$mda_start_times <- numeric(0) # (mda_start_times<-(0:1)*365.25) # 1 year No MDA
out_1_b_Tanzania <- bacteria.solve(tvec_1_b, state, parms_noMDA) # No annual
#
(parms$mda_start_times <- (0:50) * 365.25) # MDA
out_5_a_Tanzania <- bacteria.solve(tvec_5_a, state, parms)
#
parms_noMDA$mda_start_times <- numeric(0) # (mda_start_times<-(0:50)*365.25)             # No MDA
out_5_b_Tanzania <- bacteria.solve(tvec_5_b, state, parms_noMDA)
#
(parms$mda_start_times <- (0:50) * 365.25) #  MDA
out_10_a_Tanzania <- bacteria.solve(tvec_10_a, state, parms)

parms_noMDA$mda_start_times <- numeric(0) # (mda_start_times<-(0:50)*365.25)# No MDA
out_10_b_Tanzania <- bacteria.solve(tvec_10_b, state, parms_noMDA)

# ~~~~~~~~~Bi-annual MDA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parms$mda_start_times <- (0:50) * (365.25 / 2)
out_1_c_Tanzania <- bacteria.solve(tvec_1_a, state, parms) # Bi-annual(same parameters excepts mda_start)
parms$mda_start_times <- (0:50) * (365.25 / 2)
out_5_c_Tanzania <- bacteria.solve(tvec_5_a, state, parms)
parms$mda_start_times <- (0:50) * (365.25 / 2)
out_10_c_Tanzania <- bacteria.solve(tvec_10_a, state, parms)

# MDA stops_ before the end of simulation
parms$mda_start_times <- (0:5) * (365.25 / 1)
out_10_a_5_Tanzania <- bacteria.solve(tvec_10_a, state, parms)
parms$mda_start_times <- (0:6) * (365.25 / 1)
out_10_a_6_Tanzania <- bacteria.solve(tvec_10_a, state, parms)
parms$mda_start_times <- (0:7) * (365.25 / 1)
out_10_a_7_Tanzania <- bacteria.solve(tvec_10_a, state, parms)
# MDA stops_bi annual
parms$mda_start_times <- (0:5) * (365.25 / 2)
out_10_c_5_Tanzania <- bacteria.solve(tvec_10_a, state, parms)
parms$mda_start_times <- (0:6) * (365.25 / 2)
out_10_c_6_Tanzania <- bacteria.solve(tvec_10_a, state, parms)
parms$mda_start_times <- (0:7) * (365.25 / 2)
out_10_c_7_Tanzania <- bacteria.solve(tvec_10_a, state, parms)

# 50 years with 10 years of MDA
tvec_50_a <- seq(0, 20 * 365.25, 1)
tvec_50_b <- seq(0, 20 * 365.25, 1)
tvec_50_c <- seq(0, 20 * 365.25, 1)
parms_noMDA$mda_start_times <- numeric(0) # empty = no MDA
out_50_b_Tanzania <- bacteria.solve(tvec_50_b, state, parms_noMDA)
parms$mda_start_times <- (0:10) * (365.25 / 1) # annual
out_50_a_Tanzania <- bacteria.solve(tvec_50_a, state, parms)
parms$mda_start_times <- (0:10) * (365.25 / 2) # bi annual
out_50_c_Tanzania <- bacteria.solve(tvec_50_c, state, parms)
end <- Sys.time()
end - start
# Annual
results_1_a_Tanzania <- as.data.frame(out_1_a_Tanzania)
results_1_b_Tanzania <- as.data.frame(out_1_b_Tanzania)
results_1_c_Tanzania <- as.data.frame(out_1_c_Tanzania)
results_5_a_Tanzania <- as.data.frame(out_5_a_Tanzania)
results_5_b_Tanzania <- as.data.frame(out_5_b_Tanzania)
results_5_c_Tanzania <- as.data.frame(out_5_c_Tanzania)
# Bi annual
results_10_a_Tanzania <- as.data.frame(out_10_a_Tanzania)
results_10_b_Tanzania <- as.data.frame(out_10_b_Tanzania)
results_10_c_Tanzania <- as.data.frame(out_10_c_Tanzania)
# 10 years ,Stops 5-6-7
results_10_a_5_Tanzania <- as.data.frame(out_10_a_5_Tanzania)
results_10_a_6_Tanzania <- as.data.frame(out_10_a_6_Tanzania)
results_10_a_7_Tanzania <- as.data.frame(out_10_a_7_Tanzania)
results_10_c_5_Tanzania <- as.data.frame(out_10_c_5_Tanzania)
results_10_c_6_Tanzania <- as.data.frame(out_10_c_6_Tanzania)
results_10_c_7_Tanzania <- as.data.frame(out_10_c_7_Tanzania)
# 50 years,stop 10
results_50_a_Tanzania <- as.data.frame(out_50_a_Tanzania)
results_50_b_Tanzania <- as.data.frame(out_50_b_Tanzania)
results_50_c_Tanzania <- as.data.frame(out_50_c_Tanzania)

# Column names
compartment_names <- c("X", "S", "R", "Sr", "Rs", "D", "CumIncR", "AMRD")
col_names <- c("time")
for (comp in compartment_names) {
  for (age in age_groups) {
    col_names <- c(col_names, paste0(comp, "_", age))
  }
}
colnames(results_1_a_Tanzania) <- col_names
colnames(results_1_b_Tanzania) <- col_names
colnames(results_1_c_Tanzania) <- col_names
colnames(results_5_a_Tanzania) <- col_names
colnames(results_5_b_Tanzania) <- col_names
colnames(results_5_c_Tanzania) <- col_names

colnames(results_10_a_Tanzania) <- col_names
colnames(results_10_b_Tanzania) <- col_names
colnames(results_10_c_Tanzania) <- col_names
# 10 years ,stop 5-6-7
colnames(results_10_a_5_Tanzania) <- col_names
colnames(results_10_a_6_Tanzania) <- col_names
colnames(results_10_a_7_Tanzania) <- col_names
colnames(results_10_c_5_Tanzania) <- col_names
colnames(results_10_c_6_Tanzania) <- col_names
colnames(results_10_c_7_Tanzania) <- col_names
# 50 years ,stop 10
colnames(results_50_a_Tanzania) <- col_names
colnames(results_50_b_Tanzania) <- col_names
colnames(results_50_c_Tanzania) <- col_names

names(results_1_a_Tanzania)
names(results_1_b_Tanzania)
names(results_1_c_Tanzania)
names(results_5_a_Tanzania)
names(results_5_b_Tanzania)
names(results_5_c_Tanzania)
names(results_10_a_Tanzania)
names(results_10_b_Tanzania)
names(results_10_c_Tanzania)

names(results_10_a_5_Tanzania)
names(results_10_a_6_Tanzania)
names(results_10_a_7_Tanzania)
names(results_10_c_5_Tanzania)
names(results_10_c_6_Tanzania)
names(results_10_c_7_Tanzania)
# 50 Years
names(results_50_a_Tanzania)
names(results_50_b_Tanzania)
names(results_50_c_Tanzania)

# Baseline mortality
popmort <- as.data.frame(mortality %>%
  filter(Country == "United Republic of Tanzania") %>%
  filter(Year == 2023))
pyramyd_mort <- ggplot(popmort, aes(x = Age, y = Percentage)) +
  geom_col(fill = "red") +
  # scale_x_continuous(
  # limits = c(0, 100),
  # breaks = seq(0, 100, by = 10))+
  labs(
    title = "Mortality in Tanzania (2023)",
    x = "Age", y = "Deaths per 1000 pop"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
print(pyramyd_mort)
# age groups : No MDA
deaths_1_b <- (results_1_b_Tanzania[, Dindex + 1])
colnames(deaths_1_b) <- c(as.character(0:99), "100+")
deaths_1_b_last <- as.numeric(tail(deaths_1_b, 1))
deaths_1_b_prop <- round(deaths_1_b_last, 0) * 1000 / sum(round(deaths_1_b_last, 0))
barplot(deaths_1_b_prop, col = "blue")

# age groups :MDA
deaths_1_a_interv <- (results_1_a_Tanzania[, Dindex + 1])
colnames(deaths_1_a_interv) <- c(as.character(0:99), "100+")
deaths_1_a_last <- as.numeric(tail(deaths_1_a_interv, 1))
deaths_prop_1_a_interv <- round(deaths_1_a_last, 0) * 1000 / sum(round(deaths_1_a_last, 0))
barplot(deaths_prop_1_a_interv, col = "green")

death_data <- data.frame("Age" = 0:100, "NoMDA" = as.vector(round(deaths_1_b_last, 0)), "MDA" = as.vector(round(deaths_1_a_last, 0)))
# death_data <- data.frame("Age" = 0:100, "Baseline"=as.vector(round(deaths_1_b_prop, 0)), "MDA" = as.vector(round(deaths_prop_1_a_interv, 0)))
library(reshape2)
df_long <- melt(death_data, id.var = "Age")
library(ggplot2)
ggplot(df_long, aes(x = Age, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge")

#
deaths_1_c_interv <- (results_1_c_Tanzania[, Dindex + 1])
colnames(deaths_1_c_interv) <- c(as.character(0:99), "100+")
deaths_1_c_last <- as.numeric(tail(deaths_1_c_interv, 1))
deaths_prop_1_c_interv <- round(deaths_1_a_last, 0) * 1000 / sum(round(deaths_1_c_last, 0))
barplot(deaths_prop_1_c_interv, col = "purple")
#
# death_data <- data.frame("Age" = 0:100, "NoMDA"=as.vector(round(deaths_1_b_last, 0)), "MDA" = as.vector(round(deaths_1_a_last, 0)),"BiMDA" = as.vector(round(deaths_1_c_last, 0)))
death_data <- data.frame("Age" = 0:100, "Baseline" = as.vector(round(popmort[, 5], 0)), "MDA" = as.vector(round(deaths_prop_1_a_interv, 0)), "BiMDA" = as.vector(round(deaths_prop_1_c_interv, 0)))
library(reshape2)
head(death_data)
df_long <- melt(death_data, id.var = "Age")
colnames(df_long)[2:3] <- c("Scenario", "Mortality")
head(df_long)
library(ggplot2)
ggplot(df_long |>
  filter(Age < 15), aes(x = Age, y = Mortality, fill = Scenario)) +
  geom_bar(stat = "identity", position = "dodge") +
  # scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 0.1, suffix = " M")) +
  labs(y = "Deaths per 1000 population")
#
library(dplyr)
library(ggplot2)

df_plot <- df_long |>
  mutate(
    Age_group = ifelse(Age >= 101, "101+", as.character(Age)),
    Age_group = factor(Age_group, levels = c(as.character(0:100), "101+"))
  )

ggplot(
  df_plot |> filter(Age < 101 | Age >= 101),
  aes(x = Age_group, y = Mortality, fill = Scenario)
) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Age", y = "Deaths per 1000 population") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

# 6.Verification of the population over time
# Total population at each time point
names(results_1_b_Tanzania)
dim(results_1_b_Tanzania)
str(results_1_b_Tanzania)
head(results_1_b_Tanzania)
# total_pop<- rowSums(results_1_b_Tanzania[, 2:(5*n_age)])  # This will exclude time column
total_pop <- rowSums(results_1_b_Tanzania[, c(Xindex, Sindex, Rindex, Srindex, Rsindex)])
# Check conservation: Here the population is conserved if the birth=deaths
#                   : The population is not conserved as births and deaths  are not equal(we used realistic country and age specific births and deaths)
# Initial population
(initial_pop <- sum(Tanzania_pop_in_thousands[, 5]))
# Time
tim <- results_1_b_Tanzania[, 1]
# Data frame for plotting
P_t <- as.data.frame(cbind(tim, total_pop, initial_pop))
# Population variations
(change <- round(((tail(P_t$total_pop, 1) - initial_pop) / initial_pop) * 100, 2)) # Checks
# Packages for visualization
pacman::p_load(ggplot2, ggtext)
# Visualization of the initial population and total population over time
# Plot
# Population dynamic (with annotation)
population <- ggplot(P_t, aes(x = tim, y = total_pop)) +
  geom_hline(aes(yintercept = initial_pop), linetype = "dashed", color = "red") +
  geom_line(color = "blue") +
  labs(
    title = "Total population over time in Tanzania",
    subtitle = "<span style='color:red;'>Red dashed</span>: Initial population, <span style='color:blue;'>Blue line</span>: Total population over time",
    x = "Time (days)",
    y = "Total population (M)"
  ) +
  scale_y_continuous(
    labels = scales::label_number(scale = 1e-6, accuracy = 0.1, suffix = " M")
  ) +
  annotate(
    "text",
    x = max(P_t$tim) * 0.5,
    y = max(P_t$initial_pop) * 0.1,
    label = paste0("Relative change: ", change, "%"),
    size = 3,
    color = "darkgreen"
  ) +
  theme_minimal() +
  theme(
    plot.subtitle = ggtext::element_markdown(size = 11)
  )
print(population)
# Plot model
cols <- viridis(5)
par(mfrow = c(1, 1))
# Total across age groups: Here i will use Index#1e6
X_total <- rowSums(out_1_b_Tanzania[, Xindex + 1]) # Here i added +1 because first column is time
S_total <- rowSums(out_1_b_Tanzania[, Sindex + 1])
R_total <- rowSums(out_1_b_Tanzania[, Rindex + 1]) #+ rowSums(out_1_b_Tanzania[, Rsindex + 1])
Rs_total <- rowSums(out_1_b_Tanzania[, Rsindex + 1])
Sr_total <- rowSums(out_1_b_Tanzania[, Srindex + 1])
D_cum <- rowSums(out_1_b_Tanzania[, Dindex + 1])
D_daily <- c(0, diff(D_cum))
plot(D_cum)
plot(R_total)
plot(D_daily)
print(D_cum)
# D_total <- rowSums(out_1_b_Tanzania[, Dindex + 1])
prevalence <- round((R_total + Rs_total) * 100 / (S_total + R_total + Rs_total + Sr_total), 1)
summary(prevalence)
prevalence <- round(rowSums(out_1_b_Tanzania[, c(Rindex, Rsindex) + 1]) * 100 / rowSums(out_1_b_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
summary(prevalence)
mortality <- round(D_daily * 100000 / (X_total + S_total + R_total + Rs_total + Sr_total), 1)
summary(mortality)
# Simulation time in years
length(R_total) / 365.25
# Equilibrium checks
tail(R_total)
Tanzania_df_no_resisitance_mortality <- data.frame(
  time = out_1_b_Tanzania[, 1],
  S = S_total,
  X = X_total,
  R = R_total,
  Rs = Rs_total,
  Sr = Sr_total,
  D = D_daily
)
Tanzania_df_with_resisitance_mortality <- data.frame(
  time = out_1_b_Tanzania[, 1],
  S = S_total,
  X = X_total,
  R = R_total,
  Rs = Rs_total,
  Sr = Sr_total,
  D = D_daily,
  prevalence,
  mortality
)
head(Tanzania_df_with_resisitance_mortality)
baseline_prevalence <- ggplot(Tanzania_df_with_resisitance_mortality, aes(x = time, y = prevalence)) +
  geom_line(linewidth = 0.5) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10)
  ) +
  scale_x_continuous(
    breaks = seq(min(Tanzania_df_with_resisitance_mortality$time), max(Tanzania_df_with_resisitance_mortality$time), by = 3650)
  ) +
  labs(
    title = "Baseline prevalence of E.Coli infections resistant to Azithromycin",
    subtitle = "Antibiotic use in Tanzania : 23.1 DDD per 1000 inhabitants per day",
    x = "Time",
    y = "Resistance (%)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
print(baseline_prevalence)
Tanzania_df_with_resisitance_mortality$prevalence

pacman::p_load(tidyr, dplyr, ggplot2, scales)
df_b <- Tanzania_df_no_resisitance_mortality %>%
  pivot_longer(
    cols = -time,
    names_to = "Compartment",
    values_to = "Population"
  )
table(df_b$Compartment)
Resitance <- df_b %>%
  filter(Compartment %in% c("R", "Rs")) %>%
  ggplot(aes(x = time, y = Population, color = Compartment)) +
  # geom_line(linewidth = 1.1) +
  geom_point(size = 0.4) +
  scale_y_continuous(
    limits = c(0, 60000000),
    breaks = seq(0, 60000000, by = 2000000),
    labels = label_number(scale = 1e-6, suffix = " M")
  ) +
  labs(
    x = "Time",
    y = "Resistance",
    title = "Resistance over time"
  ) +
  # theme_minimal()
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    plot.title = element_text(hjust = 0, size = 14)
  )
Resitance
# Here i want to extract the last year to see if i reach equilibrium
df_b %>%
  filter(Compartment %in% c("R", "Rs")) %>%
  slice_tail(n = 10) %>% # the last  ten observations
  ggplot(aes(x = time, y = Population, color = Compartment)) +
  # geom_line(linewidth = 1.1) +
  geom_point(size = 4) +
  scale_y_continuous(
    limits = c(
      0,
      max(
        df_b %>%
          filter(Compartment %in% c("R", "Rs")) %>%
          slice_tail(n = 10) %>%
          pull(Population),
        na.rm = TRUE
      )
    ),
    labels = label_number(scale = 1e-6, suffix = " M")
  ) +
  labs(
    x = "Time",
    y = "Resistance",
    title = "Resistance over final period : Equilibrium check"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14)
  )
#
pacman::p_load(dplyr, scales, ggplot2)
# data
df_last30 <- df_b %>%
  filter(Compartment %in% c("R", "Rs")) %>%
  slice_tail(n = 10)
# plot
b <- ggplot(df_last30, aes(x = time, y = Population, color = Compartment)) +
  # Points
  geom_point(size = 4) +
  # Vertical labels below points
  geom_text(
    aes(label = label_number(scale = 1e-6, suffix = " M")(Population)),
    angle = 90,
    vjust = 1.5, # pushes text below the point
    hjust = 1.5,
    size = 3.5,
    show.legend = FALSE
  ) +
  # Y scale formatting
  scale_y_continuous(
    limits = c(
      0,
      max(df_last30$Population, na.rm = TRUE)
    ),
    labels = label_number(scale = 1e-6, suffix = " M")
  ) +
  labs(
    x = "Time",
    y = "Resistance",
    title = "Resistance over final period : Equilibrium check"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14)
  )
b
pacman::p_load(gridExtra)
# c<-grid.arrange(incidence,b,ncol=2)
# ggsave("Resistance_equilibrium_check.png",plot=c,width=10,height = 5,dpi=300)
#
dynamics <- ggplot(df_b |> dplyr::filter(Compartment != "D"), aes(x = time, y = Population, color = Compartment)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(
    limits = c(0, 50000000),
    breaks = seq(0, 50000000, by = 10000000),
    labels = label_number(scale = 1e-6, suffix = " M")
  ) +
  labs(
    x = "Time",
    y = "Population",
    title = "E.Coli Dynamics in Tanzania"
  ) +
  # theme_minimal()
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14)
  )
print(dynamics)
ggsave("Figure 1.Population and E.coli dynamics.png", plot = grid.arrange(population, dynamics, ncol = 2), width = 12, height = 5, dpi = 300)

# Keys metrics :Prevalence ,incidence,resistant cases ,mortality in under 5

# A.Prevalence  # We considered on R as in routine surveillance, resistant cases are reported
R_1_a <- round((rowSums(out_1_a_Tanzania[, c(Rindex, Rsindex) + 1])) * 100 / rowSums(out_1_a_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
R_1_b <- round((rowSums(out_1_b_Tanzania[, c(Rindex, Rsindex) + 1])) * 100 / rowSums(out_1_b_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
R_1_c <- round((rowSums(out_1_c_Tanzania[, c(Rindex, Rsindex) + 1])) * 100 / rowSums(out_1_c_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
R_5_a <- round((rowSums(out_5_a_Tanzania[, c(Rindex, Rsindex) + 1])) * 100 / rowSums(out_5_a_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
R_5_b <- round((rowSums(out_5_b_Tanzania[, c(Rindex, Rsindex) + 1])) * 100 / rowSums(out_5_b_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
R_5_c <- round((rowSums(out_5_c_Tanzania[, c(Rindex, Rsindex) + 1])) * 100 / rowSums(out_5_c_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
R_10_a <- round((rowSums(out_10_a_Tanzania[, c(Rindex, Rsindex) + 1])) * 100 / rowSums(out_10_a_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
R_10_b <- round((rowSums(out_10_b_Tanzania[, c(Rindex, Rsindex) + 1])) * 100 / rowSums(out_10_b_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
R_10_c <- round((rowSums(out_10_c_Tanzania[, c(Rindex, Rsindex) + 1])) * 100 / rowSums(out_10_c_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
R_50_a <- round((rowSums(out_50_a_Tanzania[, c(Rindex, Rsindex) + 1])) * 100 / rowSums(out_50_a_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
R_50_b <- round((rowSums(out_50_b_Tanzania[, c(Rindex, Rsindex) + 1])) * 100 / rowSums(out_50_b_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
R_50_c <- round((rowSums(out_50_c_Tanzania[, c(Rindex, Rsindex) + 1])) * 100 / rowSums(out_50_c_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
# 1-year
df_1 <- data.frame(Years = 1, Days = seq_along(R_1_a) - 1, Once = R_1_a, Baseline = R_1_b, Twice = R_1_c)
# 5-year
df_5 <- data.frame(Years = 5, Days = seq_along(R_5_a) - 1, Once = R_5_a, Baseline = R_5_b, Twice = R_5_c)
# 10-year
df_10 <- data.frame(Years = 10, Days = seq_along(R_10_a) - 1, Once = R_10_a, Baseline = R_10_b, Twice = R_10_c)
# 50 Years
df_50 <- data.frame(Years = 20, Days = seq_along(R_50_a) - 1, Once = R_50_a, Baseline = R_50_b, Twice = R_50_c)
# Combine all
df_all <- rbind(df_1, df_5, df_10, df_50)
colnames(df_all)
head(df_all)
colnames(df_all)[3:5] <- c("MDA", "No-MDA", "Bi-MDA")
library(tidyr)
#
Tanzania_df_all_long <- df_all %>%
  pivot_longer(
    cols = all_of(c("MDA", "No-MDA", "Bi-MDA")),
    names_to = "Strategy",
    values_to = "Resistance"
  )
head(Tanzania_df_all_long)
Tanzania_df_all_long <- Tanzania_df_all_long |>
  mutate(Horizon = paste0(Years, "Y"))
head(Tanzania_df_all_long)
summary(Tanzania_df_all_long$Resistance)

Tanzania_df_all_long <- Tanzania_df_all_long %>%
  mutate(Horizon = factor(Horizon, levels = c("1Y", "5Y", "10Y", "20Y")))

Tanzania_df_all_long$Resistance[2]
head(Tanzania_df_all_long)
baseline_value <- Tanzania_df_all_long %>%
  filter(Strategy == "No-MDA") %>%
  slice(1) %>%
  pull(Resistance)
colnames(Tanzania_df_all_long)[3] <- "Policy"

plot_1 <- ggplot(Tanzania_df_all_long |>
  filter(Horizon != "20Y"), aes(x = Days, y = Resistance, color = Policy)) +
  # geom_point(size=0.5)+
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_hline(
    yintercept = Tanzania_df_all_long$Resistance[2],
    color      = "gray60",
    linetype   = "dashed",
    linewidth  = 0.8
  ) +
  scale_x_continuous(
    breaks = seq(0, 5 * 3650, by = 365),
    labels = seq(0, 50) #
  ) +
  scale_y_continuous(
    limits = c(10, 35),
    breaks = seq(10, 35, by = 10), labels = seq(10, 35, by = 10)
  ) +
  scale_color_manual(values = c(
    "No-MDA" = "black",
    "MDA" = "#1f77b4", # deep blue
    "Bi-MDA" = "#d62728" # muted red
  )) +
  facet_wrap(~Horizon, scales = "free_x", ncol = 4) +
  labs(x = "Time (Years)", y = "Resistance (%)") +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  )

print(plot_1)
colnames(Tanzania_df_all_long)[3] <- "Strategy"
Tanzania_df_all_long$Resistance[2]
Tanzania_df_all_long_1 <- Tanzania_df_all_long %>%
  filter(Strategy != "No-MDA")
print(Tanzania_df_all_long_1)
plot_2 <- ggplot(Tanzania_df_all_long, aes(x = Strategy, y = Resistance, color = Strategy)) +
  geom_boxplot(color = "black", fill = "skyblue") +
  geom_hline(
    yintercept = Tanzania_df_all_long$Resistance[2], # horizontal line at 15%
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) + # y-axis 0–100 with steps of 10
  facet_wrap(~Horizon, scales = "free_x", ncol = 4)

print(plot_1)
print(plot_2)
pacman::p_load(gridExtra)
grid.arrange(plot_1, plot_2, ncol = 2)
library(tidyplots)
pacman::p_load(tidyplots, pak)
# install.packages("pak")
# pak::pak("jbengler/tidyplots")
#
Tanzania_df_all_long <- Tanzania_df_all_long %>%
  mutate(Horizon = paste0(Years, "Y"))
colnames(Tanzania_df_all_long)[3] <- "Strategy"
Tanzania_df_all_long <- Tanzania_df_all_long %>%
  mutate(Horizon = factor(Horizon, levels = c("1Y", "5Y", "10Y", "20Y")))
Tanzania_df_all_long
comp_1 <- Tanzania_df_all_long |>
  tidyplot(x = Strategy, y = Resistance, color = Strategy) |>
  adjust_size(width = 55, height = 48) |>
  add_boxplot() |> #
  add_test_pvalue(ref.group = 1) |>
  split_plot(by = Horizon)
print(comp_1)

comp_1_1 <- Tanzania_df_all_long |>
  tidyplot(x = Strategy, y = Resistance, color = Strategy) |>
  adjust_size(width = 55, height = 48) |>
  add_boxplot(outlier.alpha = 0.3) |>
  add_line(stat = "summary", fun = mean, linewidth = 1.2, color = "black") |>
  add_test_pvalue(ref.group = 1) |>
  split_plot(by = Horizon)
print(comp_1_1)
comp_1_1 <- Tanzania_df_all_long |>
  filter(Horizon != "20Y") |>
  tidyplot(x = Days, y = Resistance, color = Strategy) |>
  add_mean_line() |>
  # add_mean_dot() |>
  split_plot(by = Horizon)
print(comp_1_1)

head(Tanzania_df_all_long)
Tanzania_df_all_long_1 <- Tanzania_df_all_long |>
  filter(Horizon != "20Y")
#
plot_5 <- split_plot(
  Tanzania_df_all_long |>
    tidyplot(x = Days, y = Resistance, color = Strategy) |>
    add_areastack_absolute(),
  by = Horizon
)
print(plot_5)
head(Tanzania_df_all_long)
table(Tanzania_df_all_long$Years)

Tanzania_df_all_long |>
  dplyr::filter(Years %in% c(1, 5, 10, 20)) |>
  tidyplot(y = Resistance, color = Strategy) |>
  add_donut() |>
  adjust_size(width = 25, height = 25) |>
  split_plot(by = Years)
ggsave(
  filename = "Figure 5.Resistance_split_plot.png",
  plot = plot_5,
  width = 12,
  height = 6,
  dpi = 300
)
# Without baseline
# pacman::p_load(tidyplot,tidyverse)
# comp_2<-df_all_long_1 |>
#  tidyplot(x = Scenario, y = Resistance, color = Scenario) |>
#  adjust_size(width = 55, height = 48) |>
#  add_boxplot() |>                #
#  add_test_pvalue(ref.group = 1) |>
#  split_plot(by = Years)
# print(comp_2)

# B.Colonization level
# 1 year
c_1_a <- round(rowSums(out_1_a_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]) * 100 / rowSums(out_1_a_Tanzania[, c(Xindex, Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
c_1_b <- round(rowSums(out_1_b_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]) * 100 / rowSums(out_1_b_Tanzania[, c(Xindex, Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
c_1_c <- round(rowSums(out_1_c_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]) * 100 / rowSums(out_1_c_Tanzania[, c(Xindex, Sindex, Rsindex, Srindex, Rindex) + 1]), 1)

# 5 years
c_5_a <- round(rowSums(out_5_a_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]) * 100 / rowSums(out_5_a_Tanzania[, c(Xindex, Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
c_5_b <- round(rowSums(out_5_b_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]) * 100 / rowSums(out_5_b_Tanzania[, c(Xindex, Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
c_5_c <- round(rowSums(out_5_c_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]) * 100 / rowSums(out_5_c_Tanzania[, c(Xindex, Sindex, Rsindex, Srindex, Rindex) + 1]), 1)

# 10 years
c_10_a <- round(rowSums(out_10_a_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]) * 100 / rowSums(out_10_a_Tanzania[, c(Xindex, Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
c_10_b <- round(rowSums(out_10_b_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]) * 100 / rowSums(out_10_b_Tanzania[, c(Xindex, Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
c_10_c <- round(rowSums(out_10_c_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]) * 100 / rowSums(out_10_c_Tanzania[, c(Xindex, Sindex, Rsindex, Srindex, Rindex) + 1]), 1)

# 50 years
c_50_a <- round(rowSums(out_50_a_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]) * 100 / rowSums(out_50_a_Tanzania[, c(Xindex, Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
c_50_b <- round(rowSums(out_50_b_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]) * 100 / rowSums(out_50_b_Tanzania[, c(Xindex, Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
c_50_c <- round(rowSums(out_50_c_Tanzania[, c(Sindex, Rsindex, Srindex, Rindex) + 1]) * 100 / rowSums(out_50_c_Tanzania[, c(Xindex, Sindex, Rsindex, Srindex, Rindex) + 1]), 1)
# Data frame
# 1-year
df_c_1 <- data.frame(Years = 1, Days = seq_along(c_1_a) - 1, Once = c_1_a, Baseline = c_1_b, Twice = c_1_c)
# 5-year
df_c_5 <- data.frame(Years = 5, Days = seq_along(c_5_a) - 1, Once = c_5_a, Baseline = c_5_b, Twice = c_5_c)
# 10-year
df_c_10 <- data.frame(Years = 10, Days = seq_along(c_10_a) - 1, Once = c_10_a, Baseline = c_10_b, Twice = c_10_c)
# 50-year
df_c_50 <- data.frame(Years = 20, Days = seq_along(c_50_a) - 1, Once = c_50_a, Baseline = c_50_b, Twice = c_50_c)
# Combine all
df_c_all <- rbind(df_c_1, df_c_5, df_c_10, df_c_50)

head(df_c_all)
library(tidyr)
# Long format
df_c_all_long <- df_c_all %>%
  pivot_longer(
    cols = Once:Twice,
    names_to = "Scenario",
    values_to = "Colonisation"
  )
head(df_c_all_long)
plot_3 <- ggplot(df_c_all_long, aes(x = Days, y = Colonisation, color = Scenario)) +
  geom_line() +
  geom_area(fill = "skyblue", alpha = 0.3, position = "identity") +
  geom_hline(
    yintercept = Tanzania_df_all_long$Resistance[2],
    color      = "white",
    linetype   = "dashed",
    linewidth  = 1 # `size` is deprecated in recent ggplot2; use `linewidth`
  ) +
  scale_x_continuous(
    breaks = seq(0, 5 * 3650, by = 365),
    labels = seq(0, 50) # optionally label as years 0–10
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10)
  ) +
  facet_wrap(Scenario ~ Years, scales = "free_x", ncol = 4) +
  labs(x = "Years") +
  labs(y = "Colonisation level")
print(plot_3)
#
df_c_all_long_1 <- df_c_all_long %>%
  filter(Scenario != "Baseline")
plot_4 <- ggplot(df_c_all_long_1, aes(x = Scenario, y = Colonisation, color = Scenario)) +
  geom_point(size = 1) +
  # geom_boxplot(color="black",fill="skyblue") +
  geom_hline(
    yintercept = Tanzania_df_all_long$Resistance[2], # horizontal line at 15%
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) + # y-axis 0–100 with steps of 10
  facet_wrap(~Years, scales = "free_x", ncol = 4)
print(plot_4)
#
print(plot_3)
print(plot_4)
# print(population)
ggsave("Figure 1.Resistance.png", plot = plot_1, width = 6, height = 4)
ggsave("Figure 2.Resistance.png", plot = plot_2, width = 6, height = 4)
ggsave("Figure 3.colonisation level.png", plot = plot_3, width = 6, height = 4)
ggsave("Figure 4.colonisation level.png", plot = plot_4, width = 6, height = 4)

# Daily incidence of resistance
# ~~~ 1 ~~~
library(dplyr)
out_1_a_Tanzania <- as.data.frame(out_1_a_Tanzania)
out_1_b_Tanzania <- as.data.frame(out_1_b_Tanzania)
out_1_c_Tanzania <- as.data.frame(out_1_c_Tanzania)
out_5_a_Tanzania <- as.data.frame(out_5_a_Tanzania)
out_5_b_Tanzania <- as.data.frame(out_5_b_Tanzania)
out_5_c_Tanzania <- as.data.frame(out_5_c_Tanzania)
out_10_a_Tanzania <- as.data.frame(out_10_a_Tanzania)
out_10_b_Tanzania <- as.data.frame(out_10_b_Tanzania)
out_10_c_Tanzania <- as.data.frame(out_10_c_Tanzania)

out_50_a_Tanzania <- as.data.frame(out_50_a_Tanzania)
out_50_b_Tanzania <- as.data.frame(out_50_b_Tanzania)
out_50_c_Tanzania <- as.data.frame(out_50_c_Tanzania)
# New data frames with cumulative and daily incidence
out_1_a_Tanzania <- out_1_a_Tanzania %>%
  mutate(
    cum_inc = rowSums(across(all_of(CumIncRindex + 1))),
    Incidence = c(0, diff(cum_inc)) #
  )
out_1_b_Tanzania <- out_1_b_Tanzania %>%
  mutate(
    cum_inc = rowSums(across(all_of(CumIncRindex + 1))),
    Incidence = c(0, diff(cum_inc)) #
  )
out_1_c_Tanzania <- out_1_c_Tanzania %>%
  mutate(
    cum_inc = rowSums(across(all_of(CumIncRindex + 1))),
    Incidence = c(0, diff(cum_inc)) #
  )
# Scenario 5
out_5_a_Tanzania <- out_5_a_Tanzania %>%
  mutate(
    cum_inc = rowSums(across(all_of(CumIncRindex + 1))),
    Incidence = c(0, diff(cum_inc))
  )
out_5_b_Tanzania <- out_5_b_Tanzania %>%
  mutate(
    cum_inc = rowSums(across(all_of(CumIncRindex + 1))),
    Incidence = c(0, diff(cum_inc))
  )
out_5_c_Tanzania <- out_5_c_Tanzania %>%
  mutate(
    cum_inc = rowSums(across(all_of(CumIncRindex + 1))),
    Incidence = c(0, diff(cum_inc))
  )

# Scenario 10
out_10_a_Tanzania <- out_10_a_Tanzania %>%
  mutate(
    cum_inc = rowSums(across(all_of(CumIncRindex + 1))),
    Incidence = c(0, diff(cum_inc))
  )
out_10_b_Tanzania <- out_10_b_Tanzania %>%
  mutate(
    cum_inc = rowSums(across(all_of(CumIncRindex + 1))),
    Incidence = c(0, diff(cum_inc))
  )
out_10_c_Tanzania <- out_10_c_Tanzania %>%
  mutate(
    cum_inc = rowSums(across(all_of(CumIncRindex + 1))),
    Incidence = c(0, diff(cum_inc))
  )

# Scenario 50
out_50_a_Tanzania <- out_50_a_Tanzania %>%
  mutate(
    cum_inc = rowSums(across(all_of(CumIncRindex + 1))),
    Incidence = c(0, diff(cum_inc))
  )
out_50_b_Tanzania <- out_50_b_Tanzania %>%
  mutate(
    cum_inc = rowSums(across(all_of(CumIncRindex + 1))),
    Incidence = c(0, diff(cum_inc))
  )
out_50_c_Tanzania <- out_50_c_Tanzania %>%
  mutate(
    cum_inc = rowSums(across(all_of(CumIncRindex + 1))),
    Incidence = c(0, diff(cum_inc))
  )
# Extraction of incidence columns
incidence_1_a <- out_1_a_Tanzania[, "Incidence"]
incidence_1_b <- out_1_b_Tanzania[, "Incidence"]
incidence_1_c <- out_1_c_Tanzania[, "Incidence"]
plot(incidence_1_a)
incidence_5_a <- out_5_a_Tanzania[, "Incidence"]
incidence_5_b <- out_5_b_Tanzania[, "Incidence"]
incidence_5_c <- out_5_c_Tanzania[, "Incidence"]

incidence_10_a <- out_10_a_Tanzania[, "Incidence"]
incidence_10_b <- out_10_b_Tanzania[, "Incidence"]
incidence_10_c <- out_10_c_Tanzania[, "Incidence"]

incidence_50_a <- out_50_a_Tanzania[, "Incidence"]
incidence_50_b <- out_50_b_Tanzania[, "Incidence"]
incidence_50_c <- out_50_c_Tanzania[, "Incidence"]

df_incidence_1 <- data.frame(Years = 1, Days = seq_along(incidence_1_a) - 1, Once = incidence_1_a, Baseline = incidence_1_b, Twice = incidence_1_c)
df_incidence_5 <- data.frame(Years = 5, Days = seq_along(incidence_5_a) - 1, Once = incidence_5_a, Baseline = incidence_5_b, Twice = incidence_5_c)
df_incidence_10 <- data.frame(Years = 10, Days = seq_along(incidence_10_a) - 1, Once = incidence_10_a, Baseline = incidence_10_b, Twice = incidence_10_c)
df_incidence_50 <- data.frame(Years = 50, Days = seq_along(incidence_50_a) - 1, Once = incidence_50_a, Baseline = incidence_50_b, Twice = incidence_50_c)

# Combine all
Tanzania_df_incidence_all <- rbind(df_incidence_1, df_incidence_5, df_incidence_10, df_incidence_50)
head(Tanzania_df_incidence_all)
colnames(Tanzania_df_incidence_all)[3:5] <- c("MDA", "No-MDA", "Bi-MDA")
head(Tanzania_df_incidence_all)
library(tidyr)
# Long format
Tanzania_df_incidence_all_long <- Tanzania_df_incidence_all %>%
  pivot_longer(
    cols = c(MDA, `No-MDA`, `Bi-MDA`),
    names_to = "Strategy",
    values_to = "Incidence"
  )
# Specific Visualisation
#  plot
plot_1_a <- ggplot(out_1_a_Tanzania, aes(x = time, y = Incidence)) +
  geom_col(fill = "#ff7f0e", width = 0.8) +
  scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M")) +
  # scale_x_continuous(breaks = seq(0, max(out_1_a_Tanzania$time), by = 30)) +
  scale_x_continuous(
    breaks = seq(0, 5 * 3650, by = 365),
    labels = seq(0, 50) # optionally label as years 0–10
  ) +
  labs(
    title = "B",
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1), # rotate x labels
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0)
  )
print(plot_1_a)
# Plot B
plot_1_b <- ggplot(out_1_b_Tanzania, aes(x = time, y = Incidence)) +
  geom_col(fill = "#1f77b4", width = 0.8) + # different color for B
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(
    breaks = seq(0, 5 * 3650, by = 365),
    labels = seq(0, 50) # optionally label as years 0–10
  ) +
  # scale_x_continuous(breaks = seq(0, max(out_1_b_Tanzania$time), by = 30)) +
  labs(
    title = "A",
    x = "",
    y = "Incidence"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1), # rotate x labels
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0)
  )
print(plot_1_b)

# Plot C
plot_1_c <- ggplot(out_1_c_Tanzania, aes(x = time, y = Incidence)) +
  geom_col(fill = "#2ca02c", width = 0.8) + # different color for C
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(
    breaks = seq(0, 5 * 3650, by = 365),
    labels = seq(0, 50) # optionally label as years 0–10
  ) +
  # scale_x_continuous(breaks = seq(0, max(out_1_c_Tanzania$time), by = 30)) +
  labs(
    title = "C",
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1), # rotate x labels
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0)
  )
print(plot_1_c)
legend <- "Figure 1.Daily incidence:B.Baseline,A.MDA,C.Bi-MDA"
Y1 <- grid.arrange(plot_1_b, plot_1_a, plot_1_c, ncol = 3)
# Plot 5A
plot_5_a <- ggplot(out_5_a_Tanzania, aes(x = time, y = Incidence)) +
  geom_col(fill = "#ff7f0e", width = 0.8) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(
    breaks = seq(0, 5 * 3650, by = 365),
    labels = seq(0, 50) #
  ) +
  # scale_x_continuous(breaks = seq(0, max(out_5_a_Tanzania$time), by = 30)) +
  labs(title = "E", x = "", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0)
  )
print(plot_5_a)
# Plot 5B
plot_5_b <- ggplot(out_5_b_Tanzania, aes(x = time, y = Incidence)) +
  geom_col(fill = "#1f77b4", width = 0.8) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(
    breaks = seq(0, 5 * 3650, by = 365),
    labels = seq(0, 50) # optionally label as years 0–10
  ) +
  # scale_x_continuous(breaks = seq(0, max(out_5_b_Tanzania$time), by = 30)) +
  labs(title = "D", x = "", y = "Incidence") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0)
  )
print(plot_5_b)
# Plot 5C
plot_5_c <- ggplot(out_5_c_Tanzania, aes(x = time, y = Incidence)) +
  geom_col(fill = "#2ca02c", width = 0.8) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(
    breaks = seq(0, 5 * 3650, by = 365),
    labels = seq(0, 50) # optionally label as years 0–10
  ) +
  # scale_x_continuous(breaks = seq(0, max(out_5_c_Tanzania$time), by = 30)) +
  labs(title = "F", x = "", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0)
  )
print(plot_5_c)
# Combine 5Y plots
legend_5 <- "Figure 2. Daily incidence at 5 years: B. Baseline, A. MDA, C. Bi-MDA"
Y5 <- grid.arrange(plot_5_b, plot_5_a, plot_5_c, ncol = 3)
# Plot 10A
plot_10_a <- ggplot(out_10_a_Tanzania, aes(x = time, y = Incidence)) +
  geom_col(fill = "#ff7f0e", width = 0.8) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(
    breaks = seq(0, 5 * 3650, by = 365),
    labels = seq(0, 50) # optionally label as years 0–10
  ) +
  # scale_x_continuous(breaks = seq(0, max(out_10_a_Tanzania$time), by = 30)) +
  labs(title = "I", x = "Time(Years)", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0)
  )
print(plot_10_a)
# Plot 10B
plot_10_b <- ggplot(out_10_b_Tanzania, aes(x = time, y = Incidence)) +
  geom_col(fill = "#1f77b4", width = 0.8) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(
    breaks = seq(0, 5 * 3650, by = 365),
    labels = seq(0, 50) # optionally label as years 0–10
  ) +
  # scale_x_continuous(breaks = seq(0, max(out_10_b_Tanzania$time), by = 30)) +
  labs(title = "G", x = "Time (Years)", y = "Incidence") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0)
  )
print(plot_10_b)
# Plot 10C
plot_10_c <- ggplot(out_10_c_Tanzania, aes(x = time, y = Incidence)) +
  geom_col(fill = "#2ca02c", width = 0.8) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = "M")) +
  scale_x_continuous(
    breaks = seq(0, 5 * 3650, by = 365),
    labels = seq(0, 50) # optionally label as years 0–10
  ) +
  # scale_x_continuous(breaks = seq(0, max(out_10_c_Tanzania$time), by = 30)) +
  labs(title = "H", x = "Time (Years)", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.ticks = element_line(size = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0)
  )
print(plot_10_c)
# Combine 10Y plots
legend_10 <- "Figure 3. Daily incidence at 10 years: B. Baseline, A. MDA, C. Bi-MDA"
Y10 <- grid.arrange(plot_10_b, plot_10_a, plot_10_c, ncol = 3) # , bottom = legend_10)
Y_1_5_10 <- grid.arrange(Y1, Y5, Y10)
library(gridExtra)

Y_1_5_10 <- grid.arrange(Y1, Y5, Y10)
ggsave(
  filename = "Figure 1 Incidence poster_plot.png",
  plot = Y_1_5_10,
  width = 20, height = 10, # adjust for poster layout
  units = "in",
  dpi = 300
)

# Y_1_5_10 <- patchwork::wrap_plots(plot_10_b, plot_10_a, plot_10_c, ncol = 3) +
#  patchwork::plot_annotation(tag_levels = "A")
# Y_1_5_10

# Resistance : R_total<-R_total+Rs_total

# //a_1_a<-tail(cumsum(rowSums(out_1_a_Tanzania[, c(Rindex,Rsindex) + 1])), 1)
# //a_1_b<-tail(cumsum(rowSums(out_1_b_Tanzania[, c(Rindex,Rsindex) + 1])), 1)
# //a_1_c<-tail(cumsum(rowSums(out_1_c_Tanzania[, c(Rindex,Rsindex) + 1])), 1)
# //a_5_a<-tail(cumsum(rowSums(out_5_a_Tanzania[, c(Rindex,Rsindex) + 1])), 1)
# //a_5_b<-tail(cumsum(rowSums(out_5_b_Tanzania[, c(Rindex,Rsindex) + 1])), 1)
# //a_5_c<-tail(cumsum(rowSums(out_5_c_Tanzania[, c(Rindex,Rsindex) + 1])), 1)
# //a_10_a<-tail(cumsum(rowSums(out_10_a_Tanzania[, c(Rindex,Rsindex) + 1])), 1)
# //a_10_b<-tail(cumsum(rowSums(out_10_b_Tanzania[, c(Rindex,Rsindex) + 1])), 1)
# //a_10_c<-tail(cumsum(rowSums(out_10_c_Tanzania[, c(Rindex,Rsindex) + 1])), 1)

# Burden and  change in resistance
# I.1 Years
# Cum R : using R and Rs
cum_1_a <- rowSums(out_1_a_Tanzania[, c(Rindex, Rsindex) + 1])
cum_1_b <- rowSums(out_1_b_Tanzania[, c(Rindex, Rsindex) + 1])
cum_1_c <- rowSums(out_1_c_Tanzania[, c(Rindex, Rsindex) + 1])
#
a_1_a <- tail(cum_1_a, 1) - cum_1_a[1]
a_1_b <- tail(cum_1_b, 1) - cum_1_b[1]
a_1_c <- tail(cum_1_c, 1) - cum_1_c[1]
# Absolute changes
diff_1_a <- a_1_a - a_1_b
diff_1_c <- a_1_c - a_1_b
# Perectage change
perc_1_a <- 100 * (a_1_a - a_1_b) / cum_1_b[1]
perc_1_c <- 100 * (a_1_c - a_1_b) / cum_1_b[1]
# II.5 Years
# Cum R
cum_5_a <- rowSums(out_5_a_Tanzania[, c(Rindex, Rsindex) + 1])
cum_5_b <- rowSums(out_5_b_Tanzania[, c(Rindex, Rsindex) + 1])
cum_5_c <- rowSums(out_5_c_Tanzania[, c(Rindex, Rsindex) + 1])

# Absolute accumulation over time
a_5_a <- tail(cum_5_a, 1) - cum_5_a[1]
a_5_b <- tail(cum_5_b, 1) - cum_5_b[1]
a_5_c <- tail(cum_5_c, 1) - cum_5_c[1]

# Absolute changes (vs baseline = b)
diff_5_a <- a_5_a - a_5_b
diff_5_c <- a_5_c - a_5_b

# Percentage change
perc_5_a <- 100 * (a_5_a - a_5_b) / cum_5_b[1]
perc_5_c <- 100 * (a_5_c - a_5_b) / cum_5_b[1]

# III.3: 10Y
# Cum R
cum_10_a <- rowSums(out_10_a_Tanzania[, c(Rindex, Rsindex) + 1])
cum_10_b <- rowSums(out_10_b_Tanzania[, c(Rindex, Rsindex) + 1])
cum_10_c <- rowSums(out_10_c_Tanzania[, c(Rindex, Rsindex) + 1])

# Absolute accumulation over time
a_10_a <- tail(cum_10_a, 1) - cum_10_a[1]
a_10_b <- tail(cum_10_b, 1) - cum_10_b[1]
a_10_c <- tail(cum_10_c, 1) - cum_10_c[1]

# Absolute changes (vs baseline = b)
diff_10_a <- a_10_a - a_10_b
diff_10_c <- a_10_c - a_10_b

# Percentage change
perc_10_a <- 100 * (a_10_a - a_10_b) / cum_10_b[1]
perc_10_c <- 100 * (a_10_c - a_10_b) / cum_10_b[1]
# Data frames...........................................................................
# .A.Total resistance occured during the period

totalresistance <- data.frame(
  Scenario = c("1Y MDA", "1Y No-MDA", "1Y Bi-MDA", "5Y MDA", "5Y No-MDA", "5Y Bi-MDA", "10Y MDA", "10Y No-MDA", "10Y Bi-MDA"),
  R_final = c(a_1_a, a_1_b, a_1_c, a_5_a, a_5_b, a_5_c, a_10_a, a_10_b, a_10_c)
)

# Table
print(totalresistance)
# Strategy
totalresistance$Strategy <- c("MDA", "No-MDA", "Bi-MDA", "MDA", "No-MDA", "Bi-MDA", "MDA", "No-MDA", "Bi-MDA")
totalresistance$Horizon <- c("1Y", "1Y", "1Y", "5Y", "5Y", "5Y", "10Y", "10Y", "10Y")
# Horizon factor in desired order
totalresistance$Horizon <- factor(
  totalresistance$Horizon,
  levels = c("10Y", "5Y", "1Y")
)
# table
print(totalresistance)

# Visualization
library(ggplot2)
head(totalresistance)
p1 <- ggplot(
  totalresistance |>
    filter(Strategy != "No-MDA"),
  aes(x = Horizon, y = R_final, fill = Strategy)
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c(
    "No-MDA" = "grey40",
    "MDA" = "black",
    "Bi-MDA" = "red"
  )) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = " M")) +
  labs(
    title = "A",
    x = "Time horizon",
    y = "Excess resistant cases",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0)
  )
p1
# Change oder: 1,5,10
totalresistance$Horizon <- factor(
  totalresistance$Horizon,
  levels = c("1Y", "5Y", "10Y")
)
p2 <- ggplot(
  totalresistance,
  aes(x = Horizon, y = R_final, color = Strategy)
) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  # geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")+
  geom_line(aes(group = Strategy),
    position = position_dodge(width = 0.5),
    linewidth = 0.8
  ) +
  scale_color_manual(values = c(
    "MDA" = "black",
    "Bi-MDA" = "red", "No-MDA" = "grey"
  )) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = " M")) +
  labs(
    x = "Time horizon",
    y = "Excess resistant cases",
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal"
  )
p2
p1
p2
grid.arrange(p1, p2, ncol = 2)
# B.Excess of resistance
exces_resistant_cases_1_a <- round((a_1_a - a_1_b) * 100 / cum_1_b[1], 1)
exces_resistant_cases_1_c <- round((a_1_c - a_1_b) * 100 / cum_1_b[1], 1)
exces_resistant_cases_5_a <- round((a_5_a - a_5_b) * 100 / cum_5_b[1], 1)
exces_resistant_cases_5_c <- round((a_5_c - a_5_b) * 100 / cum_5_b[1], 1)
exces_resistant_cases_10_a <- round((a_10_a - a_10_b) * 100 / cum_10_b[1], 1)
exces_resistant_cases_10_c <- round((a_10_c - a_10_b) * 100 / cum_10_b[1], 1)

total_exces_resistant_cases <- data.frame(
  Scenario = c("1Y MDA", "1Y Bi-MDA", "5Y MDA", "5Y Bi-MDA", "10Y MDA", "10Y Bi-MDA"),
  R_final = c(exces_resistant_cases_1_a, exces_resistant_cases_1_c, exces_resistant_cases_5_a, exces_resistant_cases_5_c, exces_resistant_cases_10_a, exces_resistant_cases_10_c)
)
total_exces_resistant_cases$Horizon <- c("1Y", "1Y", "5Y", "5Y", "10Y", "10Y")
total_exces_resistant_cases$Strategy <- c("MDA", "Bi-MDA", "MDA", "Bi-MDA", "MDA", "Bi-MDA")
# order
total_exces_resistant_cases$Horizon <- factor(
  total_exces_resistant_cases$Horizon,
  levels = c("1Y", "5Y", "10Y")
)
library(ggplot2)
p3 <- ggplot(
  total_exces_resistant_cases,
  aes(x = Horizon, y = R_final, fill = Strategy)
) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_text(
    aes(label = round(R_final, 1)),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 4
  ) +
  scale_fill_manual(values = c(
    "MDA" = "black",
    "Bi-MDA" = "red"
  )) +
  labs(
    title = "B",
    x = "Time horizon",
    y = "Increase in resistance(%)",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5)
  )
p3
p4 <- ggplot(
  total_exces_resistant_cases,
  aes(x = Horizon, y = R_final, color = Strategy)
) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  # geom_hline(yintercept = 0, linetype = "dashed", color = "grey40")+
  geom_line(aes(group = Strategy),
    position = position_dodge(width = 0.5),
    linewidth = 0.8
  ) +
  scale_color_manual(values = c(
    "MDA" = "black",
    "Bi-MDA" = "red"
  )) +
  # scale_y_continuous(labels = scales::label_number(scale = 1e-9, suffix = " B")) +

  labs(
    title = "B",
    x = "Time horizon",
    y = "Increase in resistance(%)",
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0),
    legend.position = "top",
    legend.direction = "horizontal"
  )
p1
p2
p3
p4
grid.arrange(p3, p4, ncol = 2)
grid.arrange(p1, p3, ncol = 2)
#
theme_lancet <- function(base_size = 14, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      legend.position = "none"
    )
}

# Deaths : D and AMRDindex <-(7 * n_age + 1):(8 * n_age)     # Cummulative resistance
D_cum_1_a <- rowSums(out_1_a_Tanzania[, Dindex + 1])
D_cum_1_b <- rowSums(out_1_b_Tanzania[, Dindex + 1])
D_cum_1_c <- rowSums(out_1_c_Tanzania[, Dindex + 1])

D_cum_5_a <- rowSums(out_1_a_Tanzania[, Dindex + 1])
D_cum_5_b <- rowSums(out_1_b_Tanzania[, Dindex + 1])
D_cum_5_c <- rowSums(out_1_c_Tanzania[, Dindex + 1])

D_daily_1_b <- c(0, diff(D_cum_1_b))
mortality <- round(D_daily_1_b * 100000 / (X_total + S_total + R_total + Rs_total + Sr_total), 1)
summary(mortality)
plot(mortality)
# ...........................................................
# Mortality reduction due to MDA in under 5
head(out_1_a_Tanzania)
print(out_1_a_Tanzania)
# Extraction of time column and all columns for ages 0-4
ages_under_5 <- 0:4
# First, let's see what the actual column names are
colnames(results_1_a_Tanzania)
# Or see just the first 30 column names to understand the pattern
head(colnames(results_1_a_Tanzania), 30)
# Total number of columns
ncol(results_1_a_Tanzania)
# Define ages under 5 (0-4)
ages_under_5 <- 0:4
# Compartments (in the correct order)
compartments <- c("X", "S", "R", "Sr", "Rs", "D", "CumIncR", "AMRD")
# Column names for under-five across all compartments
under_five_cols <- c()
for (comp in compartments) {
  for (age in ages_under_5) {
    col_name <- paste0(comp, "_", age)
    under_five_cols <- c(under_five_cols, col_name)
  }
}
# Data for under five age groups
under_five_data <- results_1_a_Tanzania[, c("time", under_five_cols)]
head(under_five_data) #
results_1_a_Tanzania_under_five <- results_1_a_Tanzania[, c("time", under_five_cols)]
results_1_b_Tanzania_under_five <- results_1_b_Tanzania[, c("time", under_five_cols)]
results_1_c_Tanzania_under_five <- results_1_c_Tanzania[, c("time", under_five_cols)]

results_5_a_Tanzania_under_five <- results_5_a_Tanzania[, c("time", under_five_cols)]
results_5_b_Tanzania_under_five <- results_5_b_Tanzania[, c("time", under_five_cols)]
results_5_c_Tanzania_under_five <- results_5_c_Tanzania[, c("time", under_five_cols)]

results_10_a_Tanzania_under_five <- results_10_a_Tanzania[, c("time", under_five_cols)]
results_10_b_Tanzania_under_five <- results_10_b_Tanzania[, c("time", under_five_cols)]
results_10_c_Tanzania_under_five <- results_10_c_Tanzania[, c("time", under_five_cols)]

results_50_a_Tanzania_under_five <- results_50_a_Tanzania[, c("time", under_five_cols)]
results_50_b_Tanzania_under_five <- results_50_b_Tanzania[, c("time", under_five_cols)]
results_50_c_Tanzania_under_five <- results_50_c_Tanzania[, c("time", under_five_cols)]

head(results_50_c_Tanzania_under_five)

# Under five level deaths
D_cols <- grep("^D_", colnames(results_1_a_Tanzania_under_five))
AMRD_cols <- grep("^AMRD_", colnames(results_1_a_Tanzania_under_five))
D_cols <- grep("^D_", colnames(results_1_a_Tanzania_under_five))

# Total deaths
d_1_a <- sum(c(0, diff(rowSums(results_1_a_Tanzania_under_five[, D_cols]))))
d_1_b <- sum(c(0, diff(rowSums(results_1_b_Tanzania_under_five[, D_cols]))))
d_1_c <- sum(c(0, diff(rowSums(results_1_c_Tanzania_under_five[, D_cols]))))

d_5_a <- sum(c(0, diff(rowSums(results_5_a_Tanzania_under_five[, D_cols]))))
d_5_b <- sum(c(0, diff(rowSums(results_5_b_Tanzania_under_five[, D_cols]))))
d_5_c <- sum(c(0, diff(rowSums(results_5_c_Tanzania_under_five[, D_cols]))))

d_10_a <- sum(c(0, diff(rowSums(results_10_a_Tanzania_under_five[, D_cols]))))
d_10_b <- sum(c(0, diff(rowSums(results_10_b_Tanzania_under_five[, D_cols]))))
d_10_c <- sum(c(0, diff(rowSums(results_10_c_Tanzania_under_five[, D_cols]))))

d_50_a <- sum(c(0, diff(rowSums(results_50_a_Tanzania_under_five[, D_cols]))))
d_50_b <- sum(c(0, diff(rowSums(results_50_b_Tanzania_under_five[, D_cols]))))
d_50_c <- sum(c(0, diff(rowSums(results_50_c_Tanzania_under_five[, D_cols]))))

# Deaths averted (MDA effect):No MDA-MDA (minus AMR deaths )
# No MDA  MDA

# //d_1_d_a  <-  d_1_a  - d_1_a  # annual
# //d_1_d_bi  <- d_1_b  - d_1_c  # bi annual
# d_1_d_a  <-(d_1_b - AMRD_1_b) - (d_1_a - AMRD_1_a) The true version will be bellow AMRD
# d_1_d_bi  <-(d_1_b - AMRD_1_b) - (d_1_c - AMRD_1_a)
# //d_5_d_a  <- d_5_b  - d_5_a
# //d_5_d_bi <- d_5_b  - d_5_c
# d_5_d_a  <-(d_5_b - AMRD_5_b) - (d_5_a - AMRD_5_a)
# d_5_d_bi  <-(d_5_b - AMRD_5_b) - (d_5_c - AMRD_5_a)
# //d_10_d_a <- d_10_b - d_10_a
# //d_10_d_bi <- d_10_b - d_10_c
# d_10_d_a  <-(d_10_b - AMRD_10_b) - (d_10_a - AMRD_10_a)
# d_10_d_bi  <-(d_10_b - AMRD_10_b) - (d_10_c - AMRD_10_a)
# ...............................................................................
# AMR related deaths
AMRD_1_a <- sum(c(0, diff(rowSums(results_1_a_Tanzania_under_five[, AMRD_cols]))))
AMRD_1_b <- sum(c(0, diff(rowSums(results_1_b_Tanzania_under_five[, AMRD_cols]))))
AMRD_1_c <- sum(c(0, diff(rowSums(results_1_c_Tanzania_under_five[, AMRD_cols]))))

AMRD_5_a <- sum(c(0, diff(rowSums(results_5_a_Tanzania_under_five[, AMRD_cols]))))
AMRD_5_b <- sum(c(0, diff(rowSums(results_5_b_Tanzania_under_five[, AMRD_cols]))))
AMRD_5_c <- sum(c(0, diff(rowSums(results_5_c_Tanzania_under_five[, AMRD_cols]))))

AMRD_10_a <- sum(c(0, diff(rowSums(results_10_a_Tanzania_under_five[, AMRD_cols]))))
AMRD_10_b <- sum(c(0, diff(rowSums(results_10_b_Tanzania_under_five[, AMRD_cols]))))
AMRD_10_c <- sum(c(0, diff(rowSums(results_10_c_Tanzania_under_five[, AMRD_cols]))))

AMRD_50_a <- sum(c(0, diff(rowSums(results_50_a_Tanzania_under_five[, AMRD_cols]))))
AMRD_50_b <- sum(c(0, diff(rowSums(results_50_b_Tanzania_under_five[, AMRD_cols]))))
AMRD_50_c <- sum(c(0, diff(rowSums(results_50_c_Tanzania_under_five[, AMRD_cols]))))
#
d_1_d_a <- (d_1_b - AMRD_1_b) - (d_1_a - AMRD_1_a)
d_1_d_bi <- (d_1_b - AMRD_1_b) - (d_1_c - AMRD_1_a)
#
d_5_d_a <- (d_5_b - AMRD_5_b) - (d_5_a - AMRD_5_a)
d_5_d_bi <- (d_5_b - AMRD_5_b) - (d_5_c - AMRD_5_a)
#
d_10_d_a <- (d_10_b - AMRD_10_b) - (d_10_a - AMRD_10_a)
d_10_d_bi <- (d_10_b - AMRD_10_b) - (d_10_c - AMRD_10_a)

# Deaths increased by AMR (MDA effect)
AMRD_1_d_a <- AMRD_1_b - AMRD_1_a # annual
AMRD_1_d_bi <- AMRD_1_b - AMRD_1_c # bi annual
#
AMRD_5_d_a <- AMRD_5_b - AMRD_5_a
AMRD_5_d_bi <- AMRD_5_b - AMRD_5_c
#
AMRD_10_d_a <- AMRD_10_b - AMRD_10_a
AMRD_10_d_bi <- AMRD_10_b - AMRD_10_c
# ...............................................................................
# Data frames for deaths
# a.Total deaths
totalDeaths <- data.frame(
  Scenario = c("1Y No-MDA", "1Y MDA", "1Y Bi-MDA", "5Y MDA", "5Y No-MDA", "5Y Bi-MDA", "10Y MDA", "10Y No-MDA", "10Y Bi-MDA"),
  Deaths_millions = c(d_1_b, d_1_a, d_1_c, d_5_a, d_5_b, d_5_c, d_10_a, d_10_b, d_10_c)
)
print(totalDeaths)
# b.AMR deaths
AMR_Deaths <- data.frame(
  Scenario = c("1Y No-MDA", "1Y MDA", "1Y Bi-MDA", "5Y MDA", "5Y No-MDA", "5Y Bi-MDA", "10Y MDA", "10Y No-MDA", "10Y Bi-MDA"),
  Deaths_millions = c(AMRD_1_b, AMRD_1_a, AMRD_1_c, AMRD_5_a, AMRD_5_b, AMRD_5_c, AMRD_10_a, AMRD_10_b, AMRD_10_c)
)
print(AMR_Deaths)
# Data frame for deaths averted
# a.
Tanzania_Deaths_averted <- data.frame(
  Scenario = c("1Y MDA", "1Y Bi-MDA", "5Y MDA", "5Y Bi-MDA", "10Y MDA", "10Y Bi-MDA"),
  Deaths_averted_millions = c(d_1_d_a, d_1_d_bi, d_5_d_a, d_5_d_bi, d_10_d_a, d_10_d_bi)
)
print(Tanzania_Deaths_averted)
# b.
Tanzania_AMR_Deaths_averted <- data.frame(
  Scenario = c("1Y MDA", "1Y Bi-MDA", "5Y MDA", "5Y Bi-MDA", "10Y MDA", "10Y Bi-MDA"),
  AMR_Deaths_averted_millions = c(AMRD_1_d_a, AMRD_1_d_bi, AMRD_5_d_a, AMRD_5_d_bi, AMRD_10_d_a, AMRD_10_d_bi)
)
print(Tanzania_AMR_Deaths_averted)
# I.Net Benefit (Total)
# A.Calculation
# Excess AMR deaths caused by MDA (MDA scenario minus No-MDA baseline)
excess_AMRD_1_a <- AMRD_1_a - AMRD_1_b # annual MDA adds this many extra AMR deaths
excess_AMRD_1_bi <- AMRD_1_c - AMRD_1_b
excess_AMRD_5_a <- AMRD_5_a - AMRD_5_b
excess_AMRD_5_bi <- AMRD_5_c - AMRD_5_b
excess_AMRD_10_a <- AMRD_10_a - AMRD_10_b
excess_AMRD_10_bi <- AMRD_10_c - AMRD_10_b
# Net benefit  = Deaths averted by MDA − Excess AMR deaths caused by MDA
# Positive = net benefit; Negative = net harm
net_benefit_1_a <- d_1_d_a - excess_AMRD_1_a
net_benefit_1_bi <- d_1_d_bi - excess_AMRD_1_bi
net_benefit_5_a <- d_5_d_a - excess_AMRD_5_a
net_benefit_5_bi <- d_5_d_bi - excess_AMRD_5_bi
net_benefit_10_a <- d_10_d_a - excess_AMRD_10_a
net_benefit_10_bi <- d_10_d_bi - excess_AMRD_10_bi

# Summary table
net_benefit_summary <- data.frame(
  Scenario = c("1Y MDA", "1Y Bi-MDA", "5Y MDA", "5Y Bi-MDA", "10Y MDA", "10Y Bi-MDA"),
  Strategy = rep(c("MDA", "Bi-MDA"), 3),
  Horizon = rep(c("1Y", "5Y", "10Y"), each = 2),
  Deaths_averted = c(d_1_d_a, d_1_d_bi, d_5_d_a, d_5_d_bi, d_10_d_a, d_10_d_bi),
  Excess_AMR_deaths = c(
    excess_AMRD_1_a, excess_AMRD_1_bi,
    excess_AMRD_5_a, excess_AMRD_5_bi,
    excess_AMRD_10_a, excess_AMRD_10_bi
  ),
  Net_benefit = c(
    net_benefit_1_a, net_benefit_1_bi,
    net_benefit_5_a, net_benefit_5_bi,
    net_benefit_10_a, net_benefit_10_bi
  )
)

net_benefit_summary$Horizon <- factor(net_benefit_summary$Horizon,
  levels = c("1Y", "5Y", "10Y")
)

print(net_benefit_summary)
# Percentage: net benefit as % of No-MDA baseline deaths
net_benefit_summary$Net_benefit_pct_1Y <- round(net_benefit_summary$Net_benefit * 100 / d_1_b, 2)
net_benefit_summary$Net_benefit_pct_5Y <- round(net_benefit_summary$Net_benefit * 100 / d_5_b, 2)
net_benefit_summary$Net_benefit_pct_10Y <- round(net_benefit_summary$Net_benefit * 100 / d_10_b, 2)
# ...............................................................................
# B.Visualisation: Net Benefit of MDA in values
# ...............................................................................
pacman::p_load(ggplot2, tidyr, dplyr, scales, gridExtra)
# Long format : i will use stacked bar
# net_long <- net_benefit_summary %>%
# select(Scenario, Strategy, Horizon, Deaths_averted, Excess_AMR_deaths) %>%
# utate(Excess_AMR_deaths = -Excess_AMR_deaths) %>%  # flip sign: harm shown below zero
# ivot_longer(
# cols = c(Deaths_averted, Excess_AMR_deaths),
# names_to  = "Component",
# values_to = "Deaths"
# %>%
# utate(Component = factor(Component,
# levels = c("Deaths_averted", "Excess_AMR_deaths"),
# labels = c("Deaths averted (MDA benefit)", "Excess AMR deaths (MDA harm)")))


net_long <- net_benefit_summary %>%
  dplyr::select(Scenario, Strategy, Horizon, Deaths_averted, Excess_AMR_deaths) %>%
  dplyr::mutate(Excess_AMR_deaths = -Excess_AMR_deaths) %>%
  tidyr::pivot_longer(
    cols = c(Deaths_averted, Excess_AMR_deaths),
    names_to = "Component",
    values_to = "Deaths"
  ) %>%
  dplyr::mutate(Component = factor(Component,
    levels = c("Deaths_averted", "Excess_AMR_deaths"),
    labels = c("Deaths averted (MDA benefit)", "Excess AMR deaths (MDA harm)")
  ))
# Plot A:benefit up, harm down
p_net_A <- ggplot(net_long, aes(x = Horizon, y = Deaths / 1000, fill = Component)) +
  geom_col(position = "stack", width = 0.6, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  # Net benefit point
  geom_point(
    data = net_benefit_summary,
    aes(x = Horizon, y = Net_benefit / 1000),
    inherit.aes = FALSE,
    shape = 18, size = 4, color = "black"
  ) +
  scale_fill_manual(values = c(
    "Deaths averted (MDA benefit)"  = "#009E73",
    "Excess AMR deaths (MDA harm)"  = "#D55E00" ## D55E00
  )) +
  scale_y_continuous(labels = label_number(suffix = " K")) +
  facet_wrap(~Strategy, ncol = 2) +
  labs(
    # title = "A",#"Net mortality benefit of MDA in under-5s",
    # subtitle = "Green = deaths averted | Orange = excess AMR deaths | Diamond = net balance",
    x = "Time horizon",
    y = "Deaths (thousands)",
    fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )
print(p_net_A)
# Plot B: net benefit only (positive = net gain)
p_net_B <- ggplot(
  net_benefit_summary,
  aes(x = Horizon, y = Net_benefit / 1000, fill = ifelse(Net_benefit >= 0, "Net benefit", "Net harm"))
) +
  geom_col(position = position_dodge(0.7), width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(
    aes(
      label = round(Net_benefit / 1000, 1),
      vjust = ifelse(Net_benefit >= 0, -0.4, 1.2)
    ),
    position = position_dodge(0.7), size = 3.5
  ) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5, linetype = "dashed") +
  scale_fill_manual(values = c("Net benefit" = "#009E73", "Net harm" = "#D55E00")) +
  scale_y_continuous(labels = label_number(suffix = " K")) +
  facet_wrap(~Strategy, ncol = 2) +
  labs(
    title = "B",
    x = "Time horizon",
    y = "Net deaths averted (thousands)",
    fill = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
print(p_net_A)
print(p_net_B)
grid.arrange(p_net_A, p_net_B, ncol = 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II. NET BENEFIT IN MORTALITY — PERCENTAGE
# Here,all values will be  expressed as % of No-MDA baseline under-5 deaths
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reference baselines will be No-MDA under-5 deaths over each horizon)
baseline_1Y <- d_1_b
baseline_5Y <- d_5_b
baseline_10Y <- d_10_b
# Baseline per horizon
get_baseline <- function(horizon) {
  switch(horizon,
    "1Y" = baseline_1Y,
    "5Y" = baseline_5Y,
    "10Y" = baseline_10Y
  )
}

# Percentage version of net benefit summary
net_benefit_pct <- net_benefit_summary %>%
  mutate(
    baseline = case_when(
      Horizon == "1Y" ~ baseline_1Y,
      Horizon == "5Y" ~ baseline_5Y,
      Horizon == "10Y" ~ baseline_10Y
    ),
    Deaths_averted_pct = round(Deaths_averted * 100 / baseline, 2),
    Excess_AMR_deaths_pct = round(Excess_AMR_deaths * 100 / baseline, 2),
    Net_benefit_pct = round(Net_benefit * 100 / baseline, 2)
  )

print(net_benefit_pct[, c("Scenario", "Deaths_averted_pct", "Excess_AMR_deaths_pct", "Net_benefit_pct")])

# Visualisation of net benefit (%)
pacman::p_load(ggplot2, tidyr, dplyr, scales, gridExtra)
# Long format
# et_pct_long <- net_benefit_pct %>%
# select(Scenario, Strategy, Horizon, Deaths_averted_pct, Excess_AMR_deaths_pct) %>%
# mutate(Excess_AMR_deaths_pct = -Excess_AMR_deaths_pct) %>%
# ivot_longer(
# cols      = c(Deaths_averted_pct, Excess_AMR_deaths_pct),
# names_to  = "Component",
# values_to = "Pct"
# %>%
# mutate(Component = factor(Component,
# levels = c("Deaths_averted_pct", "Excess_AMR_deaths_pct"),
# labels = c("Deaths averted (MDA benefit)", "Excess AMR deaths (MDA harm)")))
# Plot A: stacked butterfly (%)
# Here is the new formula, MASS and dplyr conflict
net_pct_long <- net_benefit_pct %>%
  dplyr::select(Scenario, Strategy, Horizon, Deaths_averted_pct, Excess_AMR_deaths_pct) %>%
  dplyr::mutate(Excess_AMR_deaths_pct = -Excess_AMR_deaths_pct) %>%
  tidyr::pivot_longer(
    cols      = c(Deaths_averted_pct, Excess_AMR_deaths_pct),
    names_to  = "Component",
    values_to = "Pct"
  ) %>%
  dplyr::mutate(Component = factor(Component,
    levels = c("Deaths_averted_pct", "Excess_AMR_deaths_pct"),
    labels = c("Deaths averted (MDA benefit)", "Excess AMR deaths (MDA harm)")
  ))

p_net_pct_A <- ggplot(net_pct_long, aes(x = Horizon, y = Pct, fill = Component)) +
  geom_col(position = "stack", width = 0.6, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  # Net benefit diamond
  geom_point(
    data = net_benefit_pct,
    aes(x = Horizon, y = Net_benefit_pct),
    inherit.aes = FALSE,
    shape = 18, size = 4, color = "black"
  ) +
  scale_fill_manual(values = c(
    "Deaths averted (MDA benefit)"  = "#009E73",
    "Excess AMR deaths (MDA harm)"  = "#D55E00"
  )) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  facet_wrap(~Strategy, ncol = 2) +
  labs(
    title    = "A", # "Net mortality benefit of MDA in under-5s (%)",
    subtitle = "", # Green = deaths averted | Orange = excess AMR deaths | Diamond = net balance\nExpressed as % of No-MDA baseline deaths",
    x        = "Time horizon",
    y        = "% of baseline deaths",
    fill     = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )
print(p_net_pct_A)

# Plot B: net benefit only (%)
p_net_pct_B <- ggplot(
  net_benefit_pct,
  aes(
    x = Horizon, y = Net_benefit_pct,
    fill = ifelse(Net_benefit_pct >= 0, "Net benefit", "Net harm")
  )
) +
  geom_col(
    position = position_dodge(0.7), width = 0.6,
    color = "black", linewidth = 0.3
  ) +
  geom_text(
    aes(
      label = paste0(Net_benefit_pct, "%"),
      vjust = ifelse(Net_benefit_pct >= 0, -0.4, 1.2)
    ),
    position = position_dodge(0.7), size = 3.5
  ) +
  geom_hline(
    yintercept = 0, color = "grey40",
    linewidth = 0.5, linetype = "dashed"
  ) +
  scale_fill_manual(values = c("Net benefit" = "#009E73", "Net harm" = "#D55E00")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  facet_wrap(~Strategy, ncol = 2) +
  labs(
    title = "B",
    x     = "Time horizon",
    y     = "Net deaths averted (% of baseline)",
    fill  = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  )
print(p_net_pct_B)
# Combined
grid.arrange(p_net_pct_A, p_net_pct_B, ncol = 2)

ggsave("Figure_NetMortalityBenefit_Pct.png",
  plot = grid.arrange(p_net_pct_A, p_net_pct_B, ncol = 2),
  width = 14, height = 6, dpi = 300, bg = "white"
)


# Cut-off
#---------------------------------
# Step 1: Function to compute cut-off
# Function to compute cumulative net benefit and extract cut-off year
compute_cutoff <- function(out_mda, out_nomda, D_cols, AMRD_cols, horizon_years = 20) {
  # Daily deaths: No-MDA minus MDA (deaths averted by MDA)
  daily_deaths_averted <-
    c(0, diff(rowSums(out_nomda[, D_cols]))) -
    c(0, diff(rowSums(out_mda[, D_cols])))

  # Daily excess AMR deaths: MDA minus No-MDA
  daily_excess_AMRD <-
    c(0, diff(rowSums(out_mda[, AMRD_cols]))) -
    c(0, diff(rowSums(out_nomda[, AMRD_cols])))

  # Daily net benefit
  daily_net <- daily_deaths_averted - daily_excess_AMRD

  # Cumulative net benefit
  cum_net <- cumsum(daily_net)

  # Time in years
  time_years <- seq(0, horizon_years, length.out = length(cum_net))

  # Cut-off 1: year of peak cumulative net benefit (optimal stopping)
  peak_day <- which.max(cum_net)
  peak_year <- time_years[peak_day]

  # Cut-off 2: year cumulative net benefit crosses zero (harm threshold)
  # (NA if it never crosses zero within the horizon)
  cross_idx <- which(cum_net < 0)[1]
  cross_year <- if (!is.na(cross_idx)) time_years[cross_idx] else NA

  return(list(
    time_years   = time_years,
    cum_net      = cum_net,
    daily_net    = daily_net,
    peak_year    = peak_year,
    cross_year   = cross_year,
    peak_net     = max(cum_net)
  ))
}

# Step 2: Parameter grid and simulation loop
pacman::p_load(deSolve, dplyr, tidyr, ggplot2, scales)

# Parameter grid
mda_durations <- c(14, 30, 60) # days
mda_freqs <- c(1, 2) # per year (1=annual, 2=bi-annual)
horizon_years <- 10
tvec_long <- seq(0, horizon_years * 365.25, 1)

# Storage
results_grid <- list()

for (dur in mda_durations) {
  for (freq in mda_freqs) {
    cat("Running: duration =", dur, "| frequency =", freq, "\n")

    # MDA parameters
    parms_run <- parms # start from your baseline parms
    parms_run$mda_duration <- dur
    parms_run$r_mda <- -log(1 - parms_run$mda_cov) / dur
    parms_run$mda_start_times <- (0:100) * (365.25 / freq)

    # No-MDA parameters (same structure, MDA off)
    parms_nomda_run <- parms_noMDA
    parms_nomda_run$mda_start_times <- numeric(0)

    # Run model
    out_mda <- bacteria.solve(tvec_long, state, parms_run)
    out_nomda <- bacteria.solve(tvec_long, state, parms_nomda_run)

    # Data frames
    out_mda <- as.data.frame(out_mda)
    out_nomda <- as.data.frame(out_nomda)

    # Column indices in output (offset +1 for time column)
    D_cols_out <- Dindex + 1
    AMRD_cols_out <- AMRDindex + 1

    # Under-5 only (ages 0-4)
    under5_D_cols <- (Dindex[1:5]) + 1
    under5_AMRD_cols <- (AMRDindex[1:5]) + 1

    # Compute cut-off
    res <- compute_cutoff(
      out_mda = out_mda,
      out_nomda = out_nomda,
      D_cols = under5_D_cols,
      AMRD_cols = under5_AMRD_cols,
      horizon_years = horizon_years
    )

    # Store results
    results_grid[[paste(dur, freq, sep = "_")]] <- data.frame(
      Duration = dur,
      Frequency = freq,
      Freq_label = ifelse(freq == 1, "Annual MDA", "Bi-annual MDA"),
      Peak_year = round(res$peak_year, 2),
      Cross_year = round(res$cross_year, 2),
      Peak_net = res$peak_net,
      # Store trajectory for curve plots
      stringsAsFactors = FALSE
    )

    # Also store full trajectory for curve plot
    results_grid[[paste(dur, freq, "traj", sep = "_")]] <- data.frame(
      Duration    = dur,
      Frequency   = freq,
      Freq_label  = ifelse(freq == 1, "Annual MDA", "Bi-annual MDA"),
      Time_years  = res$time_years,
      Cum_net     = res$cum_net,
      Daily_net   = res$daily_net
    )
  }
}

# Combine summary and trajectory data frames
df_summary <- bind_rows(
  results_grid[!grepl("traj", names(results_grid))]
)

df_traj <- bind_rows(
  results_grid[grepl("traj", names(results_grid))]
)

# Factor labels
df_summary$Duration_label <- paste0(df_summary$Duration, "-day MDA")
df_traj$Duration_label <- paste0(df_traj$Duration, "-day MDA")

df_summary$Duration_label <- factor(df_summary$Duration_label,
  levels = c("14-day MDA", "30-day MDA", "60-day MDA")
)
df_traj$Duration_label <- factor(df_traj$Duration_label,
  levels = c("14-day MDA", "30-day MDA", "60-day MDA")
)

print(df_summary)

# Step 3: Plots
# Plot A — Cumulative net benefit trajectories (with cut-off marked)

# Vertical lines at peak year per scenario
peak_lines <- df_summary %>%
  select(Duration_label, Freq_label, Peak_year)

p_traj <- ggplot(
  df_traj,
  aes(
    x = Time_years, y = Cum_net / 1000,
    color = Freq_label, linetype = Freq_label
  )
) +
  geom_line(linewidth = 0.9) +
  geom_hline(
    yintercept = 0, linetype = "dashed",
    color = "grey40", linewidth = 0.5
  ) +
  # Mark peak year
  geom_vline(
    data = peak_lines,
    aes(xintercept = Peak_year, color = Freq_label),
    linetype = "dotted", linewidth = 0.7
  ) +
  # Mark zero-crossing
  geom_vline(
    data = df_summary %>% filter(!is.na(Cross_year)),
    aes(xintercept = Cross_year, color = Freq_label),
    linetype = "dashed", linewidth = 0.7
  ) +
  scale_color_manual(values = c(
    "Annual MDA"    = "#1f77b4",
    "Bi-annual MDA" = "#d62728"
  )) +
  scale_linetype_manual(values = c(
    "Annual MDA"    = "solid",
    "Bi-annual MDA" = "solid"
  )) +
  scale_y_continuous(labels = label_number(suffix = " K")) +
  scale_x_continuous(breaks = seq(0, horizon_years, by = 2)) +
  facet_wrap(~Duration_label, ncol = 3) +
  labs(
    # title    = "A ",#. Cumulative net mortality benefit over time",
    subtitle = "Dotted vertical = optimal stopping year | Dashed vertical = harm threshold",
    x = "Time (years)",
    y = "Cumulative net deaths averted",
    color = NULL, linetype = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(face = "bold", size = 12)
  )
print(p_traj)
# Plot B — Cut-off year as the outcome (sensitivity analysis heatmap)
# Pivot to long for both cut-off types
df_cutoff_long <- df_summary %>%
  pivot_longer(
    cols      = c(Peak_year, Cross_year),
    names_to  = "Cutoff_type",
    values_to = "Cutoff_year"
  ) %>%
  mutate(Cutoff_type = recode(Cutoff_type,
    "Peak_year"  = "Optimal stopping year\n(peak net benefit)",
    "Cross_year" = "Harm threshold year\n(net benefit = 0)"
  ))

p_heatmap <- ggplot(
  df_cutoff_long,
  aes(x = Freq_label, y = Duration_label, fill = Cutoff_year)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = ifelse(is.na(Cutoff_year), "Never\ncrosses",
      paste0("Year ", round(Cutoff_year, 1))
    )),
    size = 4.5, color = "black"
  ) +
  scale_fill_gradient2(
    low      = "#2166ac",
    mid      = "#f7f7f7",
    high     = "#d6604d",
    midpoint = horizon_years / 2,
    na.value = "grey85",
    name     = "Year"
  ) +
  facet_wrap(~Cutoff_type, ncol = 2) +
  labs(
    # title = "B", #Sensitivity of MDA cut-off year to intervention parameters",
    x     = "MDA frequency",
    y     = "MDA duration"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11),
    legend.position = "right"
  )
print(p_heatmap)
# grid.arrange(p_traj, p_heatmap,ncol=2)
# Plot C — Dot-and-line sensitivity plot (publication style)
p_dot <- ggplot(
  df_summary,
  aes(x = Duration_label, color = Freq_label)
) +
  # Optimal stopping
  geom_point(aes(y = Peak_year, shape = "Optimal stopping"),
    size = 4, position = position_dodge(0.4)
  ) +
  # Harm threshold
  geom_point(aes(y = Cross_year, shape = "Harm threshold"),
    size = 4, position = position_dodge(0.4)
  ) +
  geom_line(aes(y = Peak_year, group = Freq_label),
    position = position_dodge(0.4), linewidth = 0.8
  ) +
  geom_line(aes(y = Cross_year, group = Freq_label),
    linetype = "dashed",
    position = position_dodge(0.4), linewidth = 0.8
  ) +
  scale_color_manual(values = c(
    "Annual MDA"    = "#1f77b4",
    "Bi-annual MDA" = "#d62728"
  )) +
  scale_shape_manual(values = c(
    "Optimal stopping" = 16,
    "Harm threshold"   = 17
  )) +
  scale_y_continuous(
    breaks = seq(0, horizon_years, by = 2),
    limits = c(0, horizon_years)
  ) +
  labs(
    title  = "C", # Cut-off duration by MDA design parameters",
    x      = "MDA duration",
    y      = "Cut-off year",
    color  = "MDA frequency",
    shape  = "Cut-off definition"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")

print(p_dot)


library(ggplot2)

# Filter data for Cross_year (only where available)
df_cross <- subset(df_summary, !is.na(Cross_year))

p_dot_1 <- ggplot(
  df_summary,
  aes(x = Duration_label, color = Freq_label, group = Freq_label)
) +

  # Optimal stopping (Peak_year)
  geom_point(
    aes(y = Peak_year, shape = "Optimal stopping"),
    size = 4,
    position = position_dodge(width = 0.4)
  ) +
  geom_line(
    aes(y = Peak_year),
    linewidth = 0.8,
    position = position_dodge(width = 0.4)
  ) +

  # ---- Harm threshold (Cross_year) ----
  geom_point(
    data = df_cross,
    aes(y = Cross_year, shape = "Harm threshold"),
    size = 4,
    position = position_dodge(width = 0.4)
  ) +

  # (No geom_line for Cross_year → avoids misleading trends)

  # ---- Scales ----
  scale_color_manual(values = c(
    "Annual MDA"    = "#1f77b4",
    "Bi-annual MDA" = "#d62728"
  )) +
  scale_shape_manual(values = c(
    "Optimal stopping" = 16,
    "Harm threshold"   = 17
  )) +
  scale_y_continuous(
    breaks = seq(0, horizon_years, by = 5),
    limits = c(0, horizon_years)
  ) +

  # Labels
  labs(
    title = "C",
    x = "MDA duration",
    y = "Cut-off year",
    color = "MDA frequency",
    shape = "Cut-off definition"
  ) +

  #  Theme
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

print(p_dot_1)
# Combine and save
library(gridExtra)
print(p_traj)
print(p_heatmap)
print(p_dot)
print(p_dot_1)

fig_final <- grid.arrange(p_traj,
  grid.arrange(p_heatmap, p_dot, ncol = 2),
  nrow = 2, heights = c(1.2, 1)
)

ggsave("Figure_Cutoff_Sensitivity.png",
  plot = fig_final,
  width = 16, height = 12,
  dpi = 300, bg = "white"
)

#------------------------------------

# Impact assessment on mortality
change_1Y_a <- round((d_1_a - d_1_b) * 1000 / d_1_b, 2)
change_1Y_a
change_1Y_c <- round((d_1_c - d_1_b) * 1000 / d_1_b, 2)
change_1Y_c
#
change_5Y_a <- round((d_5_a - d_5_b) * 1000 / d_5_b, 2)
change_5Y_a
change_5Y_c <- round((d_5_c - d_5_b) * 1000 / d_5_b, 2)
change_5Y_c
#
change_10Y_a <- round((d_10_a - d_10_b) * 1000 / d_10_b, 2)
change_10Y_a
change_10Y_c <- round((d_10_c - d_10_b) * 1000 / d_10_b, 2)
change_10Y_c
#
pacman::p_load(ggplot2, dplyr, tidyr)
# Data frame
impact_deaths_Tanzania <- data.frame(
  Year = rep(c(1, 5, 10), each = 2),
  Scenario = rep(c("MDA", "Bi-MDA"), times = 3),
  Change = c(
    change_1Y_a, change_1Y_c,
    change_5Y_a, change_5Y_c,
    change_10Y_a, change_10Y_c
  )
)
# Visualisation
deaths_decrease <- ggplot(impact_deaths_Tanzania, aes(x = factor(Year), y = Change, fill = Scenario)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = paste0(Change, "‰")),
    position = position_dodge(width = 0.9), vjust = -0.5
  ) +
  labs(
    title = "B",
    # title = "Impact of MDA on under five mortality reduction in Tanzania",
    x = "Years",
    y = "Percentage change of deaths (‰)"
  ) +
  scale_fill_manual(values = c("skyblue", "tomato")) +
  theme_minimal()
print(deaths_decrease)
# grid.arrange(cases_increases,deaths_decrease,ncol=2)
# ggsave("Figure 2.Impact of MDA on mortality and resistance.png", plot=grid.arrange(cases_increases,deaths_decrease,ncol=2),width=12,height=5,dpi=300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#            Prevalence of AMR                                                 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Sheetal approach : Load required libraries
pacman::p_load(deSolve, doParallel, QuantPsyc)
#
# out_1_a_Tanzania <- bacteria.solve(tvec_1_a, state, parms)
results_1_a_Tanzania <- as.data.frame(out_1_a_Tanzania)
compartment_names <- c("X", "S", "R", "Sr", "Rs", "D", "CumIncR", "AMRD")
col_names <- c("time")
for (comp in compartment_names) {
  for (age in age_groups) {
    col_names <- c(col_names, paste0(comp, "_", age))
  }
}
colnames(results_1_a_Tanzania) <- col_names
names(results_1_a_Tanzania)
# Outcomes: Baseline Prevalence
R_total_1 <- rowSums(out_1_b_Tanzania[, c(Rindex, Rsindex) + 1])
R_total_5 <- rowSums(out_5_b_Tanzania[, c(Rindex, Rsindex) + 1])
R_total_10 <- rowSums(out_10_b_Tanzania[, c(Rindex, Rsindex) + 1])
R_total_50 <- rowSums(out_50_b_Tanzania[, c(Rindex, Rsindex) + 1])
# Outcomes: Cumulative Incidence
R_total_1_a <- cumsum(rowSums(out_1_a_Tanzania[, c(Rindex, Rsindex) + 1])) # annual
R_total_1_b <- cumsum(rowSums(out_1_b_Tanzania[, c(Rindex, Rsindex) + 1])) # no MDA
R_total_1_c <- cumsum(rowSums(out_1_c_Tanzania[, c(Rindex, Rsindex) + 1])) # bi-annual

R_total_5_a <- cumsum(rowSums(out_5_a_Tanzania[, c(Rindex, Rsindex) + 1]))
R_total_5_b <- cumsum(rowSums(out_5_b_Tanzania[, c(Rindex, Rsindex) + 1]))
R_total_5_c <- cumsum(rowSums(out_5_c_Tanzania[, c(Rindex, Rsindex) + 1]))

R_total_10_a <- cumsum(rowSums(out_10_a_Tanzania[, c(Rindex, Rsindex) + 1]))
R_total_10_b <- cumsum(rowSums(out_10_b_Tanzania[, c(Rindex, Rsindex) + 1]))
R_total_10_c <- cumsum(rowSums(out_10_c_Tanzania[, c(Rindex, Rsindex) + 1]))

R_total_50_a <- cumsum(rowSums(out_50_a_Tanzania[, c(Rindex, Rsindex) + 1]))
R_total_50_b <- cumsum(rowSums(out_50_b_Tanzania[, c(Rindex, Rsindex) + 1]))
R_total_50_c <- cumsum(rowSums(out_50_c_Tanzania[, c(Rindex, Rsindex) + 1]))

# Total colonized  population : We exclude x
N_1_a <- rowSums(out_1_a_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_1_b <- rowSums(out_1_b_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_1_c <- rowSums(out_1_c_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_5_a <- rowSums(out_5_a_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_5_b <- rowSums(out_5_b_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_5_c <- rowSums(out_5_c_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_10_a <- rowSums(out_10_a_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_10_b <- rowSums(out_10_b_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_10_c <- rowSums(out_10_c_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
# MDA reduction
N_10_a_5 <- rowSums(out_10_a_5_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_10_a_6 <- rowSums(out_10_a_6_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_10_a_7 <- rowSums(out_10_a_7_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_10_c_5 <- rowSums(out_10_c_5_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_10_c_6 <- rowSums(out_10_c_6_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_10_c_7 <- rowSums(out_10_c_7_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])

N_50_a <- rowSums(out_50_a_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_50_b <- rowSums(out_50_b_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
N_50_c <- rowSums(out_50_c_Tanzania[, c(Sindex, Srindex, Rindex, Rsindex) + 1])
# Proportion
prop_R_1_a <- round(rowSums(out_1_a_Tanzania[, c(Rindex, Rsindex) + 1]) / N_1_a * 100, 1)
prop_R_1_b <- round(rowSums(out_1_b_Tanzania[, c(Rindex, Rsindex) + 1]) / N_1_b * 100, 1)
prop_R_1_c <- round(rowSums(out_1_c_Tanzania[, c(Rindex, Rsindex) + 1]) / N_1_c * 100, 1)

prop_R_5_a <- round(rowSums(out_5_a_Tanzania[, c(Rindex, Rsindex) + 1]) / N_5_a * 100, 1)
prop_R_5_b <- round(rowSums(out_5_b_Tanzania[, c(Rindex, Rsindex) + 1]) / N_5_b * 100, 1)
prop_R_5_c <- round(rowSums(out_5_c_Tanzania[, c(Rindex, Rsindex) + 1]) / N_5_c * 100, 1)

prop_R_10_a <- round(rowSums(out_10_a_Tanzania[, c(Rindex, Rsindex) + 1]) / N_10_a * 100, 1)
prop_R_10_b <- round(rowSums(out_10_b_Tanzania[, c(Rindex, Rsindex) + 1]) / N_10_b * 100, 1)
prop_R_10_c <- round(rowSums(out_10_c_Tanzania[, c(Rindex, Rsindex) + 1]) / N_10_c * 100, 1)

prop_R_10_a_5 <- round(rowSums(out_10_a_5_Tanzania[, c(Rindex, Rsindex) + 1]) / N_10_a_5 * 100, 1)
prop_R_10_a_6 <- round(rowSums(out_10_a_6_Tanzania[, c(Rindex, Rsindex) + 1]) / N_10_a_6 * 100, 1)
prop_R_10_a_7 <- round(rowSums(out_10_a_7_Tanzania[, c(Rindex, Rsindex) + 1]) / N_10_a_7 * 100, 1)

prop_R_10_c_5 <- round(rowSums(out_10_c_5_Tanzania[, c(Rindex, Rsindex) + 1]) / N_10_c_5 * 100, 1)
prop_R_10_c_6 <- round(rowSums(out_10_c_6_Tanzania[, c(Rindex, Rsindex) + 1]) / N_10_c_6 * 100, 1)
prop_R_10_c_7 <- round(rowSums(out_10_c_7_Tanzania[, c(Rindex, Rsindex) + 1]) / N_10_c_7 * 100, 1)

prop_R_50_a <- round(rowSums(out_50_a_Tanzania[, c(Rindex, Rsindex) + 1]) / N_50_a * 100, 1)
prop_R_50_b <- round(rowSums(out_50_b_Tanzania[, c(Rindex, Rsindex) + 1]) / N_50_b * 100, 1)
prop_R_50_c <- round(rowSums(out_50_c_Tanzania[, c(Rindex, Rsindex) + 1]) / N_50_c * 100, 1)


# Here
# plot(time, run_d[, 7], type = "l")
dataset_1 <- cbind(time = out_1_a_Tanzania[, 1], R_total_1, R_total_1_a, R_total_1_b, R_total_1_c, N_1_a, N_1_b, N_1_c, prop_R_1_a, prop_R_1_b, prop_R_1_c)
dataset_1 <- as.data.frame(dataset_1)
#
dataset_5 <- cbind(time = out_5_a_Tanzania[, 1], R_total_5, R_total_5_a, R_total_5_b, R_total_5_c, N_5_a, N_5_b, N_5_c, prop_R_5_a, prop_R_5_b, prop_R_5_c)
dataset_5 <- as.data.frame(dataset_5)
#
dataset_10 <- cbind(time = out_10_a_Tanzania[, 1], R_total_10, R_total_10_a, R_total_10_b, R_total_10_c, N_10_a, N_10_b, N_10_c, prop_R_10_a, prop_R_10_b, prop_R_10_c)
dataset_10 <- as.data.frame(dataset_10)
#
dataset_10_567_annual <- cbind(time = out_10_a_Tanzania[, 1], prop_R_10_b, prop_R_10_a, prop_R_10_a_5, prop_R_10_a_6, prop_R_10_a_7)
dataset_10_567_annual <- as.data.frame(dataset_10_567_annual)

#
dataset_50 <- cbind(time = out_50_a_Tanzania[, 1], R_total_50, R_total_50_a, R_total_50_b, R_total_50_c, N_50_a, N_50_b, N_50_c, prop_R_50_a, prop_R_50_b, prop_R_50_c)
dataset_50 <- as.data.frame(dataset_50)
#

colnames(dataset_10_567_annual) <- c(
  "time",
  "Baseline", # prop_R_10_b
  "Ten", # prop_R_10_a
  "Five", # prop_R_10_a_5
  "Six", # prop_R_10_a_6
  "Seven" # prop_R_10_a_7
)
#
dataset_10_567_bi_annual <- cbind(time = out_10_c_Tanzania[, 1], prop_R_10_b, prop_R_10_c, prop_R_10_c_5, prop_R_10_c_6, prop_R_10_c_7)
dataset_10_567_bi_annual <- as.data.frame(dataset_10_567_bi_annual)
colnames(dataset_10_567_bi_annual) <- c(
  "time",
  "Baseline", # prop_R_10_b
  "Ten", # prop_R_10_c
  "Five", # prop_R_10_c_5
  "Six", # prop_R_10_c_6
  "Seven" # prop_R_10_c_7
)
head(dataset_10_567_annual)
head(dataset_10_567_bi_annual)
library(tidyverse)
# Convert to long format
annual_long_10_567 <- dataset_10_567_annual %>%
  pivot_longer(
    cols = -time,
    names_to = "Scenario",
    values_to = "Resistance"
  )
print(annual_long_10_567)

bi_annual_long_10_567 <- dataset_10_567_bi_annual %>%
  pivot_longer(
    cols = -time,
    names_to = "Scenario",
    values_to = "Resistance"
  )
# Re-order:annual
table(annual_long_10_567$Scenario)
annual_long_10_567$Scenario <- factor(
  annual_long_10_567$Scenario,
  levels = c("Baseline", "Five", "Six", "Seven", "Ten")
)
# Re-order: bi-annual
table(bi_annual_long_10_567$Scenario)
bi_annual_long_10_567$Scenario <- factor(
  bi_annual_long_10_567$Scenario,
  levels = c("Baseline", "Five", "Six", "Seven", "Ten")
)

table(annual_long_10_567$Scenario)
annual <- ggplot(annual_long_10_567, aes(x = time, y = Resistance, color = Scenario)) +
  geom_ribbon(
    data = subset(annual_long_10_567, Scenario == "Baseline"),
    aes(
      x = time,
      ymin = 0,
      ymax = Resistance
    ),
    inherit.aes = FALSE,
    # fill = "#C2A" ,
    fill = "grey80",
    alpha = 0.4
  ) +
  # geom_point(size = 1) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(min(annual_long_10_567$time), max(annual_long_10_567$time), by = 365)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10)
  ) +
  labs(
    title = "A.Annual MDA",
    # title = "Proportion of infection due to Escherichia coli resistant to macrolides  in Tanzania",
    # subtitle = "Age-structured mixed-carriage model integrating within-host bacterial competition,E.Coli",
    x = "Time(days)",
    y = "Resistance(%)",
    color = "Scenario"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
bi_annual <- ggplot(bi_annual_long_10_567, aes(x = time, y = Resistance, color = Scenario)) +
  geom_ribbon(
    data = subset(bi_annual_long_10_567, Scenario == "Baseline"),
    aes(
      x = time,
      ymin = 0,
      ymax = Resistance
    ),
    inherit.aes = FALSE,
    # fill = "#C2A" ,
    fill = "grey80",
    alpha = 0.4
  ) +
  # geom_point(size = 1) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(min(bi_annual_long_10_567$time), max(bi_annual_long_10_567$time), by = 365)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10)
  ) +
  labs(
    title = "B.Bi-annual",
    # title = "Proportion of infection due to Escherichia coli resistant to macrolides  in Tanzania",
    # subtitle = "Age-structured mixed-carriage model integrating within-host bacterial competition,E.Coli",
    x = "Time(days)",
    y = "Resistance(%)",
    color = "Scenario"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
# Plot
print(annual)
print(bi_annual)
grid.arrange(annual, bi_annual, ncol = 2)
combined_plot <- arrangeGrob(annual, bi_annual, ncol = 2)
# Save to a file (PNG, PDF, etc.)
ggsave("Figure 1.Effect of MDA stop during  annual MDA.png", annual, width = 12, height = 10) # adjust size
ggsave("Figure 1.Effect of MDA stop during  bi-annual MDA.png", bi_annual, width = 12, height = 10) # adjust size
ggsave("Figure 1.combined_annual_and_bi_annual_plot.png", combined_plot, width = 12, height = 6) # adjust size
# Tiyplots
pacman::p_load(pak, tidyplots) # https://jbengler.github.io/tidyplots/
# Order
table(annual_long_10_567$Scenario)
table(bi_annual_long_10_567$Scenario)
# Annual
annual_long_10_567$Scenario <- factor(
  annual_long_10_567$Scenario,
  levels = c("Baseline", "Five", "Six", "Seven", "Ten")
)
# Bi-annual
bi_annual_long_10_567$Scenario <- factor(
  annual_long_10_567$Scenario,
  levels = c("Baseline", "Five", "Six", "Seven", "Ten")
)
# 1. annual visualization
baseline_value <- annual_long_10_567 |>
  dplyr::filter(Scenario == "Baseline") |>
  dplyr::summarise(mean_res = mean(Resistance, na.rm = TRUE)) |>
  dplyr::pull(mean_res)
my_colors <- c(
  "Baseline" = "#000000", # black
  "Five"     = "#1b9e77",
  "Six"      = "#d95f02",
  "Seven"    = "#7570b3",
  "Ten"      = "#e7298a"
)
graph_1_1 <- annual_long_10_567 |>
  tidyplot(x = Scenario, y = Resistance, color = Scenario) |>
  add_mean_bar(alpha = 0.8) |>
  adjust_size(width = 125, height = 105) |>
  # add_sem_errorbar() |>
  # add_data_points_beeswarm() |>
  # view_plot(title = "Default color scheme: 'friendly'") |>
  # adjust_colors(colors_discrete_apple) |>
  adjust_colors(my_colors) |>
  add_title("A. Annual MDA")
# view_plot(title = "A.Annual MDA")
print(graph_1_1)
ggsave("Graph 1.annual_MDA.png", plot = graph_1_1, width = 8, height = 6)

# 2.Bi annual visualization
graph_2_1 <- bi_annual_long_10_567 |>
  tidyplot(x = Scenario, y = Resistance, color = Scenario) |>
  add_mean_bar(alpha = 0.8) |>
  adjust_size(width = 85, height = 105) |>
  add_sem_errorbar() |>
  # add_data_points_beeswarm() |>
  # view_plot(title = "Default color scheme: 'friendly'") |>
  # adjust_colors(colors_discrete_apple) |>
  adjust_colors(my_colors) |>
  add_title("B. Annual MDA")
# view_plot(title = "B.Bi-Annual MDA")
print(graph_2_1)
ggsave("Graph 2.annual_MDA.png", plot = graph_2_1, width = 8, height = 6)
#
graph_3 <- grid.arrange(graph_1_1, graph_2_1, ncol = 2)
ggsave("Graph 3.annual_and_bi_annual_MDA.png", plot = graph_3, width = 8, height = 6)

bi_annual_long_10_567 |>
  tidyplot(x = time, y = Resistance, color = Scenario) |>
  adjust_size(width = 105, height = 105) |>
  add_mean_line() |>
  # add_mean_dot() |>
  # add_sem_ribbon()|>
  view_plot(title = "B.Bi-Annual MDA")

# 1 year
names(dataset_1)
library(tidyverse)
# Long dataset for 1-year scenarios
dataset_1_long <- dataset_1 %>%
  dplyr::select(time, prop_R_1_a, prop_R_1_b, prop_R_1_c) %>%
  tidyr::pivot_longer(
    cols = starts_with("prop_R_1_"),
    names_to = "scenario",
    values_to = "Proportion"
  )
dataset_1_long$scenario <- factor(dataset_1_long$scenario)
table(dataset_1_long$scenario)
levels(dataset_1_long$scenario)[levels(dataset_1_long$scenario) == "prop_R_1_a"] <- "1Y MDA"
levels(dataset_1_long$scenario)[levels(dataset_1_long$scenario) == "prop_R_1_b"] <- "1Y No-MDA"
levels(dataset_1_long$scenario)[levels(dataset_1_long$scenario) == "prop_R_1_c"] <- "1Y Bi-MDA"
table(dataset_1_long$scenario)
# Long dataset for 5-year scenarios
dataset_5_long <- dataset_5 %>%
  dplyr::select(time, prop_R_5_a, prop_R_5_b, prop_R_5_c) %>%
  tidyr::pivot_longer(
    cols = starts_with("prop_R_5_"),
    names_to = "scenario",
    values_to = "Proportion"
  )
# Factor
dataset_5_long$scenario <- factor(dataset_5_long$scenario)
# Rename levels
levels(dataset_5_long$scenario)[levels(dataset_5_long$scenario) == "prop_R_5_a"] <- "5Y MDA"
levels(dataset_5_long$scenario)[levels(dataset_5_long$scenario) == "prop_R_5_b"] <- "5Y No-MDA"
levels(dataset_5_long$scenario)[levels(dataset_5_long$scenario) == "prop_R_5_c"] <- "5Y Bi-MDA"
table(dataset_5_long$scenario)
# Long dataset for 10-year scenarios
dataset_10_long <- dataset_10 %>%
  dplyr::select(time, prop_R_10_a, prop_R_10_b, prop_R_10_c) %>%
  tidyr::pivot_longer(
    cols = starts_with("prop_R_10_"),
    names_to = "scenario",
    values_to = "Proportion"
  )
# Factor
dataset_10_long$scenario <- factor(dataset_10_long$scenario)
# Rename levels
levels(dataset_10_long$scenario)[levels(dataset_10_long$scenario) == "prop_R_10_a"] <- "10Y MDA"
levels(dataset_10_long$scenario)[levels(dataset_10_long$scenario) == "prop_R_10_b"] <- "10Y No-MDA"
levels(dataset_10_long$scenario)[levels(dataset_10_long$scenario) == "prop_R_10_c"] <- "10Y Bi-MDA"
# Check
table(dataset_10_long$scenario)
# Plots
table(dataset_1_long$scenario)
dataset_1_long$scenario <- factor(
  dataset_1_long$scenario,
  levels = c("1Y No-MDA", "1Y MDA", "1Y Bi-MDA")
)

table(dataset_5_long$scenario)
dataset_1_long$scenario <- factor(
  dataset_1_long$scenario,
  levels = c("1Y No-MDA", "1Y MDA", "1Y Bi-MDA")
)

dataset_5_long$scenario <- factor(
  dataset_5_long$scenario,
  levels = c("5Y No-MDA", "5Y MDA", "5Y Bi-MDA")
)


dataset_10_long$scenario <- factor(
  dataset_10_long$scenario,
  levels = c("10Y No-MDA", "10Y MDA", "10Y Bi-MDA")
)
g_0_1 <- ggplot(dataset_1_long, aes(x = time, y = Proportion, color = scenario)) +
  geom_ribbon(
    data = subset(dataset_1_long, scenario == "1Y No-MDA"),
    aes(
      x = time,
      ymin = 0,
      ymax = Proportion
    ),
    inherit.aes = FALSE,
    fill = "#C2A5CF",
    # fill = "grey80",
    alpha = 0.4
  ) +
  # geom_point(size = 1) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(min(dataset_1_long$time), max(dataset_1_long$time), by = 30)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10)
  ) +
  labs(
    title = "Proportion of infection due to Escherichia coli resistant to macrolides  in Tanzania",
    subtitle = "Age-structured mixed-carriage model integrating within-and between host bacterial competitions",
    x = "Time(days)",
    y = "Resistance(%)",
    color = "Scenario"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.text = element_text(size = 10)
  )
print(g_0_1)
# Plot
g_0_5 <- ggplot(dataset_5_long, aes(x = time, y = Proportion, color = scenario)) +
  geom_ribbon(
    data = subset(dataset_5_long, scenario == "5Y No-MDA"),
    aes(
      x = time,
      ymin = 0,
      ymax = Proportion
    ),
    inherit.aes = FALSE,
    fill = "#C2A5CF",
    alpha = 0.4
  ) +
  # geom_point(size = 1) +
  geom_line(linewidth = 1) + #
  scale_x_continuous(
    breaks = seq(min(dataset_5_long$time), max(dataset_5_long$time), by = 365)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10)
  ) +
  labs(
    title = "Proportion of infection due to Escherichia coli resistant to macrolides  in Tanzania",
    subtitle = "Age-structured mixed-carriage model integrating within-and between host bacterial competitions",
    x = "Time (days)",
    y = "Resistance (%)",
    color = "Scenario"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
# Print
print(g_0_5)
# Plot
g_0_10 <- ggplot(dataset_10_long, aes(x = time, y = Proportion, color = scenario)) +
  geom_ribbon(
    data = subset(dataset_10_long, scenario == "10Y No-MDA"),
    aes(
      x = time,
      ymin = 0,
      ymax = Proportion
    ),
    inherit.aes = FALSE,
    fill = "#C2A5CF",
    alpha = 0.4
  ) +
  # geom_point(size = 1) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(min(dataset_10_long$time), max(dataset_10_long$time), by = 365)
  ) +
  scale_y_continuous(
    limits = c(0, 40),
    breaks = seq(0, 40, by = 10)
  ) +
  labs(
    # title = "Proportion of infection due to Escherichia coli resistant to macrolides  in Tanzania",
    # subtitle = "Age-structured mixed-carriage model integrating within-and between host bacterial competitions",
    x = "Time (days)",
    y = "Resistance (%)",
    color = "Scenario"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.title.x.top = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
# Print
head(dataset_10_long, 100)
print(g_0_10)

print(b)
#
print(g_0_1)
print(g_0_5)
print(g_0_10)
# graph<-gridExtra::grid.arrange(incidence,g_0_1,g_0_5,g_0_10,nrow=2)
ggsave("Figure 1.Scenario analysis.png", plot = graph, height = 14, width = 22, dpi = 300)
ggsave("Figure 1.One year prevalence of resistance.png", plot = g_0_1, height = 8, width = 12, dpi = 300)
ggsave("Figure 1.Five years prevalence of resistance.png", plot = g_0_5, height = 8, width = 12, dpi = 300)
ggsave("Figure 1.Ten years prevalence of resistance.png", plot = g_0_10, height = 8, width = 12, dpi = 300)

g1 <- ggplot(dataset_1, aes(x = time, y = R_total_1)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +
  geom_line(color = "blue") +
  labs(
    # title = "Resitant infections  in Tanzania",
    title = "A",
    x = "Time (days)",
    y = "Resistant infections"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.text = element_text(size = 10)
  )
names(dataset_5)
tail((dataset_5[, 2]), 1)
g2 <- ggplot(dataset_5, aes(x = time, y = R_total_5)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +
  geom_line(color = "blue") +
  labs(
    # title = "Resitant infections  in Tanzania",
    title = "B",
    x = "Time (days)",
    y = "Resistant infections"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.text = element_text(size = 10)
  )
#
names(dataset_10)
tail((dataset_10[, 2]), 1)
tail((dataset_10[, 3]), 1)
tail((dataset_10[, 4]), 1)
tail((dataset_10[, 4]), 1)
g3 <- ggplot(dataset_10, aes(x = time, y = R_total_10)) +
  geom_line(color = "blue") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +
  labs(
    title = "C",
    # title = "Resitant infections  in Tanzania",
    x = "Time (days)",
    y = "Resistant infections"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.text = element_text(size = 10)
  )
#
g4 <- ggplot(dataset_1, aes(x = time)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, accuracy = 1, suffix = " B")) +
  geom_vline(
    xintercept = 365.25 / 2,
    linetype = "dotted",
    linewidth = 0.8,
    color = "gray"
  ) +
  annotate(
    "text",
    x = 365.25 / 2,
    y = 0.95 * max(dataset_1$R_total_1_c, na.rm = TRUE),
    label = "2nd MDA",
    angle = 90,
    vjust = -0.4,
    size = 3.2
  ) +
  geom_line(
    aes(
      y = R_total_1_c,
      color = "1 year Bi-MDA",
      linetype = "1 year Bi-MDA"
    ),
    linewidth = 0.9
  ) +
  geom_line(
    aes(
      y = R_total_1_a,
      color = "1 year MDA",
      linetype = "1 year MDA"
    ),
    linewidth = 0.9
  ) +
  geom_line(
    aes(
      y = R_total_1_b,
      color = "Pre-MDA",
      linetype = "Pre-MDA"
    ),
    linewidth = 0.9
  ) +
  scale_linetype_manual(
    values = c(
      "1 year Bi-MDA" = "solid",
      "1 year MDA"    = "solid",
      "Pre-MDA"       = "dashed"
    )
  ) +
  labs(
    title = "D",
    x = "Time (days)",
    y = "Cumulative resistance",
    color = "Scenario",
    linetype = "Scenario"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
g4
g5 <- ggplot(dataset_5, aes(x = time)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, accuracy = 1, suffix = " B ")) +
  geom_line(
    aes(
      y = R_total_5_c,
      color = "5 years Bi-MDA",
      linetype = "5 years Bi-MDA"
    ),
    linewidth = 0.9
  ) +
  geom_line(
    aes(
      y = R_total_5_a,
      color = "5 years MDA",
      linetype = "5 years MDA"
    ),
    linewidth = 0.9
  ) +
  geom_line(
    aes(
      y = R_total_5_b,
      color = "5 years No-MDA",
      linetype = "5 years No-MDA"
    ),
    linewidth = 0.9
  ) +
  scale_linetype_manual(
    values = c(
      "5 years Bi-MDA" = "solid",
      "5 years MDA"    = "solid",
      "5 years No-MDA" = "dashed"
    )
  ) +
  labs(
    title = "E",
    x = "Time (days)",
    y = "Cumulative resistance",
    color = "Scenario",
    linetype = "Scenario"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
g5
theme(legend.position = "right")
g6 <- ggplot(dataset_10, aes(x = time)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, accuracy = 1, suffix = " B ")) +
  geom_line(
    aes(
      y = R_total_10_c,
      color = "10 years Bi-MDA",
      linetype = "10 years Bi-MDA"
    ),
    linewidth = 0.9
  ) +
  geom_line(
    aes(
      y = R_total_10_a,
      color = "10 years MDA",
      linetype = "10 years MDA"
    ),
    linewidth = 0.9
  ) +
  geom_line(
    aes(
      y = R_total_10_b,
      color = "10 years No-MDA",
      linetype = "10 years No-MDA"
    ),
    linewidth = 0.9
  ) +
  scale_linetype_manual(
    values = c(
      "10 years Bi-MDA" = "solid",
      "10 years MDA" = "solid",
      "10 years No-MDA" = "dashed"
    )
  ) +
  labs(
    title = "F",
    x = "Time (days)",
    y = "Cumulative resistance",
    color = "Scenario",
    linetype = "Scenario"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
g6
print(g4)
print(g5)
print(g6)
g <- gridExtra::grid.arrange(g4, g5, g6, nrow = 3)
print(g)
g8 <- gridExtra::grid.arrange(g_0_1, g4, g_0_5, g5, g_0_10, g6, nrow = 3)
g8
ggsave(
  filename = "Figure 2_Scenarios_Tanzania.png",
  plot = g8,
  width = 14,
  height = 8,
  dpi = 300
)
library(grid)
library(gridExtra)
caption_text <- paste(
  "Figure 2.",
  "A: Resistant infections(1 year) ;",
  "B: Resistant infections(5 years);",
  "C: Resistant infections(10 years);",
  "D: Cumulative resistant infections (1 year);",
  "E: Cumulative resistant infections (5 years);",
  "F: Cumulative resistant infections (10 years)."
)
wrapped_caption <- paste(strwrap(caption_text, width = 120), collapse = "\n")
g8 <- grid.arrange(
  g_0_1, g4, g_0_5, g5, g_0_10, g6,
  nrow = 3,
  top = textGrob("Figure 2. Scenarios in Tanzania",
    gp = gpar(fontsize = 14, fontface = "bold")
  ),
  bottom = textGrob(
    wrapped_caption,
    gp = gpar(fontsize = 10, fontface = "italic")
  )
)
ggsave(
  filename = "Figure_2_Scenarios_Tanzania.png",
  plot = g8,
  width = 14,
  height = 9, # slightly increase height
  dpi = 300
)
pacman::p_load(grid, gridExtra)
# TInc = run_d[length(run_d[, 7]), 7]
TInc <- tail(cumsum(rowSums(out_1_b_Tanzania[, Rindex + 1])), 1)
#-------------------------------------------------------------------------------
#              Sensitivity analysis
#-------------------------------------------------------------------------------
# 1.One-at-a-time sensitivity analysis
print(parms)
names(parms)
parms_list <- list(
  beta.S = 0.164271,
  u.S = 0.03285421,
  u.R = 0.03285421,
  u.C = 0.03285421, # example
  a = 0.16,
  a.C = 0.16,
  a.use = 0.08, # Antibiotic use (~0.01–0.05:routine use,0.05 – 0.20: High use communitie,0.20 – 0.80 (short period) :MDA  )
  a.use.eff = 0.05, # Antibiotic effect: 0.005−0.05
  k = 0.5,
  m_contact = m_contact_1y_Tanzania, # contact matrix
  mda_cycle = 365,
  mda_duration = 30,
  mda_cov = 0.9,
  theta = 0.13, # add alpha and a use
  kappa = 0,
  r_mda = 0.03054302,
  c = 0.01
)
# .........All parameters.................................................
params_to_test <- c("beta.S", "u.S", "u.R", "u.C", "a", "a.C", "a.use", "a.use.eff", "k", "mda_cycle", "mda_duration", "mda_cov", "theta", "kappa", "c")
length(params_to_test)
sens_results <- NULL
start <- Sys.time()
for (param in params_to_test) {
  base_val <- parms_list[[param]]

  values <- seq(0.8 * base_val, 1.2 * base_val, length.out = 50)

  for (val in values) {
    parms <- parms_list
    parms[[param]] <- val

    out <- bacteria.solve(tvec_1_a, state, parms)
    # TInc <- tail(cumsum(rowSums(out[, Rindex + 1])), 1)
    TInc <- cumsum(rowSums(out[, Rindex + 1]))
    sens_results <- rbind(
      sens_results,
      data.frame(
        parameter = param,
        value = val,
        cumulative_infections = TInc
      )
    )
  }
}
end <- Sys.time()
print(end - start)
library(ggplot2)
# Visualization
# 1.Histogram
sens_results <- as.data.frame(sens_results)
head(sens_results)
table(sens_results$parameter)
sens_results <- sens_results |>
  filter(!parameter %in% c("kappa", "mda_cov", "mda_cycle"))

table(sens_results$parameter)
ggplot(
  sens_results,
  aes(x = cumulative_infections / 1e9)
) +
  geom_histogram(bins = 20, fill = "steelblue", color = "black") +
  facet_wrap(~parameter, scales = "free") +
  labs(
    x = "Cumulative resitant infections",
    y = "Frequency",
    title = "Local sensitivity analysis across parameters"
  ) +
  theme_bw()
# 2.Box plot
ggplot(
  sens_results,
  aes(x = parameter, y = cumulative_infections / 1e9, fill = parameter)
) +
  geom_boxplot(fill = "steelblue") +
  facet_wrap(~parameter, scales = "free_x") +
  # facet_wrap(value~parameter$u.c, scales = "free_x") +
  labs(
    x = "Parameter",
    y = "Cumulative resitant infections (Billions)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# install.packages("tidyplots")
pacman::p_load(tidyplots)
# By value
sens_results <- as.data.frame(sens_results)
head(sens_results)
# 3.Scatter plot
ggplot(
  sens_results,
  aes(x = value, y = cumulative_infections / 1e9)
) +
  geom_point(fill = "steelblue") +
  facet_wrap(~parameter, scales = "free_x") +
  labs(
    x = "Parameter value",
    y = "Cumulative resistant infections(×10^8) " # (Billions)
  ) +
  theme_bw()
# 4.Tornado plot
library(dplyr)
sens_summary <- sens_results %>%
  group_by(parameter) %>%
  summarise(
    min_inf = min(TInc),
    max_inf = max(TInc),
    impact = max_inf - min_inf
  ) %>%
  arrange(desc(impact))
head(sens_summary)
library(dplyr)
library(ggplot2)
ggplot(
  sens_summary,
  aes(
    x = reorder(parameter, impact),
    y = impact
  )
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "Parameter",
    y = "Impact on cumulative infections (Billions)",
    title = "Sensitivity Analysis (Tornado Plot)"
  ) +
  theme_bw()
# ......................Each parameters..........................................
betasens <- NULL
# Let's say you want to vary beta.S from 50% to 150% of its value
beta_values <- seq(0.8 * parms_list$beta.S, 1.2 * parms_list$beta.S, length.out = 50)

for (beta_i in beta_values) {
  parms <- parms_list # copy the full set
  parms$beta.S <- beta_i # only change beta.S
  # Run your model
  out_1_b_Tanzania <- bacteria.solve(tvec_1_a, state, parms_noMDA)
  # Calculate infections
  # TInc <- tail(rowSums(out_1_a_Tanzania[, Rindex + 1]), 1)
  # TInc <- tail(cumsum(rowSums(out_1_b_Tanzania[, Rindex + 1])), 1)
  # Store results
  betasens <- rbind(betasens, c(beta_i, TInc))
}
# Data frame for plotting
betasens <- as.data.frame(betasens)
colnames(betasens) <- c("beta.S", "Cumulative_Infections")
# Plot histogram
h_Beta.S <- hist(betasens$Cumulative_Infections, col = "skyblue", main = "Sensitivity on beta.S", ylab = "Frequency", xlab = "Total resistant infections")
head(betasens)
library(ggplot2)
ggplot(betasens, aes(x = Cumulative_Infections / 1e9)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "black") +
  labs(
    x = "Cumulative Infections (Billions)",
    y = "Frequency"
  ) +
  theme_minimal()
ggplot(betasens, aes(x = factor(beta.S), y = Cumulative_Infections / 1e9)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(
    x = expression(beta[S]),
    y = "Cumulative Infections (Billions)"
  ) +
  theme_minimal()

library(ggplot2)
p <- ggplot(betasens, aes(x = beta.S, y = Cumulative_Infections)) +
  geom_point(alpha = 0.6) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    span = 0.75,
    linewidth = 1
  ) +
  labs(
    title = "Local sensitivity of cumulative resistant infections to beta.S",
    x = expression(beta[S]),
    y = "Resistant infections"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.text = element_text(size = 10)
  )
p
#
usens <- NULL
# Vary u.S from 50% to 150% of its baseline value
uS_values <- seq(
  0.8 * parms_list$u.S,
  1.2 * parms_list$u.S,
  length.out = 50
)
for (uS_i in uS_values) {
  parms <- parms_list # Baseline parameters
  parms$u.S <- uS_i # Here i change only u.S
  # Run the model
  out_1_b_Tanzania <- bacteria.solve(tvec_1_a, state, parms_noMDA)
  # Cumulative resistant infections
  TInc <- tail(cumsum(rowSums(out_1_b_Tanzania[, Rindex + 1])), 1)
  # Store results
  usens <- rbind(usens, c(uS_i, TInc))
}
usens <- as.data.frame(usens)
colnames(usens) <- c("u.S", "Cumulative_Infections")
h_uS <- hist(usens$Cumulative_Infections,
  col  = "skyblue",
  main = "Sensitivity on u.S",
  xlab = "Cumulative resistant infections",
  ylab = "Frequency"
)
# Plot sensitivity analysis results
par(mfrow = c(1, 2))
h_Beta.S
h_uS
print(h_uS)
ggsave(
  filename = "Figure_S1_Sensitivity_Beta.S_Tanzania_LOESS.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)
# Global sensitivity analysis using LHS
pacman::p_load(lhs, foreach, doParallel)
library(lhs)
library(foreach)
library(doParallel)
library(deSolve)
# Ensure parms_list has all parameters
parms_list <- list(
  beta.S = 0.164271,
  u.S = 0.03285421,
  u.R = 0.03285421,
  u.c = 0.03285421, # <- MUST exist
  a = 0.16,
  a.C = 0.16,
  k = 0.5,
  theta = 0.13,
  r_mda = 0.03054302,
  c = 0.01,
  m_contact = m_contact_1y_Tanzania,
  mda_cycle = 365,
  mda_duration = 30,
  mda_cov = 0.6,
  kappa = 0
)
# Parameter ranges
param_ranges <- list(
  beta.S = c(0.5 * parms_list$beta.S, 1.5 * parms_list$beta.S),
  u.S = c(0.5 * parms_list$u.S, 1.5 * parms_list$u.S),
  u.R = c(0.5 * parms_list$u.R, 1.5 * parms_list$u.R),
  a = c(0.5 * parms_list$a, 1.5 * parms_list$a),
  k = c(0.5 * parms_list$k, 1.5 * parms_list$k),
  theta = c(0.5 * parms_list$theta, 1.5 * parms_list$theta),
  r_mda = c(0.5 * parms_list$r_mda, 1.5 * parms_list$r_mda),
  c = c(0.5 * parms_list$c, 1.5 * parms_list$c)
)
param_ranges
# LHS sampling
n_sim <- 50
lhs_samples <- randomLHS(n_sim, length(param_ranges))
param_samples <- as.data.frame(sapply(seq_along(param_ranges), function(i) {
  x <- param_ranges[[i]]
  lhs_samples[, i] * (x[2] - x[1]) + x[1]
}))
colnames(param_samples) <- names(param_ranges)
# Stop any existing cluster
# Parallel setup
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
# Export everything needed to workers
results <- foreach(
  i = 1:n_sim, .combine = c,
  .packages = "deSolve"
) %dopar% {
  parms_noMDA <- parms_list
  for (p in names(param_ranges)) parms_noMDA[[p]] <- param_samples[i, p]
  out <- bacteria.solve(tvec_1_a, state, parms_noMDA)
  tail(cumsum(rowSums(out[, Rindex + 1])), 1)
}
# Stop cluster when done
stopCluster(cl)
# Stop cluster
stopCluster(cl)
# Combine parameter samples and results
data <- cbind(param_samples, Cumulative_Infections = results)
# Fit linear model to get coefficients
lm_fit <- lm(Cumulative_Infections ~ ., data = data)
summary(lm_fit)
# Rank parameters by influence
coef_df <- as.data.frame(summary(lm_fit)$coefficients)
coef_df <- coef_df[order(abs(coef_df$Estimate), decreasing = TRUE), ]
print(coef_df)
library(ggplot2)
# Intercept for clarity (optional)
coef_plot <- coef_df[rownames(coef_df) != "(Intercept)", ]
# Make a bar plot sorted by absolute effect
ggplot(coef_plot, aes(
  x = reorder(rownames(coef_plot), abs(Estimate)),
  y = Estimate,
  fill = Estimate
)) +
  geom_bar(stat = "identity") +
  coord_flip() + # horizontal bars
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Parameter", y = "Coefficient estimate", fill = "Estimate") +
  theme_minimal() +
  theme(text = element_text(size = 12))
# ...............................................................................
c_vector <- round(seq(0, 1.5, by = 0.25) * parms.orig[["c"]], 4)
k_vector <- round(seq(0, 1.5, by = 0.25) * parms.orig[["k"]], 4)
beta_vector <- seq(0, 1.5, by = 0.25) * parms.orig[["beta.S"]]
result <- matrix(0, nrow = length(c_vector), length(k_vector))
result_beta <- matrix(0, nrow = length(beta_vector), length(k_vector))
c_vector
k_vector
result
result_beta
# A.Loop over C and K values
for (i in 1:length(c_vector)) {
  for (j in 1:length(k_vector)) {
    # for (j in 1:length(beta_vector)){
    parms <- parms.orig # copy baseline
    parms[["c"]] <- c_vector[i] # update c
    parms[["k"]] <- k_vector[j] # update k
    # parms[["beta.S"]] <- beta_vector[j]  # update beta


    (mda_start_times <- (0:50) * 365.25) # Early MDA start
    out_1_a_Tanzania <- bacteria.solve(tvec_1_a, state, parms)

    # (mda_start_times<-(0:50)*365.25)             # Delay MDA start
    # out_1_b_Tanzania <- bacteria.solve(tvec_1_b, state, parms_noMDA)

    # (mda_start_times<-(0:50)*365.25/2)             # Delay MDA start
    # out_1_c_Tanzania <- bacteria.solve(tvec_1_a, state, parms)


    # (mda_start_times<-(0:50)*365.25)             # Early MDA start
    # out_5_a_Tanzania <- bacteria.solve(tvec_5_a, state, parms)

    # (mda_start_times<-(6:50)*365.25)             # Delay MDA start
    # out_5_b_Tanzania <- bacteria.solve(tvec_5_b, state, parms)

    # (mda_start_times<-(0:50)*365.25)             # Early MDA start
    # out_10_a_Tanzania <- bacteria.solve(tvec_10_a, state, parms)

    # (mda_start_times<-(51:52)*365.25)             # Delay MDA start
    # out_10_b_Tanzania <- bacteria.solve(tvec_10_b, state, parms)
    # Store R at equilibrium
    result[i, j] <- tail(((rowSums(out_1_a_Tanzania[, Rindex + 1])) * 100) / (rowSums(out_1_a_Tanzania[, Sindex + 1]) + rowSums(out_1_a_Tanzania[, Rindex + 1]) + rowSums(out_1_a_Tanzania[, Rsindex + 1]) + rowSums(out_1_a_Tanzania[, Srindex + 1]) + rowSums(out_1_a_Tanzania[, Xindex + 1])), 1)
    # result[i,j] <- tail(((rowSums(out_1_b_Tanzania[, Rindex + 1]))*100)  /(rowSums(out_1_b_Tanzania[, Sindex + 1])+rowSums(out_1_b_Tanzania[, Rindex + 1])+rowSums(out_1_b_Tanzania[, Rsindex + 1])+rowSums(out_1_b_Tanzania[, Srindex + 1])+rowSums(out_1_b_Tanzania[, Xindex + 1])), 1)
    # result[i,j] <- tail(((rowSums(out_5_a_Tanzania[, Rindex + 1]))*100) / (rowSums(out_5_a_Tanzania[, Sindex + 1]) + rowSums(out_5_a_Tanzania[, Rindex + 1]) + rowSums(out_5_a_Tanzania[, Rsindex + 1]) + rowSums(out_5_a_Tanzania[, Srindex + 1]) + rowSums(out_5_a_Tanzania[, Xindex + 1])), 1)
    # result[i,j] <- tail(((rowSums(out_5_b_Tanzania[, Rindex + 1]))*100) / (rowSums(out_5_b_Tanzania[, Sindex + 1]) + rowSums(out_5_b_Tanzania[, Rindex + 1]) + rowSums(out_5_b_Tanzania[, Rsindex + 1]) + rowSums(out_5_b_Tanzania[, Srindex + 1]) + rowSums(out_5_b_Tanzania[, Xindex + 1])), 1)
    # result[i,j] <- tail(((rowSums(out_10_a_Tanzania[, Rindex + 1]))*100)/ (rowSums(out_10_a_Tanzania[, Sindex + 1]) + rowSums(out_10_a_Tanzania[, Rindex + 1]) + rowSums(out_10_a_Tanzania[, Rsindex + 1]) + rowSums(out_10_a_Tanzania[, Srindex + 1]) + rowSums(out_10_a_Tanzania[, Xindex + 1])), 1)
    # result[i,j] <- tail(((rowSums(out_10_b_Tanzania[, Rindex + 1]))*100)/ (rowSums(out_10_b_Tanzania[, Sindex + 1]) + rowSums(out_10_b_Tanzania[, Rindex + 1]) + rowSums(out_10_b_Tanzania[, Rsindex + 1]) + rowSums(out_10_b_Tanzania[, Srindex + 1]) + rowSums(out_10_b_Tanzania[, Xindex + 1])), 1)
    # This are results for beta
    # result_beta[i,j] <- tail(rowSums(out_10_a_Tanzania[, Rindex + 1]), 1)
  }
}
# Here i will store R at equilibrium in result[I,j]
contour(result, xaxt = "n", yaxt = "n", ylab = "k", xlab = "c", nlevels = 10, col = hcl.colors(14, "Temps"), lwd = 2)
axis(1, at = 1:length(c_vector) / length(c_vector), labels = c_vector)
axis(2, at = 1:length(k_vector) / length(k_vector), labels = k_vector)
# Advanced visualization
pacman::p_load(ggplot2, reshape2, viridis)
# Convert result matrix to long format
df <- melt(result)
# df <- melt(result_beta)
names(df) <- c("k_index", "c_index", "R_equilibrium")
# names(df) <- c("beta_index", "c_index", "R_equilibrium")
# Map indices to actual parameter values
df$c <- c_vector[df$c_index]
df$k <- k_vector[df$k_index]
# df$beta <- beta_vector[df$beta_index]
# Heatmap with custom scales
ggplot(df, aes(x = c, y = k, fill = R_equilibrium)) +
  geom_tile(color = "white") + # grid
  scale_fill_viridis_c(option = "plasma", name = "Resitance(%)") +
  geom_text(aes(label = round(R_equilibrium, 2)), size = 3, color = "white") +
  scale_x_continuous(breaks = c_vector) + # custom x-axis ticks
  scale_y_continuous(breaks = k_vector) + # custom y-axis ticks
  labs(
    title = "Prevalence of E.Coli infections resistant to Azithromycin",
    subtitle = "A.One year,annual MDA",
    x = "Cost of resistance (c)",
    y = "Co-colonisation efficiency (k)"
  ) +
  theme_classic() +
  # theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
# B. Loop over C and Beta
c_vector <- seq(0.1, 1, by = 0.25) * parms.orig[["c"]]
beta_vector <- seq(0.1, 1, by = 0.25) * parms.orig[["beta.S"]]
result_beta <- matrix(0, nrow = length(beta_vector), length(k_vector))
c_vector
beta_vector
result_beta

for (i in 1:length(c_vector)) {
  for (j in 1:length(beta_vector)) {
    parms <- parms.orig # copy baseline
    parms[["c"]] <- c_vector[i] # update c
    # parms[["k"]] <- k_vector[j]  # update k
    parms[["beta.S"]] <- beta_vector[j] # update beta

    (mda_start_times <- (0:50) * 365.25) # Early MDA start
    out_1_a_Tanzania <- bacteria.solve(tvec_1_a, state, parms)

    # (mda_start_times<-(2:50)*365.25)             # Delay MDA start
    # out_1_b_Tanzania <- bacteria.solve(tvec_1_b, state, parms)

    # (mda_start_times<-(0:50)*365.25)             # Early MDA start
    # out_5_a_Tanzania <- bacteria.solve(tvec_5_a, state, parms)

    # (mda_start_times<-(6:50)*365.25)             # Delay MDA start
    # out_5_b_Tanzania <- bacteria.solve(tvec_5_b,state, parms)

    # (mda_start_times<-(0:50)*365.25)              # Early MDA start
    # out_10_a_Tanzania <- bacteria.solve(tvec_10_a, state, parms)

    # (mda_start_times<-(51:52)*365.25)             # Delay MDA start
    # out_10_b_Tanzania <- bacteria.solve(tvec_10_b, state, parms)

    # Store R at equilibrium:beta and c
    result_beta[i, j] <- tail(((rowSums(out_1_a_Tanzania[, Rindex + 1]) + rowSums(out_1_a_Tanzania[, Rsindex + 1])) * 100) / (rowSums(out_1_a_Tanzania[, Sindex + 1]) + rowSums(out_1_a_Tanzania[, Rindex + 1]) + rowSums(out_1_a_Tanzania[, Rsindex + 1]) + rowSums(out_1_a_Tanzania[, Srindex + 1]) + rowSums(out_1_a_Tanzania[, Xindex + 1])), 1)
    # result_beta[i,j] <- tail(((rowSums(out_1_b_Tanzania[, Rindex + 1])+rowSums(out_1_b_Tanzania[, Rsindex + 1]))*100)/(rowSums(out_1_b_Tanzania[, Sindex + 1])+rowSums(out_1_b_Tanzania[, Rindex + 1])+rowSums(out_1_b_Tanzania[, Rsindex + 1])+rowSums(out_1_b_Tanzania[, Srindex + 1])+rowSums(out_1_b_Tanzania[, Xindex + 1])), 1)
    # result_beta[i,j] <- tail(((rowSums(out_5_a_Tanzania[, Rindex + 1]) + rowSums(out_5_a_Tanzania[, Rsindex + 1])) * 100) / (rowSums(out_5_a_Tanzania[, Sindex + 1]) + rowSums(out_5_a_Tanzania[, Rindex + 1]) + rowSums(out_5_a_Tanzania[, Rsindex + 1]) + rowSums(out_5_a_Tanzania[, Srindex + 1]) + rowSums(out_5_a_Tanzania[, Xindex + 1])), 1)
    # result_beta[i,j] <- tail(((rowSums(out_5_b_Tanzania[, Rindex + 1]) + rowSums(out_5_b_Tanzania[, Rsindex + 1])) * 100) / (rowSums(out_5_b_Tanzania[, Sindex + 1]) + rowSums(out_5_b_Tanzania[, Rindex + 1]) + rowSums(out_5_b_Tanzania[, Rsindex + 1]) + rowSums(out_5_b_Tanzania[, Srindex + 1]) + rowSums(out_5_b_Tanzania[, Xindex + 1])), 1)
    # result_beta[i,j] <- tail(((rowSums(out_10_a_Tanzania[, Rindex + 1]) + rowSums(out_10_a_Tanzania[, Rsindex + 1])) * 100) / (rowSums(out_10_a_Tanzania[, Sindex + 1]) + rowSums(out_10_a_Tanzania[, Rindex + 1]) + rowSums(out_10_a_Tanzania[, Rsindex + 1]) + rowSums(out_10_a_Tanzania[, Srindex + 1]) + rowSums(out_10_a_Tanzania[, Xindex + 1])), 1)
    # result_beta[i,j] <- tail(((rowSums(out_10_b_Tanzania[, Rindex + 1]) + rowSums(out_10_b_Tanzania[, Rsindex + 1])) * 100) / (rowSums(out_10_b_Tanzania[, Sindex + 1]) + rowSums(out_10_b_Tanzania[, Rindex + 1]) + rowSums(out_10_b_Tanzania[, Rsindex + 1]) + rowSums(out_10_b_Tanzania[, Srindex + 1]) + rowSums(out_10_b_Tanzania[, Xindex + 1])), 1)
  }
}
# Here i wil store R at equilibrium in result[I,j]
contour(result_beta, xaxt = "n", yaxt = "n", ylab = "k", xlab = "c", nlevels = 10, col = hcl.colors(14, "Temps"), lwd = 2)
axis(1, at = 1:length(c_vector) / length(c_vector), labels = c_vector)
axis(2, at = 1:length(beta_vector) / length(beta_vector), labels = beta_vector)
# Advanced visualization
pacman::p_load(ggplot2, reshape2, viridis)
# Convert result matrix to long format
df_beta <- melt(result_beta)
names(df_beta) <- c("beta_index", "c_index", "R_beta_equilibrium")
# Map indices to actual parameter values
df_beta$c <- c_vector[df_beta$c_index]
df_beta$beta <- beta_vector[df_beta$beta_index]
# Heatmap with custom scales
ggplot(df_beta, aes(x = c, y = beta, fill = R_beta_equilibrium)) +
  geom_tile(color = "white") + # grid
  scale_fill_viridis_c(option = "viridis", name = "Resitance(%)") +
  geom_text(aes(label = round(R_beta_equilibrium, 2)), size = 3, color = "white") +
  scale_x_continuous(breaks = c_vector) + # custom x-axis ticks
  scale_y_continuous(breaks = beta_vector) + # custom y-axis ticks
  labs(
    title = "Prevalence of E.Coli infections resistant to Azithromycin ",
    subtitle = "B.One year,annual MDA",
    x = "Cost of resistance (c)",
    y = "Transmission rate (Beta.S)"
  ) +
  theme_classic() +
  # theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
# ...............................................................................
R_total_10_a <- (rowSums(out_10_a_Tanzania[, Rindex + 1]) + rowSums(out_10_a_Tanzania[, Rsindex + 1]))
R_total_10_b <- (rowSums(out_10_b_Tanzania[, Rindex + 1]) + rowSums(out_10_b_Tanzania[, Rsindex + 1]))
Rs_total_10_a <- rowSums(out_10_a_Tanzania[, Rsindex + 1])
Rs_total_10_b <- rowSums(out_10_b_Tanzania[, Rsindex + 1])
# Visualization
plot(
  out_1_b_Tanzania$time, X_total,
  type = "l",
  col = cols[1],
  las = 1,
  xaxs = "i", yaxs = "i",
  ylim = c(0, max(c(X_total, S_total, R_total, Rs_total, Sr_total))),
  bty = "n",
  lwd = 3.5,
  xlab = "Day", ylab = "Population (in Millions)",
  # xlab = "Day", ylab = "Proportion",
  main = "Pre-MDA population level Bacterial Dynamics in Tanzania"
)
lines(out_1_b_Tanzania$time, S_total, col = cols[2], lwd = 3.5)
lines(out_1_b_Tanzania$time, R_total, col = cols[3], lwd = 3.5)
lines(out_1_b_Tanzania$time, Rs_total, col = cols[4], lwd = 3.5)
lines(out_1_b_Tanzania$time, Sr_total, col = cols[5], lwd = 3.5)
legend(
  "topright",
  bty = "n",
  col = c(cols),
  # legend = c("X_total", "S_total", "R_total","Rs_total", "Sr_total"),
  legend = c("X_total", "S_total", "R_total", "Rs_total", "Sr_total"),
  lty = 1,
  lwd = 3.5,
  ncol = 1
)
# Visualization for the ten years
# a.Visualization for the 10 years R
plot(
  out_10_a_Tanzania$time, R_total_10_a,
  type = "l",
  col = cols[1],
  las = 1,
  xaxs = "i", yaxs = "i",
  ylim = c(0, max(c(R_total_10_a, R_total_10_b))),
  bty = "n",
  lwd = 3.5,
  xlab = "Day", ylab = "Population (R) in millions",
  # xlab = "Day", ylab = "Proportion (R)",
  main = "Post 10 Years MDA population level Bacterial Dynamics in Tanzania"
)
abline(v = 365 * 1, col = "red", lwd = 2, lty = 2)
# Add text above or next to the vertical line
text(
  x = 3000,
  y = max(c(R_total_10_a, R_total_10_b)) * 0.90, # vertical position (95% of ymax)
  labels = "MDA start",
  pos = 3, # 3 = above the point (can use 1,2,4)
  cex = 1,
  col = "red"
)
# lines(out_10_a_Tanzania$time, Rs_total_10_a, col = cols[2], lwd = 3.5)
lines(out_10_b_Tanzania$time, R_total_10_b, col = cols[4], lwd = 3.5)
# lines(out_10_b_Tanzania$time, Rs_total_10_b, col = cols[4], lwd = 3.5)
legend(
  "topright",
  bty = "n",
  col = c(cols[c(1, 4)]),
  # col = c("red","blue"),
  legend = c("MDA", "No MDA"),
  lty = 1,
  lwd = 3.5,
  ncol = 2
)
# b.Visualization for the 10 years Rs
plot(
  out_10_b_Tanzania$time, Rs_total_10_a,
  type = "l",
  col = cols[1],
  las = 1,
  xaxs = "i", yaxs = "i",
  ylim = c(0, 2), # max(c(Rs_total_10_a, Rs_total_10_b))),
  bty = "n",
  lwd = 3.5,
  # xlab = "Day", ylab = "Proportion(Rs)",
  xlab = "Day", ylab = "Population(Rs) in millions",
  main = "Population level Bacterial Dynamics over time"
)
abline(v = 365 * 1, col = "red", lwd = 2, lty = 2)
# lines(out_10_a_Tanzania$time, Rs_total_10_a, col = cols[2], lwd = 3.5)
lines(out_10_b_Tanzania$time, R_total_10_b, col = cols[3], lwd = 3.5)
# lines(out_10_b_Tanzania$time, Rs_total_10_b, col = cols[4], lwd = 3.5)
legend(
  "topright",
  bty = "n",
  col = c(cols[c(1, 3)]),
  legend = c("MDA", "No MDA"),
  lty = 1,
  lwd = 3.5,
  ncol = 2
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#              Additional-visualization                                  #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
pacman::p_load(data.table) # This package will allow us to reshape data set faster
results_1_Tanzania <- as.data.table(results_1_b_Tanzania)
head(results_1_Tanzania)
# results_1_Tanzania <- results_1_Tanzania[, 1:(5 * n_age + 1)]
head(results_1_Tanzania) # I will need to remove X as the resistance is calculated on colonized only
library(dplyr)
cols_keep <- grepl("^(time|R_|Rs_|S_|Sr_)", names(results_1_Tanzania))
results_1_Tanzania <- results_1_Tanzania[, ..cols_keep] # .. need to be there
head(results_1_Tanzania)

# results_1_Tanzania<- as.data.table(results_1_b_Tanzania)
# results_1_Tanzania<- as.data.table(results_1_c_Tanzania)
# Long format : Here i will be using melt to be faster
results_1_Tanzania
results_1_Tanzania_long <- melt(
  results_1_Tanzania,
  id.vars = "time",
  variable.name = "variable",
  value.name = "value"
)

head(results_1_Tanzania_long)
table(results_1_Tanzania_long$variable)
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
# Structure of the long format
names(results_1_Tanzania_long)
head(results_1_Tanzania_long)
table(results_1_Tanzania_long$compartment)
table(results_1_Tanzania_long$age_group)
# Here i will need to speed the visualization
pacman::p_load(data.table, ggplot2, scales)
# Data.table
setDT(results_1_Tanzania_long)
# Transformations in one data.table
# Proportions (in one step)
results_1_Tanzania_long[as.integer(age_group) <= 101, `:=`(
  total_by_age = sum(value),
  proportion = value / sum(value) * 100
), by = age_group]
# second checks
names(results_1_Tanzania_long)
head(results_1_Tanzania_long)
table(results_1_Tanzania_long$compartment)
table(results_1_Tanzania_long$age_group)
table(results_1_Tanzania_long$compartment)
R_RS_only <- results_1_Tanzania_long %>%
  filter(compartment %in% c("R", "Rs")) # Rs
R_RS_only$proportion
plot_a_0 <- ggplot(
  R_RS_only,
  aes(x = age_group, y = proportion, fill = compartment)
) +
  geom_col(position = "stack", col = NA) + # or dodge
  geom_hline(yintercept = 18.2, linetype = "dashed", color = "red", size = 1) +
  # or use color = "white"
  scale_y_continuous(
    limits = c(0, 60),
    breaks = seq(0, 60, by = 20),
    expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 5)) +
  labs(title = " Population level Bacterial Dynamics over time in Tanzania", x = "Age", y = "Percent ") +
  scale_fill_manual(values = c(
    "X" = "#00c4aa", "S" = "#e573f3",
    "R" = "#00b3f4", "Sr" = "#9b9602", "Rs" = "#fc726c"
  ))
print(plot_a_0)
plot_a_0 <- ggplot(
  R_RS_only,
  aes(x = proportion, y = age_group, fill = compartment)
) +
  geom_col(position = "stack", col = NA) + # or dodge
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = " Population level Bacterial Dynamics over time in Tanzania", x = "Age", y = "Percent ") +
  scale_fill_manual(values = c(
    "X" = "#00c4aa", "S" = "#e573f3",
    "R" = "#00b3f4", "Sr" = "#9b9602", "Rs" = "#fc726c"
  ))
# print(plot_a_0)
# Pre-aggregate data
Final_1_Tanzania_summary <- results_1_Tanzania_long[
  as.integer(age_group) <= 101, # Exclude 100 + age
  .(total_value = sum(value)),
  by = .(age_group, compartment)
][, proportion := total_value / sum(total_value) * 100, by = age_group]
# Plot-aggregated data
names(Final_1_Tanzania_summary)
dim(Final_1_Tanzania_summary)
table(Final_1_Tanzania_summary$compartment)
plot_a <- ggplot(
  Final_1_Tanzania_summary,
  aes(x = age_group, y = proportion, fill = compartment)
) +
  geom_col(position = "stack", col = NA) + # or dodge
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = " Population level Bacterial Dynamics over time in Tanzania", x = "Age", y = "Percent ") +
  scale_fill_manual(values = c(
    "X" = "#00c4aa", "S" = "#e573f3",
    "R" = "#00b3f4", "Sr" = "#9b9602", "Rs" = "#fc726c"
  ))
print(plot_a)
plot_b <- ggplot(
  Final_1_Tanzania_summary,
  aes(x = age_group, y = proportion, fill = compartment)
) +
  geom_col(position = "dodge", col = NA) + # or dodge
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
  labs(title = " Population level Bacterial Dynamics over time in Tanzania", x = "Age", y = "Percent ") +
  scale_fill_manual(values = c(
    "X" = "#00c4aa", "S" = "#e573f3",
    "R" = "#00b3f4", "Sr" = "#9b9602", "Rs" = "#fc726c"
  ))
print(plot_b)
table(Final_1_Tanzania_summary$compartment)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#  MDA visualization                                                           #
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
  compartment_names <- c("S", "R", "Sr", "Rs")
  col_names <- c("time")
  for (comp in compartment_names) {
    for (age in age_groups) {
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
# results_1_a_Tanzania<-results_1_a_Tanzania[,1:(5*n_age+1)]
# results_1_b_Tanzania<-results_1_b_Tanzania[,1:(5*n_age+1)]
# results_1_c_Tanzania<-results_1_c_Tanzania[,1:(5*n_age+1)]

cols_keep <- grepl("^(time|R_|Rs_|S_|Sr_)", names(results_1_a_Tanzania))
results_1_a_Tanzania <- results_1_a_Tanzania[, ..cols_keep] # .. need to be there
head(results_1_a_Tanzania)

cols_keep <- grep("^(time$|S_|Sr_|R_|Rs_)", names(results_1_a_Tanzania), value = TRUE)
results_1_a_Tanzania <- results_1_a_Tanzania[, cols_keep]
head(results_1_a_Tanzania)
results_1_b_Tanzania <- results_1_b_Tanzania[, cols_keep]
results_1_c_Tanzania <- results_1_c_Tanzania[, cols_keep]

# results_5_a_Tanzania<-results_5_a_Tanzania[,1:(5*n_age+1)]
# results_5_b_Tanzania<-results_5_b_Tanzania[,1:(5*n_age+1)]
# results_5_c_Tanzania<-results_5_c_Tanzania[,1:(5*n_age+1)]
results_5_a_Tanzania <- results_5_a_Tanzania[, cols_keep]
results_5_b_Tanzania <- results_5_b_Tanzania[, cols_keep]
results_5_c_Tanzania <- results_5_c_Tanzania[, cols_keep]

# results_10_a_Tanzania<-results_10_a_Tanzania[,1:(5*n_age+1)]
# results_10_b_Tanzania<-results_10_b_Tanzania[,1:(5*n_age+1)]
# results_10_c_Tanzania<-results_10_c_Tanzania[,1:(5*n_age+1)]
results_10_a_Tanzania <- results_10_a_Tanzania[, cols_keep]
results_10_b_Tanzania <- results_10_b_Tanzania[, cols_keep]
results_10_c_Tanzania <- results_10_c_Tanzania[, cols_keep]

scenario_1yr_no_MDA_Tanzania <- process_scenario(results_1_b_Tanzania, "1Y No-MDA")
scenario_1yr_MDA_Tanzania <- process_scenario(results_1_a_Tanzania, "1Y MDA")
scenario_1yr_Bi_MDA_Tanzania <- process_scenario(results_1_c_Tanzania, "1Y Bi-MDA")

scenario_5yr_no_MDA_Tanzania <- process_scenario(results_5_b_Tanzania, "5Y No-MDA")
scenario_5yr_MDA_Tanzania <- process_scenario(results_5_a_Tanzania, "5Y MDA")
scenario_5yr_Bi_MDA_Tanzania <- process_scenario(results_5_c_Tanzania, "5Y Bi-MDA")

scenario_10yr_no_MDA_Tanzania <- process_scenario(results_10_b_Tanzania, "10Y No-MDA")
scenario_10yr_MDA_Tanzania <- process_scenario(results_10_a_Tanzania, "10Y MDA")
scenario_10yr_Bi_MDA_Tanzania <- process_scenario(results_10_c_Tanzania, "10Y Bi-MDA")
# Scenario labels
scenario_1yr_no_MDA_Tanzania$scenario <- "1Y No-MDA"
scenario_1yr_MDA_Tanzania$scenario <- "1Y MDA"
scenario_1yr_Bi_MDA_Tanzania$scenario <- "1Y Bi-MDA"

scenario_5yr_no_MDA_Tanzania$scenario <- "5Y No-MDA"
scenario_5yr_MDA_Tanzania$scenario <- "5Y MDA"
scenario_5yr_Bi_MDA_Tanzania$scenario <- "5Y Bi-MDA"

scenario_10yr_no_MDA_Tanzania$scenario <- "10Y No-MDA"
scenario_10yr_MDA_Tanzania$scenario <- "10Y MDA"
scenario_10yr_Bi_MDA_Tanzania$scenario <- "10Y Bi-MDA"
#
pacman::p_load(dplyr)
scenario_1yr_no_MDA_Tanzania <- mutate(scenario_1yr_no_MDA_Tanzania, scenario = "1Y No-MDA")
scenario_1yr_MDA_Tanzania <- mutate(scenario_1yr_MDA_Tanzania, scenario = "1Y MDA")
scenario_1yr_Bi_MDA_Tanzania <- mutate(scenario_1yr_Bi_MDA_Tanzania, scenario = "1Y Bi-MDA")

scenario_5yr_no_MDA_Tanzania <- mutate(scenario_5yr_no_MDA_Tanzania, scenario = "5Y No-MDA")
scenario_5yr_MDA_Tanzania <- mutate(scenario_5yr_MDA_Tanzania, scenario = "5Y MDA")
scenario_5yr_Bi_MDA_Tanzania <- mutate(scenario_5yr_Bi_MDA_Tanzania, scenario = "5Y Bi-MDA")

scenario_10yr_no_MDA_Tanzania <- mutate(scenario_10yr_no_MDA_Tanzania, scenario = "10Y No-MDA")
scenario_10yr_MDA_Tanzania <- mutate(scenario_10yr_MDA_Tanzania, scenario = "10Y MDA")
scenario_10yr_Bi_MDA_Tanzania <- mutate(scenario_10yr_Bi_MDA_Tanzania, scenario = "10Y Bi-MDA")

plot_c <- ggplot(
  scenario_10yr_no_MDA_Tanzania,
  aes(x = age_group, y = proportion, fill = compartment)
) +
  geom_col(position = "dodge", col = NA) + # or dodge
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  labs(title = " Population level Bacterial Dynamics over time in Tanzania", x = "Age", y = "Percent ") +
  scale_fill_manual(values = c(
    "X" = "#00c4aa", "S" = "#e573f3",
    "R" = "#00b3f4", "Sr" = "#9b9602", "Rs" = "#fc726c"
  ))
print(plot_c)
# All scenarios
all_scenarios_Tanzania <- rbind(scenario_1yr_no_MDA_Tanzania, scenario_1yr_MDA_Tanzania, scenario_1yr_Bi_MDA_Tanzania, scenario_5yr_no_MDA_Tanzania, scenario_5yr_MDA_Tanzania, scenario_5yr_Bi_MDA_Tanzania, scenario_10yr_no_MDA_Tanzania, scenario_10yr_MDA_Tanzania, scenario_10yr_Bi_MDA_Tanzania)
table(all_scenarios_Tanzania$scenario)
# Factor levels for scenarios
# all_scenarios_Tanzania$scenario <- factor(all_scenarios_Tanzania$scenario,
#                                          levels = c("Pre-MDA", "5 Years MDA", "10 Years MDA","10 Years No MDA"))
all_scenarios_Tanzania$scenario <- factor(
  all_scenarios_Tanzania$scenario,
  levels = c(
    "1Y Bi-MDA",
    "1Y MDA",
    "1Y No-MDA",
    "5Y Bi-MDA",
    "5Y MDA",
    "5Y No-MDA",
    "10Y Bi-MDA",
    "10Y MDA",
    "10Y No-MDA"
  )
)
table(all_scenarios_Tanzania$scenario)
names(all_scenarios_Tanzania)
table(all_scenarios_Tanzania$age_group)
# I want to remove some scenarios and  compartments for better visualization
all_scenarios_Tanzania_10 <- all_scenarios_Tanzania |>
  filter(scenario %in% c("10Y Bi-MDA", "10Y MDA", "10Y No-MDA")) |>
  filter(compartment %in% c("R", "Rs"))
table(all_scenarios_Tanzania_10$scenario)
table(all_scenarios_Tanzania$scenario)

all_scenarios_Tanzania_5 <- all_scenarios_Tanzania |>
  filter(scenario %in% c("5Y Bi-MDA", "5Y MDA", "5Y No-MDA")) |>
  filter(compartment %in% c("R", "Rs"))
table(all_scenarios_Tanzania_5$scenario)

all_scenarios_Tanzania_1 <- all_scenarios_Tanzania |>
  filter(scenario %in% c("1Y Bi-MDA", "1Y MDA", "1Y No-MDA")) |>
  filter(compartment %in% c("R", "Rs"))
table(all_scenarios_Tanzania_1$scenario)
# Plot 1: Stacked bar plots comparing all three scenarios
table(all_scenarios_Tanzania$age_group, all_scenarios_Tanzania$compartment)
plot_dodged_1 <- ggplot(
  all_scenarios_Tanzania_1,
  aes(x = age_group, y = proportion, fill = compartment)
) +
  geom_col(position = "stack", col = NA) +
  geom_hline(yintercept = 23, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~scenario, ncol = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  ) +
  labs(
    # title = "Population level Bacterial Dynamics over time in Tanzania",
    title = "A",
    subtitle = "",
    x = "Age",
    y = "Percent",
    fill = "Compartment"
  ) +
  scale_fill_manual(values = c(
    "X" = "#00c4aa",
    "S" = "#e573f3",
    "R" = "#fc726c",
    "Sr" = "#9b9602",
    "Rs" = "#00b3f4"
  )) +
  scale_y_continuous(limits = c(0, 100))
print(plot_dodged_1)
#
plot_dodged_5 <- ggplot(
  all_scenarios_Tanzania_5,
  aes(x = age_group, y = proportion, fill = compartment)
) +
  geom_col(position = "stack", col = NA) +
  geom_hline(yintercept = 23, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~scenario, ncol = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 3),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  ) +
  labs(
    title = "B.",
    subtitle = "",
    x = "Age",
    y = "Percent",
    fill = "Compartment"
  ) +
  scale_fill_manual(values = c(
    "X" = "#00c4aa",
    "S" = "#e573f3",
    "R" = "#fc726c",
    "Sr" = "#9b9602",
    "Rs" = "#00b3f4"
  )) +
  scale_y_continuous(limits = c(0, 100))
print(plot_dodged_5)
# Plot 2: Dodged bar plots comparing all three scenarios
str(all_scenarios_Tanzania_10$age_group)
age_levels_leq_15 <- as.character(0:15)
plot_dodged_10 <- ggplot(
  all_scenarios_Tanzania_10 |>
    filter(age_group %in% age_levels_leq_15),
  aes(x = age_group, y = proportion, fill = compartment)
) +
  geom_col(position = "stack", col = NA) +
  geom_hline(yintercept = 18.2, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~scenario, ncol = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 10),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  ) +
  labs(
    # title = "C.",
    # subtitle = "",
    x = "Age",
    y = "Resistance(%)",
    fill = "Compartment"
  ) +
  scale_fill_manual(values = c(
    "X" = "#00c4aa",
    "S" = "#e573f3",
    "R" = "#fc726c",
    "Sr" = "#9b9602",
    "Rs" = "#00b3f4"
  )) +
  scale_y_continuous(limits = c(0, 43))
print(plot_dodged_10)
# par(mfrow=c(1,3))
print(plot_dodged_1)
print(plot_dodged_5)
print(plot_dodged_10)

ggsave("Figure_1.Prevalence of resistance 1 year.png", plot = plot_dodged_1, width = 8, height = 6, dpi = 300)
ggsave("Figure_5.Prevalence of resistance 5 year.png", plot = plot_dodged_5, width = 8, height = 6, dpi = 300)
ggsave("Figure_10.Prevalence of resistance 10 year.png", plot = plot_dodged_10, width = 8, height = 6, dpi = 300)
#
head(all_scenarios_Tanzania)
library(ggplot2)
plot <- all_scenarios_Tanzania %>%
  filter(compartment %in% c("R")) %>%
  ggplot(aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "stack") + # stack or dodge
  facet_wrap(~scenario, ncol = 3) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  ) +
  labs(
    title = "Resistant Bacterial Compartments in Tanzania",
    x = "Age group",
    y = "Percent",
    fill = "Compartment"
  ) +
  scale_fill_manual(values = c(
    "R"  = "#fc726c",
    "Rs" = "pink" # "#00b3f4"
  ))
print(plot)
#
library(dplyr)
head(all_scenarios_Tanzania_1)
dataset <- rbind(
  all_scenarios_Tanzania_10,
  all_scenarios_Tanzania_5,
  all_scenarios_Tanzania_1
)

plot_combined <- ggplot(
  dataset,
  aes(x = age_group, y = proportion, fill = compartment)
) +
  geom_col(position = "stack") +
  facet_grid(period ~ scenario) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    title = "Population-level Bacterial Dynamics in Tanzania",
    x = "Age",
    y = "Percent",
    fill = "Compartment"
  ) +
  scale_fill_manual(values = c(
    "X"  = "#00c4aa",
    "S"  = "#e573f3",
    "R"  = "#fc726c",
    "Sr" = "#9b9602",
    "Rs" = "#00b3f4"
  ))
print(plot_combined)
# grid.arrange(plot_dodged_1,plot_dodged_5,plot_dodged_10)
# Plot 3: Side-by-side comparison (stacked) - all scenarios in one row
table(all_scenarios_Tanzania$scenario)
all_scenarios_Tanzania$scenario <- factor(all_scenarios_Tanzania$scenario,
  levels = c("Pre-MDA", "1Year MDA", "1Year Bi-MDA", "5 Years Bi-MDA", "5 Years MDA", "5 Years NO MDA", "10 Years Bi-MDA", "10 Years MDA", "10 Years No MDA")
)
#
plot_horizontal <- ggplot(
  all_scenarios_Tanzania,
  aes(x = age_group, y = proportion, fill = compartment)
) +
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
# print(plot_horizontal)

# Plot 4: Focus on Resistance (R) only across scenarios

table(all_scenarios_Tanzania$compartment)
table(all_scenarios_Tanzania$scenario)
resistance_only <- all_scenarios_Tanzania[compartment == "R"]
table(all_scenarios_Tanzania$scenario)
table(resistance_only$scenario)
plot_resistance <- ggplot(
  resistance_only,
  aes(x = age_group, y = proportion, fill = scenario)
) +
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
  # scale_fill_viridis_d(option = "viridis", end = 0.8)
  scale_fill_manual(values = c(
    "10 Years MDA" = "#d95f02",
    "10 Years No MDA" = "#1b9e77"
  ))
print(plot_resistance)
#
# Plot 5: Total Resistance (R + Rs) comparison
#
total_resistance <- all_scenarios_Tanzania[compartment %in% c("R", "Rs"),
  .(proportion = sum(proportion)),
  by = .(age_group, scenario)
]

plot_total_resistance <- ggplot(
  total_resistance,
  aes(
    x = age_group, y = proportion,
    color = scenario, group = scenario
  )
) +
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
plot_by_compartment <- ggplot(
  all_scenarios_Tanzania,
  aes(x = age_group, y = proportion, fill = scenario)
) +
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
if (save.fig) {
  ggsave("comparison_stacked.png", plot_stacked, width = 12, height = 10, dpi = 300)
  ggsave("comparison_dodged.png", plot_dodged, width = 12, height = 10, dpi = 300)
  ggsave("comparison_horizontal.png", plot_horizontal, width = 16, height = 6, dpi = 300)
  ggsave("resistance_comparison.png", plot_resistance, width = 12, height = 6, dpi = 300)
  ggsave("total_resistance_line.png", plot_total_resistance, width = 12, height = 6, dpi = 300)
  ggsave("compartment_multiples.png", plot_by_compartment, width = 12, height = 10, dpi = 300)
  cat("\nPlots saved successfully!\n")
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
# ...............................................................................
