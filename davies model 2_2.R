
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
# In the above mentioned script, you need to change
# the WHO EMRO countries by WHO Afro countries to get Tanzania
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#Work envronment
getwd()
setwd("C:/Disk F/4.Oxford Modelling for Global Health/Afox_Ubuntu/Afox Placement with Ben Cooper")
getwd()
#Packages
pacman::p_load(deSolve, viridis, ggplot2, tidyr, dplyr, readr) 
#Section I. Routine data collected at global and state level over time
#Source: WHO:https://data.who.int/indicators/i/918081E/745F475

#Global distribution and incidence of multidrug resistant and ESBL producing Escherichia coli: an observational study of the ATLAS dataset
#Source : https://www.sciencedirect.com/science/article/pii/S2213398425002398

#Pfizer – Antimicrobial Testing Leadership and Surveillance (ATLAS)
#Source : https://www.amrindustryalliance.org/case-study/antimicrobial-testing-leadership-and-surveillance-atlas/
pacman::p_load(readr)
data <- read_csv("E.coli_resistance_C3_ALL_LATEST.csv")
print(data)
dim(data)
names(data)
names(data)[5] <- "Period"
names(data)[11] <- "Country"
names(data)[12] <- "Percentage"
names(data)
#Filter data
pacman::p_load(dplyr)
table(data$Country)
data$Country[data$Country == "United Republic of Tanzania"] <- "Tanzania"
data1 <- data %>%
  filter(Country %in% c("World","Nigeria","Namibia ","Rwanda","Sudan", "Tanzania", "Zambia","Uganda", "Malawi"))%>%
  filter(Country %in% c("World","Niger", "Tanzania", "Malawi"))
#Visualisation
pacman::p_load(ggplot2)
ggplot(data1, aes(x = Period, y =Percentage, color = Country, group = Country)) +
  #geom_line(linewidth = 1) +
  geom_point(size = 3) +
  #geom_smooth(method = "loess", se = TRUE, linewidth = 0.8) + 
  scale_x_continuous(breaks = seq(min(data1$Period),
    max(data1$Period), by = 1)) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20)
  )+
  labs(
    title = "Proportion of bloodstream infection due to Escherichia coli resistant to  C3 (%)",
    subtitle = "Annual estimates by country",
    #title = "Proportion of bloodstream infection due to Escherichia coli resistant to third-generation cephalosporins (%)",
    x = "Period",
    y = "Percentage(%)",
    #y = "Rate per 100,000",
    color = "Country"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.text = element_text(size = 10)
    )
#
g_0<-ggplot(data1, aes(x = Period, y = Percentage, color = Country, group = Country)) +
  geom_ribbon(
    data = subset(data1, Country == "World"),
    aes(
      x = Period,
      ymin = 0,
      ymax = Percentage
    ),
    inherit.aes = FALSE,
    fill = "#C2A5CF" ,
    #fill = "grey80",
    alpha = 0.4
  ) +
  geom_point(size = 3) +
 #geom_line(linewidth = 1) +
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
  plot=last_plot(),
  bg="white",
  width= 10,
  height=8,
  dpi= 300)
#Control parameters
save.fig <- FALSE
# 1. Demographic parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
#
pyramyd <- ggplot(Tanzania_pop_in_thousands, aes(x = Age_Category, y = Population_age)) +
  geom_col(fill = "steelblue") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +  
  labs(title = "Tanzania's population structure (2023)",
    x = "Age", y = "Population") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(pyramyd)
# Shorten the name
popstruc <- Tanzania_pop_in_thousands
#  Number of age groups 
A <- n_age <- nrow(Tanzania_pop_in_thousands)  # Fixed: use nrow instead of length
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
  labs(title = "Birth in Tanzania (2023)",
    x = "Age", y = "Births") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#
print(pyramyd_birth)
#Convert from 1000s per 1 year period to per person per day
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
  labs(title = "Mortality in Tanzania (2023)",
    x = "Age", y = "Deaths per 1000 pop") +
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
#1.h.Contact matrix  and visualization
(m_contact_1y_Tanzania <- as.matrix(read.csv("3.U.1.contact_Tanzania_1y.csv")))
dim(m_contact_1y_Tanzania)
colSums(ageing)
#
for(i in 1:n_age){
  for(j in 1:n_age){
    m_contact_1y_Tanzania[i,j]<-m_contact_1y_Tanzania[i,j]/25
  }
}
colnames(m_contact_1y_Tanzania ) <-c(as.character(0:99), "100+")
rownames(m_contact_1y_Tanzania ) <- c(as.character(0:99), "100+")
m_contact_1y_Tanzania 
#Visualization of my contact matrix
pacman::p_load(ggplot2,reshape2)
#data frame
#df <- melt(m_contact_1y_Tanzania )
df <- reshape2::melt(m_contact_1y_Tanzania)
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
Srindex <-(3*n_age+1):(4*n_age)          # drug-sensitive, treated
Rsindex <-(4*n_age+1):(5*n_age)          # drug-resistant, treated
Dindex <- (5 * n_age + 1):(6 * n_age)    # Cummulative deaths
CumRindex <-(6 * n_age + 1):(7 * n_age)  # Cummulative resistance 
#...............................................................................
#Considering how to add MDA as an intervention
#mda_start_times <- c(365, 79, 4380, 4745)
(mda_start_times<-(0:50)*365.25)
#(mda_duration <- 30)
(mda_duration <- 30)
#//rate<--log(1-(0.8))/mda_duration
#//rate
mda_active <- function(time, mda_starts, duration) {
  any(sapply(mda_starts, function(start) {
    time >= start && time < (start + duration)
  }))
 }
#...............................................................................
bacteria.odes <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Extract compartments 
    X <- state[Xindex]   # Uninfected, untreated
    S <- state[Sindex]   # drug-sensitive, untreated
    R <- state[Rindex]   # drug-resistant, untreated
    Sr <- state[Srindex] # drug-sensitive, treated
    Rs <- state[Rsindex] # drug-resistant, treated
    D <- state[Dindex]   # cumulative deaths
    CumR <- state[CumRindex] 
    # Total population
    N <- X + S + R + Sr + Rs 
    # Births
    births <- rep(0, n_age)            
    #births[1] <- sum(popbirth[,5] * N) #for dynamic population
    total_deaths <- sum(mort * N)      #for static population
    births[1] <- total_deaths          #for static population
    #
    S.tot <- S + Sr  # Susceptible co-colonised total
    R.tot <- R + Rs  # Resistance co-colonised total
    
    #lamda.S <- beta.S * (m_contact_1y_Tanzania %*% S.tot)     #Within host transmission
    #lamda.R <- beta.R * (m_contact_1y_Tanzania %*% R.tot)
    
    lamda.S <- beta.S * (m_contact_1y_Tanzania %*% (S.tot/N)) #Between host transmission
    lamda.R <- beta.R * (m_contact_1y_Tanzania %*% (R.tot/N))
    
    #Intervention : MDA implementation
    #...........................................................................
    # Considering how to add the   MDA intervention into the odes
    mda_targeted_ages <- 0:5         # Targeting 
    azt <- rep(0, n_age)             # Initialize azt vector
    azt[mda_targeted_ages] <- 0.8    # 1 # Ages targeted 
    #b <- ifelse(mda_active(t, mda_start_times, mda_duration), (1*12 / 365.25), 0)
    #//b <- ifelse(mda_active(t, mda_start_times, mda_duration), 0.16, 0)
    #//a <- b * azt                     # Apply a only on targeted ages
    #...........................................................................
    # Schedule for MDA (Pulse)
    #//time_in_cycle <- t %% mda_cycle
   #// mda_active <- (time_in_cycle < mda_duration) * 1
    
    # MDA starting after 2025 (25 years from start year 2000)
    #//rate <- ifelse(((t > 0*365) & mda_active), r_mda, 0)
    #mort <- ifelse(((t > 0*365) & mda_active), r_mda, 0)
    b <- ifelse(mda_active(t, mda_start_times, mda_duration),a+(alpha*a_use), (alpha*a_use))
    #b <- ifelse(mda_active(t, mda_start_times, mda_duration),a+u.S, u.S)
    a_t <- b * azt
    bc <- ifelse(mda_active(t, mda_start_times, mda_duration),a.c+(alpha*a_use), (alpha*a_use))
    #bc <- ifelse(mda_active(t, mda_start_times, mda_duration),a.c+u.S, u.S)
    a.c_t <- bc * azt
    #Mortality reduction during MDA
    mort_eff <- mort  # start from baseline
    if (mda_active(t, mda_start_times, mda_duration)) {
      mort_eff[mda_targeted_ages] <- mort[mda_targeted_ages] * (1 - theta)
    }
    #browser()
    #............................................................................
    # ODEs system
    dX <- births + u.S * S + u.R * R + u.c * (Sr + Rs) - (lamda.S + lamda.R ) * X +  ageing %*% X - mort_eff * X  # + MDA * S # MDA * (1-kappa) * S  # Successfully treated return to X 
    #dX <- births + u.S * S + u.R * R + u.c * (Sr + Rs)+ a * S - (lamda.S + lamda.R ) * X +  ageing %*% X - mort * X  # New
    
    dS <- lamda.S * X - u.S * S - k * lamda.R * S - a_t * S   +  ageing %*% S - mort_eff * S #ttt                     # - MDA  * # Remove all in S
    
    dR <- lamda.R * X - u.R * R - k * lamda.S * R + a.c_t * (Sr + Rs) +  ageing %*% R - mort_eff * R 
    #dR <- lamda.R * X - u.R * R - k * lamda.S * R + a * (Sr + S) +  ageing %*% R - mort * R  #New
    
    dSr <- k * lamda.R * S - Sr * u.c - a.c_t * Sr  +  ageing %*% Sr - mort_eff * Sr #ttt 
    
    dRs <- k * lamda.S * R - Rs * u.c - a.c_t * Rs  + ageing %*% Rs - mort_eff * Rs  #a=epsilon
    #dRs <- k * lamda.S * R - Rs * u.c  + ageing %*% Rs - mort * Rs   #a=epsilon #New
    dD<- mort_eff * X+mort_eff * S+mort_eff * R+mort_eff * Sr+mort_eff * Rs
    dCumR <- lamda.R * X + a_t * (Sr + Rs)

  list(c(dX, dS, dR, dSr, dRs,dD,dCumR))
  })
}
bacteria.solve <- function(t, state, parameters) {
  parameters[["beta.R"]] <- parameters[["beta.S"]] * (1 - parameters[["c"]])
  out_Tanzania <- as.data.frame(ode(state, t, bacteria.odes, parameters))
  return(out_Tanzania)
}
# Parameters
parms.orig <- list(
  beta.S = 5,              # Transmission of sensitive         : (β = 5 month−1)   
  u.S = 1,                 # Clearance sensitive (natural)     : (u = 1 month−1)
  u.R = 1,                 # Clearance resistant (natural)     : (u = 1 month−1) :lower than susceptible?
  u.c = 1,                 # Clearance co-colonised (natural)
  a = 0.16,                # Clearance sensitive (drug-induced)
  a.c = 0.16,              # Clearance co-colonised (drug-induced)
  a_use = 0.08,          # Antibiotic use (~0.01–0.05:routine use,0.05 – 0.20: High use communitie,0.20 – 0.80 (short period) :MDA  )
  alpha= 0.05 ,            # Antibiotic effect: 0.005−0.05
  k=  0.05,                # The efficiency of co-colonisation : (k= 0.25, 0.5, and 1.0)
  m_contact = m_contact_1y_Tanzania, # Social contacts per day
  # MDA parameters
  mda_cycle= 365,
  mda_duration = 30,
  mda_cov =  0.6,
  theta   =  0,    # (for static population) under-five mortality reduction
  #theta   = 0.13,#(for dynamic population) under-five mortality reduction
  kappa =    0    # 0.05  # Proportion that develop/select resistance (Assumed)
)
# MDA rate calculation: Exponential decay
(parms.orig$r_mda<--log(1-parms.orig$mda_cov)/parms.orig$mda_duration)
parms.orig
parms.orig[1:4]
# Convert daily
#[1:5]
parms.orig[1:4] <- lapply(parms.orig[1:4], function(x) x*12/365.25)# convert to daily
parms.orig[1:4]
# Adjusted clearance rates: Esther et al study
#parms.orig["u.S"]= 0.0098
#parms.orig["u.R"]= 0.0098
#parms.orig["u.C"]= 0.0098
#βR=−ln(1−p)/D where p is the daily probability of acquisition
(prob<-3.92/100) # 3.92% from Rebecca et al
(ar<- log(1/(1-prob)))
#[1] 0.03998901
#parms.orig["beta.S"]= 0.03998901  # 5/30#3.92/100 #0.0098 # (5/30: Davies) (3.92%: Rebacca Lynn et al)
# MDA parameters
#parms.orig["mda_cov"]= 0
#parms.orig["mda_cycle"] = 365        # Period between 2 mda
#parms.orig["mda_duration"] = 2*30    # mda campaign duration in days
parms.orig[["c"]] <- 0.01  # 0.1|--> cost of resistance on transmission
parms.orig[["k"]] <- 0.5   # efficiency of co-colonisation
# Initial conditions
names(Tanzania_pop_in_thousands)
head(Tanzania_pop_in_thousands)
dim(Tanzania_pop_in_thousands)
initP <- Tanzania_pop_in_thousands[,5]
sum(initP)
#initS <- Tanzania_pop_in_thousands[,5]
initX <-   0.95 * Tanzania_pop_in_thousands[,5]
initS <-   0.025 * Tanzania_pop_in_thousands[,5]
initR <-   0.025 * Tanzania_pop_in_thousands[,5]
initSr <-   0 * Tanzania_pop_in_thousands[,5]
initRs <-   0 * Tanzania_pop_in_thousands[,5]
initD <-    0 * Tanzania_pop_in_thousands[,5]
initCumR<-     0 * Tanzania_pop_in_thousands[,5]
# Combine initial states
state.orig<-c(initX,initS,initR,initSr,initRs,initD,initCumR)
#States
tvec_1_a <- seq(0,    1*365.25 ,1)      #1 Year-MDA
tvec_1_b <- seq(0,    1*365.25 ,1)      #Pre-MDA
tvec_5_a <- seq(0,    5*365.25 ,1)      #5 years of MDA
tvec_5_b <- seq(0,    5*365.25 ,1)      #5 years of MDA
tvec_10_a <- seq(0,  10*365.25 ,1)      #10 years of MDA
tvec_10_b <- seq(0,  10*365.25 ,1)      #10 years No MDA, i will need a = 10
#Run model 
#~~~~~~~~~Annual MDA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
parms <- parms.orig # Annual MDA parameters
#~~~~~~~~~No-annual MDA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
parms_noMDA <- parms  # No annual MDA  parameters
parms_noMDA$r_mda <- 0
parms_noMDA$a <- 0
parms_noMDA$a.c <- 0
parms_noMDA$theta <- 0
parms_noMDA$rho   <- 0
parms_noMDA
state <- state.orig
start <- Sys.time()
#Baseline: Equilibrium at 200 years 
tvec_0_b<- seq(0,    200*365.25 ,1) 
(mda_start_times<-(0:200)*365.25)
out_0_b_Tanzania <- bacteria.solve(tvec_0_b, state, parms_noMDA)
last10 <- tail(out_0_b_Tanzania, 3650)
#Total population stability
total_pop <- rowSums(last10[, c(Xindex, Sindex, Rindex, Srindex, Rsindex)])
#plot(total_pop, type="l")
# Max difference
max(abs(diff(total_pop)))
state <- as.numeric(out_0_b_Tanzania[nrow(out_0_b_Tanzania), -1])
#print(state )
#Modelling different scenarios
parms 
parms_noMDA
(mda_start_times<-(0:50)*365.25)             # Early MDA start
out_1_a_Tanzania <- bacteria.solve(tvec_1_a, state, parms)  # Annual 
(mda_start_times<-(0:50)*365.25)            # Delay MDA start
out_1_b_Tanzania <- bacteria.solve(tvec_1_b, state, parms_noMDA)
#
#//out_1_c_Tanzania <- bacteria.solve(tvec_1_a, state, parms_bi) # Bi-annual
#
(mda_start_times<-(0:50)*365.25)             # Early MDA start
out_5_a_Tanzania <- bacteria.solve(tvec_5_a, state, parms)
#
(mda_start_times<-(0:50)*365.25)             # Delay MDA start
out_5_b_Tanzania <- bacteria.solve(tvec_5_b, state, parms_noMDA)
#
#//out_5_c_Tanzania <- bacteria.solve(tvec_5_a, state, parms_bi)
(mda_start_times<-(0:50)*365.25)              # Early MDA start
out_10_a_Tanzania <- bacteria.solve(tvec_10_a, state, parms)
# For no MDA scenario,set a = 0 or start times very high
(mda_start_times<-(0:50)*365.25)# Delay MDA start
out_10_b_Tanzania <- bacteria.solve(tvec_10_b, state, parms_noMDA)
#
#//out_10_c_Tanzania <- bacteria.solve(tvec_10_a, state, parms_bi)
#~~~~~~~~~Bi-annual MDA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mda_start_times <- (0:50) * (365.25/2)
out_1_c_Tanzania <- bacteria.solve(tvec_1_a, state, parms) # Bi-annual(same parameters excepts mda_start)
mda_start_times <- (0:50) * (365.25/2)
out_5_c_Tanzania <- bacteria.solve(tvec_5_a, state, parms)
mda_start_times <- (0:50) * (365.25/2)
out_10_c_Tanzania <- bacteria.solve(tvec_10_a, state,parms)
end<-Sys.time()
end-start
#MDA stops_ before the end of simulation
mda_start_times <- (0:5) * (365.25/1)
out_10_a_5_Tanzania <- bacteria.solve(tvec_10_a, state,parms)
mda_start_times <- (0:6) * (365.25/1)
out_10_a_6_Tanzania <- bacteria.solve(tvec_10_a, state,parms)
mda_start_times <- (0:7) * (365.25/1)
out_10_a_7_Tanzania <- bacteria.solve(tvec_10_a, state,parms)
#MDA stops_bi annual
mda_start_times <- (0:5) * (365.25/2)
out_10_c_5_Tanzania <- bacteria.solve(tvec_10_a, state,parms)
mda_start_times <- (0:6) * (365.25/2)
out_10_c_6_Tanzania <- bacteria.solve(tvec_10_a, state,parms)
mda_start_times <- (0:7) * (365.25/2)
out_10_c_7_Tanzania <- bacteria.solve(tvec_10_a, state,parms)
#Annual
results_1_a_Tanzania <- as.data.frame(out_1_a_Tanzania)
results_1_b_Tanzania <- as.data.frame(out_1_b_Tanzania)
results_1_c_Tanzania <- as.data.frame(out_1_c_Tanzania)
results_5_a_Tanzania <- as.data.frame(out_5_a_Tanzania)
results_5_b_Tanzania <- as.data.frame(out_5_b_Tanzania)
results_5_c_Tanzania <- as.data.frame(out_5_c_Tanzania)
#Bi annual
results_10_a_Tanzania <- as.data.frame(out_10_a_Tanzania)
results_10_b_Tanzania <- as.data.frame(out_10_b_Tanzania)
results_10_c_Tanzania <- as.data.frame(out_10_c_Tanzania)
#Stops
results_10_a_5_Tanzania <- as.data.frame(out_10_a_5_Tanzania)
results_10_a_6_Tanzania <- as.data.frame(out_10_a_6_Tanzania)
results_10_a_7_Tanzania <- as.data.frame(out_10_a_7_Tanzania)
results_10_c_5_Tanzania <- as.data.frame(out_10_c_5_Tanzania)
results_10_c_6_Tanzania <- as.data.frame(out_10_c_6_Tanzania)
results_10_c_7_Tanzania <- as.data.frame(out_10_c_7_Tanzania)
# Column names 
compartment_names <- c("X", "S", "R", "Sr", "Rs","D","CumR")
col_names <- c("time")
for(comp in compartment_names) {
  for(age in age_groups) {
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

colnames(results_10_a_5_Tanzania) <- col_names
colnames(results_10_a_6_Tanzania) <- col_names
colnames(results_10_a_7_Tanzania) <- col_names
colnames(results_10_c_5_Tanzania) <- col_names
colnames(results_10_c_6_Tanzania) <- col_names
colnames(results_10_c_7_Tanzania) <- col_names

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

#6.Verification of the population over time
#Total population at each time point
names(results_1_b_Tanzania)
dim(results_1_b_Tanzania)
str(results_1_b_Tanzania)
total_pop<- rowSums(results_1_b_Tanzania[, 2:(5*n_age)])  # This will exclude time column
# Check conservation: Here the population is conserved if the birth=deaths
#                   : The population is not conserved as births and deaths  are not equal(we used realistic country and age specific births and deaths) 
#Initial population
(initial_pop <- sum(Tanzania_pop_in_thousands[,5]))
# Time 
tim<-results_1_b_Tanzania[,1]
#Data frame for plotting
P_t<-as.data.frame(cbind(tim,total_pop,initial_pop))
#Population variations
(change<-round(((tail(P_t$total_pop,1)-initial_pop)/initial_pop)*100,2)) # Checks
#Packages for visualization
pacman::p_load(ggplot2,ggtext)
#Visualization of the initial population and total population over time
#Plot 
#Population dynamic (with annotation)
population<-ggplot(P_t, aes(x = tim, y = total_pop)) +
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
#Plot model
cols <- viridis(5)
par(mfrow = c(1, 1))
#Total across age groups: Here i will use Index#1e6 
X_total <- rowSums(out_1_b_Tanzania[,  Xindex + 1])  #Here i added +1 because first column is time
S_total <- rowSums(out_1_b_Tanzania[,  Sindex + 1])
R_total <- rowSums(out_1_b_Tanzania[,  Rindex + 1])  #+ rowSums(out_1_b_Tanzania[, Rsindex + 1]) 
Rs_total <- rowSums(out_1_b_Tanzania[, Rsindex + 1])
Sr_total <- rowSums(out_1_b_Tanzania[, Srindex + 1])
D_total <- rowSums(out_1_b_Tanzania[, Dindex + 1])
prevalence<-round(R_total*100/(X_total+S_total+R_total+Rs_total+Sr_total),1)
mortality<-round(D_total*100/(X_total+S_total+R_total+Rs_total+Sr_total),1)
summary(prevalence)
summary(mortality)

#Simulation time in years 
length(R_total)/365.25
#Equilibrium checks
tail(R_total)
df_no_resisitance_mortality <- data.frame(
  time = out_1_b_Tanzania[, 1],
  S  = S_total,
  X  = X_total,
  R  = R_total,
  Rs = Rs_total,
  Sr = Sr_total,
  D  = D_total
  )
df_with_resisitance_mortality <- data.frame(
  time = out_1_b_Tanzania[, 1],
  S  = S_total,
  X  = X_total,
  R  = R_total,
  Rs = Rs_total,
  Sr = Sr_total,
  D  = D_total,
  prevalence,
  mortality
)
head(df_with_resisitance_mortality)
ggplot(df_with_resisitance_mortality, aes(x = time, y = prevalence)) +
  geom_line(size = 1.2) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10)
  )+
  scale_x_continuous(
    breaks = seq(min(df_with_resisitance_mortality$time), max(df_with_resisitance_mortality$time), by = 365)
  )  +
  labs(
    title = "Baseline prevalence of E.Coli infections resistant to Azithromycin",
    x = "Time",
    y = "Resistance (%)",
    colour = "prevalence"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
pacman::p_load(tidyr,dplyr,ggplot2,scales)
df_b <- df_no_resisitance_mortality %>%
  pivot_longer(
    cols = -time,
    names_to = "Compartment",
    values_to = "Population"
  )
table(df_b$Compartment)
incidence<-df_b %>%
  filter(Compartment %in% c("R")) %>%
  ggplot(aes(x = time, y = Population, color = Compartment)) +
  #geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_y_continuous(
    limits = c(0, 60000000),
    breaks = seq(0, 60000000, by = 2000000),
    labels = label_number(scale = 1e-6, suffix = " M")
  )+
  labs(
    x = "Time",
    y = "Resistance",
    title = "Resistance over time"
  ) +
  #theme_minimal()
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14)
  )
incidence
#Here i want to extract the last year to see if i reach equilibrium
df_b %>%
  filter(Compartment == "R") %>%
  slice_tail(n = 10) %>% #the last  ten observations
  ggplot(aes(x = time, y = Population, color = Compartment)) +
  #geom_line(linewidth = 1.1) +
  geom_point(size = 4) +
  scale_y_continuous(
    limits = c(
      0,
      max(
        df_b %>%
          filter(Compartment == "R") %>%
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
pacman::p_load(dplyr,scales,ggplot2)
# data 
df_last30 <- df_b %>%
  filter(Compartment == "R") %>%
  slice_tail(n = 10)
# plot
b<-ggplot(df_last30, aes(x = time, y = Population, color = Compartment)) +
  
  # Points
  geom_point(size = 4) +
  # Vertical labels below points
  geom_text(
    aes(label = label_number(scale = 1e-6, suffix = " M")(Population)),
    angle = 90,
    vjust = 1.5,   # pushes text below the point
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
c<-grid.arrange(incidence,b,ncol=2)
ggsave("Resistance_equilibrium_check.png",plot=c,width=10,height = 5,dpi=300)
#
dynamics<-ggplot(df_b|> dplyr::filter(Compartment != "D"), aes(x = time, y = Population, color = Compartment)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(
    limits = c(0, 50000000),
    breaks = seq(0, 50000000, by = 10000000),
    labels = label_number(scale = 1e-6, suffix = " M")
  )+
  labs(
    x = "Time",
    y = "Population",
    title = "E.Coli Dynamics in Tanzania"
  ) +
  #theme_minimal()
theme(
  axis.text.x = element_text(angle = 0, hjust = 1),
  plot.title = element_text(hjust = 0.5, size = 14)
)
print(dynamics)
ggsave("Figure 1.Population and E.coli dynamics.png",plot=grid.arrange(population,dynamics,ncol=2),width=12,height=5,dpi=300)
#R_total<-R_total+Rs_total
#tails
head(out_1_b_Tanzania)
a_1_a<-tail(cumsum(rowSums(out_1_a_Tanzania[, Rindex + 1])), 1)
a_1_b<-tail(cumsum(rowSums(out_1_b_Tanzania[, Rindex + 1])), 1)
a_1_c<-tail(cumsum(rowSums(out_1_c_Tanzania[, Rindex + 1])), 1)
a_5_a<-tail(cumsum(rowSums(out_5_a_Tanzania[, Rindex + 1])), 1)
a_5_b<-tail(cumsum(rowSums(out_5_b_Tanzania[, Rindex + 1])), 1)
a_5_c<-tail(cumsum(rowSums(out_5_c_Tanzania[, Rindex + 1])), 1)
a_10_a<-tail(cumsum(rowSums(out_10_a_Tanzania[, Rindex + 1])), 1)
a_10_b<-tail(cumsum(rowSums(out_10_b_Tanzania[, Rindex + 1])), 1)
a_10_c<-tail(cumsum(rowSums(out_10_c_Tanzania[, Rindex + 1])), 1)
a_1_a
a_1_b
a_1_c
a_1_d_a<-round((a_1_a-a_1_b)*100/a_1_b,1)  # annual
a_1_d_bi<-round((a_1_c-a_1_b)*100/a_1_b,1) # bi annual
a_1_d_a
a_1_d_bi
a_5_a
a_5_b
a_5_c
a_5_d_a<-round((a_5_a-a_5_b)*100/a_5_b,1)
a_5_d_bi<-round((a_5_c-a_5_b)*100/a_5_b,1)
a_5_d_a
a_5_d_bi
a_10_a
a_10_b
a_10_c
a_10_d_a<-round((a_10_a-a_10_b)*100/a_10_b,1)
a_10_d_bi<-round((a_10_c-a_10_b)*100/a_10_b,1)
a_10_d_a
a_10_d_bi
totalresistance <- data.frame(
  Scenario = c("Pre-MDA","1Y MDA","1Y Bi-MDA", "5Y MDA","5Y NO-MDA","5Y Bi-MDA", "10Y MDA", "10Y NO-MDA","10Y Bi-MDA"),
  R_final = c(a_1_b,a_1_a,a_1_c, a_5_a,a_5_b,a_5_c, a_10_a, a_10_b,a_10_c)
)
AMR_burden <- data.frame(
  Scenario = c("1Y MDA","1Y Bi-MDA", "5Y MDA","5Y Bi-MDA", "10Y MDA","10Y Bi-MDA"),
  R_final = c(a_1_d_a,a_1_d_bi,a_5_d_a,a_5_d_bi,a_10_d_a,a_10_d_bi)
)
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
#Plot
table(totalresistance$Scenario)
pacman::p_load(dplyr)
resistance_1Y <- totalresistance %>%
  filter(Scenario %in% c("Pre-MDA", "1Y MDA","1Y Bi-MDA"))
#Filter(Scenario %in% c( "5Y MDA","5Y NO-MDA"))
ggplot(resistance_1Y, aes(x = Scenario, y = R_final, fill = Scenario)) +
  geom_bar(stat = "identity", color = "black") +
 scale_y_continuous(labels = scales::label_number(scale = 1e-9, accuracy = 1, suffix = " M ")) +  
  scale_fill_manual(
    values = c(
      "Pre-MDA"   = "#0072B2",  # 
      "1Y MDA"    = "#009E73",  # 
      "1Y Bi-MDA" = "#D55E00"   #
    )
  )+
#scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "Effect of MDA AZT on resistance",
    x = "MDA Scenario",
    y = "Resistant population "
  ) +
   #theme_lancet() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14)
    )
#
resistance_5Y <- totalresistance %>%
  filter(Scenario %in% c("5Y MDA","5Y Bi-MDA", "5Y NO-MDA"))
ggplot(resistance_5Y, aes(x = Scenario, y = R_final, fill = Scenario)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, accuracy = 1, suffix = " M ")) +  
  #scale_y_continuous(labels = scales::label_number(scale = 1e-9, accuracy = 1, suffix = " M ")) +  
  scale_fill_manual(
    values = c(
      "5Y NO-MDA"   = "#0072B2",  # 
      "5Y MDA"    = "#009E73",  # 
      "5Y Bi-MDA" = "#D55E00"   #
    )
  )+
  #scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "Burden of AMR post-MDA",
    x = "MDA Scenario",
    y = "Resistant population"
  ) +
  #theme_lancet() +
  #theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14)
    )

resistance_10Y <- totalresistance %>%
  filter(Scenario %in% c("10Y MDA","10Y Bi-MDA", "10Y NO-MDA"))
ggplot(resistance_10Y, aes(x = Scenario, y = R_final, fill = Scenario)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, accuracy = 1, suffix = " M ")) +  
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(
    values = c(
      "10Y NO-MDA"   = "#0072B2",  # 
      "10Y MDA"    = "#009E73",  # 
      "10Y Bi-MDA" = "#D55E00"   #
    )
  )+
  #scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "Burden of AMR post-MDA",
    x = "MDA Scenario",
    y = "Resistant population "
  ) +
  #theme_lancet() +
  #theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14)
  )
#Re-order
table(totalresistance$Scenario)
totalresistance$Scenario <- factor(
  totalresistance$Scenario,
  levels = c("10Y Bi-MDA","10Y MDA","10Y NO-MDA" ,"5Y Bi-MDA","5Y MDA","5Y NO-MDA","1Y Bi-MDA","1Y MDA","Pre-MDA")
)
ggplot(totalresistance, aes(x = Scenario, y = R_final, fill = Scenario)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, accuracy = 1, suffix = " M ")) +  
  scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "Effect of MDA AZT on resistance",
    x = "MDA Scenario",
    y = "Resistant population"
  ) +
  #theme_lancet() +
  #theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14)
  )
#Burden of resistance
table(AMR_burden$Scenario)
#change the order
AMR_burden$Scenario <- factor(
  AMR_burden$Scenario,
  levels = c("10Y Bi-MDA","10Y MDA", "5Y Bi-MDA","5Y MDA","1Y Bi-MDA","1Y MDA")
)
cases_increases<-ggplot(AMR_burden, aes(x = Scenario, y = R_final, fill = Scenario)) +
  geom_bar(stat = "identity", color = "black") +
  # Add percentage labels on bars
  geom_text(aes(label = paste0(round(R_final, 1), " %")),
    vjust = -0.5, size = 4) +
  scale_y_continuous(labels = scales::label_number(scale = 1, accuracy = 1, suffix = " % ")) +  
  scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "A",
    #title = "Population-level Post-MDA Azithromycin Resistance in Escherichia coli",
    x = "MDA Scenario",
    y = "Percentage of change in resistance(%)"
  ) +
  #theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14)
  )
#print
print(cases_increases)
#Mortality reduction due to MDA
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
compartments <- c("X", "S", "R", "Sr", "Rs", "D", "CumR")
#Column names for under-five across all compartments
under_five_cols <- c()
for (comp in compartments) {
  for (age in ages_under_5) {
    col_name <- paste0(comp, "_", age)
    under_five_cols <- c(under_five_cols, col_name)
  }
}
#Extract the data for under five age groups
under_five_data <- results_1_a_Tanzania[, c("time", under_five_cols)]
#
results_1_a_Tanzania_under_five <- results_1_a_Tanzania[, c("time", under_five_cols)]
results_1_b_Tanzania_under_five <- results_1_b_Tanzania[, c("time", under_five_cols)]
results_1_c_Tanzania_under_five <- results_1_c_Tanzania[, c("time", under_five_cols)]

results_5_a_Tanzania_under_five <- results_5_a_Tanzania[, c("time", under_five_cols)]
results_5_b_Tanzania_under_five <- results_5_b_Tanzania[, c("time", under_five_cols)]
results_5_c_Tanzania_under_five <- results_5_c_Tanzania[, c("time", under_five_cols)]

results_10_a_Tanzania_under_five <- results_10_a_Tanzania[, c("time", under_five_cols)]
results_10_b_Tanzania_under_five <- results_10_b_Tanzania[, c("time", under_five_cols)]
results_10_c_Tanzania_under_five <- results_10_c_Tanzania[, c("time", under_five_cols)]
#Under five level deaths 
D_cols <- grep("^D_", colnames(results_1_a_Tanzania_under_five))
#
d_1_a  <- tail(cumsum(rowSums(results_1_a_Tanzania_under_five[, D_cols])), 1) # annual mda
d_1_b  <- tail(cumsum(rowSums(results_1_b_Tanzania_under_five[, D_cols])), 1)
d_1_c  <- tail(cumsum(rowSums(results_1_c_Tanzania_under_five[, D_cols])), 1) # bi annual mda
#
d_5_a  <- tail(cumsum(rowSums(results_5_a_Tanzania_under_five[, D_cols])), 1)
d_5_b  <- tail(cumsum(rowSums(results_5_b_Tanzania_under_five[, D_cols])), 1)
d_5_c  <- tail(cumsum(rowSums(results_5_c_Tanzania_under_five[, D_cols])), 1)
#
d_10_a <- tail(cumsum(rowSums(results_10_a_Tanzania_under_five[, D_cols])), 1)
d_10_b <- tail(cumsum(rowSums(results_10_b_Tanzania_under_five[, D_cols])), 1)
d_10_c <- tail(cumsum(rowSums(results_10_c_Tanzania_under_five[, D_cols])), 1)
#Population level deaths
#d_1_a  <- tail(cumsum(rowSums(out_1_a_Tanzania[, Dindex + 1])), 1) # annual mda
#d_1_b  <- tail(cumsum(rowSums(out_1_b_Tanzania[, Dindex + 1])), 1)
#d_1_c  <- tail(cumsum(rowSums(out_1_c_Tanzania[, Dindex + 1])), 1) # bi annual mda

#d_5_a  <- tail(cumsum(rowSums(out_5_a_Tanzania[, Dindex + 1])), 1)
#d_5_b  <- tail(cumsum(rowSums(out_5_b_Tanzania[, Dindex + 1])), 1)
#d_5_c  <- tail(cumsum(rowSums(out_5_c_Tanzania[, Dindex + 1])), 1)

#d_10_a <- tail(cumsum(rowSums(out_10_a_Tanzania[, Dindex + 1])), 1)
#d_10_b <- tail(cumsum(rowSums(out_10_b_Tanzania[, Dindex + 1])), 1)
#d_10_c <- tail(cumsum(rowSums(out_10_c_Tanzania[, Dindex + 1])), 1)
#Deaths averted (MDA effect)
d_1_d_a  <- d_1_b  - d_1_a  # annual
d_1_d_bi  <- d_1_b  - d_1_c # bi annual
#
d_5_d_a  <- d_5_b  - d_5_a
d_5_d_bi <- d_5_b  - d_5_c
#
d_10_d_a <- d_10_b - d_10_a
d_10_d_bi <- d_10_b - d_10_c
#Total deaths
totalDeaths <- data.frame(
  Scenario =        c("Pre-MDA", "1Y MDA","1Y Bi-MDA","5Y MDA","5Y NO-MDA","5Y Bi-MDA", "10Y MDA","10Y NO-MDA","10Y Bi-MDA"),
  Deaths_millions = c(d_1_b,      d_1_a,    d_1_c,       d_5_a,  d_5_b,      d_5_c,        d_10_a,  d_10_b,      d_10_c)
)
print(totalDeaths)
#
Deaths_averted <- data.frame(
  Scenario = c("1Y MDA","1Y Bi-MDA", "5Y MDA","5Y Bi-MDA", "10Y MDA","10Y Bi-MDA"),
  Deaths_averted_millions = c(d_1_d_a,d_1_d_bi, d_5_d_a,d_5_d_bi, d_10_d_a,d_10_d_bi)
)
#
deaths_1Y <- totalDeaths %>%
filter(Scenario %in% c("Pre-MDA", "1Y MDA","1Y Bi-MDA"))

ggplot(deaths_1Y, aes(x = Scenario, y = Deaths_millions, fill = Scenario)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " K ")) +  
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(
    values = c(
      "Pre-MDA"   = "#0072B2",  # 
      "1Y MDA"    = "#009E73",  # 
      "1Y Bi-MDA" = "#D55E00"   #
    )
  )+
  #scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "Effect of MDA on Mortality",
    x = "Scenario",
    y = "Deaths "
  ) +
  theme(axis.text.x = element_text(hjust = 1))
#
deaths_1Y <- totalDeaths %>%
  filter(Scenario %in% c("Pre-MDA", "1Y MDA","1Y Bi-MDA"))

ggplot(deaths_1Y, aes(x = Scenario, y = Deaths_millions, fill = Scenario)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +  
  scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "Effect of MDA on Mortality",
    x = "Scenario",
    y = "Deaths"
  ) +
  theme(axis.text.x = element_text(hjust = 1))
#
deaths_5Y <- totalDeaths %>%
  filter(Scenario %in% c("5Y NO-MDA", "5Y MDA","5Y Bi-MDA"))

ggplot(deaths_5Y, aes(x = Scenario, y = Deaths_millions, fill = Scenario)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +  
  scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "Effect of 5-Year MDA on Mortality",
    x = "Scenario",
    y = "Deaths "
  ) +
  theme(
    axis.text.x = element_text(hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
#
deaths_10Y <- totalDeaths %>%
  filter(Scenario %in% c("10Y NO-MDA", "10Y MDA","10Y Bi-MDA"))

ggplot(deaths_10Y, aes(x = Scenario, y = Deaths_millions, fill = Scenario)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +  
  scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "Effect of 10-Year MDA on Mortality",
    x = "Scenario",
    y = "Deaths "
  ) +
  theme(
    axis.text.x = element_text(hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
#
names(Deaths_averted)
table(Deaths_averted$Scenario)
#change the order
Deaths_averted$Scenario <- factor(
  Deaths_averted$Scenario,
  levels = c("10Y Bi-MDA","10Y MDA", "5Y Bi-MDA","5Y MDA", "1Y Bi-MDA","1Y MDA")
)
ggplot(Deaths_averted , aes(x = Scenario, y = Deaths_averted_millions, fill = Scenario)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +  
  scale_fill_viridis_d(option = "viridis") +
  labs(
    title = "Effect of 10-Year MDA on Mortality",
    x = "Scenario",
    y = "Deaths averted "
  ) +
  theme(
    axis.text.x = element_text(hjust = 1),
    plot.title = element_text(hjust = 0.5)
    )
#Impact assessment on mortality
change_1Y_a <-round((d_1_a-d_1_b)*100/d_1_b,2)
change_1Y_a
change_1Y_c <-round((d_1_c-d_1_b)*100/d_1_b,2)
change_1Y_c
#
change_5Y_a <-round((d_5_a-d_5_b)*100/d_5_b,2)
change_5Y_a
change_5Y_c <-round((d_5_c-d_5_b)*100/d_5_b,2)
change_5Y_c
#
change_10Y_a <-round((d_10_a-d_10_b)*100/d_10_b,2)
change_10Y_a
change_10Y_c <-round((d_10_c-d_10_b)*100/d_10_b,2)
change_10Y_c
#
pacman::p_load(ggplot2,dplyr,tidyr)
#Data frame
impact_deaths_Tanzania<- data.frame(
  Year = rep(c(1, 5, 10), each = 2),
  Scenario = rep(c("annual MDA", "bi-annual MDA"), times = 3),
  Change = c(change_1Y_a, change_1Y_c,
    change_5Y_a, change_5Y_c,
    change_10Y_a, change_10Y_c)
)
#Visualisation
deaths_decrease<-ggplot(impact_deaths_Tanzania, aes(x = factor(Year), y = Change, fill = Scenario)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = paste0(Change, "%")), 
    position = position_dodge(width = 0.9), vjust = -0.5) +
  labs(
    title = "B",
    #title = "Impact of MDA on under five mortality reduction in Tanzania",
    x = "Years",
    y = "Percentage change of deaths (%)") +
  scale_fill_manual(values = c("skyblue", "tomato")) +
  theme_minimal()
print(deaths_decrease)
grid.arrange(cases_increases,deaths_decrease,ncol=2)
ggsave("Figure 2.Impact of MDA on mortality and resistance.png", plot=grid.arrange(cases_increases,deaths_decrease,ncol=2),width=12,height=5,dpi=300)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#            Prevalence of AMR                                                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Sheetal approach : Load required libraries
pacman::p_load(deSolve,doParallel,QuantPsyc)
#
#out_1_a_Tanzania <- bacteria.solve(tvec_1_a, state, parms)
results_1_a_Tanzania <- as.data.frame(out_1_a_Tanzania)
compartment_names <- c("X", "S", "R", "Sr", "Rs","D","CumR")
col_names <- c("time")
for(comp in compartment_names) {
  for(age in age_groups) {
    col_names <- c(col_names, paste0(comp, "_", age))
  }
}
colnames(results_1_a_Tanzania) <- col_names
names(results_1_a_Tanzania)
#Outcomes: Baseline Prevalence
R_total_1 <- rowSums(out_1_b_Tanzania[, Rindex + 1])
R_total_5 <- rowSums(out_5_b_Tanzania[, Rindex + 1])
R_total_10 <- rowSums(out_10_b_Tanzania[, Rindex + 1])
#Outcomes: Cumulative Incidence
R_total_1_a <- cumsum(rowSums(out_1_a_Tanzania[, Rindex + 1]))  #annual
R_total_1_b <- cumsum(rowSums(out_1_b_Tanzania[, Rindex + 1]))  # no MDA
R_total_1_c <- cumsum(rowSums(out_1_c_Tanzania[, Rindex + 1]))  # bi-annual

R_total_5_a <- cumsum(rowSums(out_5_a_Tanzania[, Rindex + 1])) 
R_total_5_b <- cumsum(rowSums(out_5_b_Tanzania[, Rindex + 1])) 
R_total_5_c <- cumsum(rowSums(out_5_c_Tanzania[, Rindex + 1])) 

R_total_10_a <- cumsum(rowSums(out_10_a_Tanzania[, Rindex + 1])) 
R_total_10_b <- cumsum(rowSums(out_10_b_Tanzania[, Rindex + 1])) 
R_total_10_c <- cumsum(rowSums(out_10_c_Tanzania[, Rindex + 1])) 
#Total population overtime
N_1_a<-rowSums(out_1_a_Tanzania[, Sindex + 1])+rowSums(out_1_a_Tanzania[, Rindex + 1])+rowSums(out_1_a_Tanzania[, Rsindex + 1])+rowSums(out_1_a_Tanzania[, Srindex + 1])+rowSums(out_1_a_Tanzania[, Xindex + 1])
N_1_b <- rowSums(out_1_b_Tanzania[, Sindex+1]) + rowSums(out_1_b_Tanzania[, Rindex+1]) + rowSums(out_1_b_Tanzania[, Rsindex+1]) + rowSums(out_1_b_Tanzania[, Srindex+1]) + rowSums(out_1_b_Tanzania[, Xindex+1])
N_1_c <- rowSums(out_1_c_Tanzania[, Sindex+1]) + rowSums(out_1_c_Tanzania[, Rindex+1]) + rowSums(out_1_c_Tanzania[, Rsindex+1]) + rowSums(out_1_c_Tanzania[, Srindex+1]) + rowSums(out_1_c_Tanzania[, Xindex+1])
N_5_a <- rowSums(out_5_a_Tanzania[, Sindex+1]) + rowSums(out_5_a_Tanzania[, Rindex+1]) + rowSums(out_5_a_Tanzania[, Rsindex+1]) + rowSums(out_5_a_Tanzania[, Srindex+1]) + rowSums(out_5_a_Tanzania[, Xindex+1])
N_5_b <- rowSums(out_5_b_Tanzania[, Sindex+1]) + rowSums(out_5_b_Tanzania[, Rindex+1]) + rowSums(out_5_b_Tanzania[, Rsindex+1]) + rowSums(out_5_b_Tanzania[, Srindex+1]) + rowSums(out_5_b_Tanzania[, Xindex+1])
N_5_c <- rowSums(out_5_c_Tanzania[, Sindex+1]) + rowSums(out_5_c_Tanzania[, Rindex+1]) + rowSums(out_5_c_Tanzania[, Rsindex+1]) + rowSums(out_5_c_Tanzania[, Srindex+1]) + rowSums(out_5_c_Tanzania[, Xindex+1])
N_10_a <- rowSums(out_10_a_Tanzania[, Sindex+1]) + rowSums(out_10_a_Tanzania[, Rindex+1]) + rowSums(out_10_a_Tanzania[, Rsindex+1]) + rowSums(out_10_a_Tanzania[, Srindex+1]) + rowSums(out_10_a_Tanzania[, Xindex+1])
N_10_b <- rowSums(out_10_b_Tanzania[, Sindex+1]) + rowSums(out_10_b_Tanzania[, Rindex+1]) + rowSums(out_10_b_Tanzania[, Rsindex+1]) + rowSums(out_10_b_Tanzania[, Srindex+1]) + rowSums(out_10_b_Tanzania[, Xindex+1])
N_10_c <- rowSums(out_10_c_Tanzania[, Sindex+1]) + rowSums(out_10_c_Tanzania[, Rindex+1]) + rowSums(out_10_c_Tanzania[, Rsindex+1]) + rowSums(out_10_c_Tanzania[, Srindex+1]) + rowSums(out_10_c_Tanzania[, Xindex+1])
#MDA reduction
N_10_a_5 <- rowSums(out_10_a_5_Tanzania[, Sindex+1]) + rowSums(out_10_a_5_Tanzania[, Rindex+1]) + rowSums(out_10_a_5_Tanzania[, Rsindex+1]) + rowSums(out_10_a_5_Tanzania[, Srindex+1]) + rowSums(out_10_a_5_Tanzania[, Xindex+1])
N_10_a_6 <- rowSums(out_10_a_6_Tanzania[, Sindex+1]) + rowSums(out_10_a_6_Tanzania[, Rindex+1]) + rowSums(out_10_a_6_Tanzania[, Rsindex+1]) + rowSums(out_10_a_6_Tanzania[, Srindex+1]) + rowSums(out_10_a_6_Tanzania[, Xindex+1])
N_10_a_7 <- rowSums(out_10_a_7_Tanzania[, Sindex+1]) + rowSums(out_10_a_7_Tanzania[, Rindex+1]) + rowSums(out_10_a_7_Tanzania[, Rsindex+1]) + rowSums(out_10_a_7_Tanzania[, Srindex+1]) + rowSums(out_10_a_7_Tanzania[, Xindex+1])
N_10_c_5 <- rowSums(out_10_c_5_Tanzania[, Sindex+1]) + rowSums(out_10_c_5_Tanzania[, Rindex+1]) + rowSums(out_10_c_5_Tanzania[, Rsindex+1]) + rowSums(out_10_c_5_Tanzania[, Srindex+1]) + rowSums(out_10_c_5_Tanzania[, Xindex+1])
N_10_c_6 <- rowSums(out_10_c_6_Tanzania[, Sindex+1]) + rowSums(out_10_c_6_Tanzania[, Rindex+1]) + rowSums(out_10_c_6_Tanzania[, Rsindex+1]) + rowSums(out_10_c_6_Tanzania[, Srindex+1]) + rowSums(out_10_c_6_Tanzania[, Xindex+1])
N_10_c_7 <- rowSums(out_10_c_7_Tanzania[, Sindex+1]) + rowSums(out_10_c_7_Tanzania[, Rindex+1]) + rowSums(out_10_c_7_Tanzania[, Rsindex+1]) + rowSums(out_10_c_7_Tanzania[, Srindex+1]) + rowSums(out_10_c_7_Tanzania[, Xindex+1])
#Proportion
prop_R_1_a <- round(rowSums(out_1_a_Tanzania[, Rindex + 1]) / N_1_a * 100, 1)
prop_R_1_b <- round(rowSums(out_1_b_Tanzania[, Rindex + 1]) / N_1_b * 100, 1)
prop_R_1_c <- round(rowSums(out_1_c_Tanzania[, Rindex + 1]) / N_1_c * 100, 1)

prop_R_5_a  <- round(rowSums(out_5_a_Tanzania[, Rindex + 1]) / N_5_a * 100, 1)
prop_R_5_b  <- round(rowSums(out_5_b_Tanzania[, Rindex + 1]) / N_5_b * 100, 1)
prop_R_5_c  <- round(rowSums(out_5_c_Tanzania[, Rindex + 1]) / N_5_c * 100, 1)

prop_R_10_a <- round(rowSums(out_10_a_Tanzania[, Rindex + 1]) / N_10_a * 100, 1)
prop_R_10_b <- round(rowSums(out_10_b_Tanzania[, Rindex + 1]) / N_10_b * 100, 1)
prop_R_10_c <- round(rowSums(out_10_c_Tanzania[, Rindex + 1]) / N_10_c * 100, 1)

prop_R_10_a_5 <- round(rowSums(out_10_a_5_Tanzania[, Rindex + 1]) / N_10_a_5 * 100,1)
prop_R_10_a_6 <- round(rowSums(out_10_a_6_Tanzania[, Rindex + 1]) / N_10_a_6 * 100,1)
prop_R_10_a_7 <- round(rowSums(out_10_a_7_Tanzania[, Rindex + 1]) / N_10_a_7 * 100,1)

prop_R_10_c_5 <- round(rowSums(out_10_c_5_Tanzania[, Rindex + 1]) / N_10_c_5 * 100,1)
prop_R_10_c_6 <- round(rowSums(out_10_c_6_Tanzania[, Rindex + 1]) / N_10_c_6 * 100,1)
prop_R_10_c_7 <- round(rowSums(out_10_c_7_Tanzania[, Rindex + 1]) / N_10_c_7 * 100,1)
#Here
#plot(time, run_d[, 7], type = "l")
dataset_1 <- cbind(time = out_1_a_Tanzania[, 1], R_total_1, R_total_1_a,R_total_1_b,R_total_1_c,N_1_a,N_1_b,N_1_c,prop_R_1_a,prop_R_1_b,prop_R_1_c)
dataset_1 <- as.data.frame(dataset_1)
#
dataset_5 <- cbind(time = out_5_a_Tanzania[, 1], R_total_5, R_total_5_a,R_total_5_b,R_total_5_c,N_5_a,N_5_b,N_5_c,prop_R_5_a,prop_R_5_b,prop_R_5_c)
dataset_5 <- as.data.frame(dataset_5)
#
dataset_10 <- cbind(time = out_10_a_Tanzania[, 1], R_total_10, R_total_10_a,R_total_10_b,R_total_10_c,N_10_a,N_10_b,N_10_c,prop_R_10_a,prop_R_10_b,prop_R_10_c)
dataset_10 <- as.data.frame(dataset_10)
#
dataset_10_567_annual <- cbind(time = out_10_a_Tanzania[, 1],prop_R_10_b,prop_R_10_a,prop_R_10_a_5,prop_R_10_a_6,prop_R_10_a_7)
dataset_10_567_annual <- as.data.frame(dataset_10_567_annual)
colnames(dataset_10_567_annual) <- c(
  "time",
  "Baseline",#prop_R_10_b
  "Ten",     #prop_R_10_a
  "Five",    #prop_R_10_a_5
  "Six",     #prop_R_10_a_6
  "Seven"    #prop_R_10_a_7
  )
#
dataset_10_567_bi_annual <- cbind(time = out_10_c_Tanzania[, 1],prop_R_10_b,prop_R_10_c,prop_R_10_c_5,prop_R_10_c_6,prop_R_10_c_7)
dataset_10_567_bi_annual <- as.data.frame(dataset_10_567_bi_annual)
colnames(dataset_10_567_bi_annual) <- c(
  "time",
  "Baseline",#prop_R_10_b
  "Ten",     #prop_R_10_c
  "Five",    #prop_R_10_c_5
  "Six",     #prop_R_10_c_6
  "Seven"    #prop_R_10_c_7
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

#Re-order:annual
table(annual_long_10_567$Scenario)
annual_long_10_567$Scenario <- factor(
  annual_long_10_567$Scenario,
  levels = c("Baseline","Five","Six" ,"Seven","Ten")
)
#Re-order: bi-annual
table(bi_annual_long_10_567$Scenario)
bi_annual_long_10_567$Scenario <- factor(
  bi_annual_long_10_567$Scenario,
  levels = c("Baseline","Five","Six" ,"Seven","Ten")
)
table(annual_long_10_567$Scenario)
annual<-ggplot(annual_long_10_567, aes(x = time, y = Resistance, color = Scenario)) +
  geom_ribbon(
    data = subset(annual_long_10_567, Scenario == "Baseline"),
    aes(
      x = time,
      ymin = 0,
      ymax =  Resistance
    ),
    inherit.aes = FALSE,
    #fill = "#C2A" ,
   fill = "grey80",
    alpha = 0.4
  ) +
  #geom_point(size = 1) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(min(annual_long_10_567$time), max(annual_long_10_567$time), by = 365)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10))+
  labs(
    title = "A.Annual MDA",
    #title = "Proportion of infection due to Escherichia coli resistant to macrolides  in Tanzania",
    #subtitle = "Age-structured mixed-carriage model integrating within-host bacterial competition,E.Coli",
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
bi_annual<-ggplot(bi_annual_long_10_567, aes(x = time, y = Resistance, color = Scenario)) +
  geom_ribbon(
    data = subset(bi_annual_long_10_567, Scenario == "Baseline"),
    aes(
      x = time,
      ymin = 0,
      ymax =  Resistance
    ),
    inherit.aes = FALSE,
    #fill = "#C2A" ,
    fill = "grey80",
    alpha = 0.4
  ) +
  #geom_point(size = 1) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(min(bi_annual_long_10_567$time), max(bi_annual_long_10_567$time), by = 365)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10))+
  labs(
    title = "B.Bi-annual",
    #title = "Proportion of infection due to Escherichia coli resistant to macrolides  in Tanzania",
    #subtitle = "Age-structured mixed-carriage model integrating within-host bacterial competition,E.Coli",
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
grid.arrange(annual,bi_annual,ncol=2)
combined_plot <- arrangeGrob(annual, bi_annual, ncol = 2)
# Save to a file (PNG, PDF, etc.)
ggsave("combined_annual_and_bi_annual_plot.png", combined_plot, width = 12, height = 6)  # adjust size
#1 year
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
levels(dataset_1_long$scenario)[levels(dataset_1_long$scenario) == "prop_R_1_a"] <- "1 year annual-MDA"
levels(dataset_1_long$scenario)[levels(dataset_1_long$scenario) == "prop_R_1_b"] <- "1 year No-MDA"
levels(dataset_1_long$scenario)[levels(dataset_1_long$scenario) == "prop_R_1_c"] <- "1 year bi-annual-MDA"
table(dataset_1_long$scenario)
# Long dataset for 5-year scenarios
dataset_5_long <- dataset_5 %>%
  dplyr::select(time, prop_R_5_a, prop_R_5_b, prop_R_5_c) %>%
  tidyr::pivot_longer(
    cols = starts_with("prop_R_5_"),
    names_to = "scenario",
    values_to = "Proportion"
  )
#Factor 
dataset_5_long$scenario <- factor(dataset_5_long$scenario)
# Rename levels
levels(dataset_5_long$scenario)[levels(dataset_5_long$scenario) == "prop_R_5_a"] <- "5 year annual-MDA"
levels(dataset_5_long$scenario)[levels(dataset_5_long$scenario) == "prop_R_5_b"] <- "5 year No-MDA"
levels(dataset_5_long$scenario)[levels(dataset_5_long$scenario) == "prop_R_5_c"] <- "5 year bi-annual-MDA"
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
#Rename levels
levels(dataset_10_long$scenario)[levels(dataset_10_long$scenario) == "prop_R_10_a"] <- "10 year annual-MDA"
levels(dataset_10_long$scenario)[levels(dataset_10_long$scenario) == "prop_R_10_b"] <- "10 year No-MDA"
levels(dataset_10_long$scenario)[levels(dataset_10_long$scenario) == "prop_R_10_c"] <- "10 year bi-annual-MDA"
# Check
table(dataset_10_long$scenario)
#Plots
table(dataset_1_long$scenario)
dataset_1_long$scenario <- factor(
  dataset_1_long$scenario,
  levels = c("1 year No-MDA","1 year annual-MDA","1 year bi-annual-MDA")
  )
dataset_5_long$scenario <- factor(
  dataset_5_long$scenario,
  levels = c("5 year No-MDA","5 year annual-MDA","5 year bi-annual-MDA")
)
dataset_10_long$scenario <- factor(
  dataset_10_long$scenario,
  levels = c("10 year No-MDA","10 year annual-MDA","10 year bi-annual-MDA")
)
g_0_1<-ggplot(dataset_1_long, aes(x = time, y = Proportion, color = scenario)) +
  geom_ribbon(
    data = subset(dataset_1_long, scenario == "1 year No-MDA"),
    aes(
      x = time,
      ymin = 0,
      ymax =  Proportion
    ),
    inherit.aes = FALSE,
    fill = "#C2A5CF" ,
    #fill = "grey80",
    alpha = 0.4
  ) +
  #geom_point(size = 1) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = seq(min(dataset_1_long$time), max(dataset_1_long$time), by = 30)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10))+
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
    legend.text = element_text(size = 10)
    )
print(g_0_1)
# Plot
g_0_5 <- ggplot(dataset_5_long, aes(x = time, y = Proportion, color = scenario)) +
  geom_ribbon(
    data = subset(dataset_5_long, scenario == "5 year No-MDA"),
    aes(
      x = time,
      ymin = 0,
      ymax = Proportion
    ),
    inherit.aes = FALSE,
    fill = "#C2A5CF",
    alpha = 0.4
  ) +
  #geom_point(size = 1) +
  geom_line(linewidth = 1) +  # 
  scale_x_continuous(
    breaks = seq(min(dataset_5_long$time), max(dataset_5_long$time), by = 365)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10)) +
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
#Print
print(g_0_5)
# Plot
g_0_10 <- ggplot(dataset_10_long, aes(x = time, y = Proportion, color = scenario)) +
  geom_ribbon(
    data = subset(dataset_10_long, scenario == "10 year No-MDA"),
    aes(
      x = time,
      ymin = 0,
      ymax = Proportion
    ),
    inherit.aes = FALSE,
    fill = "#C2A5CF",
    alpha = 0.4
  ) +
  #geom_point(size = 1) +
  geom_line(linewidth = 1) + 
  scale_x_continuous(
    breaks = seq(min(dataset_10_long$time), max(dataset_10_long$time), by = 365)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 10)) +
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
#Print
print(g_0_10)
print(incidence)
print(b)
graph<-gridExtra::grid.arrange(incidence,g_0_1,g_0_5,g_0_10,nrow=2)
ggsave("Figure 1.Scenario analysis.png",plot = graph,height=14,width = 22,dpi=300)
g1<-ggplot(dataset_1, aes(x = time, y = R_total_1)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +  
  geom_line(color = "blue") +  labs(
    #title = "Resitant infections  in Tanzania",
    title = "A",
    x = "Time (days)",
    y = "Resistant infections"
  )+
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.text = element_text(size = 10)
  )
g1
names(dataset_5)
tail((dataset_5[,2]),1)
g2<-ggplot(dataset_5, aes(x = time, y = R_total_5)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +  
  geom_line(color = "blue") +
  labs(
    #title = "Resitant infections  in Tanzania",
    title = "B",
    x = "Time (days)",
    y = "Resistant infections"
  )+
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.text = element_text(size = 10)
  )
g2
#
names(dataset_10)
tail((dataset_10[,2]),1)
tail((dataset_10[,3]),1)
tail((dataset_10[,4]),1)
tail((dataset_10[,4]),1)
g3<-ggplot(dataset_10, aes(x = time, y = R_total_10)) +
  geom_line(color = "blue") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, accuracy = 1, suffix = " M ")) +  
  labs(
    title = "C",
    #title = "Resitant infections  in Tanzania",
    x = "Time (days)",
    y = "Resistant infections"
  )+
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11),
    #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.text = element_text(size = 10)
  )
g3
gridExtra::grid.arrange(g1,g2,g3,nrow=1)
#
library(grid)
d<-grid.arrange(
  g1, g2, g3,
  nrow = 1,top = "Figure 1: Resistants infections  over time",
  bottom = textGrob(
    "Figure 1.A:Resistants infections after 1 year no-MDA,Figure 1.B:Resistants infections after 5 years no-MDA,,Figure 1.c:Resistants infections after 10 years no-MDA",
    gp = gpar(fontsize = 10, fontface = "italic")
  )
)
ggsave(
  filename = "figure1.png",
  plot = d,
  width = 12,
  height = 6,
  dpi = 300
)
#
g4 <- ggplot(dataset_1, aes(x = time)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, accuracy = 1, suffix = " B")) +  
  geom_vline(
    xintercept = 365.25 / 2,
    linetype = "dotted",
    linewidth = 0.8,
    color = "gray"
  )+
  annotate(
    "text",
    x = 365.25 / 2,
    y = 0.95*max(dataset_1$R_total_1_c, na.rm = TRUE),
    label = "2nd MDA",
    angle = 90,
    vjust = -0.4,
    size = 3.2
  ) +
  geom_line(aes(y = R_total_1_c,
    color = "1 year Bi-MDA",
    linetype = "1 year Bi-MDA"),
    linewidth = 0.9) +
  geom_line(aes(y = R_total_1_a,
    color = "1 year MDA",
    linetype = "1 year MDA"),
    linewidth = 0.9) +
  geom_line(aes(y = R_total_1_b,
    color = "Pre-MDA",
    linetype = "Pre-MDA"),
    linewidth = 0.9) +
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
  geom_line(aes(y = R_total_5_c,
    color = "5 years Bi-MDA",
    linetype = "5 years Bi-MDA"),
    linewidth = 0.9) +
  geom_line(aes(y = R_total_5_a,
    color = "5 years MDA",
    linetype = "5 years MDA"),
    linewidth = 0.9) +
  geom_line(aes(y = R_total_5_b,
    color = "5 years No-MDA",
    linetype = "5 years No-MDA"),
    linewidth = 0.9) +
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
g6 <- ggplot(dataset_10, aes(x = time)) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-9, accuracy = 1, suffix = " B ")) +  
  geom_line(aes(y = R_total_10_c,
    color = "10 years Bi-MDA",
    linetype = "10 years Bi-MDA"),
    linewidth = 0.9) +
  geom_line(aes(y = R_total_10_a,
    color = "10 years MDA",
    linetype = "10 years MDA"),
    linewidth = 0.9) +
  geom_line(aes(y = R_total_10_b,
    color = "10 years No-MDA",
    linetype = "10 years No-MDA"),
    linewidth = 0.9) +
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
g<- gridExtra::grid.arrange(g4,g5, g6, nrow = 3)
print(g)
g8<- gridExtra::grid.arrange(g_0_1,g4, g_0_5,g5, g_0_10, g6, nrow = 3)
g8
g7<-grid.arrange(
  g1, g4, g2,g5,g3,g6,
  nrow = 3,top = "Figure 2_Scenarios_Tanzania.png",
  bottom = textGrob(
    "Figure 2.A:Resistants infections after 1 year no-MDA,Figure 2.B:Resistants infections after 5 years no-MDA,,Figure 2.c:Resistants infections after 10 years no-MDA,Figure 2.D:Cummulative resistants infections after 1 year,Figure 2.E: Cummulative resistants infections after 5 years,,Figure 2.F: Cummulative Resistants infections after 10 years",
    gp = gpar(fontsize = 10, fontface = "italic")
  )
)
g7
ggsave(
  filename = "Figure 2_Scenarios_Tanzania.png",
  plot = g7,
  width = 14,
  height = 8,
  dpi = 300
)
library(grid)
library(gridExtra)
caption_text <- paste(
  "Figure 2.",
  "A: Resistant infections after 1 year no-MDA;",
  "B: Resistant infections after 5 years no-MDA;",
  "C: Resistant infections after 10 years no-MDA;",
  "D: Cumulative resistant infections after 1 year;",
  "E: Cumulative resistant infections after 5 years;",
  "F: Cumulative resistant infections after 10 years."
)
wrapped_caption <- paste(strwrap(caption_text, width = 120), collapse = "\n")
g7 <- grid.arrange(
  g1, g4, g2, g5, g3, g6,
  nrow = 3,
  top = textGrob("Figure 2. Scenarios in Tanzania", 
    gp = gpar(fontsize = 14, fontface = "bold")),
  bottom = textGrob(
    wrapped_caption,
    gp = gpar(fontsize = 10, fontface = "italic")
  )
)
ggsave(
  filename = "Figure_2_Scenarios_Tanzania.png",
  plot = g7,
  width = 14,
  height = 9,   # slightly increase height
  dpi = 300
)
pacman::p_load(grid,gridExtra)
#TInc = run_d[length(run_d[, 7]), 7]
TInc <- tail(cumsum(rowSums(out_1_b_Tanzania[,  Rindex + 1]) ),1) 
#-------------------------------------------------------------------------------
#              Sensitivity analysis
#-------------------------------------------------------------------------------
#1.One-at-a-time sensitivity analysis
print(parms)
names(parms)
parms_list <- list(
  beta.S = 0.164271,
  u.S = 0.03285421,
  u.R = 0.03285421,
  u.c = 0.03285421,  # example
  a = 0.16,
  a.c = 0.16,
  k = 0.5,
  m_contact = m_contact_1y_Tanzania,  # contact matrix
  mda_cycle = 365,
  mda_duration = 30,
  mda_cov = 0.6,
  theta = 0.13,
  kappa = 0,
  r_mda = 0.03054302,
  c = 0.01
)
#.........All parameters.................................................
params_to_test <- c("beta.S","u.R","u.S","k","c","a","u.c","mda_duration","mda_cov")
sens_results <- NULL
start<-Sys.time()
for (param in params_to_test) {
  
  base_val <- parms_list[[param]]
  
  values <- seq(0.8 * base_val, 1.2 * base_val, length.out = 50)
  
  for (val in values) {
    
    parms <- parms_list
    parms[[param]] <- val
    
    out <- bacteria.solve(tvec_1_a, state, parms)
    
    TInc <- tail(cumsum(rowSums(out[, Rindex + 1])), 1)
    
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
end<-Sys.time()
print(end-start)
library(ggplot2)
#Visualization
#1.Histogram
ggplot(sens_results,
  aes(x = cumulative_infections/1e9)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "black") +
  facet_wrap(~parameter, scales = "free") +
  labs(
    x = "Cumulative resitant infections",
    y = "Frequency",
    title = "Local sensitivity analysis across parameters"
  ) +
  theme_bw()
#2.Box plot
ggplot(sens_results,
  aes(x = parameter, y = cumulative_infections/1e9, fill = parameter)) +
  geom_boxplot(fill = "steelblue") +
  #facet_wrap(~parameter, scales = "free_x") +
  facet_wrap(value~parameter, scales = "free_x") +
  labs(
    x = "Parameter",
    y = "Cumulative resitant infections (Billions)"
  ) +
  theme_bw() +
  theme(legend.position = "none")
#3.Scatter plot
ggplot(sens_results,
  aes(x = value, y = cumulative_infections/1e9)) +
  geom_point(fill = "steelblue") +
  facet_wrap(~parameter, scales = "free_x") +
  labs(
    x = "Parameter value",
    y = "Cumulative resistant infections "#(Billions)
  ) +
  theme_bw()
#4.Tornado plot
library(dplyr)
sens_summary <- sens_results %>%
  group_by(parameter) %>%
  summarise(
    min_inf = min(TInc),
    max_inf = max(TInc),
    impact = max_inf - min_inf
  ) %>%
  arrange(desc(impact))
library(dplyr)
library(ggplot2)
ggplot(sens_summary,
  aes(x = reorder(parameter, impact),
    y = impact/1e9)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "Parameter",
    y = "Impact on cumulative infections (Billions)",
    title = "Sensitivity Analysis (Tornado Plot)"
  ) +
  theme_bw()
#......................Each parameters..........................................
betasens <- NULL
# Let's say you want to vary beta.S from 50% to 150% of its value
beta_values <- seq(0.8 * parms_list$beta.S, 1.2 * parms_list$beta.S, length.out = 50)

for (beta_i in beta_values) {
  parms <- parms_list           # copy the full set
  parms$beta.S <- beta_i        # only change beta.S
  # Run your model
  out_1_b_Tanzania <- bacteria.solve(tvec_1_a, state, parms_noMDA)
  # Calculate infections
  #TInc <- tail(rowSums(out_1_a_Tanzania[, Rindex + 1]), 1)
  TInc <- tail(cumsum(rowSums(out_1_b_Tanzania[, Rindex + 1])), 1)
  
  # Store results
  betasens <- rbind(betasens, c(beta_i, TInc))
}
# Data frame for plotting
betasens <- as.data.frame(betasens)
colnames(betasens) <- c("beta.S", "Cumulative_Infections")
# Plot histogram
h_Beta.S<-hist(betasens$Cumulative_Infections, col = "skyblue", main = "Sensitivity on beta.S",ylab="Frequency",xlab="Total resistant infections")
head(betasens)
library(ggplot2)
ggplot(betasens, aes(x = Cumulative_Infections/1e9)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "black") +
  labs(
    x = "Cumulative Infections (Billions)",
    y = "Frequency"
  ) +
  theme_minimal()
ggplot(betasens, aes(x = factor(beta.S), y = Cumulative_Infections/1e9)) +
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
    #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
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
  parms <- parms_list      # Baseline parameters
  parms$u.S <- uS_i        # Here i change only u.S
  # Run the model
  out_1_b_Tanzania <- bacteria.solve(tvec_1_a, state, parms_noMDA)
  # Cumulative resistant infections
  TInc <- tail(cumsum(rowSums(out_1_b_Tanzania[, Rindex + 1]) ), 1)
  # Store results
  usens <- rbind(usens, c(uS_i, TInc))
}
usens <- as.data.frame(usens)
colnames(usens) <- c("u.S", "Cumulative_Infections")
h_uS<-hist(usens$Cumulative_Infections,
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
pacman::p_load(lhs,foreach,doParallel)
library(lhs)
library(foreach)
library(doParallel)
library(deSolve)
# Ensure parms_list has all parameters
parms_list <- list(
  beta.S = 0.164271,
  u.S = 0.03285421,
  u.R = 0.03285421,
  u.c = 0.03285421,  # <- MUST exist
  a = 0.16,
  a.c = 0.16,
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
  beta.S = c(0.5*parms_list$beta.S, 1.5*parms_list$beta.S),
  u.S = c(0.5*parms_list$u.S, 1.5*parms_list$u.S),
  u.R = c(0.5*parms_list$u.R, 1.5*parms_list$u.R),
  a = c(0.5*parms_list$a, 1.5*parms_list$a),
  k = c(0.5*parms_list$k, 1.5*parms_list$k),
  theta = c(0.5*parms_list$theta, 1.5*parms_list$theta),
  r_mda = c(0.5*parms_list$r_mda, 1.5*parms_list$r_mda),
  c = c(0.5*parms_list$c, 1.5*parms_list$c)
)
param_ranges
# LHS sampling
n_sim <- 50
lhs_samples <- randomLHS(n_sim, length(param_ranges))
param_samples <- as.data.frame(sapply(seq_along(param_ranges), function(i){
  x <- param_ranges[[i]]
  lhs_samples[, i] * (x[2] - x[1]) + x[1]
}))
colnames(param_samples) <- names(param_ranges)
# Stop any existing cluster
# Parallel setup
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)
#Export everything needed to workers
results <- foreach(i = 1:n_sim, .combine = c,
  .packages = "deSolve") %dopar% {
    parms_noMDA <- parms_list
    for(p in names(param_ranges)) parms_noMDA[[p]] <- param_samples[i, p]
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
#Rank parameters by influence
coef_df <- as.data.frame(summary(lm_fit)$coefficients)
coef_df <- coef_df[order(abs(coef_df$Estimate), decreasing = TRUE), ]
print(coef_df)
library(ggplot2)
#Intercept for clarity (optional)
coef_plot <- coef_df[rownames(coef_df) != "(Intercept)", ]
# Make a bar plot sorted by absolute effect
ggplot(coef_plot, aes(x = reorder(rownames(coef_plot), abs(Estimate)),
  y = Estimate,
  fill = Estimate)) +
  geom_bar(stat = "identity") +
  coord_flip() +                                            # horizontal bars
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Parameter", y = "Coefficient estimate", fill = "Estimate") +
  theme_minimal() +
  theme(text = element_text(size = 12))
#...............................................................................
c_vector<-round(seq(0,1.5,by=0.25)*parms.orig[["c"]],4)
k_vector<-round(seq(0,1.5,by=0.25)*parms.orig[["k"]],4)
beta_vector<-seq(0,1.5,by=0.25)*parms.orig[["beta.S"]]
result<-matrix(0,nrow = length(c_vector),length(k_vector))
result_beta<-matrix(0,nrow = length(beta_vector),length(k_vector))
c_vector
k_vector
result
result_beta
# A.Loop over C and K values
for (i in 1:length(c_vector)){
  for (j in 1:length(k_vector)){
  #for (j in 1:length(beta_vector)){
    parms <- parms.orig          # copy baseline
    parms[["c"]] <- c_vector[i]  # update c
    parms[["k"]] <- k_vector[j]  # update k
    #parms[["beta.S"]] <- beta_vector[j]  # update beta
   
    
    #(mda_start_times<-(0:50)*365.25)             # Early MDA start
    #out_1_a_Tanzania <- bacteria.solve(tvec_1_a, state, parms)
    
    (mda_start_times<-(0:50)*365.25)             # Delay MDA start
    out_1_b_Tanzania <- bacteria.solve(tvec_1_b, state, parms_noMDA)
    
    #(mda_start_times<-(0:50)*365.25/2)             # Delay MDA start
    #out_1_c_Tanzania <- bacteria.solve(tvec_1_a, state, parms)
    
    
    #(mda_start_times<-(0:50)*365.25)             # Early MDA start
    #out_5_a_Tanzania <- bacteria.solve(tvec_5_a, state, parms)
    
    #(mda_start_times<-(6:50)*365.25)             # Delay MDA start
    #out_5_b_Tanzania <- bacteria.solve(tvec_5_b, state, parms)
   
    #(mda_start_times<-(0:50)*365.25)             # Early MDA start
    #out_10_a_Tanzania <- bacteria.solve(tvec_10_a, state, parms)
    
    #(mda_start_times<-(51:52)*365.25)             # Delay MDA start
    #out_10_b_Tanzania <- bacteria.solve(tvec_10_b, state, parms)
    # Store R at equilibrium
    #result[i,j] <-tail(((rowSums(out_1_a_Tanzania[, Rindex + 1]))*100)  /(rowSums(out_1_a_Tanzania[, Sindex + 1])+rowSums(out_1_a_Tanzania[, Rindex + 1])+rowSums(out_1_a_Tanzania[, Rsindex + 1])+rowSums(out_1_a_Tanzania[, Srindex + 1])+rowSums(out_1_a_Tanzania[, Xindex + 1])), 1)
    result[i,j] <- tail(((rowSums(out_1_b_Tanzania[, Rindex + 1]))*100)  /(rowSums(out_1_b_Tanzania[, Sindex + 1])+rowSums(out_1_b_Tanzania[, Rindex + 1])+rowSums(out_1_b_Tanzania[, Rsindex + 1])+rowSums(out_1_b_Tanzania[, Srindex + 1])+rowSums(out_1_b_Tanzania[, Xindex + 1])), 1)
    #result[i,j] <- tail(((rowSums(out_5_a_Tanzania[, Rindex + 1]))*100) / (rowSums(out_5_a_Tanzania[, Sindex + 1]) + rowSums(out_5_a_Tanzania[, Rindex + 1]) + rowSums(out_5_a_Tanzania[, Rsindex + 1]) + rowSums(out_5_a_Tanzania[, Srindex + 1]) + rowSums(out_5_a_Tanzania[, Xindex + 1])), 1)
    #result[i,j] <- tail(((rowSums(out_5_b_Tanzania[, Rindex + 1]))*100) / (rowSums(out_5_b_Tanzania[, Sindex + 1]) + rowSums(out_5_b_Tanzania[, Rindex + 1]) + rowSums(out_5_b_Tanzania[, Rsindex + 1]) + rowSums(out_5_b_Tanzania[, Srindex + 1]) + rowSums(out_5_b_Tanzania[, Xindex + 1])), 1)
    #result[i,j] <- tail(((rowSums(out_10_a_Tanzania[, Rindex + 1]))*100)/ (rowSums(out_10_a_Tanzania[, Sindex + 1]) + rowSums(out_10_a_Tanzania[, Rindex + 1]) + rowSums(out_10_a_Tanzania[, Rsindex + 1]) + rowSums(out_10_a_Tanzania[, Srindex + 1]) + rowSums(out_10_a_Tanzania[, Xindex + 1])), 1)
    #result[i,j] <- tail(((rowSums(out_10_b_Tanzania[, Rindex + 1]))*100)/ (rowSums(out_10_b_Tanzania[, Sindex + 1]) + rowSums(out_10_b_Tanzania[, Rindex + 1]) + rowSums(out_10_b_Tanzania[, Rsindex + 1]) + rowSums(out_10_b_Tanzania[, Srindex + 1]) + rowSums(out_10_b_Tanzania[, Xindex + 1])), 1)
    # This are results for beta
    #result_beta[i,j] <- tail(rowSums(out_10_a_Tanzania[, Rindex + 1]), 1)
  }
}
#Here i will store R at equilibrium in result[I,j]
contour(result, xaxt = "n", yaxt = "n",ylab="k",xlab="c",nlevels = 10,col = hcl.colors(14, "Temps"),lwd =2)
axis(1, at=1:length(c_vector)/length(c_vector), labels=c_vector)
axis(2, at=1:length(k_vector)/length(k_vector), labels=k_vector)
#Advanced visualization 
pacman::p_load(ggplot2,reshape2,viridis)
# Convert result matrix to long format
df <- melt(result)
#df <- melt(result_beta)
names(df) <- c("k_index", "c_index", "R_equilibrium")
#names(df) <- c("beta_index", "c_index", "R_equilibrium")
# Map indices to actual parameter values
df$c <- c_vector[df$c_index]
df$k <- k_vector[df$k_index]
#df$beta <- beta_vector[df$beta_index]
# Heatmap with custom scales
ggplot(df, aes(x = c, y = k, fill = R_equilibrium)) +
  geom_tile(color = "white") +  # grid
  scale_fill_viridis_c(option = "plasma", name = "Resitance(in%)") +
  geom_text(aes(label = round(R_equilibrium, 2)), size = 3, color = "white") +
  scale_x_continuous(breaks = c_vector) +  # custom x-axis ticks
  scale_y_continuous(breaks = k_vector) +  # custom y-axis ticks
  labs(
    title = "Prevalence of E.Coli's resistance at Equilibrium after 5 Years of MDA",
    x = "Cost of resistance (c)",
    y = "Co-colonisation efficiency (k)"
  ) +
  theme_classic()+
  #theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
#B. Loop over C and Beta
c_vector<-seq(0.1,1,by=0.25)*parms.orig[["c"]]
beta_vector<-seq(0.1,1,by=0.25)*parms.orig[["beta.S"]]
result_beta<-matrix(0,nrow = length(beta_vector),length(k_vector))
c_vector
beta_vector
result_beta

for (i in 1:length(c_vector)){
  for (j in 1:length(beta_vector)){
    parms <- parms.orig          # copy baseline
    parms[["c"]] <- c_vector[i]  # update c
    #parms[["k"]] <- k_vector[j]  # update k
    parms[["beta.S"]] <- beta_vector[j]  # update beta
    
    #(mda_start_times<-(0:50)*365.25)             # Early MDA start
    #out_1_a_Tanzania <- bacteria.solve(tvec_1_a, state, parms)
    
    #(mda_start_times<-(2:50)*365.25)             # Delay MDA start
    #out_1_b_Tanzania <- bacteria.solve(tvec_1_b, state, parms)
    
    (mda_start_times<-(0:50)*365.25)             # Early MDA start
    out_5_a_Tanzania <- bacteria.solve(tvec_5_a, state, parms)
    
    #(mda_start_times<-(6:50)*365.25)             # Delay MDA start
    # out_5_b_Tanzania <- bacteria.solve(tvec_5_b,state, parms)
    
    #(mda_start_times<-(0:50)*365.25)              # Early MDA start
    # out_10_a_Tanzania <- bacteria.solve(tvec_10_a, state, parms)
    
    #(mda_start_times<-(51:52)*365.25)             # Delay MDA start
    #out_10_b_Tanzania <- bacteria.solve(tvec_10_b, state, parms)
    
    # Store R at equilibrium:beta and c
    #result_beta[i,j] <- tail(((rowSums(out_1_a_Tanzania[, Rindex + 1])+rowSums(out_1_a_Tanzania[, Rsindex + 1]))*100)/(rowSums(out_1_a_Tanzania[, Sindex + 1])+rowSums(out_1_a_Tanzania[, Rindex + 1])+rowSums(out_1_a_Tanzania[, Rsindex + 1])+rowSums(out_1_a_Tanzania[, Srindex + 1])+rowSums(out_1_a_Tanzania[, Xindex + 1])), 1)
    #result_beta[i,j] <- tail(((rowSums(out_1_b_Tanzania[, Rindex + 1])+rowSums(out_1_b_Tanzania[, Rsindex + 1]))*100)/(rowSums(out_1_b_Tanzania[, Sindex + 1])+rowSums(out_1_b_Tanzania[, Rindex + 1])+rowSums(out_1_b_Tanzania[, Rsindex + 1])+rowSums(out_1_b_Tanzania[, Srindex + 1])+rowSums(out_1_b_Tanzania[, Xindex + 1])), 1)
    result_beta[i,j] <- tail(((rowSums(out_5_a_Tanzania[, Rindex + 1]) + rowSums(out_5_a_Tanzania[, Rsindex + 1])) * 100) / (rowSums(out_5_a_Tanzania[, Sindex + 1]) + rowSums(out_5_a_Tanzania[, Rindex + 1]) + rowSums(out_5_a_Tanzania[, Rsindex + 1]) + rowSums(out_5_a_Tanzania[, Srindex + 1]) + rowSums(out_5_a_Tanzania[, Xindex + 1])), 1)
    #result_beta[i,j] <- tail(((rowSums(out_5_b_Tanzania[, Rindex + 1]) + rowSums(out_5_b_Tanzania[, Rsindex + 1])) * 100) / (rowSums(out_5_b_Tanzania[, Sindex + 1]) + rowSums(out_5_b_Tanzania[, Rindex + 1]) + rowSums(out_5_b_Tanzania[, Rsindex + 1]) + rowSums(out_5_b_Tanzania[, Srindex + 1]) + rowSums(out_5_b_Tanzania[, Xindex + 1])), 1)
    #result_beta[i,j] <- tail(((rowSums(out_10_a_Tanzania[, Rindex + 1]) + rowSums(out_10_a_Tanzania[, Rsindex + 1])) * 100) / (rowSums(out_10_a_Tanzania[, Sindex + 1]) + rowSums(out_10_a_Tanzania[, Rindex + 1]) + rowSums(out_10_a_Tanzania[, Rsindex + 1]) + rowSums(out_10_a_Tanzania[, Srindex + 1]) + rowSums(out_10_a_Tanzania[, Xindex + 1])), 1)
    # result_beta[i,j] <- tail(((rowSums(out_10_b_Tanzania[, Rindex + 1]) + rowSums(out_10_b_Tanzania[, Rsindex + 1])) * 100) / (rowSums(out_10_b_Tanzania[, Sindex + 1]) + rowSums(out_10_b_Tanzania[, Rindex + 1]) + rowSums(out_10_b_Tanzania[, Rsindex + 1]) + rowSums(out_10_b_Tanzania[, Srindex + 1]) + rowSums(out_10_b_Tanzania[, Xindex + 1])), 1)
    
  }
}
#Here i wil store R at equilibrium in result[I,j]
contour(result_beta, xaxt = "n", yaxt = "n",ylab="k",xlab="c",nlevels = 10,col = hcl.colors(14, "Temps"),lwd =2)
axis(1, at=1:length(c_vector)/length(c_vector), labels=c_vector)
axis(2, at=1:length(beta_vector)/length(beta_vector), labels=beta_vector)
#Advanced visualization 
pacman::p_load(ggplot2,reshape2,viridis)
# Convert result matrix to long format
df_beta <- melt(result_beta)
names(df_beta) <- c("beta_index", "c_index", "R_beta_equilibrium")
# Map indices to actual parameter values
df_beta$c <- c_vector[df_beta$c_index]
df_beta$beta <- beta_vector[df_beta$beta_index]
# Heatmap with custom scales
ggplot(df_beta, aes(x = c, y = beta, fill = R_beta_equilibrium)) +
  geom_tile(color = "white") +  # grid
  scale_fill_viridis_c(option = "viridis", name = "Resitance(in%)") +
  geom_text(aes(label = round(R_beta_equilibrium, 2)), size = 3, color = "white") +
  scale_x_continuous(breaks = c_vector) +  # custom x-axis ticks
  scale_y_continuous(breaks = beta_vector) +  # custom y-axis ticks
  labs(
    title = "Prevalence of E.Coli's resistance at Equilibrium after 50Y No MDA",
    x = "Cost of resistance (c)",
    y = "Transmission rate (Beta.S)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
#...............................................................................
R_total_10_a <- (rowSums(out_10_a_Tanzania[,  Rindex + 1]) +rowSums(out_10_a_Tanzania[,  Rsindex + 1]))  
R_total_10_b <- (rowSums(out_10_b_Tanzania[,  Rindex + 1]) +rowSums(out_10_b_Tanzania[,  Rsindex + 1])) 
Rs_total_10_a <- rowSums(out_10_a_Tanzania[,  Rsindex + 1]) 
Rs_total_10_b <- rowSums(out_10_b_Tanzania[,  Rsindex + 1])
#Visualization
plot(
  out_1_b_Tanzania$time, X_total,
  type = "l",
  col = cols[1],
  las = 1,
  xaxs = "i", yaxs = "i",
  ylim = c(0,max(c(X_total, S_total, R_total,Rs_total, Sr_total))),
  bty = "n",
  lwd = 3.5,
  xlab = "Day", ylab = "Population (in Millions)",
  #xlab = "Day", ylab = "Proportion",
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
  #legend = c("X_total", "S_total", "R_total","Rs_total", "Sr_total"),
  legend = c("X_total", "S_total","R_total", "Rs_total", "Sr_total"),
  lty = 1,
  lwd = 3.5,
  ncol = 1
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
  xlab = "Day", ylab = "Population (R) in millions",
  #xlab = "Day", ylab = "Proportion (R)",
  main = "Post 10 Years MDA population level Bacterial Dynamics in Tanzania"
)
abline(v = 365*1, col = "red", lwd = 2, lty = 2)
# Add text above or next to the vertical line
text(
  x = 3000,
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
  ylim = c(0,2),#max(c(Rs_total_10_a, Rs_total_10_b))),
  bty = "n",
  lwd = 3.5,
  #xlab = "Day", ylab = "Proportion(Rs)",
  xlab = "Day", ylab = "Population(Rs) in millions",
  main = "Population level Bacterial Dynamics over time"
)
abline(v = 365*1, col = "red", lwd = 2, lty = 2)
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
#              Additional-visualization                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
pacman::p_load(data.table)   # This package will allow us to reshape data set faster
results_1_Tanzania<- as.data.table(results_1_b_Tanzania)
results_1_Tanzania <- results_1_Tanzania[, 1:(5 * n_age + 1)]
#results_1_Tanzania<- as.data.table(results_1_b_Tanzania)
#results_1_Tanzania<- as.data.table(results_1_c_Tanzania)
#Long format : Here i will be using melt to be faster
results_1_Tanzania_long <- melt(
  results_1_Tanzania,
  id.vars = "time",
  variable.name = "variable",
  value.name = "value"
)
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
#Structure of the long format
names(results_1_Tanzania_long)
head(results_1_Tanzania_long)
table(results_1_Tanzania_long$compartment)
table(results_1_Tanzania_long$age_group)
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
#second checks
names(results_1_Tanzania_long)
head(results_1_Tanzania_long)
table(results_1_Tanzania_long$compartment)
table(results_1_Tanzania_long$age_group)
table(results_1_Tanzania_long$compartment)
R_RS_only <- results_1_Tanzania_long %>%
  filter(compartment %in% c("R")) #Rs
plot_a_0 <- ggplot(R_RS_only, 
                 aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "stack",col= NA) +#or dodge
  scale_y_continuous(
    limits = c(0, 60),
    breaks = seq(0, 60, by = 20),
    expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = " Population level Bacterial Dynamics over time in Tanzania", x = "Age", y = "Percent ") +
  scale_fill_manual(values = c("X" = "#00c4aa", "S" = "#e573f3", 
                               "R" = "#00b3f4", "Sr" = "#9b9602", "Rs" = "#fc726c")
    )
print(plot_a_0)
plot_a_0 <- ggplot(R_RS_only, 
  aes(x = proportion , y = age_group, fill = compartment)) +
  geom_col(position = "stack",col= NA) +#or dodge
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = " Population level Bacterial Dynamics over time in Tanzania", x = "Age", y = "Percent ") +
  scale_fill_manual(values = c("X" = "#00c4aa", "S" = "#e573f3", 
    "R" = "#00b3f4", "Sr" = "#9b9602", "Rs" = "#fc726c"))
print(plot_a_0)
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
results_1_a_Tanzania<-results_1_a_Tanzania[,1:(5*n_age+1)]
results_1_b_Tanzania<-results_1_b_Tanzania[,1:(5*n_age+1)]
results_1_c_Tanzania<-results_1_c_Tanzania[,1:(5*n_age+1)]

results_5_a_Tanzania<-results_5_a_Tanzania[,1:(5*n_age+1)]
results_5_b_Tanzania<-results_5_b_Tanzania[,1:(5*n_age+1)]
results_5_c_Tanzania<-results_5_c_Tanzania[,1:(5*n_age+1)]

results_10_a_Tanzania<-results_10_a_Tanzania[,1:(5*n_age+1)]
results_10_b_Tanzania<-results_10_b_Tanzania[,1:(5*n_age+1)]
results_10_c_Tanzania<-results_10_c_Tanzania[,1:(5*n_age+1)]

scenario_1yr_no_MDA_Tanzania <- process_scenario(results_1_b_Tanzania, "Pre-MDA")
scenario_1yr_MDA_Tanzania <- process_scenario(results_1_a_Tanzania, "1 Year MDA")
scenario_1yr_Bi_MDA_Tanzania <- process_scenario(results_1_c_Tanzania, "1 Year Bi-MDA")

scenario_5yr_no_MDA_Tanzania <- process_scenario(results_5_b_Tanzania, "5 Years No MDA")
scenario_5yr_MDA_Tanzania <- process_scenario(results_5_a_Tanzania, "5 Years MDA")
scenario_5yr_Bi_MDA_Tanzania <- process_scenario(results_5_c_Tanzania, "5 Years Bi-MDA")

scenario_10yr_no_MDA_Tanzania <- process_scenario(results_10_b_Tanzania, "10 Years No MDA")
scenario_10yr_MDA_Tanzania <- process_scenario(results_10_a_Tanzania, "10 Years MDA")
scenario_10yr_Bi_MDA_Tanzania <- process_scenario(results_10_c_Tanzania, "10 Years Bi-MDA")
#Scenario labels 
scenario_1yr_no_MDA_Tanzania$scenario<-"Pre-MDA"
scenario_1yr_MDA_Tanzania$scenario<-"1 Year MDA"
scenario_1yr_Bi_MDA_Tanzania$scenario<-"1 Year Bi-MDA"

scenario_5yr_no_MDA_Tanzania$scenario<-"5 Years No MDA"
scenario_5yr_MDA_Tanzania$scenario<-"5 Years MDA"
scenario_5yr_Bi_MDA_Tanzania$scenario<-"5 Years Bi-MDA"

scenario_10yr_no_MDA_Tanzania$scenario<-"10 Years No MDA"
scenario_10yr_MDA_Tanzania$scenario<-"10 Years MDA"
scenario_10yr_Bi_MDA_Tanzania$scenario<-"10 Years Bi-MDA"
#
pacman::p_load(dplyr)
scenario_1yr_no_MDA_Tanzania  <- mutate(scenario_1yr_no_MDA_Tanzania,  scenario = "Pre-MDA")
scenario_1yr_MDA_Tanzania <- mutate(scenario_1yr_MDA_Tanzania,  scenario = "1 Year MDA")
scenario_1yr_Bi_MDA_Tanzania <- mutate(scenario_1yr_Bi_MDA_Tanzania,  scenario = "1 Year Bi-MDA")

scenario_5yr_no_MDA_Tanzania  <- mutate(scenario_5yr_no_MDA_Tanzania,  scenario = "5 Years No MDA")
scenario_5yr_MDA_Tanzania <- mutate(scenario_5yr_MDA_Tanzania,         scenario = "5 Years MDA")
scenario_5yr_Bi_MDA_Tanzania <- mutate(scenario_5yr_Bi_MDA_Tanzania,         scenario = "5 Years Bi-MDA")

scenario_10yr_no_MDA_Tanzania <- mutate(scenario_10yr_no_MDA_Tanzania, scenario = "10 Years No MDA")
scenario_10yr_MDA_Tanzania <- mutate(scenario_10yr_MDA_Tanzania,        scenario = "10 Years MDA")
scenario_10yr_Bi_MDA_Tanzania <- mutate(scenario_10yr_Bi_MDA_Tanzania,        scenario = "10 Years Bi-MDA")

plot_c <- ggplot(scenario_10yr_no_MDA_Tanzania, 
                 aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "dodge",col= NA) +#or dodge
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = " Population level Bacterial Dynamics over time in Tanzania", x = "Age", y = "Percent ") +
  scale_fill_manual(values = c("X" = "#00c4aa", "S" = "#e573f3", 
                               "R" = "#00b3f4", "Sr" = "#9b9602", "Rs" = "#fc726c"))
print(plot_c)
# All scenarios
all_scenarios_Tanzania <- rbind(scenario_1yr_no_MDA_Tanzania,scenario_1yr_MDA_Tanzania,scenario_1yr_Bi_MDA_Tanzania, scenario_5yr_no_MDA_Tanzania,scenario_5yr_MDA_Tanzania,scenario_5yr_Bi_MDA_Tanzania,scenario_10yr_no_MDA_Tanzania,scenario_10yr_MDA_Tanzania,scenario_10yr_Bi_MDA_Tanzania)
table(all_scenarios_Tanzania$scenario)
# Factor levels for scenarios
#all_scenarios_Tanzania$scenario <- factor(all_scenarios_Tanzania$scenario, 
#                                          levels = c("Pre-MDA", "5 Years MDA", "10 Years MDA","10 Years No MDA"))
all_scenarios_Tanzania$scenario <- factor(
  all_scenarios_Tanzania$scenario,
  levels = c(
    "1 Year Bi-MDA",
    "1 Year MDA",
    "Pre-MDA",
    "5 Years Bi-MDA",
    "5 Years MDA",
    "5 Years No MDA",
    "10 Years Bi-MDA",
    "10 Years MDA",
    "10 Years No MDA"
    )
  )
table(all_scenarios_Tanzania$scenario)
names(all_scenarios_Tanzania)
table(all_scenarios_Tanzania$age_group)
#I want to remove some scenarios and  compartments for better visualization
all_scenarios_Tanzania_10 <- all_scenarios_Tanzania |>
  filter(scenario %in% c("10 Years Bi-MDA","10 Years MDA","10 Years No MDA")) |>
  filter(compartment %in% c("R"))
table(all_scenarios_Tanzania_10$scenario)
table(all_scenarios_Tanzania$scenario)

all_scenarios_Tanzania_5 <- all_scenarios_Tanzania |>
  filter(scenario %in% c("5 Years Bi-MDA","5 Years MDA","5 Years No MDA")) |>
  filter(compartment %in% c("R"))
table(all_scenarios_Tanzania_5$scenario)

all_scenarios_Tanzania_1 <- all_scenarios_Tanzania |>
  filter(scenario %in% c("1 Year Bi-MDA","1 Year MDA","Pre-MDA")) |>
  filter(compartment %in% c("R"))
table(all_scenarios_Tanzania_1$scenario)
#Plot 1: Stacked bar plots comparing all three scenarios
table(all_scenarios_Tanzania$age_group,all_scenarios_Tanzania$compartment)
plot_dodged_1 <- ggplot(all_scenarios_Tanzania_1, 
  aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "stack", col = NA) +
  facet_wrap(~scenario, ncol = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
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
    "R" = "#fc726c", 
    "Sr" = "#9b9602", 
    "Rs" = "#00b3f4"
  ))+
  scale_y_continuous(limits = c(0, 100))
print(plot_dodged_1)
#
plot_dodged_5 <- ggplot(all_scenarios_Tanzania_5, 
                       aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "stack", col = NA) +
  facet_wrap(~scenario, ncol = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
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
    "R" = "#fc726c", 
    "Sr" = "#9b9602", 
    "Rs" = "#00b3f4"
  ))+
  scale_y_continuous(limits = c(0, 100))
print(plot_dodged_5)
# Plot 2: Dodged bar plots comparing all three scenarios
plot_dodged_10 <- ggplot(all_scenarios_Tanzania_10, 
                      aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "stack", col = NA) +
  facet_wrap(~scenario, ncol = 3) +
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
    "R" = "#fc726c", 
    "Sr" = "#9b9602", 
    "Rs" = "#00b3f4"
  ))+
  scale_y_continuous(limits = c(0, 100))
print(plot_dodged_10)
#par(mfrow=c(1,3))
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
  filter(compartment %in% c("R", "Rs")) %>%
  ggplot(aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "stack") + #stack or dodge
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
    "Rs" = "pink" #"#00b3f4"
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

plot_combined <- ggplot(dataset, 
  aes(x = age_group, y = proportion, fill = compartment)) +
  geom_col(position = "stack") +
  facet_grid(period ~ scenario) +
  scale_y_continuous(limits = c(0,100)) +
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
#grid.arrange(plot_dodged_1,plot_dodged_5,plot_dodged_10)
# Plot 3: Side-by-side comparison (stacked) - all scenarios in one row
table(all_scenarios_Tanzania$scenario)
all_scenarios_Tanzania$scenario <- factor(all_scenarios_Tanzania$scenario,
  levels = c("Pre-MDA","1Year MDA","1Year Bi-MDA","5 Years Bi-MDA","5 Years MDA","5 Years NO MDA","10 Years Bi-MDA", "10 Years MDA","10 Years No MDA"))
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
#print(plot_horizontal)

# Plot 4: Focus on Resistance (R) only across scenarios

table(all_scenarios_Tanzania$compartment)
table(all_scenarios_Tanzania$scenario)
resistance_only <- all_scenarios_Tanzania[compartment == "R"]
table(all_scenarios_Tanzania$scenario)
table(resistance_only$scenario)
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
  #scale_fill_viridis_d(option = "viridis", end = 0.8)
  scale_fill_manual(values = c(
    "10 Years MDA" = "#d95f02",
    "10 Years No MDA" = "#1b9e77"
  ))
print(plot_resistance)
#
#Plot 5: Total Resistance (R + Rs) comparison
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
#...............................................................................

