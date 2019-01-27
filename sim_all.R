# title: Sim all 
# date: 1/2/2019
# updated: 
# purpose: Run all simulation settings 

# required packages
library(dplyr)
library(CVXR)

# source functions 
source('fairness_functions_sim.R')
source('run_simulation.R') # sim function

# read in simulation data
simdata<-read.csv('simulation_data.csv')

# 2 different sample sizes 
n1 = 1000
n2 = 10000

# Number of simulations
nsims=500

#### Scenario 1 - Y is a function of all variables (Y1)
keepvars<-c('Y1','X1','X2','X3','X5','X6','X7','X8','X9','A1','A2')
sim1data<-simdata[,keepvars]

# Run simulation sim(scenario, df, simn, nsims, yvar, avar, avar2)
sim1_1000<-sim(1,sim1data, n1, nsims, 'Y1', 'A1', 'A2')
sim1_10000<-sim(1,sim1data, n2, nsims, 'Y1', 'A1', 'A2')

#### Scenario 2 - Y is a function of some variables (Y2)
keepvars<-c('Y2','X1','X2','X3','X4','X5','X6','X7','X8','X9','A1','A2')
sim2data<-simdata[,keepvars]

sim2_1000<-sim(2,sim2data, n1, nsims, 'Y2', 'A1', 'A2')
sim2_10000<-sim(2,sim2data, n2, nsims, 'Y2', 'A1', 'A2')

#### Scenario 3 - Y is a function of some variables (Y3) but not all available for prediction
keepvars<-c('Y2','X1','X4','X6','X7','X8','A1','A2')
sim3data<-simdata[,keepvars]

sim3_1000<-sim(3,sim3data, n1, nsims, 'Y2', 'A1', 'A2')
sim3_10000<-sim(3,sim3data, n2, nsims, 'Y2', 'A1', 'A2')

sim1_1000$id<-'sim 1 n=1000'
sim1_10000$id<-'sim 1 n=10000'
sim2_1000$id<-'sim 2 n=1000'
sim2_10000$id<-'sim 2 n=10000'
sim3_1000$id<-'sim 3 n=1000'
sim3_10000$id<-'sim 3 n=10000'

simout<-rbind(sim1_1000,sim1_10000,sim2_1000,sim2_10000,sim3_1000,sim3_10000)
write.csv(simout, 'simresults.csv')







