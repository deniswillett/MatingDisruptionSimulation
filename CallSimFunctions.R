library(ggplot2)
library(plyr)
library(FNN)
library(foreach)
library(doSNOW)
library(rlecuyer)

############# Make Sure WD Is set ############


################# Source Requisite Functions for Simulation ############

# Placement Functions
source('PlacementFunctions.R')

# Different Disruption Functions
source('DisruptionFunctions.R')

# Run Simulation Function
source('SimulationFunction.R')

################### Set Simulation Parameters #############

# Number of simulation replications
nRep=1000

# Range of dispensers
rDisp=c(seq(2,100, by=6), 100)

# Range of radii of effectiveness
rRad=c(0.01, seq(0.05,0.3, by=0.05))

# Range of Replications
rRep=1:nRep

# Range of Males, Females, Traps, Mating Events
rMales  <- 100
rFemales <- 100
rTraps <- 3
rME <- 5

# All combinations of relevant variables for simulation
input <- expand.grid(Males = rMales, Females = rFemales, Traps = rTraps, 
                     Disp = rDisp, Rep = rRep, MateEvent = rME, Rad = rRad)

# Add ID for plyr functions
input.id <- cbind(ID = 1:nrow(input), input)


nCores <- 4

cl <- makeCluster(nCores, type='SOCK')
registerDoSNOW(cl)
clusterEvalQ(cl, library(FNN))
clusterEvalQ(cl, 
             source("PlacementFunctions.R"))
clusterEvalQ(cl, 
             source("DisruptionFunctions.R"))
clusterEvalQ(cl, 
             source("SimulationFunction.R"))
clusterSetupRNG(cl, seed = rep(24, 6))



tm.p <- system.time(
        out.p <- ddply(input.id, .(ID), .fun = runSimulation, .parallel = TRUE)
        )


write.csv(out.p, 'out1test.csv')

stopCluster(cl)


tm.np <- system.time(
        out <- ddply(input.id, .(ID), runSimulation, .progress = 'text')
        )

