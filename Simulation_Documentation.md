# Mating Disruption Simulation
Denis Willett  
August 12, 2014  





## Background
&nbsp; &nbsp; &nbsp; &nbsp; Miller and Gut have used a probabilistic approach to model the effects of competitive and noncompetitive disruption with pheromone dispensers used in mating disruption (Miller et al. 2006a, 2006b, 2010).  While trap catch under these two scenarios was never compared directly in their papers (they only compared male-female visitation rates, a metric difficult to measure in the field), slight modification of their formulas based on their discussion results in __Figure 1.__  

>__Figure 1:__  Trap catch as a function of dispenser number under competitive and noncompetitive disruption scenarios.  Competitive disruption was modeled as per (Miller et al. 2006a) with
>
>$$Trap Catch= \frac{Trap Density*Male Density}{Female Density*Trap Density*Dispensor Density}$$
>
>Similarly, noncompetitive disruption was modeled as 
>
>$$Trap Catch= \frac{Trap Density*Male Density}{Female Density*Trap Density} * Area Not Camouflaged $$
>
>Densities are considered equivalent to absolute number.  For both scenarios, there were three traps, 100 males, 100 females, and varied dispenser number.  Dispensers and traps were considered to be single female equivalents.  For noncompetitive disruption, the plume area covered by a given dispenser was assumed to be 1/1000 of total area; the remaining area not covered by plumes was considered the area not camouflaged.  All equations and discussion of their formulation can be found in (Miller et al. 2006a).  
>
>
![plot of chunk unnamed-chunk-2](./Simulation_Documentation_files/figure-html/unnamed-chunk-2.png) 
        
&nbsp; &nbsp; &nbsp; &nbsp; The primary insights offered through this analysis is a differentiation of the two scenarios based on field data.  Competitive disruption could be distinguished from noncompetitive disruption by curvature.  In addition, the fitting of field data to these probabilistic curves would yield parameters useful for determining population levels.  These probabilistic models imply that competitive disruption would outperform noncompetitive disruption for control of a given pest, _ceteris paribus_.  
	
&nbsp; &nbsp; &nbsp; &nbsp; In addition, John Byers (2007, 2008, 2012) has used simulations to explore the effects of an effective attraction radius on trap catch and mating disruption.  

&nbsp; &nbsp; &nbsp; &nbsp; The below model is similar in many respects to Byers 2007 model, but simplifies a number of inputs.  

## Methods
&nbsp; &nbsp; &nbsp; &nbsp; This simulation approach consists of a Monte Carlo simulation of an individual agent-based spatially explicit model similar in some aspects to a continuous two-dimensional cellular automata model.  For both disruption scenarios, points representing a given number of males, females, traps, and dispensers were placed in a 1 by 1 box.  Males and females were always randomly placed using a uniform distribution.  Trap and dispenser positions were randomly assigned using circle packing algorithms to ensure uniformity of converage.  Traps, females, and dispensers were assigned a radius of effectiveness within which males could be affected.  

### Assumptions
&nbsp; &nbsp; &nbsp; &nbsp; This model makes a number of assumptions.  Namely: 

1. Females, Traps and Dispensers all have an equivalent radius of effectiveness.  

2.  Males can mate more than once.  

3.  Females are stationary.  

4.  Attractors are perfect, if a male is within the radius of effectiveness of an attractor, it is certain he will be affected by it.  


### Competitive Disruption
&nbsp; &nbsp; &nbsp; &nbsp; For each time step under competitive disruption, males were assigned to their nearest attractor (whether a female, a trap, or a dispenser) if their position was within an attractor’s radius of effectiveness.  Any males assigned to a trap were removed from the simulation (_i.e._ they were captured).  Any females to which a male had been assigned were removed from the simulation (_i.e._ they had mated).  Immediately prior to the end of the time step, any remaining males exhibited random movement.  At the end of each time step, the total number of males caught in traps and the total number of mated females were recorded.  

### Noncompetitive Disruption
&nbsp; &nbsp; &nbsp; &nbsp; For each time step under noncompetitive disruption, any males within the radius of effectiveness of a dispenser were temporarily removed from that time step (simulating camouflage or ‘jamming’).  After that initial step, the simulation proceeded similarly to competitive disruption with available males assigned to the nearest attractor, males assigned to traps permanently removed, and females to which males have been assigned considered ‘mated’ and permanently removed.  As in the model for competitive disruption, remaining males in a given time step exhibited random movement with the number of trapped males and mated females recorded.  

### Simulation

&nbsp; &nbsp; &nbsp; &nbsp; The above models for competitive and noncompetitive disruption were run with 1000 replications over five time steps to simulate multiple mating events over time.  For each replication, total trap catch (the number of all males caught over all time steps) and total mated females (the number of females that mated over all time steps) were reported.  


## Results  

> __Figure 1:__ Competitive and noncompetitive disruption scenarios can yield lines of different curvature depending on the radius of effectiveness (RoE).  
>
>
![plot of chunk unnamed-chunk-3](./Simulation_Documentation_files/figure-html/unnamed-chunk-3.png) 

---------------------------------------


> __Figure 2:__ Response surfaces for competitive (__A__) and noncompetitive (__B__) scenarios.  
>
>
![plot of chunk unnamed-chunk-4](./Simulation_Documentation_files/figure-html/unnamed-chunk-4.png) 


## Code

### Mating Disruption Functions


```r
runCompetitive <- function(Males, Females, Traps, Disp, ti, rad, rep) {
    nMales <- nrow(Males)
    nFemales <- nrow(Females)
    nTraps <- nrow(Traps)
    nDisp <- nrow(Disp)
    
    attractors <- rbind(Traps, Females, Disp)
    
    MatedM <- MatedF <- TrapC <- numeric(ti)
    
    # Loop over mating event time steps
    for (i in 1:ti) {
        
        # reset number of available females,
        
        attractors <- rbind(Traps, subset(attractors, Type == "Females"), Disp)
        
        if (nrow(Males) > 0) {
            # If there are still males left
            
            # return nearest attractor to each Male
            MtoA <- get.knnx(attractors[, 1:2], Males, k = 1)
            
            # If male is within attraction radius, set attractor to NA
            attractors[MtoA$nn.index[which(MtoA$nn.dist <= rad)], 1:2] <- NA
            
            # Count those females with a male in their attraction radius
            MatedF[i] <- length(which(is.na(subset(attractors, Type == "Females")$X) == 
                TRUE))
            
            # Assume the number of females mated is equivalent to the number of males
            MatedM[i] <- MatedF[i]
            
            # Label those males that have been trapped for eventual removal
            Males[which(MtoA$nn.index >= range(which(attractors$Type == "Trap"))[1] & 
                MtoA$nn.index <= range(which(attractors$Type == "Trap"))[2]), 
                ] <- NA
            
            # Count Trapped Males
            TrapC[i] <- length(which(is.na(Males$X) == TRUE))
            
            # Remove Trapped Males
            Males <- na.omit(Males)
            
            # Remove Mated Females
            attractors <- na.omit(attractors)
            
            
            # Randomly move remaining males
            Males$X <- Males$X + runif(1, min = -0.12, max = 0.12)
            Males$Y <- Males$Y + runif(1, min = -0.12, max = 0.12)
            
            # Make sure randomly moved males stay in box
            Males$X[which(Males$X > 1)] <- 1
            Males$Y[which(Males$Y > 1)] <- 1
            Males$X[which(Males$X < 0)] <- 0
            Males$Y[which(Males$Y < 0)] <- 0
        } else {
            TrapC[i] <- MatedF[i] <- MatedM[i] <- 0
        }
        
    }
    # Sum stats over mating events
    TC = sum(TrapC)
    MF = sum(MatedF)
    MM = sum(MatedM)
    out = data.frame(Males = nMales, Females = nFemales, Disp = nDisp, Traps = nTraps, 
        Rad = rad, Rep = rep, MateEvent = ti, TC = TC, MF = MF, MM = MM)
    return(out)
    
}

runNonCompetitive <- function(Males, Females, Traps, Disp, ti, rad, rep) {
    nMales <- nrow(Males)
    nFemales <- nrow(Females)
    nTraps <- nrow(Traps)
    nDisp <- nrow(Disp)
    
    
    # Combine Traps and Females as attractors
    attractors = rbind(Traps, Females)
    
    MatedM <- MatedF <- TrapC <- numeric(ti)
    
    # Loop over mating event time steps
    for (i in 1:ti) {
        # identify jammed males within radius of Disp
        
        if (nrow(Males) > 0) {
            # Provided there are still Males left
            
            # Get nearest dispenser for each male
            DtoM <- get.knnx(Disp[, 1:2], Males, k = 1)
            
            # Determine which Males can be influenced by Dispenser ie are within radius
            # of effectiveness
            jamM <- Males[which(DtoM$nn.dist <= rad), ]
            
            # Determine which Males are still in play
            Males <- Males[which(DtoM$nn.dist > rad), ]
            
            # Reset attractors
            attractors = rbind(Traps, subset(attractors, Type == "Females"))
            
            if (nrow(Males) > 0) {
                # Find nearest attractor to each Male
                MtoA <- get.knnx(attractors[, 1:2], Males, k = 1)
                
                # Set those attractors with a male within its radius of Effectiveness to NA
                attractors[MtoA$nn.index[which(MtoA$nn.dist <= rad)], 1:2] <- NA
                
                # Count those females with a male in radius of attractiveness ie those that
                # are mated
                MatedF[i] <- length(which(is.na(subset(attractors, Type == "Females")$X) == 
                  TRUE))
                
                # Assume mated females = mated males
                MatedM[i] <- MatedF[i]
                
                # Set those males within trap range to NA
                Males[which(MtoA$nn.index >= range(which(attractors$Type == 
                  "Trap"))[1] & MtoA$nn.index <= range(which(attractors$Type == 
                  "Trap"))[2]), ] <- NA
                
                # Count those males within trap range
                TrapC[i] <- length(which(is.na(Males$X) == TRUE))
                
                # Remove trapped males
                Males = na.omit(Males)
                
                # Remove mated females
                attractors = na.omit(attractors)
                
                # Restore jammed males to full population for next iteration
                Males = rbind(Males, jamM)
            } else {
                TrapC[i] <- MatedF[i] <- MatedM[i] <- 0
            }
        } else {
            TrapC[i] <- MatedF[i] <- MatedM[i] <- 0
        }
        # Randomly move remaining males
        Males$X <- Males$X + runif(1, min = -0.12, max = 0.12)
        Males$Y <- Males$Y + runif(1, min = -0.12, max = 0.12)
        
        # Make sure randomly moved males stay in box
        Males$X[which(Males$X > 1)] <- 1
        Males$Y[which(Males$Y > 1)] <- 1
        Males$X[which(Males$X < 0)] <- 0
        Males$Y[which(Males$Y < 0)] <- 0
    }
    # Sum stats over mating events
    TC <- sum(TrapC)
    MF <- sum(MatedF)
    MM <- sum(MatedM)
    out <- data.frame(Males = nMales, Females = nFemales, Disp = nDisp, Traps = nTraps, 
        Rad = rad, Rep = rep, MateEvent = ti, TC = TC, MF = MF, MM = MM)
    return(out)
}

runControl <- function(Males, Females, Traps, ti, rad, rep) {
    # Note: This function is identical to the runCompetitive function except for
    # the lack of dispensers.
    
    nMales <- nrow(Males)
    nFemales <- nrow(Females)
    nTraps <- nrow(Traps)
    nDisp <- 0
    attractors <- rbind(Traps, Females)
    
    MatedM <- MatedF <- TrapC <- numeric(ti)
    
    # Loop over mating event time steps
    for (i in 1:ti) {
        
        # reset number of available females,
        
        attractors <- rbind(Traps, subset(attractors, Type == "Females"))
        
        if (nrow(Males) > 0) {
            # If there are still males left
            
            # return nearest attractor to each Male
            MtoA <- get.knnx(attractors[, 1:2], Males, k = 1)
            
            # If male is within attraction radius, set attractor to NA
            attractors[MtoA$nn.index[which(MtoA$nn.dist <= rad)], 1:2] <- NA
            
            # Count those females with a male in their attraction radius
            MatedF[i] <- length(which(is.na(subset(attractors, Type == "Females")$X) == 
                TRUE))
            
            # Assume the number of females mated is equivalent to the number of males
            MatedM[i] <- MatedF[i]
            
            # Label those males that have been trapped for eventual removal
            Males[which(MtoA$nn.index >= range(which(attractors$Type == "Trap"))[1] & 
                MtoA$nn.index <= range(which(attractors$Type == "Trap"))[2]), 
                ] <- NA
            
            # Count Trapped Males
            TrapC[i] <- length(which(is.na(Males$X) == TRUE))
            
            # Remove Trapped Males
            Males <- na.omit(Males)
            
            # Remove Mated Females
            attractors <- na.omit(attractors)
            
            
            # Randomly move remaining males
            Males$X <- Males$X + runif(1, min = -0.12, max = 0.12)
            Males$Y <- Males$Y + runif(1, min = -0.12, max = 0.12)
            
            # Make sure randomly moved males stay in box
            Males$X[which(Males$X > 1)] <- 1
            Males$Y[which(Males$Y > 1)] <- 1
            Males$X[which(Males$X < 0)] <- 0
            Males$Y[which(Males$Y < 0)] <- 0
        } else {
            TrapC[i] <- MatedF[i] <- MatedM[i] <- 0
        }
        
    }
    # Sum stats over mating events
    TC <- sum(TrapC)
    MF <- sum(MatedF)
    MM <- sum(MatedM)
    out = data.frame(Males = nMales, Females = nFemales, Disp = nDisp, Traps = nTraps, 
        Rad = rad, Rep = rep, MateEvent = ti, TC = TC, MF = MF, MM = MM)
    return(out)
    
}
```

### Simulation Function


```r
runSimulation <- function(run.Param) {
    # runParam should be a data frame with colnames: Males, Females, Traps,
    # Disp, Rep, MateEvent, Radius
    
    # Divy up appropriate variables for clarity
    nMales <- run.Param$Males  #number of males
    nFemales <- run.Param$Females  #number of females
    nTraps <- run.Param$Traps  # number of traps
    nDisp <- run.Param$Disp  # number of Dispensers
    ti <- run.Param$MateEvent  # number of mating events/time steps 
    rad <- run.Param$Rad  # radius of effectiveness
    rep <- run.Param$Rep  # number of replications
    
    # Set Positions Males and Females get randomly placed
    Males <- Rand_Place(n = nMales, Type = "Males")[, 1:2]
    Females <- Rand_Place(n = nFemales, Type = "Females")
    
    # Traps and Dispensers are randomly uniformly spaced using pack circles
    # function Set up Trap positions
    TrapPos = pack.circles(matrix(c(rad, nTraps), ncol = 2), size = c(1, 1))[, 
        1:2]
    Traps = data.frame(TrapPos, Type = rep("Trap", nrow(TrapPos)))
    
    # Set up Dispenser Positions
    DispPos = pack.circles(matrix(c(rad, nDisp), ncol = 2), size = c(1, 1))[, 
        1:2]
    Disp = data.frame(DispPos, Type = rep("Disp", nrow(DispPos)))
    
    Comp <- cbind(runCompetitive(Males, Females, Traps, Disp, ti, rad, rep), 
        Type = "Comp")
    NonComp <- cbind(runNonCompetitive(Males, Females, Traps, Disp, ti, rad, 
        rep), Type = "NonComp")
    Control <- cbind(runControl(Males, Females, Traps, ti, rad, rep), Type = "Control")
    
    out <- rbind(Comp, NonComp, Control)
    return(out)
}
```

### Run Simulation

This uses parallel architecture built into the _foreach_, _snow_, and _plyr_ packages and takes approximately 1.5 hours to run on the current computer with four cores (see session info for more details on the computing environment).  


```r
library(ggplot2)
library(plyr)
library(FNN)
library(foreach)
library(doSNOW)
library(rlecuyer)

############# Make Sure WD Is set ############


################# Source Requisite Functions for Simulation ############

# Placement Functions
source("PlacementFunctions.R")

# Different Disruption Functions
source("DisruptionFunctions.R")

# Run Simulation Function
source("SimulationFunction.R")

################### Set Simulation Parameters #############

# Number of simulation replications
nRep = 1000

# Range of dispensers
rDisp = c(seq(2, 100, by = 6), 100)

# Range of radii of effectiveness
rRad = c(0.01, seq(0.05, 0.3, by = 0.05))

# Range of Replications
rRep = 1:nRep

# Range of Males, Females, Traps, Mating Events
rMales <- 100
rFemales <- 100
rTraps <- 3
rME <- 5

# All combinations of relevant variables for simulation
input <- expand.grid(Males = rMales, Females = rFemales, Traps = rTraps, Disp = rDisp, 
    Rep = rRep, MateEvent = rME, Rad = rRad)

# Add ID for plyr functions
input.id <- cbind(ID = 1:nrow(input), input)


nCores <- 4

cl <- makeCluster(nCores, type = "SOCK")
registerDoSNOW(cl)
clusterEvalQ(cl, library(FNN))
clusterEvalQ(cl, source("PlacementFunctions.R"))
clusterEvalQ(cl, source("DisruptionFunctions.R"))
clusterEvalQ(cl, source("SimulationFunction.R"))
clusterSetupRNG(cl, seed = rep(24, 6))



tm.p <- system.time(out.p <- ddply(input.id, .(ID), .fun = runSimulation, .parallel = TRUE))


write.csv(out.p, "out1test.csv")

stopCluster(cl)
```


### Session Info


```
## R version 3.1.1 (2014-07-10)
## Platform: x86_64-apple-darwin10.8.0 (64-bit)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] rlecuyer_0.3-3  doSNOW_1.0.12   snow_0.3-13     iterators_1.0.7
##  [5] foreach_1.4.2   FNN_1.1         plot3D_1.0-1    mgcv_1.8-2     
##  [9] nlme_3.1-117    dplyr_0.2       reshape2_1.4    ggplot2_1.0.0  
## [13] knitr_1.6      
## 
## loaded via a namespace (and not attached):
##  [1] assertthat_0.1   codetools_0.2-8  colorspace_1.2-4 digest_0.6.4    
##  [5] evaluate_0.5.5   formatR_0.10     grid_3.1.1       gtable_0.1.2    
##  [9] htmltools_0.2.4  labeling_0.2     lattice_0.20-29  magrittr_1.0.1  
## [13] MASS_7.3-33      Matrix_1.1-4     misc3d_0.8-4     munsell_0.4.2   
## [17] parallel_3.1.1   plyr_1.8.1       proto_0.3-10     Rcpp_0.11.2     
## [21] rmarkdown_0.2.49 scales_0.2.4     stringr_0.6.2    tools_3.1.1     
## [25] yaml_2.1.13
```


