runCompetitive <- function(Males, Females, Traps, Disp, ti, rad, rep) {
        nMales <- nrow(Males)
        nFemales <- nrow(Females)
        nTraps <- nrow(Traps)
        nDisp  <- nrow(Disp)
        
        attractors <- rbind(Traps, Females, Disp)  
        
        MatedM <- MatedF <- TrapC <- numeric(ti)
        
        #Loop over mating event time steps
        for (i in 1:ti) {
                
                #reset number of available females,
                
                attractors <- rbind(Traps, subset(attractors, Type=='Females'), Disp)
                
                if(nrow(Males) > 0) { # If there are still males left
                        
                        # return nearest attractor to each Male
                        MtoA <- get.knnx(attractors[,1:2], Males, k=1) 
                        
                        # If male is within attraction radius, set attractor to NA
                        attractors[MtoA$nn.index[which(MtoA$nn.dist <=rad)],1:2] <- NA 
                        
                        # Count those females with a male in their attraction radius
                        MatedF[i] <- length(which(is.na(subset(attractors, Type == 'Females')$X) == TRUE))
                        
                        # Assume the number of females mated is equivalent to the number of males
                        MatedM[i] <- MatedF[i]
                        
                        # Label those males that have been trapped for eventual removal
                        Males[which(MtoA$nn.index >=range(which(attractors$Type == 'Trap'))[1] & 
                                            MtoA$nn.index <=range(which(attractors$Type == 'Trap'))[2]),] <- NA
                        
                        # Count Trapped Males
                        TrapC[i] <- length(which(is.na(Males$X)==TRUE))
                        
                        # Remove Trapped Males
                        Males <- na.omit(Males) 
                        
                        #Remove Mated Females
                        attractors <- na.omit(attractors)
                        

                #Randomly move remaining males
                Males$X <- Males$X+runif(1, min=-0.12, max=0.12)
                Males$Y <- Males$Y+runif(1, min=-0.12, max=0.12)
                
                #Make sure randomly moved males stay in box
                Males$X[which(Males$X > 1)] <- 1
                Males$Y[which(Males$Y > 1)] <- 1 
                Males$X[which(Males$X < 0)] <- 0
                Males$Y[which(Males$Y < 0)] <- 0 
        } else {
                TrapC[i] <- MatedF[i] <- MatedM[i] <- 0
        }
        
        }
        #Sum stats over mating events
        TC=sum(TrapC)
        MF=sum(MatedF)
        MM=sum(MatedM)
        out=data.frame(Males=nMales, Females=nFemales, Disp=nDisp, 
                       Traps=nTraps, Rad=rad, Rep=rep, MateEvent=ti, 
                       'TC'=TC, 'MF'=MF,'MM'= MM)
        return(out)
        
}

runNonCompetitive <- function(Males, Females, Traps, Disp, ti, rad, rep) {
        nMales <- nrow(Males)
        nFemales <- nrow(Females)
        nTraps <- nrow(Traps)
        nDisp  <- nrow(Disp)
        
        
        #Combine Traps and Females as attractors
        attractors=rbind(Traps, Females)
        
        MatedM <- MatedF <- TrapC <- numeric(ti)
        
        #Loop over mating event time steps
        for (i in 1:ti) {
                #identify jammed males within radius of Disp
                
                if(nrow(Males)>0) { # Provided there are still Males left
                        
                        # Get nearest dispenser for each male
                        DtoM <- get.knnx(Disp[,1:2], Males, k=1)
                        
                        # Determine which Males can be influenced by Dispenser
                        # ie are within radius of effectiveness
                        jamM <- Males[which(DtoM$nn.dist <= rad),]
                        
                        # Determine which Males are still in play
                        Males <- Males[which(DtoM$nn.dist > rad),]
                        
                        # Reset attractors
                        attractors=rbind(Traps, subset(attractors, Type=='Females'))
                        
                        if (nrow(Males) > 0) {
                                # Find nearest attractor to each Male
                                MtoA <- get.knnx(attractors[,1:2], Males, k=1)
                                
                                # Set those attractors with a male within its radius of 
                                # Effectiveness to NA
                                attractors[MtoA$nn.index[which(MtoA$nn.dist <=rad)],1:2] <- NA
                                
                                # Count those females with a male in radius of attractiveness
                                # ie those that are mated
                                MatedF[i] <- length(which(is.na(subset(attractors, Type == 'Females')$X) == TRUE))
                                
                                # Assume mated females = mated males
                                MatedM[i] <- MatedF[i]
                                
                                # Set those males within trap range to NA
                                Males[which(MtoA$nn.index >=range(which(attractors$Type == 'Trap'))[1] & 
                                                    MtoA$nn.index <=range(which(attractors$Type == 'Trap'))[2]),] <- NA
                                
                                # Count those males within trap range
                                TrapC[i] <- length(which(is.na(Males$X)==TRUE))
                                
                                # Remove trapped males
                                Males=na.omit(Males)
                                
                                # Remove mated females
                                attractors=na.omit(attractors)
                                
                                # Restore jammed males to full population for next iteration
                                Males=rbind(Males, jamM)        
                        } else {
                                TrapC[i] <- MatedF[i] <- MatedM[i] <- 0
                        }
                } else {
                        TrapC[i] <- MatedF[i] <- MatedM[i] <- 0
                }
                #Randomly move remaining males
                Males$X <- Males$X+runif(1, min=-0.12, max=0.12)
                Males$Y <- Males$Y+runif(1, min=-0.12, max=0.12)
                
                #Make sure randomly moved males stay in box
                Males$X[which(Males$X > 1)] <- 1
                Males$Y[which(Males$Y > 1)] <- 1 
                Males$X[which(Males$X < 0)] <- 0
                Males$Y[which(Males$Y < 0)] <- 0 
        }
        #Sum stats over mating events
        TC <- sum(TrapC)
        MF <- sum(MatedF)
        MM <- sum(MatedM)
        out <- data.frame(Males=nMales, Females=nFemales, Disp=nDisp, 
                       Traps=nTraps, Rad=rad, Rep=rep, MateEvent=ti, 
                       'TC'=TC, 'MF'=MF,'MM'= MM)
        return(out)     
}

runControl <- function(Males, Females, Traps, ti, rad, rep) {
        # Note: This function is identical to the runCompetitive function
        # except for the lack of dispensers. 
        
        nMales <- nrow(Males)
        nFemales <- nrow(Females)
        nTraps <- nrow(Traps)
        nDisp  <- 0
        attractors <- rbind(Traps, Females)  
        
        MatedM <- MatedF <- TrapC <- numeric(ti)
        
        #Loop over mating event time steps
        for (i in 1:ti) {
                
                #reset number of available females,
                
                attractors <- rbind(Traps, subset(attractors, Type=='Females'))
                
                if(nrow(Males) > 0) { # If there are still males left
                        
                        # return nearest attractor to each Male
                        MtoA <- get.knnx(attractors[,1:2], Males, k=1) 
                        
                        # If male is within attraction radius, set attractor to NA
                        attractors[MtoA$nn.index[which(MtoA$nn.dist <=rad)],1:2] <- NA 
                        
                        # Count those females with a male in their attraction radius
                        MatedF[i] <- length(which(is.na(subset(attractors, Type == 'Females')$X) == TRUE))
                        
                        # Assume the number of females mated is equivalent to the number of males
                        MatedM[i] <- MatedF[i]
                        
                        # Label those males that have been trapped for eventual removal
                        Males[which(MtoA$nn.index >=range(which(attractors$Type == 'Trap'))[1] & 
                                            MtoA$nn.index <=range(which(attractors$Type == 'Trap'))[2]),] <- NA
                        
                        # Count Trapped Males
                        TrapC[i] <- length(which(is.na(Males$X)==TRUE))
                        
                        # Remove Trapped Males
                        Males <- na.omit(Males) 
                        
                        #Remove Mated Females
                        attractors <- na.omit(attractors)
                        
                        
                        #Randomly move remaining males
                        Males$X <- Males$X+runif(1, min=-0.12, max=0.12)
                        Males$Y <- Males$Y+runif(1, min=-0.12, max=0.12)
                        
                        #Make sure randomly moved males stay in box
                        Males$X[which(Males$X > 1)] <- 1
                        Males$Y[which(Males$Y > 1)] <- 1 
                        Males$X[which(Males$X < 0)] <- 0
                        Males$Y[which(Males$Y < 0)] <- 0 
                } else {
                        TrapC[i] <- MatedF[i] <- MatedM[i] <- 0
                }
                
        }
        #Sum stats over mating events
        TC <- sum(TrapC)
        MF <- sum(MatedF)
        MM <- sum(MatedM)
        out=data.frame(Males=nMales, Females=nFemales, Disp=nDisp, 
                       Traps=nTraps, Rad=rad, Rep=rep, MateEvent=ti, 
                       'TC'=TC, 'MF'=MF,'MM'= MM)
        return(out)
        
}