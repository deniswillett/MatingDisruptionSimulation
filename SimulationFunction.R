runSimulation <- function(run.Param) {
        #runParam should be a data frame with colnames:
        #Males, Females, Traps, Disp, Rep, MateEvent, Radius
        
        #Divy up appropriate variables for clarity
        nMales <- run.Param$Males #number of males
        nFemales <- run.Param$Females #number of females
        nTraps <- run.Param$Traps # number of traps
        nDisp <- run.Param$Disp # number of Dispensers
        ti <- run.Param$MateEvent # number of mating events/time steps 
        rad <- run.Param$Rad # radius of effectiveness
        rep <- run.Param$Rep # number of replications
        
        # Set Positions
        # Males and Females get randomly placed
        Males <- Rand_Place(n=nMales, Type='Males')[,1:2]
        Females <- Rand_Place(n=nFemales, Type='Females')
        
        # Traps and Dispensers are randomly uniformly spaced using pack circles 
        # function
        #Set up Trap positions
        TrapPos=pack.circles(matrix(c(rad,nTraps), ncol=2), size=c(1,1))[,1:2]
        Traps=data.frame(TrapPos, Type=rep('Trap', nrow(TrapPos)))
        
        #Set up Dispenser Positions
        DispPos=pack.circles(matrix(c(rad,nDisp), ncol=2), size=c(1,1))[,1:2]
        Disp=data.frame(DispPos, Type=rep('Disp', nrow(DispPos)))
        
        Comp <- cbind(runCompetitive(Males, Females, Traps, Disp, ti, rad, rep),
                      Type = 'Comp')
        NonComp <- cbind(runNonCompetitive(Males, Females, Traps, Disp, ti, rad, rep),
                         Type = 'NonComp')
        Control <- cbind(runControl(Males, Females, Traps, ti, rad, rep), Type = 'Control')
        
        out <- rbind(Comp, NonComp, Control)
        return(out)
}