getPPD <- function(pd,delta){
    # ----------------------------------------
    # This function is to calculate posterior predictive distribution (delta interval) of each of the proteins for each patient
    # Input:
    #    pd : the output of pred.density(bms) function
    #    delta : the tuning parameter for protein activation status
    # Output : various quantities
    n = length(pd$fit)
    
    probmat = matrix(0,nrow=n,ncol=3)
    colnames(probmat) = c("positive","negative","neutral")
    
    pdden = pd$densities()
    for(i in 1:n){
        xx = pdden[[i]]$x
        yy = pdden[[i]]$y
        syy = sum(yy)
        probmat[i,1] = sum(yy[xx>delta])/syy
        probmat[i,2] = sum(yy[xx < -1*delta])/syy
    }
    probmat[,3] = 1 - probmat[,1] - probmat[,2]
    w.max = apply(probmat,1,which.max)
    active = vector(length=n)
    active[w.max==1] = 1
    active[w.max==2] = -1
    active[w.max==3] = 0
    
    results = data.frame(names(pd$fit),pd$fit,pd$std.err,probmat,active)
    colnames(results) = c("samplename","FittedValue","sd","PositiveProb","NegativeProb","NeutralProb","ActivationStatus")
    return(list(results=results,pd.fit=pd))
}

