fitPC <- function(dat,alpha,stable=FALSE,alpha.array,B,labels,verbose=F) {
    # ----------------------------------------
    # Statbility selection for CPDAG using PC algorithm
    # Input :
    #     - dat : dataset, which colnames are nodenames
    #     - alpha : tuning parameter for PC algorithm
    #     - stable : if TRUE, selection probabilities are calculated from subsampling
    #     - alpha.array : the set of alpha values for the stability selction
    #     - B : number of subsampling
    #     - labels : The node names
    # Output:
    #     - fit: pcAlgo object for CPDAG
    #     - v: set of vertices for CPDAG
    # ----------------------------------------
    
    n = nrow(dat)
    p = ncol(dat)
    w.upper = which(upper.tri(diag(p)),arr.ind=T)
    w.upper.inv = cbind(w.upper[,2],w.upper[,1])
    
    fit = pc(suffStat = list(C = cor(dat),n = n), indepTest = gaussCItest, alpha = alpha,labels=labels,verbose = verbose)
    
    if (stable) {
        if (verbose) {cat("calculating selection probabilities.....\n")}
        n.sub = floor(n/2)
        
        sub.pc <- function(w,alpha) {
            fit = pc(suffStat = list(C = cor(dat[w,]),n = n.sub), indepTest = gaussCItest, alpha = alpha,p=p,verbose = F)
            A = as(fit@graph,"matrix")
            AA = A + t(A)
            A[AA==2] = 0
            tmp1 = cbind(A[w.upper],A[w.upper.inv],as.numeric(AA[w.upper]==2))
            tmp2 = cbind(as.numeric(rowSums(tmp1)==0),tmp1)
            return(apply(tmp2,1,function(x) unlist(which(x==1)))-1)
        }
        
        pcA = array(0,dim=c(nrow(w.upper),length(alpha.array),B)) # B lists of edges x alphas matrix
        for (b in 1:B) {
            w = sample(1:n,n.sub)
            pcA[,,b] = sapply(alpha.array,function(alpha)sub.pc(w,alpha))
        }
        
        # Calculate selection probabilities #
        pcSel = array(0,dim=c(nrow(w.upper),length(alpha.array),3)) # 3 lists of edges x alphas matrices
        tmp1 = (pcA==1)  ## P(a->b)
        pcSel[,,1] = apply(tmp1,c(1,2),sum)/B
        tmp1 = (pcA==2)  ## P(a<-b)
        pcSel[,,2] = apply(tmp1,c(1,2),sum)/B
        tmp1 = (pcA==3)  ## P(a<->b)
        pcSel[,,3] = apply(tmp1,c(1,2),sum)/B
        
        # Calculate max selection probabilities #
        pcmaxSel = matrix(0,nrow(w.upper),3)
        pcmaxSel[,1] = apply(pcSel[,,1],1,max)
        pcmaxSel[,2] = apply(pcSel[,,2],1,max)
        pcmaxSel[,3] = apply(pcSel[,,3],1,max)
        colnames(pcmaxSel) = c("a->b","a<-b","a<->b")
        
        return(list(fit=fit,pcA=pcA,pcSel=pcSel,pcmaxSel=pcmaxSel,rownames=cbind(labels[w.upper[,1]],labels[w.upper[,2]])))
    }else {return(list(fit=fit))}
}

