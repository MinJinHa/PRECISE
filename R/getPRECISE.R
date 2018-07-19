getPRECISE <- function(outlist,pdlist,nodes,delta) {
    # ----------------------------------------
    # This function is to calculate PRECISE network for all patients in pd
    # Input:
    #     outlist: output file for bms function
    #     pdlist: pred.density(out) function where out is bms output
    #     nodes : the names of nodes (the order is the same as outlist and pdlist)
    # delta :  the tuning parameter for protein activation status
    
    p = length(nodes)
    n = length(pdlist[[1]]$fit)
    statusmat = pmat.pos = pmat.neg = pmat.neu = matrix(ncol=p,nrow=n)
    samplename = names(pdlist[[1]]$fit)
    for (v in 1:p) {
        cat("node",v,"\n")
        if (!all(is.na(pdlist[[v]]$fit))) {
            pd = getPPD(pdlist[[v]],delta)$results
            w = match(as.character(pd[,1]),samplename)
            statusmat[w,v] = pd[,7]
            pmat.pos[w,v] = pd[,4]
            pmat.neg[w,v] = pd[,5]
            pmat.neu[w,v] = pd[,6]
        }
    }
    numch = rowSums(getPosteriors(outlist,nodes)$G>0.5)
    Netscore.pos = apply(t(t(pmat.pos) * (numch+1)),1,sum,na.rm=T)
    Netscore.neg = apply(t(t(pmat.neg) * (numch+1)),1,sum,na.rm=T)
    Netscore.neu = apply(t(t(pmat.neu) * (numch+1)),1,sum,na.rm=T)
    Netscore.mat = data.frame(Netscore.pos,Netscore.neg,Netscore.neu)
    tmp1 = unlist(apply(Netscore.mat,1,which.max))
    score =rep(0,length(tmp1))
    score[tmp1==1]=1
    score[tmp1==2]=-1
    return(list(net=statusmat,net.status=score,score.mat=Netscore.mat,samplename=samplename))
}


