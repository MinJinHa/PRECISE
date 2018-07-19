
getBMS <- function(dat,Gmat,predictiveDensity=TRUE) {
    ####################################################
    # This function performs bms for all nodes
    # Input :
    #    dat : a list of data includes response vectors and design matices for all nodes
    #    Gmat : is an andjacency matrix for protein regulators
    #    pd : if TRUE, pred.density function is performed on the bms output
    ####################################################
    n = length(dat$ylist[[1]])
    p = length(dat$ylist)
    outlist= vector("list",p)
    unigenes = names(dat$ylist)
    for (i in 1:p) {
        y = dat$ylist[[i]]
        X = dat$Xlist[[i]]
        typevar = matrix(unlist(strsplit(colnames(dat$Xlist[[i]]),split="_")),ncol=2,byrow=T)
        paname = typevar[which(typevar[,1] == "pa"),2]
        w = match(paname,unigenes)
        mprior.size = rep(0.5,ncol(X))
        mprior.size[typevar[,1]=="pa"] = Gmat[w,i]
        
        w = sort(unique(c(which(apply(is.nan(X),2,sum)==0 & sapply(1:ncol(X),function(x) length(unique(X[,x])))>n*0.1))))
        X = X[,w,drop=F]
        mprior.size = mprior.size[w]
        
        # Exclude correlated collumns #
        rr= ncol(X)
        w.upper = which(upper.tri(diag(rr)),arr.ind=T)
        m.addr = matrix(w.upper[abs(cor(X,use="pairwise.complete.obs")[w.upper])>0.9,],ncol=2)
        if (nrow(m.addr)>0) {
            X = X[,-m.addr[,1]]
            mprior.size = mprior.size[-m.addr[,1]]
        }
        
        if (ncol(X)>1) {
            nmodel = 500
            isnt.working = TRUE
            while (isnt.working & nmodel>=20) {
                cat("nmodel=",nmodel,"\n")
                out = bms(X.data = cbind(y,X), mprior="pip",mprior.size =mprior.size,burn = 1000, iter = 2000,nmodel=nmodel)
                pd = pred.density(out,newdata = X[1,,drop=FALSE])
                isnt.working = is.na(pd$fit)
                nmodel = nmodel - 10
            }
            outlist[[i]]= out
        }else {
            outlist[[i]]=zlm(y~X)
        }
    }
    
    pdlist = vector("list",p)
    if (predictiveDensity) {
        for (i in 1:p) {
            out = outlist[[i]]
            X = dat$Xlist[[i]]
            if(ncol(X)>1) {X = X[,match(rownames(coef(out,order.by.pip=F)),colnames(X))]}
            y = dat$ylist[[i]]
            pdlist[[i]] =  pred.density(out,newdata=X)
        }
    }
    
    return(list(outlist=outlist,pdlist=pdlist))
}

