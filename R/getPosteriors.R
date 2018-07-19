getPosteriors <- function(outlist,nodes) {
    ###################################################
    ## This function calculates posterior graph using median probability model
    ## Input :
    ##    outlist: a list for bms results
    ##    dat : a list which includes responses and design matrices
    #####################################################
    p = length(nodes)
    G = matrix(0,p,p) # adjacency matrix for row -> column
    
    intGlist=newXlist= vector("list",p)
    for (i in 1:p) {
        out = outlist[[i]]
        if (!is(out, "zlm")) {
            out.coef = coef(out,order.by.pip=F)
            covars = matrix(unlist(strsplit(rownames(out.coef),split="_")),byrow=T,ncol=2)
            pip = out.coef[,1]
            
            # Set regulators #
            w = which(covars[,1]=="pa")
            G[match(covars[w,2],nodes),i] = pip[w]
            intG = cbind(covars[-w,],pip[-w])
            colnames(intG) = c("type","gene","pip")
            intGlist[[i]] =intG
        }else{
            covars = matrix(unlist(strsplit(colnames(out$model[,2]),split="_")),byrow=T,ncol=2)
            coef = out$coefficients[2]
            sd = sqrt(out$coef2moments[2]-out$coefficients[2]^2)
            if (coef-2*sd<0 & coef+2*sd>0) {pip=0
            }else{pip=1}
            w = which(covars[,1]=="pa")
            G[match(covars[w,2],nodes),i] = pip[w]
            intG = cbind(covars[-w,],pip[-w])
            colnames(intG) = c("type","gene","pip")
            intGlist[[i]] =intG
        }
    }
    colnames(G)=rownames(G) = nodes
    names(intGlist) = nodes
    return(list(G=G,intGlist=intGlist))
}

