###################
###smooth method###
###################


###zhengli

rt=read.table("case.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
rt<-rt[rowMeans(rt)>0,]
rt = as.data.frame(t(rt))


##scale function
rt = x
center <- colMeans(x, na.rm=TRUE)
x <- sweep(x, 2L, center, check.margin=FALSE)
f <- function(v){
  v <- v[!is.na(v)]
  sqrt(sum(v^2) / max(1, length(v) - 1L))
}

scale <- apply(x, 2L, f)


###test data

dat = rt[1:3000]
vars <- colnames(dat)
## covariate
id <- 1:nrow(dat)


####loess smooth improve verssion
minifun <- function(x){
  m = x
  y = c("y")
  calcSSE <<- function(y){
    loessMod <- loess(formula = paste(m, "id", sep = "~"), data=dat, span=y)
    res <- try(loessMod$residuals, silent=T)
    sse <- sum(res^2)
    return(sse)

  }    # Run optim to find span that gives min SSE, starting at 0.5
  opt = optimize(calcSSE,c(0.25,0.75))
  minimum = opt$minimum
  #return(minimum)
  filt = loess(formula = paste(m, "id", sep = "~"), data=dat, span=minimum)$fitted
  return(filt)
}

new.dat <- as.data.frame(lapply(vars, minifun),col.names = colnames(dat))

###time test: 30000 column genes cost 110min




