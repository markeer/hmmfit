###

prcomp.default <-
  function(x, retx = TRUE, center = F, scale. = FALSE, tol = NULL, ...)
  {
    chkDots(...)
    x <- as.matrix(x)
    #x <- scale(x, center = center, scale = scale.)
    cen <- attr(x, "scaled:center")
    sc <- attr(x, "scaled:scale")
    if(any(sc == 0))
      stop("cannot rescale a constant/zero column to unit variance")
    s <- svd(x, nu = 0)
    s$d <- s$d / sqrt(max(1, nrow(x) - 1))
    if (!is.null(tol)) {
      ## we get rank at least one even for a 0 matrix.
      rank <- sum(s$d > (s$d[1L]*tol))
      if (rank < ncol(x)) {
        s$v <- s$v[, 1L:rank, drop = FALSE]
        s$d <- s$d[1L:rank]
      }
    }
    dimnames(s$v) <-
      list(colnames(x), paste0("PC", seq_len(ncol(s$v))))
    r <- list(sdev = s$d, rotation = s$v,
              center = if(is.null(cen)) FALSE else cen,
              scale = if(is.null(sc)) FALSE else sc)
    if (retx) r$x <- x %*% s$v
    class(r) <- "prcomp"
    r
  }

data=read.table("input.txt",header=T,sep="\t",row.names=1)   #读取表格
data=t(as.matrix(data))   #矩阵转置
data.class <- rownames(data)

###step1 PCA scale method
PCA_scale <- function(x,pccut=1e-5,){
  center <- colMeans(x, na.rm=TRUE)
  x <- sweep(x, 2L, center, check.margin=FALSE)
  f <- function(v){
    v <- v[!is.na(v)]
    sqrt(sum(v^2) / max(1, length(v) - 1L))
  }
  scale <- apply(x, 2L, f)
  x <- sweep(x, 2L, scale, "/", check.margin=FALSE)
  data.pca <- prcomp.default(x, scale. = F)

  sig = pca.sum$importance[1,]
  m =c()
  for (i in c(1:11)){
    m = c(m,sig[i])
    tmp = rep(0, times=12)
    m = c(m,tmp)

  }
  m=c(m,sig[12])
  sig1 = matrix(m/100,12,12)

}

center <- colMeans(x, na.rm=TRUE)
x <- sweep(x, 2L, center, check.margin=FALSE)
f <- function(v){
  v <- v[!is.na(v)]
  sqrt(sum(v^2) / max(1, length(v) - 1L))
}

scale <- apply(x, 2L, f)

x <- sweep(x, 2L, scale, "/", check.margin=FALSE)


data.pca <- prcomp.default(x, scale. = F)

sig = pca.sum$importance[1,][1:12]/100

sig1 =diag(sig)


#for pca scale
t(pca.sum$x[,1:12] %*% (sig1) %*%t(pca.sum$rotation[,1:12]))

pcamat=pca.sum$x[,1:12] %*% (sig1) %*%t(pca.sum$rotation[,1:12])

pcamat <- sweep(pcamat, 2L, scale, "*", check.margin=FALSE)
pcamat <- sweep(pcamat, 2L, center, "+", check.margin=FALSE)


