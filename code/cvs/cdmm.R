cdmm <- function(y, x, lam, nlam=100, rlam=1/nlam, mu=1, std=TRUE, maxv=0.4*length(y), maxit=c(20, 50), tol=c(1e-4, 1e-7)) {
	if (std) {
		y <- scale(y, scale=FALSE)
		x <- scale(x, scale=apply(x, 2, sd)*sqrt(nrow(x)-1))
		fac <- 1/attr(x, "scaled:scale")
	} else
		fac <- rep(1, ncol(x))
	if (missing(lam)) {
		lam.max <- max(abs(crossprod(x, y)))
		lam <- lam.max*exp(seq(0, log(rlam), length=nlam))
	}
	res <- .Call("cdmm", y, x, fac, lam, mu, maxv, as.integer(maxit), tol)
	names(res) <- c("sol", "lam")
	res$sol <- res$sol*fac
	if (std) res$int <- attr(y, "scaled:center") - drop(crossprod(res$sol, attr(x, "scaled:center")))
	res
}

gic.cdmm <- function(y, x, lam, type="bic", constr=TRUE) {
	if (constr)
		res <- cdmm(y, x, lam)
	else
		res <- cdmm(y, x, lam, mu=0, maxit=c(50, 1))
	n <- length(y)
	fit <- log(colMeans((y - matrix(res$int, n, length(res$lam), byrow=TRUE) - x %*% res$sol)^2))
	a <- switch(type, bic=log(n), aic=2, ft=log(log(n))*log(max(ncol(x), n)))/n
	if (constr)
		gics <- fit + a*(colSums(res$sol != 0) - 1)
	else
		gics <- fit + a*colSums(res$sol != 0)
	ilam <- which.min(gics)
	list(bet=res$sol[, ilam], lam=res$lam[ilam], int=res$int[ilam])
}

cv.cdmm <- function(y, x, lam, foldid, nfold=10, refit=FALSE, type="min", constr=TRUE) {
	if (constr)
		res <- cdmm(y, x, lam)
	else
		res <- cdmm(y, x, lam, mu=0, maxit=c(50, 1))
	if (missing(foldid)) foldid <- sample(rep(1:nfold, length=length(y)))
	pred <- matrix(, nfold, length(res$lam))
	for (i in 1:nfold) {
		yt <- y[foldid != i]; xt <- x[foldid != i, ]
		yv <- y[foldid == i]; xv <- x[foldid == i, ]
		if (constr) {
			fit <- cdmm(yt, xt, res$lam, maxv=Inf)
			if (refit) for (j in 1:length(res$lam)) {
				supp <- fit$sol[, j] != 0
				if (any(supp)) {
					ans <- cdmm(yt, as.matrix(xt[, supp]), 0, maxv=Inf)
					fit$sol[supp, j] <- ans$sol; fit$int[j] <- ans$int
				}
			}
		} else {
			fit <- cdmm(yt, xt, res$lam, mu=0, maxv=Inf, maxit=c(50, 1))
 			if (refit) for (j in 1:length(res$lam)) {
				supp <- fit$sol[, j] != 0
				if (any(supp)) {
					ans <- cdmm(yt, as.matrix(xt[, supp]), 0, mu=0, maxv=Inf, maxit=c(50, 1))
					fit$sol[supp, j] <- ans$sol; fit$int[j] <- ans$int
				}
			}
		}
		pred[i, ] <- colSums((yv - matrix(fit$int, length(yv), length(res$lam), byrow=TRUE) - xv %*% fit$sol)^2)
	}
	cvm <- colSums(pred)/length(y)
	cvse <- sqrt((colSums(pred^2) - length(y)*cvm^2)/(length(y) - 1))
	imin <- which.min(cvm)
	ilam <- switch(type, min=imin, "1se"=match(TRUE, cvm <= cvm[imin] + cvse[imin]))
	list(bet=res$sol[, ilam], lam=res$lam[ilam], int=res$int[ilam], foldid=foldid)
}

stab.cdmm <- function(y, x, lam, nlam=100, rlam=1/nlam, nsample=100, nsub=floor(0.5*length(y)), constr=TRUE, seed=213) {
	set.seed(seed)
	if (missing(lam)) {
		y <- scale(y, scale=FALSE)
		x <- scale(x, scale=apply(x, 2, sd)*sqrt(nrow(x)-1))
		lam.max <- max(abs(crossprod(x, y)))
		lam <- lam.max*exp(seq(0, log(rlam), length=nlam))
	}
	path <- matrix(0, ncol(x), nlam)
	for (i in 1:nsample) {
		isub <- sample(1:length(y), nsub)
		ysub <- y[isub]; xsub <- x[isub, ]
		if (constr)
			sol <- cdmm(ysub, xsub, lam, maxv=Inf)$sol
		else
			sol <- cdmm(ysub, xsub, lam, mu=0, maxv=Inf, maxit=c(50, 1))$sol
		path <- path + (sol != 0)
	}
	path <- path/nsample
	prob <- apply(path, 1, max)
	list(path=path, prob=prob)
}

cdmm_loop <- function(y, z, nrep = 100, nboot = 100, seed = 23){
  # set.seed(seed1)
  # n <- length(y)
  # pe <- numeric(nrep)
  # for (i in 1:nrep) {
  #   itrn <- sample(n, ntrn)
  #   itst <- setdiff(1:n, itrn)
  #   ans <- cv.cdmm(y[itrn], z[itrn, ], refit=TRUE)
  #   bet <- ans$bet; int <- ans$int
  #   pe[i] <- mean((y[itst] - int - z[itst, ] %*% bet)^2)
  #   r2[i] <- 1-pe[i]/mean((y[itst]-mean(y[itst]))^2)
  # }
  set.seed(seed)
  p <- ncol(z)
  bet.bcv <- matrix(NA, p, nboot)
  for (i in 1:nboot) {
    bootid <- sample(1:length(y), replace=TRUE)
    bet.bcv[, i] <- cv.cdmm(y[bootid], z[bootid, ], refit=TRUE)$bet
  }
  stab.prob <- stab.cdmm(y, z)$prob
  bcv.prob <- rowMeans(bet.bcv != 0)
  isel <- bcv.prob >= 0.7
  return(data.frame(genus=colnames(z)[isel], bcv.prob=bcv.prob[isel], stab.prob=stab.prob[isel], 
                    refitted.coef=rowMeans(matrix(bet.bcv[isel,], nrow = sum(isel)))))
}