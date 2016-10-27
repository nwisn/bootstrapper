#======================================================================================================
# twogropu: bootstrap methods for comparing two groups
# author: Nicholas Wisniewski
# date: October 27, 2016
#======================================================================================================

# faster mean and standard deviation
mean.fast <- function(x) sum(x) / length(x)
sd.fast <- function(x){
    mu <- mean.fast(x)
    sqrt(sum((x - mu)^2) / (length(x)-1))
}


# function that computes test statistics and effect sizes
tstat <- function(Xbar, Ybar, Xsd, Ysd, Xn, Yn, type="unequal.var"){
    d <- Xbar - Ybar
    if(type == "equal.var"){ # Student's statistic
        s <- sqrt( ((Xn-1)*Xsd^2 + (Yn-1)*Ysd^2) / (Xn + Yn - 2) ) * sqrt((1/Xn) + (1/Yn))
    } else if(type == "unequal.var"){  # Welch's statistic
        s <- sqrt( ((Xsd^2)/Xn) + ((Ysd^2)/Yn) )  
    } else if(type == "Cohen") {  # Cohen's d
        s <- sqrt( ((Xn-1)*Xsd^2 + (Yn-1)*Ysd^2) / (Xn + Yn - 2) ) 
    } else if(type == "delta") {  # simple difference in means
        s <- 1
    }
    t <- d/s
    return(t)
}

# function to bootstrap pvalues. can do any bootstrap test
boot_stats <- function(xy, X, Y, middle=mean.fast, spread=sd.fast, type="unequal.var", rank=F){
    nx <- length(X)
    ny <- length(Y)
    x.mean <- apply(xy[, (1:nx)], 1, middle)
    y.mean <- apply(xy[,-(1:nx)], 1, middle)
    x.sd <- apply(xy[, (1:nx)], 1, spread)
    y.sd <- apply(xy[,-(1:nx)], 1, spread)
    t <- tstat(x.mean, y.mean, x.sd, y.sd, nx, ny, type=type)
    if(rank==T){xyr <- rank(c(X,Y)); X <- xyr[1:nx];  Y <- xyr[-(1:ny)] }
    t0 <- tstat(middle(X), middle(Y), spread(X), spread(Y), nx, ny, type=type)
    p <- sum(abs(t) >= abs(t0))/dim(xy)[1]
    return( list(p = p, t0 = t0, t = t, X.mu = x.mean, Y.mu = y.mean, X.s = x.sd, Y.s = y.sd))
}

# boot.ci for BCa confidence intervals
# library(boot)
# diff.means <- function(longform, w){
#     ix <- which(longform[,"group"]=="X")
#     m1 <- sum(longform[ix,"value"] * w[ix])
#     m2 <- sum(longform[-ix, "value"] * w[-ix])
#     m1-m2
# }


#' Two Group Comparison
#' 
#' Computes statistical significance using parametric, non-parametric, and bootstrap methods
#' @param X The first group
#' @param Y The second group
#' @param nboot The number of resamples
#' @param replace Whether to sample with replacement
#' @param location Function that measures central location. Defaults to mean
#' @param scale Function that measures spread. Defaults to standard deviation
#' @param dither Whether to add dithering noise. Defaults to no
#' @param amp Amplitude of dithering noise. 
#' @param use Statistics to compute. Can use multiple. Defaults to all.
#' @return A `twogroup` object
#' @export
twogroup <- function(X, Y, nboot=10000, replace=T, location=mean.fast, scale=sd.fast, dither=F, amp=.7,
                     use=c("equal.var","unequal.var","rank","Student","Welch","Wilcox","ES")){
    if(dither==T){
        noise <- matrix(rnorm(nboot*length(c(X,Y)), 0, amp*min(sd(X),sd(Y))/length(c(X,Y))^(1/2)), nboot, length(c(X,Y)))
    } else noise <- 0
    
    if("equal.var" %in% use){
        # one-box
        xy <- matrix(sample(c(X,Y), size = nboot * length(c(X,Y)), replace = replace), nboot, length(c(X,Y))) 
        onebox <- boot_stats(xy + noise, X, Y, middle=location, spread=scale, type="equal.var")
    }
    if("unequal.var" %in% use){
        # two-box
        x <- matrix(sample(X-location(X), size = nboot * length(X), replace = replace), nboot, length(X))  # x residuals
        y <- matrix(sample(Y-location(Y), size = nboot * length(Y), replace = replace), nboot, length(Y))  # y residuals
        twobox <- boot_stats(cbind(x,y) + noise, X, Y, middle=location, spread=scale, type="unequal.var")
    }
    if("rank" %in% use){
        # rank test
        xy <- matrix(sample(c(X,Y), size = nboot * length(c(X,Y)), replace = replace), nboot, length(c(X,Y)))  
        xy.rank <- t(apply(xy, 1, rank)) # rank
        onebox.rank <- boot_stats(xy.rank, X, Y, middle=location, spread=scale, type="delta", rank=T)  
    }
    if("Student" %in% use){
        student <- t.test(X,Y, var.equal=T)
    }
    if("Welch" %in% use){
        welch <- t.test(X,Y, var.equal=F)
    }
    if("Wilcox" %in% use){
        wilcox <- wilcox.test(X,Y, conf.int=T)
    }
    if("ES" %in% use){
#         longform <- data.frame(group=c(rep("X",length(X)), rep("Y",length(Y))), value=c(X,Y))
#         result <- boot(longform, diff.means, R = nboot, stype = "w", strata = longform$group)
#         bca <- boot.ci(result, type="all")
        
        # effect size with confidence interval
        x <- matrix(sample(X, size = nboot * length(X), replace = replace), nboot, length(X))  # x
        y <- matrix(sample(Y, size = nboot * length(Y), replace = replace), nboot, length(Y))  # y
        effectsize <- boot_stats(cbind(x,y) + noise, X, Y, middle=location, spread=scale, type="delta")
        effectsize.d <- boot_stats(cbind(x,y) + noise, X, Y, middle=location, spread=scale, type="Cohen")
        
    }
    
    which.return <- match(use, c("equal.var", "unequal.var", "rank", "Student", "Welch", "Wilcox", "ES"))
    output <- list(equal.var=if(exists("onebox")) onebox else NULL, 
                   unequal.var=if(exists("twobox")) twobox else NULL, 
                   rank=if(exists("onebox.rank")) onebox.rank else NULL,
                   Student=if(exists("student")) student else NULL,
                   Welch=if(exists("welch")) welch else NULL,
                   Wilcox=if(exists("wilcox")) wilcox else NULL,
                   ES=if(exists("effectsize")) list(delta=effectsize, Cohen=effectsize.d) else NULL)
    #selected <- output[which.return]
    selected <- output
    class(selected) <- "twogroup"
    return(selected)
}


#' Tables of Statistics
#' 
#' Creates a table of p-values across methods, as well as confidence intervals, and degrees of freedom
#' @param tg A `twogroup` object
#' @return A list of 3 tables
#' @export
stattable <- function(tg){
    stopifnot(inherits(tg, "twogroup"))
    p1 <- tg$equal.var$p 
    p2 <- tg$unequal.var$p 
    p3 <- tg$rank$p 
    p4 <- tg$Student$p.value 
    p5 <- tg$Welch$p.value 
    p6 <- tg$Wilcox$p.value 
    MD <-  tg$ES$delta$t0 
    MD.CI <-  as.vector(quantile(tg$ES$delta$t,c(.025,.975))) 
    CD <-  tg$ES$Cohen$t0 
    CD.CI <-  as.vector(quantile(tg$ES$Cohen$t, c(.025,.975))) 
    pmat <- rbind(p1,p2,p3,p4,p5,p6)
    df.student <- tg$Student$parameter
    df.welch <- tg$Welch$parameter
    dfmat <- rbind(df.student, df.welch)
    
    if(!is.null(pmat)){
        whichuse <- match(rownames(pmat), c("p1","p2","p3","p4","p5","p6"))
        rownames(pmat) <- c("equal.var","unequal.var","rank","Student","Welch","Wilcox")[whichuse]
        colnames(pmat) <- "p.value"
    }
    
    if(!is.null(dfmat)){
        whichuse <- match(rownames(dfmat), c("df.student", "df.welch"))
        rownames(dfmat) <- c("Student","Welch")[whichuse]
        colnames(dfmat) <- "df"
    }
    
    cimat <- rbind(c(MD, MD.CI[1], MD.CI[2]), c(CD, CD.CI[1], CD.CI[2]))
    rownames(cimat) <- c("Location Difference", "Cohen's d")
    colnames(cimat) <- c("d","95% LCI", "95% UCI") [1:ncol(cimat)]
    return(list(p=pmat, ci=cimat, df=dfmat))
}


#' Summary
#' 
#' Prints a summary of statistical significance and confidence intervals
#' @param tg A `twogroup` object
#' @return none
#' @export
summary.twogroup <- function(tg){
    stopifnot(inherits(tg, "twogroup"))
    tables <- stattable(tg)
    if(!is.null(tables$p)) print(tables$p)
    if(!is.null(tg$ES$delta)) print(tables$ci)
    if(!is.null(tables$df)) print(tables$df)
}

#' Plot
#' 
#' Plots the resampling distribution for different hypotheses
#' @param tg A `twogroup` object
#' @return none
#' @export
plot.twogroup <- function(tg, use=NULL, alpha=.05){
    stopifnot(inherits(tg, "twogroup"))
    uses <- c("equal.var","unequal.var","rank","Student","Welch","Wilcox", "ES")
    tables <- stattable(tg)
    whichp <- match(rownames(tables$p), uses)
    if(is.null(use)){
        if(min(whichp) %in% c(1,2,3,7)) use <- uses[min(whichp)] else stop("No plot methods found for the specified tests")
    }
    if(use %in% c("Student","Welch","Wilcox")) stop("No plot methods found for the specified tests")
    if(use != "ES") {t<-tg[[use]]$t; t0<-tg[[use]]$t0} else {t<-tg[[use]]$delta$t; t0<-tg[[use]]$delta$t0}
        
    par(mfrow=c(1,2)) 
    
    hist(t, breaks="FD", col="grey", border = "white", xlab="test statistic", main=paste(use, "resamples", sep=" "), probability = T, xlim=c(mean(t)-5*sd(t), mean(t)+5*sd(t)))
    legend("topright", c("Accept", "Reject"), col=c("green", "red"), lwd=5, cex=.5)
    d <- density(t)
    #lines(d, lwd=2)
    crits <- as.vector(quantile(t, c(alpha/2,1-alpha/2)))
    #abline(v=c(crits, t0), col=c("black","black","red"), lwd=2)
    abline(v=c(t0), col=c("black"), lwd=2, lty="dashed")
    lower.ix <- which(d$x<=crits[1])
    upper.ix <- which(d$x>=crits[2])
    polygon(c(d$x[lower.ix],tail(d$x[lower.ix],n=1)), c(d$y[lower.ix],0), col = "red", border = "black", density=50)
    polygon(c(head(d$x[upper.ix],n=1), d$x[upper.ix]), c(0, d$y[upper.ix]), col = "red", border = "black", density=50)
    polygon(d$x[c(tail(lower.ix,1),tail(lower.ix,1):head(upper.ix,1),head(upper.ix,1)) ], c(0,d$y[tail(lower.ix,1):head(upper.ix,1)],0), col = "green", border = "black", density=50)
    
#     m<-mean(t)
#     std<-sqrt(var(t))
#     curve(dnorm(x, mean=m, sd=std), 
#           col="blue", lwd=1, add=TRUE, yaxt="n")
    
#    qqnorm(t, pch=16, cex=.3, main="normal q-q plot")
    if(is.null(tables$df)) df <- Inf else{
        if(use=="equal.var") df <- tables$df[which(rownames(tables$df)=="Student"),]
        if(use=="unequal.var") df <- tables$df[which(rownames(tables$df)=="Welch"),]
        if(use=="rank") df <- Inf
        if(use=="ES") df <- Inf
    }

    require(limma)
    qqt(t, df = df, ylim = range(t), main = "Student's t Q-Q plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", plot.it = TRUE, pch=16, cex=.5) 
    qqline(t, distribution=function(x) qt(x,df=df))
    abline(h=c(t0), lty=c("dashed"), col=c("black"), lwd=2)
    polygon(c(-100,100,100,-100,-100),c(quantile(t,c(alpha/2,1-alpha/2))[1],quantile(t,c(alpha/2,1-alpha/2))[1],quantile(t,c(alpha/2,1-alpha/2))[2],quantile(t,c(alpha/2,1-alpha/2))[2],quantile(t,c(alpha/2,1-alpha/2))[1]) , col = "green", border = "green", density=10)
    polygon(c(-100,100,100,-100,-100),c(quantile(t,c(alpha/2,1-alpha/2))[1],quantile(t,c(alpha/2,1-alpha/2))[1],-100,-100,quantile(t,c(alpha/2,1-alpha/2))[1]) , col = "red", border = "red", density=10)
    polygon(c(-100,100,100,-100,-100),c(quantile(t,c(alpha/2,1-alpha/2))[2],quantile(t,c(alpha/2,1-alpha/2))[2],100,100,quantile(t,c(alpha/2,1-alpha/2))[2]) , col = "red", border = "red", density=10)
    par(mfrow=c(1,1)) 
}

