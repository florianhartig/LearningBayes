maxLikCoda <- function (x) 
{
    xx <- as.matrix(x)
    pars.maxp <- rep(NA, ncol(xx))
    for (i in 1:nvar(x)) {
        y <- xx[, i, drop = TRUE]
        bwf <- function(x) {
          x <- x[!is.na(as.vector(x))]
          return(1.06 * min(sd(x), IQR(x)/1.34) * length(x)^-0.2)
        }
        bw <- bwf(y)
        width <- 4 * bw
            scale <- "open"
            if (max(y) <= 1 && 1 - max(y) < 2 * bw) {
                if (min(y) >= 0 && min(y) < 2 * bw) {
                  scale <- "proportion"
                  y <- c(y, -y, 2 - y)
                }
            }
            else if (min(y) >= 0 && min(y) < 2 * bw) {
                scale <- "positive"
                y <- c(y, -y)
            }
            else scale <- "open"
            dens <- density(y, width = width)
            if (scale == "proportion") {
                dens$y <- 3 * dens$y[dens$x >= 0 & dens$x <= 
                  1]
                dens$x <- dens$x[dens$x >= 0 & dens$x <= 1]
            }
            else if (scale == "positive") {
                dens$y <- 2 * dens$y[dens$x >= 0]
                dens$x <- dens$x[dens$x >= 0]
            }
        ii <- which.max(dens$y)
        pars.maxp[i] <- dens$x[ii]
        
    }##for
    return(pars.maxp)
}##maxLikCoda
