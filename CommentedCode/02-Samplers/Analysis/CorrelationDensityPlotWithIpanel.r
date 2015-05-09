# Plot for correlation density, good for posterior correlation of MCMC chains
# by Florian Hartig http://florianhartig.wordpress.com/
# used, e.g., in http://www.biogeosciences.net/11/1261/2014/bg-11-1261-2014.html

library(IDPmisc)
 
panel.hist <- function(x, ...)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(usr[1:2], 0, 1.5) )
h <- hist(x, plot = FALSE)
breaks <- h$breaks; nB <- length(breaks)
y <- h$counts; y <- y/max(y)
rect(breaks[-nB], 0, breaks[-1], y, col="blue4", ...)
}
 
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
usr <- par("usr"); on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
r <- abs(cor(x, y, method = "spearman"))
txt <- format(c(r, 0.123456789), digits=digits)[1]
txt <- paste(prefix, txt, sep="")
if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
text(0.5, 0.5, txt, cex = cex * r)
}
 
betterPairs <- function(YourData){
return(pairs(YourData, lower.panel=function(...) {par(new=TRUE);ipanel.smooth(...)}, diag.panel=panel.hist, upper.panel=panel.cor))
}
 
# Example
 
#x = rnorm(10000)
#betterPairs(data.frame(A = x, B = 0.6 * x + 0.3 * rnorm(10000), C = rnorm(10000)))

##################################################################

# Scatter plot with histograms 



scatterhist = function(x, y, xlab="", ylab="", smooth = T){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  if (smooth == T) smoothScatter(x,y, colramp = colorRampPalette(c("white", "darkorange")))
  else plot(x,y)
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
        at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
        at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}

# Example

# 
# 
# banana=function(A,B,C1,C2,N,keep=10,init=10)
# {
#   R=init*keep+N*keep
#   x1=x2=0
#   bimat=matrix(double(2*N),ncol=2)
#   for (r in 1:R) {
#     x1=rnorm(1,mean=(B*x2+C1)/(A*(x2^2)+1),sd=sqrt(1/(A*(x2^2)+1)))
#     x2=rnorm(1,mean=(B*x2+C2)/(A*(x1^2)+1),sd=sqrt(1/(A*(x1^2)+1)))
#     if (r>init*keep && r%%keep==0) {
#       mkeep=r/keep; bimat[mkeep-init,]=c(x1,x2)
#     }
#   }
#   
#   return(bimat)
# }
# A=0.5; B=0; C1=C2=3
# sample=banana(A=A,B=B,C1=C1,C2=C2,1000)
# scatterhist(sample[,1], sample[,2])
