
# From a talk by Duncan Murdoch given in June 11, 2010, Teaching statistical computing using 3D graphics in R
# Vienna university of economics and business
# http://statmath.wu.ac.at/research/talks/resources/tscu3dg.R
# To obtain rgl on linux, it is best to do sudo apt-get install r-cran-rgl from the linux prompt. This ensures that all the dependencies with OpenGL are met.
#

# Downloaded 9.5.15 by FH from http://www.ibm.com/developerworks/library/ba-optimR-john-nash/index.html
# See copyright notice on this website

library(rgl)
r3dDefaults$windowRect <- c(0,50, 700, 700)

library(tkrgl)

# Demo 4a

showsimplex <- function(x, f, col="blue") {
	n <- nrow(x)
	z <- numeric(n)
	for (i in 1:n) z[i] <- f(x[i,])
	xyz <- cbind(x, z)
	
	# This is tricky:  
	
	# 1. draw all lines, taking vertices two at a time:
	segments3d(xyz[as.numeric(combn(n, 2)),])
	# 2. draw all faces, taking vertices three at a time:
	triangles3d(xyz[as.numeric(combn(n, 3)),], col=col, alpha=0.3)
}

neldermead <- function(x, f) {
	n <- nrow(x)
	p <- ncol(x)
	
	if (n != p + 1) stop(paste('Need', p + 1, 'starting points'))
	
	fx <- rep(NA, n)
	for (i in 1:n) fx[i] <- f(x[i,])
	
	o <- order(fx)
	fx <- fx[o]
	x <- x[o,]
	xmid <- apply(x[1:p,], 2, mean)
	z1 <- xmid - (x[n,] - xmid)
	fz1 <- f(z1)
	
	if (fz1 < fx[1]) {
		z2 <- xmid - 2*(x[n,] - xmid)
		fz2 <- f(z2)
		if (fz2 < fz1) {
			cat('Accepted reflection and expansion, f(z2)=',fz2,'\n')
			x[n,] <- z2
		} else {
			cat('Accepted good reflection, f(z1)=',fz1,'\n')
			x[n,] <- z1
		}
	} else if (fz1 < fx[p]) {
		cat('Accepted okay reflection, f(z1)=',fz1,'\n')
		x[n,] <- z1
	} else {
		if (fz1 < fx[n]) {
			x[n,] <- z1
			fx[n] <- fz1
		}
		z3 <- xmid + (x[n,] - xmid)/2
		fz3 <- f(z3)
		if (fz3 < fx[n]) {
			cat('Accepted contraction 1, f(z3)=',fz3,'\n')
			x[n,] <- z3
		} else {
			cat('Accepted contraction 2,')
			for (i in 2:n) {
				x[i,] <- x[1,] + (x[i,] - x[1,])/2
				cat(' f(z', i+2, ') = ', f(x[i,]), sep='')
			}
			cat('\n')
		}
	}
	return(x)
}


library(misc3d)

# Example taken from ?contour3d
#Example 2: Nested contours of mixture of three tri-variate normal densities
nmix3 <- function(x, y, z, m, s) {
	0.4 * dnorm(x, m, s) * dnorm(y, m, s) * dnorm(z, m, s) +
			0.3 * dnorm(x, -m, s) * dnorm(y, -m, s) * dnorm(z, -m, s) +
			0.3 * dnorm(x, m, s) * dnorm(y, -1.5 * m, s) * dnorm(z, m, s)
}
f <- function(x,y,z) nmix3(x,y,z,.5,.5)
g <- function(n = 40, k = 5, alo = 0.1, ahi = 0.5, cmap = heat.colors) {
	th <- seq(0.05, 0.2, len = k)
	col <- rev(cmap(length(th)))
	al <- seq(alo, ahi, len = length(th))
	x <- seq(-2, 2, len=n)
	contour3d(f,th,x,x,x,color=col,alpha=al)
	bg3d(col="white")
}

f3 <- function(x) -f(x[1], x[2], x[3])

set.seed(3)

open3d(zoom=0.74,windowRect = c(0,0, 800, 800))
g(40,5)                
xyz <- matrix(rnorm(12, sd=0.1) + rep(rnorm(3,sd=2), each=4), 4, 3)
showsimplex(xyz, f3)

# Demo 4b

for (i in 1:30) {
	xyz <- neldermead(xyz,f3); showsimplex(xyz, f3, "blue")
	Sys.sleep(1)
}

