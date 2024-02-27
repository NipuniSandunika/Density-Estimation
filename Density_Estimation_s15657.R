### Density Estimation ###
      ## S15657 ##

## Frequency polygon

library(MASS)
waiting = geyser$waiting #in MASS
n = length(waiting)
# freq poly bin width using normal ref rule
h = 2.15 * sqrt(var(waiting)) * n^(-1/5)
# calculate the sequence of breaks and histogram
br = pretty(waiting, diff(range(waiting)) / h)
brplus = c(min(br)-h, max(br+h))
histg = hist(waiting, breaks = br, freq = FALSE,
              main = "Frequency polygon",col = "skyblue", xlim = brplus)
vx = histg$mids
vy = histg$density
#density est at vertices of polygon
delta = diff(vx)[1] # h after pretty is applied
k = length(vx)
vx = vx + delta
# the bins on the ends
vx = c(vx[1]- 2 * delta, vx[1]- delta, vx)
vy = c(0, vy, 0)
# add the polygon to the histogram
polygon(vx, vy)


## ASH density estimate

library(MASS) 
waiting = geyser$waiting 
n = length(waiting) 
m = 20 
a = min(waiting) - .5 
b = max(waiting) + .5 
h = 5.5 
delta = h / m 

#get the bin counts on the delta-width mesh. 
br = seq(a - delta*m, b + 2*delta*m, delta) 
histg = hist(waiting, breaks = br, plot = FALSE) 
nk = histg$counts 
K = abs((1-m):(m-1)) 

fhat = function(x) { 
  # locate the leftmost interval containing x 
  i = max(which(x > br)) 
  k = (i - m + 1):(i + m - 1) 
  # get the 2m-1 bin counts centered at x 
  vk = nk[k] 
  sum((1 - K / m) * vk) / (n * h)   #f.hat 
} 

# density can be computed at any points in range of data 
z = as.matrix(seq(a, b + h, .1)) 
f.ash = apply(z, 1, fhat)   #density estimates at midpts 

# plot ASH density estimate over histogram 
br2 = seq(a, b + h, h) 
hist(waiting, breaks = br2, freq = FALSE, main = "ASH density estimate", col = "green",
     ylim = c(0, max(f.ash))) 
lines(z, f.ash, xlab = "waiting") 

## Kernel density estimate of geyser waiting time

data(geyser)
names(geyser)

hist(geyser$waiting,freq = F , col = "purple", main = "Kernel density estimate of geyser waiting time") #histogram
lines(density(geyser$waiting), lwd = 2, col = 'black')

#density estimators of kernel functions 
par(mfrow=c(3,2))
plot(density(geyser$waiting,kernel = "rectangular"))
plot(density(geyser$waiting,kernel = "gaussian"))
plot(density(geyser$waiting,kernel = "triangular"))
plot(density(geyser$waiting,kernel = "epanechnikov"))
plot(density(geyser$waiting,kernel = "biweight"))
plot(density(geyser$waiting,kernel = "cosine"))
par(mfrow=c(1,1))

# Some density estimates are produced by different bandwidth
library(MASS)
truehist(geyser$waiting ,nbin=25,col="lightgrey")
lines(density(geyser$waiting ),lwd = 2)
lines(density(geyser$waiting ,bw="SJ"), col="red",lwd = 2) # bandwidth selected by Sheather-Jones method
lines(density(geyser$waiting ,bw="SJ-dpi"), col="blue",lwd = 2) # Sheather-Jones bandwidth with a data-dependent adjustment
legend("topright", c("bw=default","h=SJ","h=SJ-dpi"), box.lty = 0,
       lty = 1, col = c("black","red","blue"), lwd = c(2, 2, 2))

#Kernel density estimates using a Gaussian kernel with different band widths
#applying different h

library(MASS)
plot(density(geyser$waiting,bw =2, kernel = "gaussian"),col ="purple", lwd=2,main = "Kernel density estimates using a Gaussian kernel with different band 
widths")
lines(density(geyser$waiting,bw =5, kernel = "gaussian"),col ="red", lwd=2)
lines(density(geyser$waiting,bw =10, kernel = "gaussian"),col ="green", lwd=2)
legend("topright", c("h=2","h=5","h=10"), box.lty = 0,
       lty = 1, col = c("purple","red","green"), lwd = c(2, 2, 2))


## Exponential density
x = rexp(1000, 1) 
plot(density(x), xlim = c(-1, 6), ylim = c(0, 1), main="") 
abline(v = 0) 
# add the true density to compare 
y = seq(.001, 6, .01) 
lines(y, dexp(y, 1), lty = 2) 

#Reflection boundary technique
xx = c(x, -x) 
g = density(xx, bw = bw.nrd0(x)) 
a = seq(0, 6, .01) 

ghat = approx(g$x, g$y, xout = a) 
fhat = 2 * ghat$y       # density estimate along a 

bw <- paste("Bandwidth = ", round(g$bw, 5)) 
plot(a, fhat, type="l", xlim=c(-1, 6), ylim=c(0, 1), 
     main = "", xlab = bw, ylab = "Density") 
abline(v = 0) 
# add the true density to compare 
y <- seq(.001, 6, .01) 
lines(y, dexp(y, 1), lty = 2) 

## Bivariate density polygon

# Generate standard bivariate normal random sample
n <- 1000
d <- 2
x <- matrix(rnorm(n * d), n, d)

# Clear any existing 3D plots
clear3d()

# 3D scatter plot using rgl
library(rgl)

plot3d(x[, 1], x[, 2], col="lightblue", size=2,
       xlab="X", ylab="Y", zlab="Frequency",
       main="Bivariate Normal Distribution")

# Bivariate frequency table: bind2d
bin2d <- 
  function(x, breaks1 = "Sturges", breaks2 = "Sturges"){ 
    # Data matrix x is n by 2 
    # breaks1, breaks2: any valid breaks for hist function 
    # using same defaults as hist 
    histg1 <- hist(x[,1], breaks = breaks1, plot = FALSE) 
    histg2 <- hist(x[,2], breaks = breaks2, plot = FALSE) 
    brx <- histg1$breaks 
    bry <- histg2$breaks 
    
    # bin frequencies 
    freq <- table(cut(x[,1], brx),  cut(x[,2], bry)) 
    
    return(list(call = match.call(), freq = freq, 
                breaks1 = brx, breaks2 = bry, 
                mids1 = histg1$mids, mids2 = histg2$mids)) 
  }
bin2d(iris[1:50,1:2]) 
#generate standard bivariate normal random sample 
n <- 2000;   d <- 2 
x <- matrix(rnorm(n*d), n, d) 

# compute the frequency table and density estimates 
# using bin2d function from the previous example 
b <- bin2d(x) 
h1 <- diff(b$breaks1) 
h2 <- diff(b$breaks2) 

# matrix h contains the areas of the bins in b 
h <- outer(h1, h2, "*") 

Z <- b$freq / (n * h)  # the density estimate 
persp(x=b$mids1, y=b$mids2, z=Z, shade=TRUE, 
      xlab="X", ylab="Y", main="", 
      theta=45, phi=30, ltheta=60) 


## Bivariate ASH density estimate

library(ash)  # for bivariate ASH density est. 
# generate N_2(0,Sigma) data 
n <- 2000 
d <- 2 
nbin <- c(30, 30) # number of bins 
m <- c(5, 5)      # smoothing parameters 

# First example with positive correlation 
Sigma <- matrix(c(1, .5, .5, 1), 2, 2) 
set.seed(345)

library(MASS) #for mvrnorm() 
x <- mvrnorm(n, c(0, 0), Sigma=Sigma) 
b <- bin2(x, nbin = nbin) 
# kopt is the kernel type, here triangular 
est <- ash2(b, m = m, kopt = c(1,0)) 
par(mfrow= c(2,2)) 

persp(x = est$x, y = est$y, z = est$z, shade=TRUE, 
      xlab = "X", ylab = "Y", zlab = "", main="", 
      theta = 30, phi = 75, ltheta = 30, box = FALSE) 
contour(x = est$x, y = est$y, z = est$z, main="") 

# Second example with negative correlation 
Sigma <- matrix(c(1, -.5, -.5, 1), 2, 2) 
set.seed(345) 
#x <- rmvn.eigen(n, c(0, 0), Sigma=Sigma) 
x <- mvrnorm(n, c(0, 0), Sigma=Sigma) 
b <- bin2(x, nbin = nbin) 
est <- ash2(b, m = m, kopt = c(1,0)) 

persp(x = est$x, y = est$y, z = est$z, shade=TRUE, 
      xlab = "X", ylab = "Y", zlab = "", main="", 
      theta = 30, phi = 75, ltheta = 30, box = FALSE) 
contour(x = est$x, y = est$y, z = est$z, main="") 

par(mfrow= c(1,1))





## Product kernel estimate of a bivariate normal mixture

library(MASS)  #for mvrnorm and kde2d 
#generate the normal mixture data 
n <- 2000 
p <- c(.2, .3, .5) 
mu <- matrix(c(0, 1, 4, 0, 3, -1), 3, 2) 
Sigma <- diag(2) 
i <- sample(1:3, replace = TRUE, prob = p, size = n) 
k <- table(i) 

x1 <- mvrnorm(k[1], mu = mu[1,], Sigma) 
x2 <- mvrnorm(k[2], mu = mu[2,], Sigma) 
x3 <- mvrnorm(k[3], mu = mu[3,], Sigma) 
X <-  rbind(x1, x2, x3)   #the mixture data 
x <- X[,1] 
y <- X[,2] 

print(c(bandwidth.nrd(x), bandwidth.nrd(y))) 


# accepting the default normal reference bandwidth 
fhat <- kde2d(x, y) 
contour(fhat) 

persp(fhat, phi = 30, theta = 20, d = 5, xlab = "x") 
# select bandwidth by unbiased cross-validation 
h = c(ucv(x), ucv(y)) 

h 

fhat <- kde2d(x, y, h = h) 
contour(fhat) 
persp(fhat, phi = 30, theta = 20, d = 5, xlab = "x") 



