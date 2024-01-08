"""
Created on Sun Jan  7 16:12:29 2024

@author: wyattpetryshen

MidPointFM2D algorithm from Suape 1988

Saupe, D. (1988). Algorithms for random fractals.
In: Peitgen, HO., Saupe, D. (eds) The Science of Fractal Images. Springer, New York, NY.
https://doi.org/10.1007/978-1-4612-3784-6_2
"""
#R implementation
f3 <- function(delta,x0,x1,x2){return((x0+x1+x2)/3 + delta *rnorm(1,0,1))}
f4 <- function(delta,x0,x1,x2,x3){return((x0+x1+x2+x3)/4 + delta *rnorm(1,0,1))}
MidPointFM2D <- function(maxlevel, sigma, H, seed, addition){
  '''
  X_arr: doubly indexed real array of size (N + 1)^2
  maxlevel: max number of recursions, N = 2^maxlevel
  sigma: initial standard deviation
  H: parameter H determines fractal dimension D = 3 - H
  seed: random seed
  '''
  if(2^maxlevel > (2^31-1)){stop("Array N greater than memory 2^N greater than 2^31-1")}
  set.seed(seed)
  N = 2^maxlevel
  delta = sigma
  X <- array(dim = c((N+1),(N+1)))
  X[,] <- 0
  #set initial random corners
  X[1,1] <- rnorm(1,0,delta)
  X[1,N+1] <- rnorm(1,0,delta)
  X[N+1,1] <- rnorm(1,0,delta)
  X[N+1,N+1] <- rnorm(1,0,delta)
  D <- N
  d <- N/2
  for(stage in seq(2,maxlevel+1,1)){
    #going from grid type I to type II
    delta = delta * 0.5^(0.5*H)
    #interpolate and offset points
    for(x in seq(d+1, to = N-d+1, by = D)){
      for(y in seq(d+1, to = N-d+1, by = D)){
        X[x,y] <- f4(delta, X[x+d,y+d], X[x+d,y-d],X[x-d,y+d],X[x-d,y-d])
      }
    }
    #Displace other points also if needed
    if(addition == TRUE){
      for(x in seq(0+1,N+1,D)){
        for(y in seq(0+1,N+1,D)){
          X[x,y] = X[x,y] + rnorm(1,0,delta)
        }
      }
    }
    #going from grid type II to type I
    delta = delta * 0.5^(0.5*H)
    #interpolate and offset boundary grid points
    for(x in seq(d+1,N-d+1,D)){
      X[x,1] <- f3(delta, X[x+d,1],X[x-d,1],X[x,d])
      X[x,N+1] <- f3(delta, X[x+d,N+1],X[x-d,N+1],X[x,N-d])
      X[1,x] <- f3(delta, X[1,x+d],X[1,x-d],X[d,x])
      X[N+1,x] <- f3(delta, X[N+1,x+d],X[N+1,x-d],X[N+1-d,x])
    }
    #interpolate and offset interior grid points
    for(x in seq(d+1,N-d+1,D)){
      tryCatch(
      for(y in seq(D+1,N-d+1,D)){
        X[x,y] = f4(delta,X[x,y+d], X[x,y-d],X[x+d,y],X[x-d,y])
      }
      , error = function(e){})
      }
    tryCatch(
    for(x in seq(D+1,N-d+1,D)){
      for(y in seq(d+1,N-d+1,D)){
        X[x,y] = f4(delta, X[x,y+d], X[x,y-d],X[x+d,y],X[x-d,y])
        }
      }, error = function(e){})
    #displace other points also if needed
    if(addition == TRUE){
      for(x in seq(0+1,N+1,D)){
        for(y in seq(0+1,N+1,D)){
          X[x,y] = X[x,y] + rnorm(1,0,delta)
        }
      }
      for(x in seq(d+1,N-d+1,D)){
        for(y in seq(d+1,N-d+1,D)){
          X[x,y] = X[x,y] + rnorm(1,0,delta)
        }
      }
    }
    D <- D/2
    d <- d/2
  }
  return(X)
}
