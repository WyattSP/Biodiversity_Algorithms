# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 16:12:29 2024

@author: wyattpetryshen

MidPointFM2D algorithm from Suape 1988

Saupe, D. (1988). Algorithms for random fractals.
In: Peitgen, HO., Saupe, D. (eds) The Science of Fractal Images. Springer, New York, NY.
https://doi.org/10.1007/978-1-4612-3784-6_2
"""
#Python implementation
import numpy as np
import random

def f3(delta,x0,x1,x2): return((x0+x1+x2)/3 + delta * random.gauss(0,1))
def f4(delta,x0,x1,x2,x3): return((x0+x1+x2+x3)/4 + delta * random.gauss(0,1))
def MidPointFM2D(maxlevel, sigma, H, seed, addition):
    np.random.seed(seed)
    N = 2**maxlevel
    delta = sigma
    #Create array of size (N+1) by (N+1)
    X = np.zeros((N+1,N+1))

    #Set initial corners
    X[0][0] = delta * random.gauss(0,1)
    X[0][N] = delta * random.gauss(0,1)
    X[N][0] = delta * random.gauss(0,1)
    X[N][N] = delta * random.gauss(0,1)
    D = int(N)
    d = int(N/2)

    for stage in np.arange(1,maxlevel+1):
        delta = delta * 0.5**(0.5*H)

        #Interpolate and offset points
        for x in range(d,N-d+1,D):
            for y in range(d,N-d+1,D):
                X[x][y] = f4(delta,X[x+d][y+d],X[x+d][y-d],X[x-d][y+d],X[x-d][y-d])
        #Displace other points also if nedded
        if addition == True:
            for x in range(0,N+1,D):
                for y in range(0,N+1,D):
                    X[x][y] = X[x][y] + delta * random.gauss(0,1)
        #Going from grid type II to type I
        delta = delta * 0.5**(0.5*H)
        #Interpolate and offset boundary grid points
        for x in range(d,N-d+1,D):
            X[x][0] = f3(delta,X[x+d][0],X[x-d][0],X[x][d])
            X[x][N] = f3(delta,X[x+d][N],X[x-d][N],X[x][N-d])
            X[0][x] = f3(delta,X[0][x+d],X[0][x-d],X[d][x])
            X[N][x] = f3(delta,X[N][x+d],X[N][x-d],X[N-d][x])
        #interpolate and offset interior gird points
        for x in range(d,N-d+1,D):
            for y in range(D,N-d+1,D):
                X[x][y] = f4(delta,X[x][y+d],X[x][y-d],X[x+d][y],X[x-d][y])
        for x in range(D,N-d+1,D):
            for y in range(d,N-d+1,D):
                X[x][y]= f4(delta,X[x][y+d],X[x][y-d],X[x+d][y],X[x-d][y])
        if addition == True:
            for x in np.arange(0,N+1,D):
                for y in np.arange(0,N+1,D):
                    X[x][y] = X[x][y] + delta * random.gauss(0,1)
            for x in np.arange(d,N-d+1,D):
                for y in np.arange(d,N-d+1,D):
                    X[x][y] = X[x][y] + delta * random.gauss(0,1)
        D = int(D/2)
        d = int(d/2)
    return(X)
