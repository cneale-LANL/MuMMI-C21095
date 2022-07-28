#=================================================================================
# determine subsampling points 
# Use simulated anealing to get close to uniform distribution in weights
# Only allow selection from a subset
#.................................................................................

#setwd("/home/tyler/LANL/pilot2/Nick_R_code")
setwd("~/LANL/Current/NCI-P2/UQ_Feb21")
rm( list=ls() )
gc()

library(parallel)
MaxCores <- detectCores() - 1

library(Rfast)
library(dplyr)
library(gam)
library(readr)
library(RANN)
library(reticulate)
use_python("/usr/bin/python3")
#use_python("/Users/treddy/python_venvs/python_376_venv/bin/python")
source_python("lib.py")

#source("functions.r")

#========================================================================
#  functions
#  kernel for W_2^\infty

kp <- function(x,y){
  v <- (x-y)
  kk <- rep(0,length(v))
  idx.nzero <- (abs(v) > 1e-8 )
  kk[idx.nzero] <- 2*(sin(v[idx.nzero])-v[idx.nzero]*cos(v[idx.nzero]))/(pi*v[idx.nzero]^3)
  kk[!idx.nzero] <- 4/(6*pi)
  prod(kk)
}

#  table function with predermined label list, allows zero counts
Table <- function(x,label){
  z <- rep(0, length(label))
  names(z) <- label
  y <- table(x)
  z[names(y)] <- y
  return(z)
}


#=======================================================================
#  read in the data

# this is a very fast reader, which places the data
# in a tibble, which we then convert to a full data
# frame
patch <- as.data.frame(read_delim("list.latentspace.patches.txt", delim = " ", col_names = FALSE))

#.......................................................................
# id of allowed points  --- should be read from a file.
# for now, is a random subset of size 30000
PTS.label   <- sample(patch[,1], 30000,replace = FALSE)
PTS         <- patch[ patch[,1] %in% PTS.label, ]


#======================================================================
#  define data matrices, center and scale

sscale <- 5

#......................................................................
# full dataset
XX      <- as.matrix( patch[,-c(1,2)] )
mu.X    <- apply(XX, 2, mean)
XX      <- sscale * (XX - matrix( mu.X, dim(XX)[1], dim(XX)[2], byrow = TRUE) )
rownames(XX) <- patch[,1]
XX.list <- as.data.frame(t(XX))

#......................................................................
# available points to be selected
ZZ           <- as.matrix( PTS[,-c(1,2)] )
ZZ           <- sscale * (ZZ - matrix(mu.X, dim(ZZ)[1], dim(ZZ)[2], byrow = TRUE ) )
rownames(ZZ) <- PTS[,1]
ZZ.list      <- as.data.frame( t(ZZ) )

#......................................................................
# data parameters

nn.tot   <- dim(XX)[1]              # total number of points
nn       <- dim(ZZ)[1]              # number of allowable points 
mm.tot   <- 500                     # total size of subsample
pdim     <- dim(XX)[2]              # dimension of points

#......................................................................
#  initialize sample by selecting points that are as uniform as possible
#  Use min-max algorithm:  maximize minimal distance of cloud to anchors
#  run with small batches to speed up calculations

# parameters
m.batch     <- 5
n.iter.0    <- ceiling( mm.tot/m.batch )
bound.dist <- floor(0.05*nn)   # 2.5% of points nearby used for density estimate
k.star     <- 2                # points in high dimensions are all close  

# initialize points
idx.Z     <- 1:nn
idx.u     <- sample( idx.Z, m.batch )
UU        <- as.matrix( ZZ[idx.u, ] )

# apply min-max algorithm in batch mode
for ( k in 2:n.iter.0 ){
  dd      <- dista(ZZ,UU)
  dm      <- rowMins(dd, value = TRUE)
  dO      <- Order(dm, descending = TRUE)
  idx.new <- dO[1:m.batch] 
  idx.u   <- c( idx.u, idx.new )
  idx.X   <- idx.Z[ !(idx.Z %in% idx.new) ]
  UU      <- as.matrix( ZZ[idx.u, ])
  print(k)
}

# set rownames for UU
rownames(UU) <- PTS.label[ idx.u ]
UU.list   <- as.data.frame( t(UU) )

#......................................................................
# estimate of the density

vv   <- 3    # variance (chosen by looking at the data)
fhat <- rep(0,mm.tot)
for ( k in 1:mm.tot ){
  iidx    <- sort.list( abs(UU[k,]))[1:3]
  idx.c1  <- abs(XX[,iidx[1]]-UU[k,iidx[1]]) < 0.7
  idx.c2  <- abs(XX[,iidx[2]]-UU[k,iidx[2]]) < 0.7
  idx.c3  <- abs(XX[,iidx[3]]-UU[k,iidx[3]]) < 0.7
  idx.cc  <- idx.c1 & idx.c2 & idx.c3
  dd      <- dista(matrix(UU[k,],1,pdim),XX[idx.cc,])
  fhat[k] <- mean( exp( -0.5*dd/4 - pdim*log(2*pi)/2 - pdim*log(vv)/2 ) ) 
  print(k)
}

#......................................................................
# calculate the list of nearest neighbors (from admissible set)
# patch.names are used instead of indice

dd.idx  <- rownames(UU)[ dista( ZZ, UU, k = 1, index = TRUE) ]
NN.list <- split( rownames(ZZ), dd.idx ) 
ww        <- sapply(NN.list, length)/nn


#......................................................................
# calculate KK (quadratic kernel)
KK   <- matrix(0,mm.tot,mm.tot)
for ( k in 1:mm.tot ){
  KK[k,] <- sapply( UU.list, FUN=kp, y=UU.list[[k]] )
}

# loss function
L0 <- -2*sum(fhat*ww) + ww %*% KK %*% ww


#======================================================================
#  initialize MCMC
#  number of iteration
niter  <- 1000

#  initialize variables used in MCMC
LL       <- rep(0, niter+1 )               # Loss function
LL[1]    <- L0
Smp      <- matrix("", niter+1, mm.tot )   # matrix of subsample names
Smp[1, ] <- names(ww)
WW       <- matrix(0, niter+1, mm.tot )    # matrix of weights
WW[1, ]  <- ww


#......................................................................
#  simulated annealing cooling sequence

C1 <- 1000      # scale. depends on size of the loss function
Ca <- 2         # lower value of cooling sequence
Cb <- 20        # upper value of cooling sequence
cooling <- C1*log( seq( Ca, Cb,length=niter) )

#......................................................................
#  selection process:  select in index, and replace the "center" by an index in its neighborhood
#  as defined by NN.list
#  everything is done using patch.names to ensure that there are no indice mismatch

for ( k in 1:niter ){
  # select which element from the subsample to replace by a randomly select point from the data cloud
  smp1 <- sample(1:mm.tot,1)             # select one point in subsample
  smp2 <- sample( NN.list[[smp1]], 1)    
  
  UUnew  <- UU
  uu.new <- ZZ[ match( smp2, rownames(ZZ) ), ]
  
  ####
  
  UUnew[ rownames(UUnew) %in% names(NN.list)[smp1] , ] <- ZZ[ match( smp2, rownames(ZZ) ), ]
  
  smp2 <- sample(1:nn,1)      # select one point in the sample
  iidx.new <- Smp[k,]
  iidx.new[ smp1 ] <- smp2
  
  # update subsample and weights
  AA.new          <- AA
  AA.new[ smp1, ] <- XX0[ smp2, ]

  # XX0 and AA are both matrix arrays of type double;
  # for each point in XX0, find the index of the closest
  # neighbor in AA.new using a kd-tree for faster search
  # than brute force
  VV.new          <- nearest_neighbor_inds_kd_tree(AA=AA.new, XX0=XX0)
  VV.idx.new      <- iidx.new[ VV.new ]
  ww.new          <- Table( VV.idx.new, iidx.new )/nn
  
  # update K-matrix
  KK.new                <- KK
  AA.list.new           <- AA.list
  AA.list.new[[ smp1 ]] <- XX0[smp2, ]
  kk.new                <- sapply( AA.list.new, FUN=kp, y=AA.list.new[[smp1]], ss=data.scale )
  KK.new[smp1, ]        <- kk.new
  KK.new[, smp1]        <- kk.new
  
  #  calculate loss 
  LL.new <- t(ww.new) %*% KK.new %*% ww.new - lambda * sum(ww.new*ww.new)
  RR     <- exp( cooling[k]*( LL.new - LL[k] ) )
  U      <- runif(1)
  if ( RR < U ){
    # reject move
    
    Smp[k+1,] <- Smp[k,]
    LL[k+1]   <- LL[k]
    WW[k+1,]  <- WW[k,]
    
    print(paste("iteration = ",k,"  reject rr = ",round(RR,5)))
  } else {
    # accept move
    
    Smp[k+1,] <- iidx.new
    LL[k+1]   <- LL.new
    WW[k+1,]  <- ww.new
    
    # update K matrix and AA
    AA        <- AA.new
    AA.list   <- AA.list.new
    KK        <- KK.new
    
    print(paste("iteration = ",k,"  accept rr = ",round(RR,5)))
  }
}


#  make a few plots
plot( 1:(niter+1), LL )


