## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----results="asis",tidy=FALSE,eval=TRUE,echo=FALSE----------------------
cat("[Simulating a single cluster](#test1)  \n")
cat("[Simulating four clusters with equal variances](#test2)  \n") 
cat("[Simulating four clusters with unequal variances](#test3)  \n")
cat("[Simulating four clusters with one cluster pushed to the outside](#test4)  \n")
cat("[Simulating four clusters with one small cluster](#test5)  \n")
cat("[Simulated five clusters with one central cluster](#test6)  \n")
cat("[Keeping track of cluster allocations](#test7)  \n")
cat("[Generating more complex multi ringed structures](#test8)  \n")
cat("[Afterthoughts](#test9)  \n")

## ----fig.width=6,fig.height=6--------------------------------------------
library(clusterlab)
synthetic <- clusterlab(centers=1,numbervec=100)

## ----fig.width=6,fig.height=6--------------------------------------------
library(clusterlab)
synthetic <- clusterlab(centers=4,r=8,sdvec=c(2.5,2.5,2.5,2.5),   
                        alphas=c(1,1,1,1),centralcluster=FALSE,   
                        numbervec=c(50,50,50,50))

## ----fig.width=6,fig.height=6--------------------------------------------
library(clusterlab)
synthetic <- clusterlab(centers=4,r=8,sdvec=c(1,1,2.5,2.5),   
                        alphas=c(1,1,1,1),centralcluster=FALSE,   
                        numbervec=c(50,50,50,50))

## ----fig.width=6,fig.height=6--------------------------------------------
library(clusterlab)
synthetic <- clusterlab(centers=4,r=8,sdvec=c(2.5,2.5,2.5,2.5),   
                        alphas=c(1,2,1,1),centralcluster=FALSE,   
                        numbervec=c(50,50,50,50))

## ----fig.width=6,fig.height=6--------------------------------------------
library(clusterlab)
synthetic <- clusterlab(centers=4,r=8,sdvec=c(2.5,2.5,2.5,2.5),   
                        alphas=c(1,1,1,1),centralcluster=FALSE,   
                        numbervec=c(15,50,50,50))

## ----fig.width=6,fig.height=6--------------------------------------------
library(clusterlab)
synthetic <- clusterlab(centers=5,r=8,sdvec=c(2.5,2.5,2.5,2.5,2.5),   
                        alphas=c(2,2,2,2,2),centralcluster=TRUE,   
                        numbervec=c(50,50,50,50,50))

## ----fig.width=6,fig.height=6--------------------------------------------
head(synthetic$identity_matrix)

## ----fig.width=6,fig.height=6--------------------------------------------
library(clusterlab)
synthetic <- clusterlab(centers=5,r=7,sdvec=c(6,6,6,6,6),   
                        alphas=c(2,2,2,2,2),centralcluster=FALSE,   
                        numbervec=c(50,50,50,50),rings=5,ringalphas=c(2,4,6,8,10,12), 
                        ringthetas = c(30,90,180,0,0,0), seed=123) # for a six cluster solution)

