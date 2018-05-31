## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

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

