#' clusterlab
#'
#' This function runs clusterlab which is a simulator for Gaussian clusters. It is simple yet highly flexible allowing
#' the user to choose parameters such as the variance and size of each cluster individually. Clusterlab allows 
#' precise control over the distance between the N clusters.
#' 
#' @param centers Numerical value: the number of clusters to simulate (N)
#' @param r Numerical value: the number of units of the radius of the circle on which the clusters are generated
#' @param sdvec Numerical vector: standard deviation of each cluster, N values are required
#' @param alphas Numerical vector: how many units to push each cluster away from the initial placement, N values are required
#' @param centralcluster Logical flag: whether to place a cluster in the middle of the rest
#' @param numbervec Numerical vector: the number of samples in each cluster, N values are required
#' @param features Numerical value: the number of features for the data
#' @param seed Numerical value: fixes the seed if you want to repeat results, set the seed to 123 for example here
#' @param rings Numerical value: the number of concentric rings to generate (previous settings apply to all ring clusters)
#' @param ringalphas Numerical vector: a vector of numbers to push each ring out by, must equal number of rings
#' @param ringthetas Numerical vector: a vector of angles to rotate each ring by, must equal number of rings
#'
#' @return A list, containing: 
#' 1) the synthetic data
#' 2) cluster membership matrix
#' @export
#'
#' @examples
#' synthetic <- clusterlab(centers=4,r=8,sdvec=c(2.5,2.5,2.5,2.5),   
#' alphas=c(1,1,1,1),centralcluster=FALSE,   
#' numbervec=c(50,50,50,50)) # for a six cluster solution)   
#' 

clusterlab <- function(centers=1,r=8,sdvec=NULL,alphas=NULL,centralcluster=FALSE,
                       numbervec=NULL,features=500,seed=NULL,rings=NULL,ringalphas=NULL,
                       ringthetas=NULL){
  
  message('running clusterlab...')
  
  if (is.null(seed) == TRUE){ # user does not give extra annotation data
    set.seed(seed)
  }
  if (is.null(sdvec) == TRUE){ # user does not give extra annotation data
    message('user has not set standard deviation of clusters, setting automatically...')
    sdvec <- rep(1,centers)
  }
  if (is.null(alphas) == TRUE){ # user does not give extra annotation data
    message('user has not set alphas of clusters, setting automatically...')
    alphas <- rep(1,centers)
  }
  if (length(numbervec) != centers){
    message('user has not set length of numbervec equal to number of clusters, setting automatically...')
    numbervec <- rep(100,centers)
  }
  if (length(sdvec) != centers){
    message('user has not set length of sdvec equal to number of clusters, setting automatically...')
    sdvec <- rep(1,centers)
  }
  if (length(alphas) != centers){
    message('user has not set length of alphas equal to number of clusters, setting automatically...')
    alphas <- rep(1,centers)
  }
  if (centers != 1){
    if (is.null(alphas) == TRUE){ # user does not give extra annotation data
      message('user has not set alphas of clusters, setting automatically...')
      alphas <- rep(1,centers)
    }
  }
  if (abs(max(sdvec) - min(sdvec)) != 0 & is.null(rings) == FALSE){
    message('multiple ring method does not currently allow individual cluster variances, setting equal to first element')
    sdvec <- rep(sdvec[1],centers)
  }
  if (abs(max(alphas) - min(alphas)) != 0 & is.null(rings) == FALSE){
    message('multiple ring method does not currently allow individual alphas, setting equal to first element')
    message('please use the ringalphas parameter to seperate the rings by a specified degree...')
    alphas <- rep(alphas[1],centers)
  }
  if (abs(max(numbervec) - min(numbervec)) != 0 & is.null(rings) == FALSE){
    message('multiple ring method does not currently allow individual cluster sizes, setting equal to first element')
    numbervec <- rep(numbervec[1],centers)
  }
  if (is.null(rings) == FALSE & centralcluster == TRUE){
    message('ring method does not currently allow a central cluster to be generated, skipping')
  }
  if (is.null(rings) == FALSE & is.null(ringalphas) == TRUE){
    message('ring alphas not set, setting...')
    ringalphas <- seq(2,rings*2,2)
  }
  if (is.null(rings) == FALSE & is.null(ringthetas) == TRUE){
    message('ring thetas not set, setting...')
    ringthetas <- rep(0,rings)
  }
  
  if (is.null(rings) == TRUE){ # doing without rings 
    if (centers != 1){ # we are generating more than one cluster
      # create N points on a circle in 2D space that are evenly spaced which are samples
      if (centralcluster==TRUE){ # if we are having a cluster at 0,0
        matrix <- matrix(nrow=centers-1,ncol=2)
        n = centers-1
        i = 1
        for (x in seq(0,n-1)){ # edited to N instead of N-1
          x1 <- cos(2*pi/n*x)*r
          y1 <- sin(2*pi/n*x)*r
          matrix[i,1] <- x1
          matrix[i,2] <- y1
          i = i + 1
        }
        matrix <- rbind(matrix,c(0,0)) # add a center at the 0,0 co ordinate 
      }else{ # no central cluster
        matrix <- matrix(nrow=centers,ncol=2)
        n = centers
        i = 1
        for (x in seq(0,n-1)){ # edited to N instead of N-1
          x1 <- cos(2*pi/n*x)*r
          y1 <- sin(2*pi/n*x)*r
          matrix[i,1] <- x1
          matrix[i,2] <- y1
          i = i + 1
        }
      }
      # use scalar multiplication to pull the centers further apart by a parameter, alpha
      matrix2 <- matrix(nrow=nrow(matrix),ncol=2)
      for (row in seq(1,nrow(matrix))){
        matrix2[row,] <- matrix[row,]*alphas[row] # scalar multiplication
      }
      # using rnorm create new points based on each center with SD=X,mean=0 (shifting point then saving)
      newsamplematrix <- matrix(ncol=2,nrow=0)
      identitymatrix <- matrix(ncol=2,nrow=sum(numbervec)) # hold the cluster identity (ground truth)
      c <- 1
      for (j in seq(1,nrow(matrix2))){ # for each center
        z <- numbervec[j] # select the number of samples for each cluster
        sd <- sdvec[j]
        for (sample in seq(1,z)){ # loop to create new samples by adding noise
          newsample <- matrix2[j,]+rnorm(2,mean=0,sd=sd) # add the noise to the cluster center
          newsamplematrix <- rbind(newsamplematrix,newsample) # add the newsample to the new sample matrix
          identitymatrix[c,1] <- paste('c',j,'s',sample,sep='')
          identitymatrix[c,2] <- j
          c = c + 1
        }
      }
      identitymatrix <- data.frame(identitymatrix)
      colnames(identitymatrix) <- c('sampleID','cluster')
    }else{ # this is just for one cluster
      matrix2 <- matrix(nrow=1,ncol=2)
      matrix2[1,1] <- 0
      matrix2[1,2] <- 0
      # using rnorm create new points based on each center with SD=X,mean=0 (shifting point then saving)
      newsamplematrix <- matrix(ncol=2,nrow=0)
      identitymatrix <- matrix(ncol=2,nrow=sum(numbervec)) # hold the cluster identity (ground truth)
      c <- 1
      for (j in seq(1,nrow(matrix2))){ # for each center
        z <- numbervec[j] # select the number of samples for each cluster
        sd <- sdvec[j]
        for (sample in seq(1,z)){ # loop to create new samples by adding noise
          newsample <- matrix2[j,]+rnorm(2,mean=0,sd=sd) # add the noise to the cluster center
          newsamplematrix <- rbind(newsamplematrix,newsample) # add the newsample to the new sample matrix
          identitymatrix[c,1] <- paste('c',j,'s',sample,sep='')
          identitymatrix[c,2] <- j
          c = c + 1
        }
      }
      identitymatrix <- data.frame(identitymatrix)
      colnames(identitymatrix) <- c('sampleID','cluster')
    }
  }else{ # otherwise we are generating rings
    message('we are generating clusters arranged in rings...')
    # make a new alpha vector for each ring change alpha
    alphalist <- list()
    alphasoriginal <- alphas
    for (ring in seq(1,rings)){
      #alphas <- c(alphas,alphasoriginal*rings) # choose the ring pull apart factor here
      alphalist[[ring]] <- alphasoriginal*ringalphas[ring]
    }
    ringmatrix <- matrix(nrow=0,ncol=2) # to hold all the co ordinates of all points
    for (ring in seq(1,rings)){ # for every ring create N points
      matrix <- matrix(nrow=centers,ncol=2) # rings does not currently support central cluster
      n = centers
      i = 1
      for (x in seq(0,n-1)){ # edited to N instead of N-1
        x1 <- cos(2*pi/n*x)*r
        y1 <- sin(2*pi/n*x)*r
        matrix[i,1] <- x1
        matrix[i,2] <- y1
        # if we are rotating
        if (is.null(ringthetas)==FALSE){
          matrix[i,1] <- x1*cos(ringthetas[ring])-y1*sin(ringthetas[ring]) # x co-ordinate
          matrix[i,2] <- x1*sin(ringthetas[ring])+y1*cos(ringthetas[ring]) # y co-ordinate
        }
        i = i + 1
      }
      # find the centroid i.e. origin (approx 0,0)
      #print(mean(matrix[,1]))
      #print(mean(matrix[,2]))
      
      # use scalar multiplication to pull the centers further apart by a parameter, alpha
      # need to create new vector of alphas because of the rings
      alphas <- alphalist[[ring]]
      matrix2 <- matrix(nrow=nrow(matrix),ncol=2)
      for (row in seq(1,nrow(matrix))){
        matrix2[row,] <- matrix[row,]*alphas[row] # scalar multiplication
      }
      ringmatrix <- rbind(ringmatrix, matrix2)
    }
    
    matrix2 <- ringmatrix # over writing matrix2
    numbervec <- rep(numbervec,rings) # increase numbervec because of rings
    sdvec <- rep(sdvec,rings) # increase sdvec because of rings

    # using rnorm create new points based on each center with SD=X,mean=0 (shifting point then saving)
    newsamplematrix <- matrix(ncol=2,nrow=0)
    identitymatrix <- matrix(ncol=2,nrow=sum(numbervec)) # hold the cluster identity (ground truth)
    c <- 1
    for (j in seq(1,nrow(matrix2))){ # for each center
      z <- numbervec[j] # select the number of samples for each cluster
      sd <- sdvec[j]
      for (sample in seq(1,z)){ # loop to create new samples by adding noise
        newsample <- matrix2[j,]+rnorm(2,mean=0,sd=sd) # add the noise to the cluster center
        newsamplematrix <- rbind(newsamplematrix,newsample) # add the newsample to the new sample matrix
        identitymatrix[c,1] <- paste('c',j,'s',sample,sep='')
        identitymatrix[c,2] <- j
        c = c + 1
      }
    }
    #print(nrow(identitymatrix))
    #print(nrow(newsamplematrix)) # this is correct
    identitymatrix <- data.frame(identitymatrix)
    colnames(identitymatrix) <- c('sampleID','cluster')
  }
  
  # reverse PCA to generate N dimensional synthetic dataset
  n2 <- features
  x <- rnorm(n2, mean = 0, sd = 0.1) # fixed eigenvector1
  y <- rnorm(n2, mean = 0, sd = 0.1) # fixed eigenvector2
  matrix2 <- newsamplematrix
  res = matrix(nrow = nrow(matrix2), ncol = n2)
  for (i in seq(1,nrow(matrix2),1)){
    a <- matrix2[i,1] # get fixed co ordinate1 for entire row
    b <- matrix2[i,2] # get fixed co ordinate2 for entire row
    xk <- x
    yk <- y
    answer <- a*xk + b*yk
    res[i,] <- answer
  }
  mydata <- as.data.frame(res)
  pca1 = prcomp(mydata)
  scores <- data.frame(pca1$x) # PC score matrix
  # plot example
  p <- ggplot(data = scores, aes_string(x = 'PC1', y = 'PC2', colour = identitymatrix$cluster) ) + geom_point() +
    theme_bw() + 
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 18, colour = 'black'),
          axis.text.x = element_text(size = 18, colour = 'black'),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))
  print(p)
  mydata <- data.frame(t(mydata))
  colnames(mydata) <- identitymatrix$sampleID
  
  message('finished.')
  
  newlist <- list('synthetic_data' = mydata, 'identity_matrix' = identitymatrix)
  return(newlist)
}

