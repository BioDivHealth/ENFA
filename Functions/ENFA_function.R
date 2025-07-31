#
#~~~\>|>|>|>\\\~~~~~~~~~?/\/\/\\/\/\/|/\/\/\/~~~\>|>|>|>\\\~~~~~~~~~?/\/\/\\/\/\/|/\/\/\/
#
#                       Main Function to run the ENFA analysis
#
#~~~~~~~~~~~~?/\/\/\\/\/\/|/\/\/\/~~~\>|>|>|>\\\~~~~~~~~~?/\/\/\\/\/\/|/\/\/\/
#
# Adapted from the ENFA ade4 package function - see https://github.com/cran/adehabitatHS/blob/master/R/enfa.r for the original function
#
ENFA_function<-function(data, # Data.frame containing the environmental information with no NAs
                        presence_index, # an index of length equal to nrow(data) showing the presence values for each row (0=absence; 1-n= number of presence)
                        n_speciation_axes=2, # number of speciation axes to preserve, it should be in the range of 1 to ncol(data) (default=1)
                        row.weights= rep(1, nrow(data))/nrow(data), # weights for the rows of the data.frame
                        col.weights= rep(1, ncol(data))# weight for the variables of the data.frame
                        ){
  
  # 0. Load the needed packages
  list.of.packages<-c("tidyr","terra","sf","data.table","dplyr")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  ## Adapted from the ENFA ade4 package function - see https://github.com/cran/adehabitatHS/blob/master/R/enfa.r for the original function
  # Presence vector (0,1)
  prb <- presence_index
  presence_index <- presence_index/sum(presence_index) # Proportion of presences from the total in each cell
  
  # Proportion of total row weigths; if lw=rep(1, nrow(rast_vals$tab[,-1]))/nrow(rast_vals$tab[,-1]) their sum =1
  row.weights <- row.weights/sum(row.weights)

  # Convert the data.frame into matrix and subtract the sum of weighted rows from each variable 
  Z <- as.matrix(data)
  n <- nrow(Z)
  
  f1 <- function(v) sum(v * row.weights)
  
  center <- apply(Z, 2, f1) 
  Z <- sweep(Z, 2, center)
  Ze <- sweep(Z, 2, sqrt(col.weights), "*") # multiply by the square route of the col.weights
  
  # Ze is the environmental matrix that contains the data for the whole environment or study area
  
  ## Build the data matrices for the species presence pixels/points 
  # Inertia matrices S and G
  DpZ <- apply(Ze, 2, function(x) x*presence_index) # weighted species or inertia matrix
  
  ## Marginality computation
  # Marginality vector coordinates
  mar <- apply(Z,2,function(x) sum(x*presence_index)) # Calculate the marginality over the Row 
                                                      # weighted presence matrix. If weights are equal 
                                                      # to 1 is the same as with the inertia matrix
  # Weighted marginality vector coordinates
  me <- mar*sqrt(col.weights)
  
  # Crossproduct of the environmental and inertia matrizes
  Se <- crossprod(Ze, DpZ) 
  Ge <- crossprod(Ze, apply(Ze,2,function(x) x*row.weights))
  
  ## Computation of S^(-1/2)
  eS <- eigen(Se)
  kee <- (eS$values > 1e-9)     ## keep only the positive values
  S12 <- eS$vectors[,kee] %*% diag(eS$values[kee]^(-0.5)) %*% t(eS$vectors[,kee])
  
  ## Passage to the third problem
  W <- S12 %*% Ge %*% S12
  x <- S12 %*% me
  b <- x / sqrt(sum(x^2))
  
  ## Eigenstructure of H
  H <- (diag(ncol(Ze)) - b%*%t(b)) %*% W %*% (diag(ncol(Ze)) - b%*%t(b))
  s <- eigen(H)$values[-ncol(Z)]
  # barplot(s)
  
  ## Number of eigenvalues
  
  if (n_speciation_axes <= 0 | n_speciation_axes > (ncol(Ze) - 1)){
    n_speciation_axes <- 1}
  
  ## coordinates of the columns on the specialization axes
  co <- matrix(nrow = ncol(Z), ncol = n_speciation_axes + 1)
  tt <- data.frame((S12 %*% eigen(H)$vectors)[, 1:n_speciation_axes])
  tt <- apply(tt,2,as.numeric)
  ww <- apply(tt, 2, function(x) x/sqrt(col.weights))
  norw <- sqrt(diag(t(as.matrix(tt))%*%as.matrix(tt)))
  co[, 2:(n_speciation_axes + 1)] <- sweep(ww, 2, norw, "/")
  
  ## coordinates of the columns on the marginality axis
  m <- me/sqrt(col.weights)
  co[, 1] <- m/sqrt(sum(m^2))
  
  ## marginality
  m <- sum(m^2)
  
  ## Coordinates of the rows on these axes (values of marginality and specificity for the environmental matrix)
  li <- Z %*% apply(co, 2, function(x) x*col.weights)
  
  ## Output
  co <- as.data.frame(co) # vector coordinates
  li <- as.data.frame(li) # row.values
  names(co) <- c("Marginality", paste("Specialization", (1:n_speciation_axes), sep = ""))
  row.names(co) <- dimnames(data)[[2]]
  names(li) <- c("Marginality", paste("Specialization", (1:n_speciation_axes), sep = ""))

  ## Calculate the environmental distance to niche centroid (Mahalanobis distance)
  Zli <- li[, 1:(n_speciation_axes + 1)]
  f1 <- function(x) rep(x, prb)
  
  Sli <- apply(Zli, 2, f1)
  
  m_cord <- apply(Sli, 2, mean) # Specificity and marginality coordinates
  cov <- t(as.matrix(Sli)) %*% as.matrix(Sli)/nrow(Sli)
  maha <- data.frame(Mahalanobis.Dist = mahalanobis(Zli, center = m_cord, cov = cov))
  
  # Suitability map (check this)
  m_dist<-mahalanobis(Zli, center = m_cord, cov = cov)
  Suitability <- 1 - pchisq(m_dist, df=ncol(Zli))
  Suitability <- data.frame(Suitability = Suitability)
  
  # Return the results
  return(list(marginality=m,
              row.weights=row.weights,
              col.weights=col.weights,
              presence_index=presence_index,
              marginality_coordiantes=me,
              coordinates_axis=co,
              marginality_specificity_vals=li,
              niche_centroid_coordinates=m_cord,
              prediction=maha,
              Suitability=Suitability))
  }
