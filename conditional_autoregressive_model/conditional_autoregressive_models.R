### Conditional Autoregressive Models
### Spatial Statistics and Statistical Tools for Environmental Data

# I am reproducing the paper of Crassie and Chan 
# (SPATIAL MODELING OF REGIONAL VARIABLES, '68), where they analyze 
# analyze a data set of Sudden Infant Deaths, 1974 - 1978, in the counties 
# of North Carolina, using a Markovrandom-field approach to spatial modeling. 

### Libraries
library(sp)
library(maptools)
library(spdep)
library(stats)
require(maps)  

### Loading Dataset
data(nc.sids)				
names(nc.sids)

### Build the transformed rates variable 
sids.ft <- sqrt(1000*(nc.sids$SID74/nc.sids$BIR74))+	  
  sqrt(1000*(nc.sids$SID74+1)/nc.sids$BIR74)   

nwbir.ft <- sqrt(1000*(nc.sids$NWBIR74/nc.sids$BIR74))+	   
  sqrt(1000*(nc.sids$NWBIR74+1)/nc.sids$BIR74)	

### Creating dataframe
sids <- data.frame(long=nc.sids$lon, lat=nc.sids$lat,
                   sids.ft=sids.ft, nwbir.ft=nwbir.ft, nc.sids)

# Visualization af the data on the map. 
# Colors are related to the intensity of the phenomenon.
sids.ord <- c(27,41,2,85,1,22,57,28,96,100, 
              53,43,69,34,7,95,11,52,48,81,21,90,64,98, 
              91,82,4,4,4,56,42,40,88,30,33,25,24,76,8, 
              73,13,59,26,16,63,55,72,6,86,87,39,66,54,
              83,60,74,65,46,78,38,36,68,32,70,67,31,99,
              5,93,29,80,17,97,20,14,51,77,47,89,94,12,
              50,61,79,92,71,10,3,58,75,45,84,15,37,9,44,
              19,62,18,49,23,35)

palette(c("lightyellow","darkgoldenrod",    
          "orange","brown"))
breaks.sids <- c(0,2.2,3.0,3.5,7)       
sids.ygrp <- cut(sids$sids.ft, breaks.sids)  
sids.y <- rep(0,length(sids.ygrp)) 
sids.y[sids.ygrp=="(0,2.2]"] <- 1      
sids.y[sids.ygrp=="(2.2,3]"] <- 2           
sids.y[sids.ygrp=="(3,3.5]"] <- 3           
sids.y[sids.ygrp=="(3.5,7]"] <- 4           
sids.y <- c(sids.y[1:84],NA,sids.y[85:99])  
sids.y2 <- sids.y[sids.ord] 

map("county","north carolina",fill=T,col=sids.y2)  
title("SIDS Rates",cex=1.5,col=1)
legend(-82.5,35.05,legend=c("< 2.2","2.2-3.0",     
                            "3.0-3.5","> 3.5"),fill=1:4,cex=0.5)             

### Deleting the outlier
sids <- sids[sids$CNTY.ID!=2096,]
n = length(sids$sids.ft) # 99 counties

### Converting the dataframe in a SpatialPixelsDataFrame
coordinates(sids) = ~long+lat			
coords = coordinates(sids)	

### Creating a neighborhood structure
sids.nb <- dnearneigh(coords, 0, 48.28, longlat=T) 

### Visualization of the graph 
plot(coords)
plot(sids.nb, coords, add = T)

### Information about the neighborhood
sids.nb

# The average number of link is 3.86 and as we can see in the previous graph,
# some counties have no link with other counties. This is confirmed by 
# the informations provided by sids.nb.

### Converting the neighborhood structure to a matrix
sids.W.mat = nb2mat(sids.nb, style="B", zero.policy=T)

### Checking if the matrix is symmetric
isSymmetric(sids.W.mat)

# The matrix is not symmtric. We will build a symmetric adjacency
# matrix with weighted links.

# Computing the distances between counties in the neighborhood list
sids.dist <- nbdists(sids.nb, coords,longlat=T)     
dij <- unlist(sids.dist)                           

# Elements of the weighting system
term1 <- min(dij)/dij              

# Get id of neighboring counties by row and columns
row.id <- row(sids.W.mat)[sids.W.mat==1]      
col.id <- col(sids.W.mat)[sids.W.mat==1]          

# Numbers of new borns in 1974
numb <- sids$BIR74                               

# Second term of the weigthing system
term2 <- sqrt(numb[row.id]/numb[col.id])          

# Weights
wgts <- term1*term2                               

# Weights matrix
wgtmat <- matrix(0,nrow=n,ncol=n)                  

for (i in 1:length(wgts)){
  wgtmat[col.id[i],row.id[i]] <- wgts[i]                            
}

# The matrix must be converted into a weights list to be used 
# in some functions of spdep
sids.W <- mat2listw(wgtmat) 

isSymmetric(wgtmat)                 
# The matrix is still not symmetric.

# We build a diagonal matrix with the number of new borns per county. 
D<-diag(numb,99,99)
wmat <- D %*% wgtmat
sids.WD <- mat2listw(wmat)
isSymmetric(wmat)

# Now the matrix is symmetric. We will use this matrix forse the analysis
# as the adjacency matrix

### Pseudo-likelihood -----------
# Now I will find the regression coefficient following the pseudo-likelihood
# approach. Using the pseudo-likelihood is equivalent to performing a
# linear regression using the response variable and the weights as covariate,
# without using the intercept.

z = wmat %*% sids$sids.ft

### Null model of the regression
mod.null = lm(sids$sids.ft ~ z - 1)
summary(mod.null)

### Adding covariates
mod.cov = lm(sids$sids.ft ~ z + sids$nwbir.ft - 1)
summary(mod.cov)

# There an improvement in the second model.
# The null model has an R-squared equal to 0.6055, while the second one
# is equal to 0.8972.
# Let's compute the RMSE.

rmse.null = sqrt(mean((sids$sids.ft - mod.null$fitted.values)^2))
rmse.cov = sqrt(mean((sids$sids.ft - mod.cov$fitted.values)^2))

cbind(null_model=rmse.null, cov_model=rmse.cov)

# The rmse is less in the model with covariates.
# with respect to the null model. I will use this one to compare
# the results with the model obtained using the spautolm function.

### Plot to compare the predictions and the real values
plot(sids$sids.ft, mod.null$fitted.values, xlab="observed",     
     ylab="fitted", main="Null model")         
abline(0,1)

plot(sids$sids.ft, mod.cov$fitted.values, xlab="observed",     
     ylab="fitted", cex.lab=1.5,
     main="Model with covariates")         
abline(0,1)

### Checking the residuals of the second model
hist(mod.cov$residuals,xlab="Residuals SIDS", 
     ylab="Frequency",cex.lab=1.6,density=12, 
     main="Residuals histogram",cex.main=1.4)
qqnorm(scale(mod.cov$residuals), cex.lab=1.5,
       xlab="Normal ", ylab="Residuals SIDS",
       cex.main=1.5)        
qqline(scale(mod.cov$residuals))

# The residuals have a normal behaviour, meaning that the model that
# we assumed, a linear regression and a gaussian error, is good.

### Looking at the residuals
summary(residuals(mod.cov))

# The residuals have a maximum of 3.14, meaning that we can have an error
# of maximum 3.14 in the sids prediction, but the median and the mean 
# are around 0.2.

### Spatial CAR using spautolm ------------------------

# Now I will use the neighborhood structure of the previous analysis
# in the spautolm function.
sids.raceslm <- spautolm(sids.ft ~ nwbir.ft, data=sids, 
                         family="CAR", listw=sids.WD, weights=(BIR74))
summary(sids.raceslm)   

# I will compare the RMSE of this model with the one obtained using the
# pseudo-likelihood.

rmse.spautolm = sqrt(mean((sids$sids.ft - sids.raceslm$fit$fitted.values)^2))
cbind(cov_model=rmse.cov, spautolm=rmse.spautolm)

# As expected, the model that uses the likelihood approach is better than
# the one that uses the pseudo-likelihood.

# Residual analysis
hist(sids.raceslm$fit$residuals,xlab="Residuals SIDS", 
     ylab="Frequency",cex.lab=1.6,density=12, 
     main="Residuals histogram",cex.main=1.4)
qqnorm(scale(sids.raceslm$fit$residuals), cex.lab=1.5,
       xlab="Normal ", ylab="Residuals SIDS",
       cex.main=1.5)        
qqline(scale(sids.raceslm$fit$residuals))
# Also the residuals have a better normal behaviour using the likelihood.

### CAR model estimate using 1st order neighborhood -------------------------
# I will perform the same analysis as before using a different
# neighborhood structure.

### Creating a 1st order neighborhood structure
sids.nb = knearneigh(coords, k=4, longlat=T)

### Visualization of the graph 
plot(coords)
plot(knn2nb(sids.nb), coords, add = T)

### Information about the neighborhood
knn2nb(sids.nb)

# The average number of link is 4 and as is confirmed by 
# the informations provided by sids.nb.

### Converting the neighborhood structure to a matrix
sids.W.mat = nb2mat(knn2nb(sids.nb, sym=T), style="B", zero.policy=T)

### Checking if the matrix is symmetric
isSymmetric(sids.W.mat)

# The matrix is not symmetric. We will build a symmetric adjacency
# matrix with weighted links.

# Computing the distances between counties in the neighborhood list
sids.dist <- nbdists(knn2nb(sids.nb, sym=T), coords, longlat=T)     
dij <- unlist(sids.dist)                           

# Elements of the weighting system
term1 <- min(dij)/dij              

# Get id of neighboring counties by row and columns
row.id <- row(sids.W.mat)[sids.W.mat==1]      
col.id <- col(sids.W.mat)[sids.W.mat==1]          

# Numbers of new borns in 1974
numb <- sids$BIR74                               

# Second term of the weigthing system
term2 <- sqrt(numb[row.id]/numb[col.id])          

# Weights
wgts <- term1*term2                               

# Weights matrix
wgtmat <- matrix(0, nrow=n, ncol=n) 

for (i in 1:length(wgts)){
  wgtmat[col.id[i], row.id[i]] <- wgts[i]                            
}

# The matrix must be converted into a weights list to be used 
# in some functions of spdep
sids.W <- mat2listw(wgtmat) 

isSymmetric(wgtmat)                 
# The matrix is still not symmetric.

# We build a diagonal matrix with the number of new borns per county. 
D<-diag(numb, 99, 99)
wmat <- D %*% wgtmat
sids.WD <- mat2listw(wmat)
isSymmetric(wmat)

# Now I will use the neighborhood structure of the previous analysis
# in the spautolm function.
sids.raceslm.2 <- spautolm(sids.ft ~ nwbir.ft, data=sids, 
                         family="CAR", listw=sids.WD, weights=(BIR74))
summary(sids.raceslm.2)   

### Predicted values vs real values
plot(sids$sids.ft, sids.raceslm.2$fit$fitted.values, xlab="observed",     
     ylab="fitted", main="Null model")         
abline(0,1)

# I will compare the RMSE of the previous models
rmse.spautolm.2 = sqrt(mean((sids$sids.ft - sids.raceslm.2$fit$fitted.values)^2))
cbind(cov_model=rmse.cov, spautolm=rmse.spautolm, spautolm2=rmse.spautolm.2)

# The model spautolm2 has a slightly better RMSE with respect to the spautolm.
# Let's compare the residuals in order to choose the best one and decide which
# neighborhood structure can give us better predictions.

# Residual analysis spautolm vs spautolm.2 with 1st order neigh.
par(mfrow=c(2,2))
hist(sids.raceslm$fit$residuals,xlab="Residuals SIDS", 
     ylab="Frequency", main="Residuals histogram model 1")
qqnorm(scale(sids.raceslm$fit$residuals), cex.lab=1.5,
       xlab="Normal ", ylab="Residuals SIDS",
       cex.main=1.5)        
qqline(scale(sids.raceslm$fit$residuals))

hist(sids.raceslm.2$fit$residuals,xlab="Residuals SIDS", 
     ylab="Frequency", main="Residuals histogram model 2")
qqnorm(scale(sids.raceslm.2$fit$residuals), cex.lab=1.5,
       xlab="Normal ", ylab="Residuals SIDS",
       cex.main=1.5)        
qqline(scale(sids.raceslm.2$fit$residuals))

# There's no much difference between the plots.
# We can check the sum of residuals for both models.

cbind(sum(spautolm=sids.raceslm$fit$residuals),
      spautolm2=sum(sids.raceslm.2$fit$residuals))

# The RMSE and the sum of residuals of the second model is slightly
# less than the one of the first model. In the end I'll choose the
# second model that uses a 1st order neighborhood system.




