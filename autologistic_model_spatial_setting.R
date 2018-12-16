### Libraries
require(spdep)

# Data comes from a study in agriculture. The primary aim is
# the study of the spatial structure of the Carex Arenaria in a field
# with different species of plants. We also want to assess if
# the abundance of other plants helps in predicting the presence
# of Carex Arenaria in a certain location.

### Dataset
carex = read.table("carex.txt", header=TRUE)
colnames(carex) = c("x", "y", "carex_presence", "abund")
str(carex)
summary(carex)
attach(carex)

# aim: 

### Data visualization
par(mfrow=c(1,1), pty="m")
plot(x,y,type="n", lab=c(20,20,7),	xlab="X", ylab="Y")
text(x,y, labels=as.character(round(abund,0)))
for (i in 1:576){						   
  if (carex_presence[i]==1) 					    
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
             y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0)	   
} 
title("Abundance of the same plant", cex.main=1.5)
mtext(side=3, line=0.3, "boxes highlight the Carex arenaria presence",
      cex=1.0)

### Converting the dataframe in a SpatialPixelsDataFrame
coordinates(carex) = ~x+y	 
gridded(carex) = T
class(carex)

### 1st order neighborhood system ----
# I use the dnearneigh function in order to build a 1st order neighborhood.
par(mfrow=c(1,1))
carex.dn = dnearneigh(coordinates(carex),0,1)
plot(carex.dn, coords=coordinates(carex))
title("1st order neighborhood system")

# Informations about the neighbour list object.
carex.dn

# The average number of links is 3.833, because the nodes at the borders
# have less links with respect to the ones in the centre.

# This problem can be solved using a weighting system for each link.

# Converting the carex.dn variable in a matrix
carmat = nb2mat(carex.dn, style="B")	

# Creating 2 new variables needed for the weight.
n_sick = rep(0,nrow(carex))			
n_nbr = rep(0,nrow(carex))

# n_sick is the sum of neighbours with Carex Arenaria.
# n_nbr is the sum of neighbours.
for (i in 1:nrow(carex)){					
  n_sick[i] = sum(carmat[,i]*carex$carex_presence)
  n_nbr[i] = sum(carmat[,i])			
}					

# re-weighting the informations about the neighbors, taking into account
# which of them have Carex Arenaria and normalizing in order to avoid
# problems with the borders
ni = 4*n_sick/n_nbr	
summary(ni)

### Autologistic regression: Null model that uses only the spatial dependence
car1.null = glm(carex_presence~ni, data=carex, family=binomial)
summary(car1.null)

### AL(1): spatial dependenc and abund as covariate
car1.abund = glm(carex_presence~ni+abund, data=carex, family=binomial)
summary(car1.abund)

# Null model prediction
# Now we can check the performance of the two model by looking at the
# confusion matrix. We will compare the correct observations.
# I'm using different values of probabilities to classify the location as 1.
pr = predict(car1.null, type="response")
summary(pr)

# I use the range of the possible probabilities for the logistic regression.
# I'll choose the one that gives the highest accuracy in the prediction.
pp = seq(0.19,0.6,0.01)
fit1 = matrix(0, nrow(carex), length(pp))
accuracy = rep(0, length(pp))
for (i in 1:length(pp)) {
  fit1[,i] = ifelse(pr > pp[i], 1, 0)
}
for (i in 1:length(pp)) {
  t = table(carex$carex_presence, fit1[,i])
  accuracy[i] = (t[1,1]+t[2,2])/nrow(carex)*100
}
best_value = pp[which.max(accuracy)]
fit1 = ifelse(pr > best_value, 1, 0)
table(carex$carex_presence, fit1)
# 40 false positive
# 133 false negative
best_value # 0.38
round(max(accuracy),2) 
# The accuracy is 69.97% and there are 40 exact predictions of 
# the Carex Arenaria presence. I will use 038 as threshold also in the
# following models in order to compare the results.

# Model 2 prediction
pr = predict(car1.abund, type="response")
summary(pr)
pp = seq(0.1,0.7,0.01)
fit1.abund = matrix(0, nrow(carex), length(pp))
accuracy = rep(0, length(pp))
for (i in 1:length(pp)) {
  fit1.abund[,i] = ifelse(pr > pp[i], 1, 0)
}
for (i in 1:length(pp)) {
  t = table(carex$carex_presence, fit1.abund[,i])
  accuracy[i] = (t[1,1]+t[2,2])/nrow(carex)*100
}
best_value = pp[which.max(accuracy)]
fit1.abund = ifelse(pr > best_value, 1, 0)
table(carex$carex_presence, fit1.abund)
# 14 false positive
# 161 false negative
best_value # 0.55
round(max(accuracy),2) # 69.62%
# The accuracy is 69.62% but there are only 14 exact predictions of
# the Carex Arenaria presence. I'll use the same value of the null model
# for the probability in order to increase the exact predictions.

fit1.abund = ifelse(pr > 0.38, 1, 0)
table(carex$carex_presence, fit1.abund)
# 72 true positive
# 82 false positive
# 103 false negative
# 67.88% accuracy

### Predictive map
# In the following plots the green boxes highlight the correct predictions,
# the red ones are highlight the false negative and the black ones
# highlight the false positive.
par(mfrow=c(1,2)) 
plot(x,y,type="n", lab=c(20,20,7),	xlab="X", ylab="Y")
text(x,y, labels=as.character(round(abund,0)))
for (i in 1:576){						  
  if (carex_presence[i]==1 && fit1[i]==1) 					     
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
             y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0, 
             col="green", lwd=3)	
  if (carex_presence[i]==1 && fit1[i]==0) 					     
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
             y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0, 
             col="red", lwd=3)	
  if (carex_presence[i]==0 && fit1[i]==1) 					     
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
             y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0, 
             col="black", lwd=3)
}
title("1st Order Neigh. - no Covariate", cex.main=1.5)
mtext(side=3, line=0.3, 
      "Green=correct - Red=false negative - Black=false positive",
      cex=1.0)
plot(x,y,type="n", lab=c(20,20,7),	xlab="X", ylab="Y")
text(x,y, labels=as.character(round(abund,0)))
for (i in 1:576){						  
  if (carex_presence[i]==1 && fit1.abund[i]==1) 					     
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
             y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0, 
             col="green", lwd=3)	
  if (carex_presence[i]==1 && fit1.abund[i]==0) 					     
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
             y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0, 
             col="red", lwd=3)	
  if (carex_presence[i]==0 && fit1.abund[i]==1) 					     
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
             y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0, 
             col="black", lwd=3)
}
title("1st Order Neigh. - with Covariate", cex.main=1.5)
mtext(side=3, line=0.3,
      "Green=correct - Red=false negative - Black=false positive",
      cex=1.0)

# The second algorithm has less accuracy, but is capable to dected more
# exact locations having the Carex Arenaria with respect to the first 
# algorithm, meaning that the covariate is useful to improve the prediction.
# However the false positive cases are more in the second model.

### 2nd order neighborhood system -----
# Now I use 2nd nrighborhood system to see if I can improve the prediction.
par(mfrow=c(1,1))
carex.dn = dnearneigh(coordinates(carex), 0, 1.5)
plot(carex.dn, coords=coordinates(carex))
title("2nd order neighborhood system")

# As I did before, I use a weighting system to correct the borders effect.
carmat = nb2mat(carex.dn, style="B")	

# Information about the neighborhood
carex.dn

n_sick = rep(0,nrow(carex))			
n_nbr = rep(0,nrow(carex))

for (i in 1:nrow(carex)){					
  n_sick[i] = sum(carmat[,i]*carex$carex_presence)
  n_nbr[i] = sum(carmat[,i])			
}					
# weights for the neighborhood structure
ni = 8*n_sick/n_nbr	

### Null model: spatial dependence only
car2.null = glm(carex_presence~ni, data=carex, family=binomial)
summary(car2.null)

### model 2: spatial dependence and abundance as covariate 
car2.abund = glm(carex_presence~ni+abund, data=carex, family=binomial)
summary(car2.abund)

### null model prediction: only spatial dependence
# I will use 0.38 ad threshold probability to decide if a locations has
# value 1.
pr = predict(car2.null, type="response")
fit2 = ifelse(pr > 0.38, 1, 0)
table(carex$carex_presence, fit2)
# 77 false positive
# 106 false negative
# Accuracy: 68.23%

### model 2 prediction: spatial dependence and covariate
pr = predict(car2.abund, type="response")
fit2.abund = ifelse(pr > 0.38, 1, 0)
table(carex$carex_presence, fit2.abund)
# 87 false positive
# 94 false negative
# Accuracy: 68.58%

### Predictive map
# In the following plots the green boxes highlight the correct predictions,
# the red ones are highlight the false negative and the black ones
# highlight the false positive.
par(mfrow=c(1,2)) 
plot(x,y,type="n", lab=c(20,20,7),	xlab="X", ylab="Y")
text(x,y, labels=as.character(round(abund,0)))
for (i in 1:576){						  
  if (carex_presence[i]==1 && fit2[i]==1) 					     
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
             y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0, 
             col="green", lwd=3)	
  if (carex_presence[i]==1 && fit2[i]==0) 					     
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
             y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0, 
             col="red", lwd=3)	
  if (carex_presence[i]==0 && fit2[i]==1) 					     
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
             y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0, 
             col="black", lwd=3)
}
title("2nd Order Neigh. - no Covariate", cex.main=1.5)
mtext(side=3, line=0.3, 
      "Green=correct - Red=false negative - Black=false positive",
      cex=1.0)
plot(x,y,type="n", lab=c(20,20,7),	xlab="X", ylab="Y")
text(x,y, labels=as.character(round(abund,0)))
for (i in 1:576){						  
  if (carex_presence[i]==1 && fit2.abund[i]==1) 					     
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
             y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0, 
             col="green", lwd=3)	
  if (carex_presence[i]==1 && fit2.abund[i]==0) 					     
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
             y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0, 
             col="red", lwd=3)	
  if (carex_presence[i]==0 && fit2.abund[i]==1) 					     
    polygon (x[i]+c(-0.4, -0.4, 0.4, 0.4), 	   
             y[i]+c(-0.45, 0.4, 0.4, -0.45), density=0, 
             col="black", lwd=3)
}
title("2nd Order Neigh. - with Covariate", cex.main=1.5)
mtext(side=3, line=0.3, 
      "Green=correct - Red=false negative - Black=false positive",
      cex=1.0)

# The second model has an accuracy equal to 68.58%, slightly better 
# with respect to the first one with 68.23%. The second one also
# has less false negative predictions. 

### Final Choice of the model -----

# In order to choose the "best" model (we have no specification and the
# choice can be different with repect to what we need), it's useful to
# analize the confusion matrix for each one. I will compare the models 
# that use the covariate because their superior performance has already
# been stated in the previous analysis.

table(carex$carex_presence, fit1.abund)
# 319 correct negative
# 72 correct positive
## Accuracy: 67.88%
# 82 false positive
# 103 false negative

table(carex$carex_presence, fit2.abund)
# 314 correct negative
# 81 correct positive
## Accuracy: 68.58%
# 87 false positive
# 94 false negative

# The algorithm that uses a 2nd order neighborhood system and a covariate
# has a slightly better accuracy, with only 5 more false positive case but
# 9 less false negative with respect to the one using a 1st order 
# neighborhood system. Overall the performance of the second algorithm is
# better and without specification, meaning that we don't know if we
# need to minimize false positive or false negative, I would choose this one.





