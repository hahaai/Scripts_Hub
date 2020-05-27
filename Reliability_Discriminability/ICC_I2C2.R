

#### ICC
# the input is data which is in long format. It  has three columns: iccdata, sess, and sub. 
# Each row is a unique combine of subject and session.
# The following script does ICC over "sess".
library('lme4')
lm=tryCatch(summary(lmer(data ~ sess + (1 | sub),data=iccdata)), error=function(e) 'NA')

Var_between=as.numeric(attributes(lm$varcor$sub)$stddev)^2
Var_within=as.numeric(attributes(lm$varcor)$sc)^2
ICC=(Var_between)/(Var_between+Var_within)


#### I2C2
# from trace_only_unba.R
# An n by p data matrix containing functional responses. Each row contains measurements 
# from a function for one observation at a set of grid points, and each column contains measurements
# of all functions at a particular grid point.

# the visit/session and subject Id should match in the data. for example, in the data,
# row1 is sub1-visit1
# row2 is sub1-visit2
# row3 is sub1-visit3
# row4 is sub2-visit1
# row5 is sub2-visit2
# row6 is sub2-visit3
# then the I2C2_visit should be 1 2 3 1 2 3
#      the I2C2_id should be 1 1 1 2 2 2
i2c2= I2C2(data,I = length(I2C2_id),J = length(I2C2_visit),visit = I2C2_visit, id =I2C2_id,demean = TRUE) 


#### discriminability
# from reliability.R
# where D is the pairwise distance matrix, ids should be a string vector corresponds to the subject id of rows of D.
rd <- rdf(D,ids)
discri <- mean(rd)
