set.seed(12345)
rm(list=ls())
library(Matching)
library(haven)
library(dplyr)
library(parallel)
library(survey)
library(lattice)
library(ebal)

# Downloading and preparing the data
# Source: http://users.nber.org/~rdehejia/data/nswdata2.html

cps <- read_dta("http://www.nber.org/~rdehejia/data/cps_controls.dta")
nsw <- read_dta("http://www.nber.org/~rdehejia/data/nsw_dw.dta") %>%
  dplyr::filter(treat == 1)
foo <- rbind(nsw, cps) %>%
  dplyr::select(-data_id)

###################################
## Genetic Matching: Replication ##
###################################

# This is a slightly edited version of the replication files
# made available by (). Specifically, the relevant file 
# is "replicate.GenMatch.cps1.cluster1.popsize5000.wg20.R"
# I rename the variables to use the replication code as closely as possible.

foo$educ <- foo$education
foo$hispan <- foo$hispanic
foo <- foo %>% dplyr::select(-education, -hispanic)

X <-   as.matrix(cbind( foo$age,
                        foo$educ,
                        foo$black,
                        foo$hispan,
                        foo$married,
                        foo$nodegree,
                        I(foo$re74/1000),
                        #I(as.real(foo$re74 == 0)),
                        I(foo$re75/1000)))# ,
                        #I(as.real(foo$re75 == 0))))

BalanceMat <- as.matrix(cbind(
                      foo$age,
                      foo$educ,
                      foo$black,
                      foo$hispan,
                      foo$married,
                      foo$nodegree,
                      I(foo$re74/1000),
                      I(foo$re75/1000),
                      I((foo$re74/1000)^2),
                      I((foo$re75/1000)^2),
                      
                      I((foo$age)^2),
                      I((foo$educ)^2),
                      I((foo$black)^2),
                      I((foo$hispan)^2),
                      I((foo$married)^2),
                      I((foo$nodegree)^2),
                      
                      I(foo$age*foo$educ),
                      I(foo$age*foo$black),
                      I(foo$age*foo$hispan),
                      I(foo$age*foo$married),    
                      I(foo$age*foo$nodegree),
                      I(foo$age*(foo$re74/1000)),
                      I(foo$age*(foo$re75/1000)),
                      
                      I(foo$educ*foo$black),
                      I(foo$educ*foo$hispan),
                      I(foo$educ*foo$married),    
                      I(foo$educ*foo$nodegree),
                      I(foo$educ*(foo$re74/1000)),
                      I(foo$educ*(foo$re75/1000)),
                      
                      I(foo$black*foo$hispan),
                      I(foo$black*foo$married),
                      I(foo$black*foo$nodegree),
                      I(foo$black*(foo$re74/1000)),
                      I(foo$black*(foo$re75/1000)),
                      
                      I(foo$hispan*foo$married),
                      I(foo$hispan*foo$nodegree),
                      I(foo$hispan*(foo$re74/1000)),
                      I(foo$hispan*(foo$re75/1000)),
                      
                      I(foo$married*foo$nodegree),
                      I(foo$married*(foo$re74/1000)),
                      I(foo$married*(foo$re75/1000)),
                      
                      I(foo$nodegree*(foo$re74/1000)),
                      I(foo$nodegree*(foo$re75/1000)),                        
                      
                      I((foo$re74/1000)*(foo$re75/1000))))

### ORTHOGONALIZE TO PROPENSITY SCORE

# From Review of Econ and Statistics
model.A = foo$treat~I(foo$age) + I(foo$age^2) + I(foo$age^3) + I(foo$educ) + I(foo$educ^2) + I(foo$married) + I(foo$nodegree) + I(foo$black) + I(foo$hispan) + I(foo$re74) + I(foo$re75) + I(foo$re74 == 0) + I(foo$re75 == 0) + I(I(foo$re74)*I(foo$educ))

pscores.A <- glm(model.A, family = binomial)

X2 <- X
for (i in 1:ncol(X2)) 
  {
   lm1 = lm(X2[,i]~pscores.A$linear.pred)
   X2[,i] = lm1$residual
   }


# VARIABLES WE MATCH ON BEGIN WITH PROPENSITY SCORE
orthoX2.plus.pscore = cbind(pscores.A$linear.pred, X2)

orthoX2.plus.pscore[,1] <- orthoX2.plus.pscore[,1] -
                           mean(orthoX2.plus.pscore[,1])


# NORMALIZE ALL COVARS BY STANDARD DEVIATION
for (i in 1:(dim(orthoX2.plus.pscore)[2])) 
  {
    orthoX2.plus.pscore[,i] <-
      orthoX2.plus.pscore[,i]/sqrt(var(orthoX2.plus.pscore[,i]))
  }

# SET BOUNDS ON GENMATCH DOMAIN (OPTIONAL)
lower <- c(95, rep(0, dim(orthoX2.plus.pscore)[2] - 1))
upper <- c(105, rep(1000, dim(orthoX2.plus.pscore)[2] - 1))

Domains <- as.matrix(cbind(lower, upper))

set.seed(12345)

# I am using the line above instead of 'source("~/R/cluster1.R")' of the original replication files.
cl <- makeCluster(detectCores() - 1)

GM.out <-
  GenMatch(Tr = foo$treat,
           X = orthoX2.plus.pscore,
           BalanceMatrix = BalanceMat,
           pop.size = 5000,
           #           pop.size = 50,           
           max.generations = 25,
           wait.generations = 20,
           print.level = 3,
           cluster=cl)

#benchmark=$1794
stopCluster(cl)

mout <- Match(Y = foo$re78,
            Tr = foo$treat,
            X = orthoX2.plus.pscore,
            Weight.matrix = GM.out)
summary(mout)


########### I STOPPED HERE 
MatchBalance(treat ~ age + black + married + nodegree + re74 + re75 + educ + hispan, data = foo, match.out = mout) 

#############################################

########### ENTROPY BALANCING ###############

#############################################

dat <- foo

# unemployed
dat$u74 <- as.numeric(dat$re74==0)
dat$u75 <- as.numeric(dat$re75==0)

# covars
covars <- c("age","educ","black","hispan","married","nodegree","re74","re75","u74","u75")

# compute all interactions
X <- as.matrix(dat[,covars])
XX <- matrixmaker(X)

# prepare data: exclude non-sensical and co-linear vars, code treatment and outcome, change names
out <- c("black.black","nodegree.nodegree","nodegree.educ","married.married","hispan.hispan","hispan.black",
         "u75.u75","u74.u74","u74.re74","u75.re75","u74.re74","u75.re75","re74.re74","re75.re75","re75.re74")

XX     <- XX[,(colnames(XX) %in% out)==F]
dat    <- data.frame(dat[,(names(dat) %in% covars)==F],XX)
covars <- names(dat)[-which((names(dat) %in% c("treat","re78")))]
X     <- dat[,covars]
X     <- as.matrix(X)
colnames(X) <-gsub(".age","*Age",colnames(X))
colnames(X) <-gsub(".educ*","*Schooling",colnames(X))
colnames(X) <-gsub(".black","*Black",colnames(X))
colnames(X) <-gsub(".hispan","*Hispanic",colnames(X))
colnames(X) <-gsub(".married","*Married",colnames(X))
colnames(X) <-gsub(".re74","*Earnings 1974",colnames(X))
colnames(X) <-gsub(".re75","*Earnings 1975",colnames(X))
colnames(X) <-gsub(".u74","*Unemployed 1974",colnames(X))
colnames(X) <-gsub(".u75","*Unemployed 1975",colnames(X))
colnames(X) <-gsub(".nodegree","*HS Dropout",colnames(X))
colnames(X) <-gsub("age","Age",colnames(X))
colnames(X) <-gsub("educ","Schooling",colnames(X))
colnames(X) <-gsub("black","Black",colnames(X))
colnames(X) <-gsub("hispan","Hispanic",colnames(X))
colnames(X) <-gsub("married","Married",colnames(X))
colnames(X) <-gsub("re74","Earnings 1974",colnames(X))
colnames(X) <-gsub("re75","Earnings 1975",colnames(X))
colnames(X) <-gsub("u74","Unemployed 1974",colnames(X))
colnames(X) <-gsub("u75","Unemployed 1975",colnames(X))
colnames(X) <-gsub("nodegree","HS Dropout",colnames(X))
dat <- data.frame(dat,X)
Y <- dat$re78
W <- dat$tr

# Entropy Balancing
out.eb <- ebalance(
  Treatment=W,
  X=X
)

bout.eb <- MatchBalance(W~X,weights=c(W[W==1],out.eb$w),ks=FALSE)
bal.eb  <- baltest.collect(matchbal.out=bout.eb,var.names=colnames(X),after=F)
round(bal.eb,2)

# # The code below prouces an error
# Entropy Balancing (with trimmed weights)
out.ebtr <- ebal::ebalance.trim(out.eb, print.level = 3)

bout.ebtr <- MatchBalance(W~X,weights=c(W[W==1],out.ebtr$w),ks=FALSE)
bal.ebtr  <- baltest.collect(matchbal.out=bout.eb,var.names=colnames(X),after=FALSE)
round(bal.ebtr,2)

# Effect Estimates

dat$w   <-  c(W[W==1],out.eb$w)
des1 <- svydesign(id=~1,weights=~w, data= dat)
mod1 <- svyglm(re78 ~ treat, design = des1)
summary(mod1)

#############################################

## GENETIC MATHCING + ENTROPY BALANCING #####

#############################################

summary(mout)

dat.gen <- rbind(foo[mout$index.treated, ], foo[mout$index.control, ])

# We run the same process again...

# unemployed
dat.gen$u74 <- as.numeric(dat.gen$re74==0)
dat.gen$u75 <- as.numeric(dat.gen$re75==0)

# covars
covars <- c("age","educ","black","hispan","married","nodegree","re74","re75","u74","u75")

# compute all interactions
X <- as.matrix(dat.gen[,covars])
XX <- matrixmaker(X)

# prepare dat.gena: exclude non-sensical and co-linear vars, code treatment and outcome, change names
out <- c("black.black","nodegree.nodegree","nodegree.educ","married.married","hispan.hispan","hispan.black",
         "u75.u75","u74.u74","u74.re74","u75.re75","u74.re74","u75.re75","re74.re74","re75.re75","re75.re74")

XX     <- XX[,(colnames(XX) %in% out)==F]
dat.gen    <- data.frame(dat.gen[,(names(dat.gen) %in% covars)==F],XX)
covars <- names(dat.gen)[-which((names(dat.gen) %in% c("treat","re78")))]
X     <- dat.gen[,covars]
X     <- as.matrix(X)
colnames(X) <-gsub(".age","*Age",colnames(X))
colnames(X) <-gsub(".educ*","*Schooling",colnames(X))
colnames(X) <-gsub(".black","*Black",colnames(X))
colnames(X) <-gsub(".hispan","*Hispanic",colnames(X))
colnames(X) <-gsub(".married","*Married",colnames(X))
colnames(X) <-gsub(".re74","*Earnings 1974",colnames(X))
colnames(X) <-gsub(".re75","*Earnings 1975",colnames(X))
colnames(X) <-gsub(".u74","*Unemployed 1974",colnames(X))
colnames(X) <-gsub(".u75","*Unemployed 1975",colnames(X))
colnames(X) <-gsub(".nodegree","*HS Dropout",colnames(X))
colnames(X) <-gsub("age","Age",colnames(X))
colnames(X) <-gsub("educ","Schooling",colnames(X))
colnames(X) <-gsub("black","Black",colnames(X))
colnames(X) <-gsub("hispan","Hispanic",colnames(X))
colnames(X) <-gsub("married","Married",colnames(X))
colnames(X) <-gsub("re74","Earnings 1974",colnames(X))
colnames(X) <-gsub("re75","Earnings 1975",colnames(X))
colnames(X) <-gsub("u74","Unemployed 1974",colnames(X))
colnames(X) <-gsub("u75","Unemployed 1975",colnames(X))
colnames(X) <-gsub("nodegree","HS Dropout",colnames(X))
dat.gen <- data.frame(dat.gen,X)
Y <- dat.gen$re78
W <- dat.gen$tr

# Entropy Balancing
out.eb <- ebalance(
  Treatment = W,
  X=X
)

bout.eb <- MatchBalance(W ~ age + black + married + nodegree + re74 + re75 + educ + hispan, data = dat.gen, 
                        weights=c(dat.gen$treat[dat.gen$treat==1], out.eb$w),ks=FALSE)

dat.gen$w   <-  c(dat.gen$treat[dat.gen$treat==1],out.eb$w)
des1 <- svydesign(id=~1,weights=~w, data= dat.gen)
mod1 <- svyglm(re78 ~ treat, design = des1)
summary(mod1)

#############################################

### ENTROPY BALACING + GENETIC MATCHING #####

#############################################

dat.ebgen <- foo

# unemployed
dat.ebgen$u74 <- as.numeric(dat.ebgen$re74==0)
dat.ebgen$u75 <- as.numeric(dat.ebgen$re75==0)

# covars
covars <- c("age","educ","black","hispan","married","nodegree","re74","re75","u74","u75")

# compute all interactions
X <- as.matrix(dat.ebgen[,covars])
XX <- matrixmaker(X)

# prepare data: exclude non-sensical and co-linear vars, code treatment and outcome, change names
out <- c("black.black","nodegree.nodegree","nodegree.educ","married.married","hispan.hispan","hispan.black",
         "u75.u75","u74.u74","u74.re74","u75.re75","u74.re74","u75.re75","re74.re74","re75.re75","re75.re74")

XX     <- XX[,(colnames(XX) %in% out)==F]
dat.ebgen    <- data.frame(dat.ebgen[,(names(dat.ebgen) %in% covars)==F],XX)
covars <- names(dat.ebgen)[-which((names(dat.ebgen) %in% c("treat","re78")))]
X     <- dat.ebgen[,covars]
X     <- as.matrix(X)
colnames(X) <-gsub(".age","*Age",colnames(X))
colnames(X) <-gsub(".educ*","*Schooling",colnames(X))
colnames(X) <-gsub(".black","*Black",colnames(X))
colnames(X) <-gsub(".hispan","*Hispanic",colnames(X))
colnames(X) <-gsub(".married","*Married",colnames(X))
colnames(X) <-gsub(".re74","*Earnings 1974",colnames(X))
colnames(X) <-gsub(".re75","*Earnings 1975",colnames(X))
colnames(X) <-gsub(".u74","*Unemployed 1974",colnames(X))
colnames(X) <-gsub(".u75","*Unemployed 1975",colnames(X))
colnames(X) <-gsub(".nodegree","*HS Dropout",colnames(X))
colnames(X) <-gsub("age","Age",colnames(X))
colnames(X) <-gsub("educ","Schooling",colnames(X))
colnames(X) <-gsub("black","Black",colnames(X))
colnames(X) <-gsub("hispan","Hispanic",colnames(X))
colnames(X) <-gsub("married","Married",colnames(X))
colnames(X) <-gsub("re74","Earnings 1974",colnames(X))
colnames(X) <-gsub("re75","Earnings 1975",colnames(X))
colnames(X) <-gsub("u74","Unemployed 1974",colnames(X))
colnames(X) <-gsub("u75","Unemployed 1975",colnames(X))
colnames(X) <-gsub("nodegree","HS Dropout",colnames(X))
dat.ebgen <- data.frame(dat.ebgen,X)
Y <- dat.ebgen$re78
W <- dat.ebgen$tr

# Entropy Balancing
out.ebgen <- ebalance(
  Treatment=W,
  X=X
)

bout.eb <- MatchBalance(W~X,weights=c(W[W==1],out.ebgen$w),ks=FALSE)

# Effect Estimates

dat.ebgen$w   <-  c(W[W==1],out.ebgen$w)
des1 <- svydesign(id=~1,weights=~w, data= dat.ebgen)
mod1 <- svyglm(re78 ~ treat, design = des1)
summary(mod1)

### Now, GenMatch()

X <-   as.matrix(cbind( foo$age,
                        foo$educ,
                        foo$black,
                        foo$hispan,
                        foo$married,
                        foo$nodegree,
                        I(foo$re74/1000),
                        I(foo$re75/1000)))

BalanceMat <- as.matrix(cbind(
  foo$age,
  foo$educ,
  foo$black,
  foo$hispan,
  foo$married,
  foo$nodegree,
  I(foo$re74/1000),
  I(foo$re75/1000),
  I((foo$re74/1000)^2),
  I((foo$re75/1000)^2), 
  
  I((foo$age)^2),
  I((foo$educ)^2),
  I((foo$black)^2), 
  I((foo$hispan)^2), 
  I((foo$married)^2), 
  I((foo$nodegree)^2),
  
  I(foo$age*foo$educ),
  I(foo$age*foo$black),
  I(foo$age*foo$hispan),
  I(foo$age*foo$married),    
  I(foo$age*foo$nodegree),
  I(foo$age*(foo$re74/1000)),
  I(foo$age*(foo$re75/1000)),
  
  I(foo$educ*foo$black),
  I(foo$educ*foo$hispan),
  I(foo$educ*foo$married),    
  I(foo$educ*foo$nodegree), 
  I(foo$educ*(foo$re74/1000)),
  I(foo$educ*(foo$re75/1000)),
  
  I(foo$black*foo$hispan),
  I(foo$black*foo$married),
  I(foo$black*foo$nodegree),
  I(foo$black*(foo$re74/1000)),
  I(foo$black*(foo$re75/1000)),
  
  I(foo$hispan*foo$married),
  I(foo$hispan*foo$nodegree),
  I(foo$hispan*(foo$re74/1000)),
  I(foo$hispan*(foo$re75/1000)),
  
  I(foo$married*foo$nodegree),
  I(foo$married*(foo$re74/1000)),
  I(foo$married*(foo$re75/1000)),
  
  I(foo$nodegree*(foo$re74/1000)), 
  I(foo$nodegree*(foo$re75/1000)),                        
  
  I((foo$re74/1000)*(foo$re75/1000))                                                 
) )


### ORTHOGONALIZE TO PROPENSITY SCORE

# From Review of Econ and Statistics
model.A = foo$treat~I(foo$age) + I(foo$age^2) + I(foo$age^3) + I(foo$educ) + I(foo$educ^2) + I(foo$married) + I(foo$nodegree) + I(foo$black) + I(foo$hispan) + I(foo$re74) + I(foo$re75) + I(foo$re74 == 0) + I(foo$re75 == 0) + I(I(foo$re74)*I(foo$educ))

pscores.A <- glm(model.A, family = binomial)

X2 <- X
for (i in 1:ncol(X2)) {
  lm1 = lm(X2[,i]~pscores.A$linear.pred)
  X2[,i] = lm1$residual
}


# VARIABLES WE MATCH ON BEGIN WITH PROPENSITY SCORE
orthoX2.plus.pscore = cbind(pscores.A$linear.pred, X2)

orthoX2.plus.pscore[,1] <- orthoX2.plus.pscore[,1] -
  mean(orthoX2.plus.pscore[,1])


# NORMALIZE ALL COVARS BY STANDARD DEVIATION
for (i in 1:(dim(orthoX2.plus.pscore)[2])) {
  orthoX2.plus.pscore[,i] <-
    orthoX2.plus.pscore[,i]/sqrt(var(orthoX2.plus.pscore[,i]))
}

# SET BOUNDS ON GENMATCH DOMAIN (OPTIONAL)
lower <- c(95, rep(0, dim(orthoX2.plus.pscore)[2] - 1))
upper <- c(105, rep(1000, dim(orthoX2.plus.pscore)[2] - 1))

Domains <- as.matrix(cbind(lower, upper))

set.seed(12345)

cl <- makeCluster(detectCores())

GM.out <- GenMatch(Tr = foo$treat,
                   X = orthoX2.plus.pscore,
                   BalanceMatrix = BalanceMat,
                   weights = c(foo$treat[foo$treat==1],out.ebgen$w), ## Use the weights outputted by ebal
                   pop.size = 100,           
                   max.generations = 10,
                   wait.generations = 5,
                   print.level = 1,
                   cluster = cl)

stopCluster(cl)

mout <- Match(Y = foo$re78,
              Tr = foo$treat,
              weights = c(foo$treat[foo$treat==1],out.ebgen$w),
              X = orthoX2.plus.pscore,
              Weight.matrix = GM.out)

summary(mout)
MatchBalance(treat ~ age + black + married + nodegree + re74 + re75 + educ + hispan, data = foo, match.out = mout) 

# # # Test characteristics of the weight vector

sorted.weights <- sort(out.ebgen$w, decreasing = T)
sum.1perc <- sum(sorted.weights[1:16])

sorted.weights <- sort(out.ebgen$w)
sum <- 0

for (i in 1:length(sorted.weights)) {
  if (sum > sum.1perc) {
    print(i)
    break
  }
  sum <- sum + sorted.weights[i]
}

# # # Trim dataset to the 2 top percentiles

sorted.weights <- sort(out.ebgen$w, decreasing = T)
cutoff2perc <- sorted.weights[320]

# # # GenMatch again...

# Select observations that should be kept
kept <- c(rep(T, sum(foo$treat == 1)), out.ebgen$w >= cutoff2perc)
foo.ebal <- foo[kept, ]

X <-   as.matrix(cbind( foo.ebal$age,
                        foo.ebal$educ,
                        foo.ebal$black,
                        foo.ebal$hispan,
                        foo.ebal$married,
                        foo.ebal$nodegree,
                        I(foo.ebal$re74/1000),
                        I(foo.ebal$re75/1000)))


BalanceMat <- as.matrix(cbind(
  foo.ebal$age,
  foo.ebal$educ,
  foo.ebal$black,
  foo.ebal$hispan,
  foo.ebal$married,
  foo.ebal$nodegree,
  I(foo.ebal$re74/1000),
  I(foo.ebal$re75/1000),
  I((foo.ebal$re74/1000)^2),
  I((foo.ebal$re75/1000)^2), 
  
  I((foo.ebal$age)^2),
  I((foo.ebal$educ)^2),
  I((foo.ebal$black)^2), 
  I((foo.ebal$hispan)^2), 
  I((foo.ebal$married)^2), 
  I((foo.ebal$nodegree)^2),
  
  I(foo.ebal$age*foo.ebal$educ),
  I(foo.ebal$age*foo.ebal$black),
  I(foo.ebal$age*foo.ebal$hispan),
  I(foo.ebal$age*foo.ebal$married),    
  I(foo.ebal$age*foo.ebal$nodegree),
  I(foo.ebal$age*(foo.ebal$re74/1000)),
  I(foo.ebal$age*(foo.ebal$re75/1000)),
  
  I(foo.ebal$educ*foo.ebal$black),
  I(foo.ebal$educ*foo.ebal$hispan),
  I(foo.ebal$educ*foo.ebal$married),    
  I(foo.ebal$educ*foo.ebal$nodegree), 
  I(foo.ebal$educ*(foo.ebal$re74/1000)),
  I(foo.ebal$educ*(foo.ebal$re75/1000)),
  
  I(foo.ebal$black*foo.ebal$hispan),
  I(foo.ebal$black*foo.ebal$married),
  I(foo.ebal$black*foo.ebal$nodegree),
  I(foo.ebal$black*(foo.ebal$re74/1000)),
  I(foo.ebal$black*(foo.ebal$re75/1000)),
  
  I(foo.ebal$hispan*foo.ebal$married),
  I(foo.ebal$hispan*foo.ebal$nodegree),
  I(foo.ebal$hispan*(foo.ebal$re74/1000)),
  I(foo.ebal$hispan*(foo.ebal$re75/1000)),
  
  I(foo.ebal$married*foo.ebal$nodegree),
  I(foo.ebal$married*(foo.ebal$re74/1000)),
  I(foo.ebal$married*(foo.ebal$re75/1000)),
  
  I(foo.ebal$nodegree*(foo.ebal$re74/1000)), 
  I(foo.ebal$nodegree*(foo.ebal$re75/1000)),                        
  
  I((foo.ebal$re74/1000)*(foo.ebal$re75/1000))                                                 
) )


### ORTHOGONALIZE TO PROPENSITY SCORE

# From Review of Econ and Statistics
model.A = foo.ebal$treat~I(foo.ebal$age) + I(foo.ebal$age^2) + I(foo.ebal$age^3) + I(foo.ebal$educ) + I(foo.ebal$educ^2) + I(foo.ebal$married) + I(foo.ebal$nodegree) + I(foo.ebal$black) + I(foo.ebal$hispan) + I(foo.ebal$re74) + I(foo.ebal$re75) + I(foo.ebal$re74 == 0) + I(foo.ebal$re75 == 0) + I(I(foo.ebal$re74)*I(foo.ebal$educ))

pscores.A <- glm(model.A, family = binomial)

X2 <- X
for (i in 1:ncol(X2)) {
  lm1 = lm(X2[,i]~pscores.A$linear.pred)
  X2[,i] = lm1$residual
}


# VARIABLES WE MATCH ON BEGIN WITH PROPENSITY SCORE
orthoX2.plus.pscore = cbind(pscores.A$linear.pred, X2)

orthoX2.plus.pscore[,1] <- orthoX2.plus.pscore[,1] -
  mean(orthoX2.plus.pscore[,1])


# NORMALIZE ALL COVARS BY STANDARD DEVIATION
for (i in 1:(dim(orthoX2.plus.pscore)[2])) {
  orthoX2.plus.pscore[,i] <-
    orthoX2.plus.pscore[,i]/sqrt(var(orthoX2.plus.pscore[,i]))
}

# SET BOUNDS ON GENMATCH DOMAIN (OPTIONAL)
lower <- c(95, rep(0, dim(orthoX2.plus.pscore)[2] - 1))
upper <- c(105, rep(1000, dim(orthoX2.plus.pscore)[2] - 1))

Domains <- as.matrix(cbind(lower, upper))

set.seed(12345)

cl <- makeCluster(detectCores())

GM.out <- GenMatch(Tr = foo.ebal$treat,
                   X = orthoX2.plus.pscore,
                   BalanceMatrix = BalanceMat,
                   pop.size = 1000,           
                   max.generations = 10,
                   wait.generations = 5,
                   print.level = 1,
                   cluster = cl)

stopCluster(cl)

mout <- Match(Y = foo.ebal$re78,
              Tr = foo.ebal$treat,
              X = orthoX2.plus.pscore,
              Weight.matrix = GM.out)

sum(GM.out$value < 0.1)
sum(GM.out$value < 0.2)


summary(mout)
MatchBalance(treat ~ age + black + married + nodegree + re74 + re75 + educ + hispan, data = foo.ebal, match.out = mout) 

