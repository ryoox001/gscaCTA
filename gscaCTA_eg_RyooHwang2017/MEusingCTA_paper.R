##########################################
### Model Evaluation using CTA in GSCA ###
##########################################

### Install and load the package of 'gesca'
#install.packages("gesca")
library(gesca)
#install.packages("gdata")
library(gdata)

### Upload ECLS-K:2011 data
setwd("C:\\Users\\jr3gv\\Documents\\ProfRyoo\\Projects\\ProjectA5_GSCA\\MEusingCTA\\gscaCTA\\gscaCTA_eg_RyooHwang2017")
getwd()

### Social skill 

social <- read.csv("ECLS_CTA_SocialSkill.csv", header=T); dim(social); names(social)
colnames(social) <- c("id", "p4con", "p4int", "p4sad", "p4imp", "t4con", "t4int", "t4extb", "t4intb")
head(social); attach(social); summary(social)

# Reverse code: hihger value --> high social skills
social$p4sad <- 5-social$p4sad
social$p4imp <- 5-social$p4imp
social$t4extb <- 5-social$t4extb
social$t4intb <- 5-social$t4intb

# Likewise Mplus, "gesca" deals missing data internally
social[is.na(social)] <- -9999; head(social); summary(social)

#####
##### Step 1: Fitting three models
#####

# Base model
myModel0 <- "
  # Measurement model
    ss_p =~ p4con + p4int + p4sad + p4imp 
    ss_t =~ t4con + t4int + t4extb + t4intb
  # Structural model
    #ss_p ~ ss_t
"
fit0 <- gesca.run(myModel0, social, nbt=100, moption=3, missingvalue = -9999)
summary(fit0); latentmeasures(fit0); fit0$rho

# Structural model
myModel1 <- "
  # Measurement model
    ss_p =~ p4con + p4int + p4sad + p4imp 
    ss_t =~ t4con + t4int + t4extb + t4intb
  # Structural model
    ss_t ~ ss_p
"
fit1 <- gesca.run(myModel1, social, nbt=100, moption=3, missingvalue = -9999)
summary(fit1); latentmeasures(fit1)

# Higher order model - not-identified at CFA
myModel2 <- "
  # Measurement model
    ss_p =~ p4con + p4int + p4sad + p4imp 
    ss_t =~ t4con + t4int + t4extb + t4intb
  # 2nd order 
    ss =: ss_p + ss_t
  # Structural model
    #ss_t ~ ss_p
"
fit2 <- gesca.run(myModel2, social, nbt=100, moption=3, missingvalue = -9999)
summary(fit2); latentmeasures(fit2)

# four-factor model
myModel3 <- "
  # Measurement model
    ss_p1 =~ p4con + p4int 
    ss_p2 =~ p4sad + p4imp 
    ss_t1 =~ t4con + t4int 
    ss_t2 =~ t4extb + t4intb
"
fit3 <- gesca.run(myModel3, social, nbt=100, moption=3, missingvalue = -9999)
summary(fit3); latentmeasures(fit3); fit3$rho





#####
##### Step 2: Obtaining model-implied covariance matrices
#####

summary(fit0); fit0$corr_corres
summary(fit1); fit1$corr_corres
summary(fit2); fit2$corr_corres
summary(fit3); fit3$corr_corres

fit0$LV_VAR
fit1$LV_VAR
fit2$LV_VAR
fit3$LV_VAR

### construct correlation residual matrices

# In GSCA, indicators and latent variables are standardized

#
# Model 0
#

varcov0 <- matrix(rep(0,64), ncol=8); varcov0
upperTriangle(varcov0, diag=TRUE) <- upperTriangle(fit0$corr_corres, diag=TRUE)
lowerTriangle(varcov0) <- lowerTriangle(t(fit0$corr_corres))
for(i in 1:8){
  varcov0[i,i] <- 1
}

mivc0 <- matrix(rep(0,64), ncol=8)
# Loading
loading0 <- fit0$mC; loading0
# Latent correlation
lcorr0 <- fit0$latentcorr[1,2]; lcorr0

# diagonal entries
for(i in 1:8){
  mivc0[i,i] <- loading0[i]^2+1
}

# off-diag
mivc0[1,2] <- loading0[1]*loading0[2]
mivc0[1,3] <- loading0[1]*loading0[3]
mivc0[1,4] <- loading0[1]*loading0[4]
mivc0[1,5] <- loading0[1]*loading0[5]*lcorr0
mivc0[1,6] <- loading0[1]*loading0[6]*lcorr0
mivc0[1,7] <- loading0[1]*loading0[7]*lcorr0
mivc0[1,8] <- loading0[1]*loading0[8]*lcorr0

mivc0[2,3] <- loading0[2]*loading0[3]
mivc0[2,4] <- loading0[2]*loading0[4]
mivc0[2,5] <- loading0[2]*loading0[5]*lcorr0
mivc0[2,6] <- loading0[2]*loading0[6]*lcorr0
mivc0[2,7] <- loading0[2]*loading0[7]*lcorr0
mivc0[2,8] <- loading0[2]*loading0[8]*lcorr0

mivc0[3,4] <- loading0[3]*loading0[4]
mivc0[3,5] <- loading0[3]*loading0[5]*lcorr0
mivc0[3,6] <- loading0[3]*loading0[6]*lcorr0
mivc0[3,7] <- loading0[3]*loading0[7]*lcorr0
mivc0[3,8] <- loading0[3]*loading0[8]*lcorr0

mivc0[4,5] <- loading0[4]*loading0[5]*lcorr0
mivc0[4,6] <- loading0[4]*loading0[6]*lcorr0
mivc0[4,7] <- loading0[4]*loading0[7]*lcorr0
mivc0[4,8] <- loading0[4]*loading0[8]*lcorr0

mivc0[5,6] <- loading0[5]*loading0[6]
mivc0[5,7] <- loading0[5]*loading0[7]
mivc0[5,8] <- loading0[5]*loading0[8]

mivc0[6,7] <- loading0[6]*loading0[7]
mivc0[6,8] <- loading0[6]*loading0[8]

mivc0[7,8] <- loading0[7]*loading0[8]

lowerTriangle(mivc0) <- lowerTriangle(t(mivc0))
mivc0

write.csv(mivc0, file="mivc0.csv")

#
# Model 1
#

varcov1 <- matrix(rep(0,64), ncol=8); varcov1
upperTriangle(varcov1, diag=TRUE) <- upperTriangle(fit1$corr_corres, diag=TRUE)
lowerTriangle(varcov1) <- lowerTriangle(t(fit1$corr_corres))
for(i in 1:8){
  varcov1[i,i] <- 1
}

mivc1 <- matrix(rep(0,64), ncol=8)
# Loading
loading1 <- fit1$mC; loading1
# Latent correlation & r2
b1 <- fit1$latentcorr[1,2]; b1
zeta1 <- 1-fit1$R2[1,2]; zeta1


# diagonal entries
for(i in 1:4){
  mivc1[i,i] <- loading1[i]^2+1
}
for(i in 5:8){
  mivc1[i,i] <- loading1[i]^2*b1^2+zeta1+1
}

# off-diag
mivc1[1,2] <- loading1[1]*loading1[2]
mivc1[1,3] <- loading1[1]*loading1[3]
mivc1[1,4] <- loading1[1]*loading1[4]
mivc1[1,5] <- loading1[1]*loading1[5]*b1
mivc1[1,6] <- loading1[1]*loading1[6]*b1
mivc1[1,7] <- loading1[1]*loading1[7]*b1
mivc1[1,8] <- loading1[1]*loading1[8]*b1

mivc1[2,3] <- loading1[2]*loading1[3]
mivc1[2,4] <- loading1[2]*loading1[4]
mivc1[2,5] <- loading1[2]*loading1[5]*b1
mivc1[2,6] <- loading1[2]*loading1[6]*b1
mivc1[2,7] <- loading1[2]*loading1[7]*b1
mivc1[2,8] <- loading1[2]*loading1[8]*b1

mivc1[3,4] <- loading1[3]*loading1[4]
mivc1[3,5] <- loading1[3]*loading1[5]*b1
mivc1[3,6] <- loading1[3]*loading1[6]*b1
mivc1[3,7] <- loading1[3]*loading1[7]*b1
mivc1[3,8] <- loading1[3]*loading1[8]*b1

mivc1[4,5] <- loading1[4]*loading1[5]*b1
mivc1[4,6] <- loading1[4]*loading1[6]*b1
mivc1[4,7] <- loading1[4]*loading1[7]*b1
mivc1[4,8] <- loading1[4]*loading1[8]*b1

mivc1[5,6] <- loading1[5]*loading1[6]*(b1^2+zeta1)
mivc1[5,7] <- loading1[5]*loading1[7]*(b1^2+zeta1)
mivc1[5,8] <- loading1[5]*loading1[8]*(b1^2+zeta1)

mivc1[6,7] <- loading1[6]*loading1[7]*(b1^2+zeta1)
mivc1[6,8] <- loading1[6]*loading1[8]*(b1^2+zeta1)

mivc1[7,8] <- loading1[7]*loading1[8]*(b1^2+zeta1)

lowerTriangle(mivc1) <- lowerTriangle(t(mivc1))
mivc1

write.csv(mivc1, file="mivc1.csv")


#
# Model 2
#

varcov2 <- matrix(rep(0,64), ncol=8); varcov2
upperTriangle(varcov2, diag=TRUE) <- upperTriangle(fit2$corr_corres, diag=TRUE)
lowerTriangle(varcov2) <- lowerTriangle(t(fit2$corr_corres))
for(i in 1:8){
  varcov2[i,i] <- 1
}

mivc2 <- matrix(rep(0,64), ncol=8)
# Loading
loading2 <- fit2$mC; loading2
# Latent correlation & r2
b2 <- fit2$latentcorr[1,3]; b2
b3 <- fit2$latentcorr[2,3]; b3
zeta2 <- 1-fit2$R2[1,1]; zeta2
zeta3 <- 1-fit2$R2[1,2]; zeta3

# diagonal entries
for(i in 1:4){
  mivc2[i,i] <- loading2[i]^2*b2^2+loading2[i]^2*zeta2+1
}
for(i in 5:8){
  mivc2[i,i] <- loading2[i]^2*b3^2+loading2[i]^2*zeta3+1
}

# off-diag
mivc2[1,2] <- loading2[1]*loading2[2]*(b2^2+zeta2)
mivc2[1,3] <- loading2[1]*loading2[3]*(b2^2+zeta2)
mivc2[1,4] <- loading2[1]*loading2[4]*(b2^2+zeta2)
mivc2[1,5] <- loading2[1]*loading2[5]*b2*b3
mivc2[1,6] <- loading2[1]*loading2[6]*b2*b3
mivc2[1,7] <- loading2[1]*loading2[7]*b2*b3
mivc2[1,8] <- loading2[1]*loading2[8]*b2*b3

mivc2[2,3] <- loading2[2]*loading2[3]*(b2^2+zeta2)
mivc2[2,4] <- loading2[2]*loading2[4]*(b2^2+zeta2)
mivc2[2,5] <- loading2[2]*loading2[5]*b2*b3
mivc2[2,6] <- loading2[2]*loading2[6]*b2*b3
mivc2[2,7] <- loading2[2]*loading2[7]*b2*b3
mivc2[2,8] <- loading2[2]*loading2[8]*b2*b3

mivc2[3,4] <- loading2[3]*loading2[4]*(b2^2+zeta2)
mivc2[3,5] <- loading2[3]*loading2[5]*b2*b3
mivc2[3,6] <- loading2[3]*loading2[6]*b2*b3
mivc2[3,7] <- loading2[3]*loading2[7]*b2*b3
mivc2[3,8] <- loading2[3]*loading2[8]*b2*b3

mivc2[4,5] <- loading2[4]*loading2[5]*b2*b3
mivc2[4,6] <- loading2[4]*loading2[6]*b2*b3
mivc2[4,7] <- loading2[4]*loading2[7]*b2*b3
mivc2[4,8] <- loading2[4]*loading2[8]*b2*b3

mivc2[5,6] <- loading2[5]*loading2[6]*(b3^2+zeta3)
mivc2[5,7] <- loading2[5]*loading2[7]*(b3^2+zeta3)
mivc2[5,8] <- loading2[5]*loading2[8]*(b3^2+zeta3)

mivc2[6,7] <- loading2[6]*loading2[7]*(b3^2+zeta3)
mivc2[6,8] <- loading2[6]*loading2[8]*(b3^2+zeta3)

mivc2[7,8] <- loading2[7]*loading2[8]*(b3^2+zeta3)

lowerTriangle(mivc2) <- lowerTriangle(t(mivc2))
mivc2

write.csv(mivc2, file="mivc2.csv")

#
# Model 3
#

varcov3 <- matrix(rep(0,64), ncol=8); varcov3
upperTriangle(varcov3, diag=TRUE) <- upperTriangle(fit3$corr_corres, diag=TRUE)
lowerTriangle(varcov3) <- lowerTriangle(t(fit3$corr_corres))
for(i in 1:8){
  varcov3[i,i] <- 1
}

mivc3 <- matrix(rep(0,64), ncol=8)
# Loading
loading3 <- fit3$mC; loading3
# Latent correlation & r2
lcorr12 <- fit3$latentcorr[1,2]; lcorr12
lcorr13 <- fit3$latentcorr[1,3]; lcorr13
lcorr14 <- fit3$latentcorr[1,4]; lcorr14
lcorr23 <- fit3$latentcorr[2,3]; lcorr23
lcorr24 <- fit3$latentcorr[2,4]; lcorr24
lcorr34 <- fit3$latentcorr[3,4]; lcorr34

# diagonal entries
for(i in 1:8){
  mivc3[i,i] <- loading3[i]^2+1
}


# off-diag
mivc3[1,2] <- loading3[1]*loading3[2]
mivc3[1,3] <- loading3[1]*loading3[3]*lcorr12
mivc3[1,4] <- loading3[1]*loading3[4]*lcorr12
mivc3[1,5] <- loading3[1]*loading3[5]*lcorr13
mivc3[1,6] <- loading3[1]*loading3[6]*lcorr13
mivc3[1,7] <- loading3[1]*loading3[7]*lcorr14
mivc3[1,8] <- loading3[1]*loading3[8]*lcorr14

mivc3[2,3] <- loading3[2]*loading3[3]*lcorr12
mivc3[2,4] <- loading3[2]*loading3[4]*lcorr12
mivc3[2,5] <- loading3[2]*loading3[5]*lcorr13
mivc3[2,6] <- loading3[2]*loading3[6]*lcorr13
mivc3[2,7] <- loading3[2]*loading3[7]*lcorr14
mivc3[2,8] <- loading3[2]*loading3[8]*lcorr14

mivc3[3,4] <- loading3[3]*loading3[4]
mivc3[3,5] <- loading3[3]*loading3[5]*lcorr23
mivc3[3,6] <- loading3[3]*loading3[6]*lcorr23
mivc3[3,7] <- loading3[3]*loading3[7]*lcorr24
mivc3[3,8] <- loading3[3]*loading3[8]*lcorr24

mivc3[4,5] <- loading3[4]*loading3[5]*lcorr23
mivc3[4,6] <- loading3[4]*loading3[6]*lcorr23
mivc3[4,7] <- loading3[4]*loading3[7]*lcorr24
mivc3[4,8] <- loading3[4]*loading3[8]*lcorr24

mivc3[5,6] <- loading3[5]*loading3[6]
mivc3[5,7] <- loading3[5]*loading3[7]*lcorr34
mivc3[5,8] <- loading3[5]*loading3[8]*lcorr34

mivc3[6,7] <- loading3[6]*loading3[7]*lcorr34
mivc3[6,8] <- loading3[6]*loading3[8]*lcorr34

mivc3[7,8] <- loading3[7]*loading3[8]

lowerTriangle(mivc3) <- lowerTriangle(t(mivc3))
mivc3

write.csv(mivc3, file="mivc3.csv")