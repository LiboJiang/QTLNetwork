








#load package and function
source("sup.R")
library(splines)
library(orthogonalsplinebasis)
library(MASS)
library(grplasso)
library(parallel)


source("ODE.smooth.R")
source("ODE.varsel.R")
source("ODE.varsel1.R")


#load data
FM <- read.csv("F-marker.csv")
FM2014 <- read.csv("F-2014_sample-sig_snp.csv",header = F)[,-1]
FM201401 <- FM[as.character(unlist(c(FM2014[1,]))),]

FM2015 <- read.csv("F-2015_sample-sig_snp.csv",header = F)[,-1]
FM201501 <- FM[as.character(unlist(c(FM2015[1,]))),]

FM02 <- rbind(FM201401,FM201501)


LM <- read.csv("L-marker-1.csv")
L01 <- read.csv("L-2015_sample-sig_snp.csv",header = F)[,-1]
L02 <- LM[as.character(unlist(c(L01[1,]))),]

YM <- read.csv("Y-marker.csv")
Y01 <- read.csv("Y-2015_sample-sig_snp.csv",header = F)[,-1]
Y02 <- YM[as.character(unlist(c(Y01[1,]))),]


Fp1 <- read.csv("F_2014_-HT-DIA.csv",header = F)
Fp2 <- read.csv("F_2015-HT-DIA.csv",header = F)

Fpheno <- dat_F (marker=FM02,pheno=Fp1,pi1=1:8,pi2=9:16)
Fpheno_e <- varf(p1=Fpheno$H,p2=Fpheno$D,marker=Fpheno$Marker)

Fpheno2015 <- dat_F (marker=FM02,pheno=Fp2,pi1=1:11,pi2=12:22)
Fpheno2015_e <- varf(p1=Fpheno2015$H,p2=Fpheno2015$D,marker=Fpheno2015$Marker)

Lp1 <- read.csv("L-HT-DIA-2015.csv",header = F)
Lpheno <- dat_F (marker=L02,pheno=Lp1,pi1=1:11,pi2=12:22)
Lpheno_e <- varf(p1=Lpheno$H,p2=Lpheno$D,marker=Lpheno$Marker)

Yp1 <- read.csv("HD_Y-2015_-HT-DIA.csv",header = F)
Ypheno <- dat_F (marker=Y02,pheno=Yp1,pi1=1:11,pi2=12:22)
Ypheno_e <- varf(p1=Ypheno$H,p2=Ypheno$D,marker=Ypheno$Marker)


#F2014_H
stage2_F2014_H <- smooth.optim(times=1:8,para=rep(.1,6),y=t(Fpheno_e$eh),nt=seq(1,8,length=30))
stage3_F2014_H <- varsel(X=t(stage2_F2014_H$smooth.d),Y=t(stage2_F2014_H$dsmooth.d),tt=seq(1,8,length=120))
F2014_H.odee <- optim.parallel(connect=stage3_F2014_H$connect,effect=t(stage2_F2014_H$smooth.d),
                               n.cores=4,proc=ode.optim,order=6,times=seq(1,8,length=30),nstep=29)
F2014_H.res <- interType(con=stage3_F2014_H$connect,alle=F2014_H.odee,sme=stage2_F2014_H$smooth.d)
aaaF2014_H <- regasso(connect1=stage3_F2014_H$connect,gene=stage2_F2014_H$smooth.d,interaction=F2014_H.odee)


#F2014_D
stage2_F2014_D <- smooth.optim(times=1:8,para=rep(.1,6),y=t(Fpheno_e$ed)*100,nt=seq(1,8,length=30))
stage3_F2014_D <- varsel(X=t(stage2_F2014_D$smooth.d),Y=t(stage2_F2014_D$dsmooth.d),tt=seq(1,8,length=120))
F2014_D.odee <- optim.parallel(connect=stage3_F2014_D$connect,effect=t(stage2_F2014_D$smooth.d),
                               n.cores=4,proc=ode.optim,order=6,times=seq(1,8,length=30),nstep=29)
F2014_D.res <- interType(con=stage3_F2014_D$connect,alle=F2014_D.odee,sme=stage2_F2014_D$smooth.d)
aaaF2014_D <- regasso(connect1=stage3_F2014_D$connect,gene=stage2_F2014_D$smooth.d,interaction=F2014_D.odee)


#F2015_H
stage2_F2015_H <- smooth.optim(times=1:11,para=rep(.1,6),y=t(Fpheno2015_e$eh),nt=seq(1,11,length=30))
stage3_F2015_H <- varsel(X=t(stage2_F2015_H$smooth.d),Y=t(stage2_F2015_H$dsmooth.d),tt=seq(1,11,length=120))
F2015_H.odee <- optim.parallel(connect=stage3_F2015_H$connect,effect=t(stage2_F2015_H$smooth.d),
                               n.cores=1,proc=ode.optim,order=6,times=seq(1,11,length=30),nstep=29)
F2015_H.res <- interType(con=stage3_F2015_H$connect,alle=F2015_H.odee,sme=stage2_F2015_H$smooth.d)
aaaF2015_H <- regasso(connect1=stage3_F2015_H$connect,gene=stage2_F2015_H$smooth.d,interaction=F2015_H.odee)
#F2015_D
stage2_F2015_D <- smooth.optim(times=1:11,para=rep(.1,6),y=t(Fpheno2015_e$ed)*100,nt=seq(1,11,length=30))
stage3_F2015_D <- varsel(X=t(stage2_F2015_D$smooth.d),Y=t(stage2_F2015_D$dsmooth.d),tt=seq(1,11,length=120))
F2015_D.odee <- optim.parallel(connect=stage3_F2015_D$connect,effect=t(stage2_F2015_D$smooth.d),
                               n.cores=4,proc=ode.optim,order=6,times=seq(1,11,length=30),nstep=29)
F2015_D.res <- interType(con=stage3_F2015_D$connect,alle=F2015_D.odee,sme=stage2_F2015_D$smooth.d)

aaaF2015_D <- regasso(connect1=stage3_F2015_D$connect,gene=stage2_F2015_D$smooth.d,interaction=F2015_D.odee)

#L_H
stage2_L_H <- smooth.optim(times=1:11,para=rep(.1,6),y=t(Lpheno_e$eh),nt=seq(1,11,length=30))
stage3_L_H <- varsel(X=t(stage2_L_H$smooth.d),Y=t(stage2_L_H$dsmooth.d),tt=seq(1,11,length=120))
L_H.odee <- optim.parallel(connect=stage3_L_H$connect,effect=t(stage2_L_H$smooth.d),
                           n.cores=4,proc=ode.optim,order=6,times=seq(1,11,length=30),nstep=29)
L_H.res <- interType(con=stage3_L_H$connect,alle=L_H.odee,sme=stage2_L_H$smooth.d)
aaaL_H <- regasso(connect1=stage3_L_H$connect,gene=stage2_L_H$smooth.d,interaction=L_H.odee)

#L_D
stage2_L_D <- smooth.optim(times=1:11,para=rep(.1,6),y=t(Lpheno_e$ed)*100,nt=seq(1,11,length=30))
stage3_L_D <- varsel(X=t(stage2_L_D$smooth.d),Y=t(stage2_L_D$dsmooth.d),tt=seq(1,11,length=120))
L_D.odee <- optim.parallel(connect=stage3_L_D$connect,effect=t(stage2_L_D$smooth.d),
                           n.cores=4,proc=ode.optim,order=6,times=seq(1,11,length=30),nstep=29)
L_D.res <- interType(con=stage3_L_D$connect,alle=L_D.odee,sme=stage2_L_D$smooth.d)
aaaL_D <- regasso(connect1=stage3_L_D$connect,gene=stage2_L_D$smooth.d,interaction=L_D.odee)

#Y_H
stage2_Y_H <- smooth.optim(times=1:11,para=rep(.1,6),y=t(Ypheno_e$eh),nt=seq(1,11,length=30))
stage3_Y_H <- varsel(X=t(stage2_Y_H$smooth.d),Y=t(stage2_Y_H$dsmooth.d),tt=seq(1,11,length=120))
Y_H.odee <- optim.parallel(connect=stage3_Y_H$connect,effect=t(stage2_Y_H$smooth.d),
                           n.cores=4,proc=ode.optim,order=6,times=seq(1,11,length=30),nstep=29)
Y_H.res <- interType(con=stage3_Y_H$connect,alle=Y_H.odee,sme=stage2_Y_H$smooth.d)
aaaY_H <- regasso(connect1=stage3_Y_H$connect,gene=stage2_Y_H$smooth.d,interaction=Y_H.odee)

#Y_D
stage2_Y_D <- smooth.optim(times=1:11,para=rep(.1,6),y=t(Ypheno_e$ed)*100,nt=seq(1,11,length=30))
stage3_Y_D <- varsel(X=t(stage2_Y_D$smooth.d),Y=t(stage2_Y_D$dsmooth.d),tt=seq(1,11,length=120))
Y_D.odee <- optim.parallel(connect=stage3_Y_D$connect,effect=t(stage2_Y_D$smooth.d),
                           n.cores=4,proc=ode.optim,order=6,times=seq(1,11,length=30),nstep=29)
Y_D.res <- interType(con=stage3_Y_D$connect,alle=Y_D.odee,sme=stage2_Y_D$smooth.d)
aaaY_D <- regasso(connect1=stage3_Y_D$connect,gene=stage2_Y_D$smooth.d,interaction=Y_D.odee)


