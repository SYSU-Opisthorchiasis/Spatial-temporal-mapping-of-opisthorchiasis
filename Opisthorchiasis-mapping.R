library(foreign)
library(raster)
library(xts)
library(plotKML)
library(maptools)
library(spatstat)
library(rgdal)
library(sp)
library(rgeos)
library(readxl)

####1 read data####
setwd('F:\\data_path\\')#data path

ADM1<-read.dta("data_adm1.dta")
ADM1<-as.data.frame(ADM1)
dataset_p<-ADM1
dataset_p$year_ORG<-round((dataset_p$YearStart+dataset_p$YearFinish)/2)#survey year
dataset_p$Diag1<-NA
dataset_p$Diag2<-NA
dataset_p[which(dataset_p$Diagnostic1==2),"Diag1"]<-1
dataset_p[which(is.na(dataset_p$Diag1)==TRUE),"Diag1"]<-0
dataset_p[which(dataset_p$Diagnostic1==3),"Diag2"]<-1
dataset_p[which(is.na(dataset_p$Diag2)==TRUE),"Diag2"]<-0#dummy variables for diagnostic methods (Diag1:FECT;Diag2:Others;Reference group:Kato-Katz)
ADM1<-dataset_p

ADM2<-read.dta("data_adm2.dta")
ADM2<-as.data.frame(ADM2)
dataset_p<-ADM2
dataset_p$year_ORG<-round((dataset_p$YearStart+dataset_p$YearFinish)/2)
dataset_p$Diag1<-NA
dataset_p$Diag2<-NA
dataset_p[which(dataset_p$Diagnostic1==2),"Diag1"]<-1
dataset_p[which(is.na(dataset_p$Diag1)==TRUE),"Diag1"]<-0
dataset_p[which(dataset_p$Diagnostic1==3),"Diag2"]<-1
dataset_p[which(is.na(dataset_p$Diag2)==TRUE),"Diag2"]<-0
ADM2<-dataset_p

ADM3<-read.dta("data_adm3.dta")
ADM3<-as.data.frame(ADM3)
dataset_p<-ADM3
dataset_p$year_ORG<-round((dataset_p$YearStart+dataset_p$YearFinish)/2)
dataset_p$Diag1<-NA
dataset_p$Diag2<-NA
dataset_p[which(dataset_p$Diagnostic1==2),"Diag1"]<-1
dataset_p[which(is.na(dataset_p$Diag1)==TRUE),"Diag1"]<-0
dataset_p[which(dataset_p$Diagnostic1==3),"Diag2"]<-1
dataset_p[which(is.na(dataset_p$Diag2)==TRUE),"Diag2"]<-0
ADM3<-dataset_p

point<-read.dta("data_point.dta")
point<-as.data.frame(point)
dataset_p<-point
dataset_p$year_ORG<-round((dataset_p$YearStart+dataset_p$YearFinish)/2)
dataset_p$Diag1<-NA
dataset_p$Diag2<-NA
dataset_p[which(dataset_p$Diagnostic1==2),"Diag1"]<-1
dataset_p[which(is.na(dataset_p$Diag1)==TRUE),"Diag1"]<-0
dataset_p[which(dataset_p$Diagnostic1==3),"Diag2"]<-1
dataset_p[which(is.na(dataset_p$Diag2)==TRUE),"Diag2"]<-0
point<-dataset_p

data.whole<-read.dta("data_whole.dta")#prediction dataset
data.whole<-as.data.frame(data.whole)


####2 model setting####
library(INLA)
library(maps)
library(ggplot2)
library(MASS)
library(MBA)
library(fields)
library(combinat)
library(spdep)
library(gtools)
library(matrixStats)

#2.1 setting####
data_p<-point
data_p$Intercept=1
data_a1<-ADM1
data_a1$Intercept=1
data_a2<-ADM2
data_a2$Intercept=1
data_a3<-ADM3
data_a3$Intercept=1
data.whole<-data.whole
data.whole$Intercept=1

#2.2 mesh####
a<-shapefile("southeast_adm0.shp")#ADM0 shapefile of southeast asia
plot(a)
IDs <- a$ID_0
IDs[IDs %in% c(min(IDs):max(IDs))] <- "1"
IDs[IDs %in% c(max(IDs):min(IDs))] <- "1"
outer_border <- unionSpatialPolygons(a, IDs)
plot(outer_border)
boundaries<-fortify(outer_border)
boundaries<-boundaries[,1:2]
colnames(boundaries)<-c("x","y")
plot(boundaries)

#coordinates for making mesh
data_a3_thai<-data_a3[which(data_a3$ID0==228),]
data.whole_thai<-data.whole[which(data.whole$ID0==228),]
data_add3_samp<-NULL
for(i in 1:length(table(data_a3_thai$ID3))){
  print(i)
  samp0<-data.whole_thai[which(data.whole_thai$ID3==names(table(data_a3_thai$ID3))[i]),]
  if(nrow(samp0)<=5){
    data.whole_samp<-samp0
  }else{data.whole_samp<-samp0[sample(c(1:nrow(samp0)),5),]}
  data_add3_samp<-rbind(data_add3_samp,data.whole_samp)
}

data_a2_thai<-data_a2[which(data_a2$ID0==228),]
data.whole_thai<-data.whole[which(data.whole$ID0==228),]
data_add2_samp_thai<-NULL
for(i in 1:length(table(data_a2_thai$ID2))){
  print(i)
  samp0<-data.whole_thai[which(data.whole_thai$ID2==names(table(data_a2_thai$ID2))[i]),]
  if(nrow(samp0)<=5){
    data.whole_samp<-samp0
  }else{data.whole_samp<-samp0[sample(c(1:nrow(samp0)),5),]}
  data_add2_samp_thai<-rbind(data_add2_samp_thai,data.whole_samp)
}

data_a2_laos<-data_a2[which(data_a2$ID0==123),]
data.whole_laos<-data.whole[which(data.whole$ID0==123),]
data_add2_samp_laos<-NULL
for(i in 1:length(table(data_a2_laos$ID2))){
  print(i)
  samp0<-data.whole_laos[which(data.whole_laos$ID2==names(table(data_a2_laos$ID2))[i]),]
  if(nrow(samp0)<=5){
    data.whole_samp<-samp0
  }else{data.whole_samp<-samp0[sample(c(1:nrow(samp0)),5),]}
  data_add2_samp_laos<-rbind(data_add2_samp_laos,data.whole_samp)
}

data_a2_cam<-data_a2[which(data_a2$ID0==40),]
data.whole_cam<-data.whole[which(data.whole$ID0==40),]
data_add2_samp_cam<-NULL
for(i in 1:length(table(data_a2_cam$ID2))){
  print(i)
  samp0<-data.whole_cam[which(data.whole_cam$ID2==names(table(data_a2_cam$ID2))[i]),]
  if(nrow(samp0)<=5){
    data.whole_samp<-samp0
  }else{data.whole_samp<-samp0[sample(c(1:nrow(samp0)),5),]}
  data_add2_samp_cam<-rbind(data_add2_samp_cam,data.whole_samp)
}

data_a2_viet<-data_a2[which(data_a2$ID0==250),]
data.whole_viet<-data.whole[which(data.whole$ID0==250),]
data_add2_samp_viet<-NULL
for(i in 1:length(table(data_a2_viet$ID2))){
  print(i)
  samp0<-data.whole_viet[which(data.whole_viet$ID2==names(table(data_a2_viet$ID2))[i]),]
  if(nrow(samp0)<=5){
    data.whole_samp<-samp0
  }else{data.whole_samp<-samp0[sample(c(1:nrow(samp0)),5),]}
  data_add2_samp_viet<-rbind(data_add2_samp_viet,data.whole_samp)
}

data_add2_samp<-rbind(data_add2_samp_thai,data_add2_samp_laos,data_add2_samp_cam,data_add2_samp_viet)

loc1<-data_a1[,c("Longitude","Latitude")]
loc2<-data_add2_samp[,c("Longitude","Latitude")]
loc3<-data_add3_samp[,c("Longitude","Latitude")]
loc<-rbind(loc1,loc2,loc3)
loc_final<-unique(loc)

bnd <- inla.nonconvex.hull(coords.pred,convex=-0.06)
mesh <- inla.mesh.2d(loc=loc, boundary = bnd, offset=c(0.1, 0.4),
                     max.edge=c(0.35,0.5))
plot(mesh)
points(boundaries,lwd=3)
points(as.matrix(coords),pch=20,cex=1,col=2)


#2.3 SPDE model####
nu <- 1 #Matern smoothness parameter
set.seed(1)
dist<-dist(as.matrix(coords.pred), method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
ran.pr <- median(dist)
kap.pr <- sqrt(8*nu)/ran.pr
spde <- inla.spde2.matern(mesh=mesh, alpha=2,
                          B.tau=matrix(c(0, 1, 0),nrow=1,ncol=3),
                          B.kappa=matrix(c(0, 0, 1),nrow=1,ncol=3),
                          theta.prior.mean=c(0, log(kap.pr)),
                          theta.prior.prec=c(1, 1)
)


#2.4 Time Knots####
gtime<-seq(1978,2018,by=5)
mesh.t<-inla.mesh.1d(gtime,degree = 2)
idx.st<-inla.spde.make.index('i',n.spde = spde$n.spde, n.group = mesh.t$m)


#2.5 Point observations####
data_p_bino<-data_p[which(data_p$type=="binominal"),] #For data have both NumPositiv and NumExamine, use binominal-likelihood.
coords_p_bino <- cbind(data_p_bino$Longitude,data_p_bino$Latitude)
colnames(coords_p_bino)<-c("x","y")

data_p_beta<-data_p[which(data_p$type=="beta"),] #For data only have prevalence, use beta-likelihood.
data_p_beta$pos<-data_p_beta$NumPositiv/data_p_beta$NumExamine #prevalence
data_p_beta[which(data_p_beta$pos==0),"pos"]<-0.00000000001 #Since prevalence=0 is invalid in beta-likelihood, the data with a prevalence of 0 is replaced with a very small value.
coords_p_beta <- cbind(data_p_beta$Longitude,data_p_beta$Latitude)
colnames(coords_p_beta)<-c("x","y")

Ap_bino = inla.spde.make.A(mesh, loc=as.matrix(coords_p_bino),group.mesh=mesh.t,group=data_p_bino$year_ORG)
stack.p_bino = inla.stack(tag="point_bino",
                          data = list(y=cbind(data_p_bino$NumPositiv,NA),Ntrials=data_p_bino$NumExamine,link="logit"),
                          A = list(Ap_bino, 1),
                          effects = list(c(idx.st,list(Intercept=1)), list(
                            data.frame(SurveyType=data_p_bino$SurveyType,
                                       slst_night11=data_p_bino$slst_night11, slst_night12=data_p_bino$slst_night12,
                                       shii=data_p_bino$shii, 
                                       swb_new=data_p_bino$swb_new,
                                       selevation=data_p_bino$selevation,
                                       scity=data_p_bino$scity,
                                       Diag1=data_p_bino$Diag1,Diag2=data_p_bino$Diag2))))

Ap_beta = inla.spde.make.A(mesh, loc=as.matrix(coords_p_beta),group.mesh=mesh.t,group=data_p_beta$year_ORG)
stack.p_beta = inla.stack(tag="point_beta",
                          data = list(y=cbind(NA,data_p_beta$pos),link="logit"),
                          A = list(Ap_beta, 1),
                          effects = list(c(idx.st,list(Intercept=1)), list(
                            data.frame(SurveyType=data_p_beta$SurveyType,
                                       slst_night11=data_p_beta$slst_night11, slst_night12=data_p_beta$slst_night12,
                                       shii=data_p_beta$shii, 
                                       swb_new=data_p_beta$swb_new,
                                       selevation=data_p_beta$selevation,
                                       scity=data_p_beta$scity,
                                       Diag1=data_p_beta$Diag1,Diag2=data_p_beta$Diag2))))


#2.6 ADM1 observations####
area.poly1 <- shapefile("southeast_adm1.shp")#ADM1 shapefile of southeast asia
mesh.coord.in1 <- NULL
for(i in 1:length(area.poly1)){
  mesh.coord1 <- SpatialPoints(mesh$loc,proj4string=CRS(proj4string(area.poly1)))
  x <- as.vector(which(!is.na(over(mesh.coord1, area.poly1[i,]))))
  x <- x[which(x<=length(mesh.coord1))]
  mesh.coord.in1 <- rbind(mesh.coord.in1 , mesh$loc[x,])
}

mesh.coord.in1 <- mesh.coord.in1[,1:2]
block1 <- rep(0, nrow(mesh.coord.in1))
print(length(block1))
for(i in 1:length(area.poly1)){
  locin1<-SpatialPoints(mesh.coord.in1,proj4string=CRS(proj4string(area.poly1)))
  x<-as.vector(which(!is.na((over(locin1,area.poly1[i,])))))
  if(length(x)!=0){
    for(j in 1:length(x)){
      if(x[j]<=nrow(mesh.coord.in1)){
        block1[x[j]] <- i
      }
    }
  }
}

data_a1<-ADM1
data_a1_bino<-data_a1[which(data_a1$type=="binominal"),]
data_a1_bino$surveyid<-seq(1,nrow(data_a1_bino))
data_a1_bino<-data_a1_bino[with(data_a1_bino,order(data_a1_bino$ID,data_a1_bino$surveyid)),]
data_a1_bino$bs_id<-cumsum(!duplicated(data_a1_bino[,c("ID","surveyid")]))

locint1<-cbind(data.frame(locin1),block1)
colnames(locint1)<-c("long_loc","lat_loc","ID")
merg1_bino<-merge(locint1,data_a1_bino)
merg1_bino<-merg1_bino[with(merg1_bino,order(merg1_bino$bs_id)),]

Aa_bino <- inla.spde.make.A(mesh=mesh,loc=as.matrix(data.frame(merg1_bino[,c("long_loc","lat_loc")])),
                            group.mesh=mesh.t,group=merg1_bino$year_ORG,
                            block=merg1_bino$bs_id, block.rescale="sum")

stack.a1_bino <- inla.stack(tag='areal1_bino',
                            data = list(y=cbind(data_a1_bino$NumPositiv,NA),Ntrials=data_a1_bino$NumExamine,link="logit"),
                            A=list(Aa_bino, 1),
                            effects = list(c(idx.st,list(Intercept=1)),
                                           list(data.frame(SurveyType=data_a1_bino$SurveyType,
                                                           slst_night11=data_a1_bino$slst_night11, slst_night12=data_a1_bino$slst_night12,
                                                           shii=data_a1_bino$shii, 
                                                           swb_new=data_a1_bino$swb_new,
                                                           selevation=data_a1_bino$selevation,
                                                           scity=data_a1_bino$scity,
                                                           Diag1=data_a1_bino$Diag1,Diag2=data_a1_bino$Diag2))))

data_a1_beta<-data_a1[which(data_a1$type=="bino"),]
data_a1_beta$pos<-data_a1_beta$NumPositiv/data_a1_beta$NumExamine
data_a1_beta[which(data_a1_beta$pos==0),"pos"]<-0.00000000001

data_a1_beta$surveyid<-seq(1,nrow(data_a1_beta))
data_a1_beta<-data_a1_beta[with(data_a1_beta,order(data_a1_beta$ID,data_a1_beta$surveyid)),]
data_a1_beta$bs_id<-cumsum(!duplicated(data_a1_beta[,c("ID","surveyid")]))

locint1<-cbind(data.frame(locin1),block1)
colnames(locint1)<-c("long_loc","lat_loc","ID")
merg1_beta<-merge(locint1,data_a1_beta)
merg1_beta<-merg1_beta[with(merg1_beta,order(merg1_beta$bs_id)),]

Aa_beta <- inla.spde.make.A(mesh=mesh,loc=as.matrix(data.frame(merg1_beta[,c("long_loc","lat_loc")])),
                            group.mesh=mesh.t,group=merg1_beta$year_ORG,
                            block=merg1_beta$bs_id, block.rescale="sum")

stack.a1_beta <- inla.stack(tag='areal1_beta',
                            data = list(y=cbind(NA,data_a1_beta$pos),link="logit"),
                            A=list(Aa_beta, 1),
                            effects = list(c(idx.st,list(Intercept=1)),
                                           list(data.frame(SurveyType=data_a1_beta$SurveyType,
                                                           slst_night11=data_a1_beta$slst_night11, slst_night12=data_a1_beta$slst_night12,
                                                           shii=data_a1_beta$shii, 
                                                           swb_new=data_a1_beta$swb_new,
                                                           selevation=data_a1_beta$selevation,
                                                           scity=data_a1_beta$scity,
                                                           Diag1=data_a1_beta$Diag1,Diag2=data_a1_beta$Diag2))))


#2.7 ADM2 observations####
area.poly2 <- shapefile("southeast_adm2.shp")#ADM2 shapefile of southeast asia
mesh.coord.in2 <- NULL
for(i in 1:length(area.poly2)){
  mesh.coord2 <- SpatialPoints(mesh$loc,proj4string=CRS(proj4string(area.poly2)))
  x <- as.vector(which(!is.na(over(mesh.coord2, area.poly2[i,]))))
  x <- x[which(x<=length(mesh.coord2))]
  mesh.coord.in2 <- rbind(mesh.coord.in2 , mesh$loc[x,])
}

mesh.coord.in2 <- mesh.coord.in2[,1:2]
block2 <- rep(0, nrow(mesh.coord.in2))
print(length(block2))
for(i in 1:length(area.poly2)){
  locin2<-SpatialPoints(mesh.coord.in2,proj4string=CRS(proj4string(area.poly2)))
  x<-as.vector(which(!is.na((over(locin2,area.poly2[i,])))))
  if(length(x)!=0){
    for(j in 1:length(x)){
      if(x[j]<=nrow(mesh.coord.in2)){
        block2[x[j]] <- i
      }
    }
  }
}

data_a2<-ADM2
data_a2_bino<-data_a2[which(data_a2$type=="binominal"),]
data_a2_bino$surveyid<-seq(1,nrow(data_a2_bino))
data_a2_bino<-data_a2_bino[with(data_a2_bino,order(data_a2_bino$ID,data_a2_bino$surveyid)),]
data_a2_bino$bs_id<-cumsum(!duplicated(data_a2_bino[,c("ID","surveyid")]))

locint2<-cbind(data.frame(locin2),block2)
colnames(locint2)<-c("long_loc","lat_loc","ID")
merg2_bino<-merge(locint2,data_a2_bino)
merg2_bino<-merg2_bino[with(merg2_bino,order(merg2_bino$bs_id)),]

Aa2_bino <- inla.spde.make.A(mesh=mesh,loc=as.matrix(data.frame(merg2_bino[,c("long_loc","lat_loc")])),
                             group.mesh=mesh.t,group=merg2_bino$year_ORG,
                             block=merg2_bino$bs_id, block.rescale="sum")

stack.a2_bino <- inla.stack(tag='areal2_bino',
                            data = list(y=cbind(data_a2_bino$NumPositiv,NA),Ntrials=data_a2_bino$NumExamine,link="logit"),
                            A=list(Aa2_bino, 1),
                            effects = list(c(idx.st,list(Intercept=1)),
                                           list(data.frame(SurveyType=data_a2_bino$SurveyType,
                                                           slst_night11=data_a2_bino$slst_night11, slst_night12=data_a2_bino$slst_night12,
                                                           shii=data_a2_bino$shii, 
                                                           swb_new=data_a2_bino$swb_new,
                                                           selevation=data_a2_bino$selevation,
                                                           scity=data_a2_bino$scity,
                                                           Diag1=data_a2_bino$Diag1,Diag2=data_a2_bino$Diag2))))


data_a2_beta<-data_a2[which(data_a2$type=="beta"),]
data_a2_beta$pos<-data_a2_beta$NumPositiv/data_a2_beta$NumExamine
data_a2_beta[which(data_a2_beta$pos==0),"pos"]<-0.00000000001

data_a2_beta$surveyid<-seq(1,nrow(data_a2_beta))
data_a2_beta<-data_a2_beta[with(data_a2_beta,order(data_a2_beta$ID,data_a2_beta$surveyid)),]
data_a2_beta$bs_id<-cumsum(!duplicated(data_a2_beta[,c("ID","surveyid")]))

locint2<-cbind(data.frame(locin2),block2)
colnames(locint2)<-c("long_loc","lat_loc","ID")
merg2_beta<-merge(locint2,data_a2_beta)
merg2_beta<-merg2_beta[with(merg2_beta,order(merg2_beta$bs_id)),]

Aa2_beta <- inla.spde.make.A(mesh=mesh,loc=as.matrix(data.frame(merg2_beta[,c("long_loc","lat_loc")])),
                             group.mesh=mesh.t,group=merg2_beta$year_ORG,
                             block=merg2_beta$bs_id, block.rescale="sum")

stack.a2_beta <- inla.stack(tag='areal2_beta',
                            data = list(y=cbind(NA,data_a2_beta$pos),link="logit"),
                            A=list(Aa2_beta, 1),
                            effects = list(c(idx.st,list(Intercept=1)),
                                           list(data.frame(SurveyType=data_a2_beta$SurveyType,
                                                           slst_night11=data_a2_beta$slst_night11, slst_night12=data_a2_beta$slst_night12,
                                                           shii=data_a2_beta$shii, 
                                                           swb_new=data_a2_beta$swb_new,
                                                           selevation=data_a2_beta$selevation,
                                                           scity=data_a2_beta$scity,
                                                           Diag1=data_a2_beta$Diag1,Diag2=data_a2_beta$Diag2))))


#2.8 ADM3 observations####
area.poly3 <- shapefile("southeast_adm3.shp")#ADM3 shapefile of southeast asia
mesh.coord.in3 <- NULL
for(i in 1:length(area.poly3)){
  mesh.coord3 <- SpatialPoints(mesh$loc,proj4string=CRS(proj4string(area.poly3)))
  x <- as.vector(which(!is.na(over(mesh.coord3, area.poly3[i,]))))
  x <- x[which(x<=length(mesh.coord3))]
  mesh.coord.in3 <- rbind(mesh.coord.in3 , mesh$loc[x,])
}

mesh.coord.in3 <- mesh.coord.in3[,1:2]
block3 <- rep(0, nrow(mesh.coord.in3))
print(length(block3))
for(i in 1:length(area.poly3)){
  locin3<-SpatialPoints(mesh.coord.in3,proj4string=CRS(proj4string(area.poly3)))
  x<-as.vector(which(!is.na((over(locin3,area.poly3[i,])))))
  if(length(x)!=0){
    for(j in 1:length(x)){
      if(x[j]<=nrow(mesh.coord.in3)){
        block3[x[j]] <- i
      }
    }
  }
}

data_a3<-ADM3
data_a3_bino<-data_a3 #All the data in ADM3 have both NumPositiv and NumExamine. So all the data in ADM3 use binominal-likelihood.
data_a3_bino$surveyid<-seq(1,nrow(data_a3_bino))
data_a3_bino<-data_a3_bino[with(data_a3_bino,order(data_a3_bino$ID,data_a3_bino$surveyid)),]
data_a3_bino$bs_id<-cumsum(!duplicated(data_a3_bino[,c("ID","surveyid")]))

locint3<-cbind(data.frame(locin3),block3)
colnames(locint3)<-c("long_loc","lat_loc","ID")
merg3_bino<-merge(locint3,data_a3_bino)
merg3_bino<-merg3_bino[with(merg3_bino,order(merg3_bino$bs_id)),]

Aa3_bino <- inla.spde.make.A(mesh=mesh,loc=as.matrix(data.frame(merg3_bino[,c("long_loc","lat_loc")])),
                             group.mesh=mesh.t,group=merg3_bino$year_ORG,
                             block=merg3_bino$bs_id, block.rescale="sum")

stack.a3_bino <- inla.stack(tag='areal3_bino',
                            data = list(y=cbind(data_a3_bino$NumPositiv,NA),Ntrials=data_a3_bino$NumExamine,link="logit"),
                            A=list(Aa3_bino, 1),
                            effects = list(c(idx.st,list(Intercept=1)),
                                           list(data.frame(SurveyType=data_a3_bino$SurveyType,
                                                           slst_night11=data_a3_bino$slst_night11, slst_night12=data_a3_bino$slst_night12,
                                                           shii=data_a3_bino$shii, 
                                                           swb_new=data_a3_bino$swb_new,
                                                           selevation=data_a3_bino$selevation,
                                                           scity=data_a3_bino$scity,
                                                           Diag1=data_a3_bino$Diag1,Diag2=data_a3_bino$Diag2))))



####3 model fitting####
for1<-c("-1","Intercept","f(i, model=spde,group=i.group, control.group=list(model='ar1'))")
para2<-c("SurveyType","slst_day11+slst_day12","shii","swb_new","selevation","scity","Diag1+Diag2")

formula <- as.formula(paste("y ~ ", paste(paste(for1, collapse= "+"),paste(para2, collapse= "+"),sep="+")))
stack.bino <- inla.stack(stack.p_bino,stack.a1_bino,stack.a2_bino,stack.a3_bino)
stack.beta <- inla.stack(stack.p_beta,stack.a1_beta,stack.a2_beta)
stack.full <- inla.stack(stack.bino,stack.beta)

link<-rep(NA,length(stack.full$data$data$y.2))
link[which(is.na(link))]<-1

result<- inla(formula, data=inla.stack.data(stack.full), family=c("binomial","beta"),
              Ntrials=inla.stack.data(stack.full)$Ntrials,
              control.predictor=list(compute=TRUE, A=inla.stack.A(stack.full), quantiles=c(0.025,0.5,0.975), link=link),
              control.compute = list(dic=TRUE,config=TRUE,cpo=TRUE), 
              control.inla = list(int.strategy = "eb", reordering = "metis"),
              control.mode = list(theta=c(1.9132225,-3.750117,1.224630,1.572346)),#set initial for parameters
              verbose=TRUE)

res.s <- inla.spde.result(result,"i",spde,do.transform = TRUE)
#eg for the range
inla.qmarginal(c(0.025,0.5,0.975),res.s$marginals.range.nominal[[1]])
#for point spatial variance
inla.qmarginal(c(0.025,0.5,0.975), res.s$marginals.variance.nominal[[1]])
summary(result)
log.score<--mean(log(result$cpo$cpo),na.rm = TRUE)
log.score
dic<-result$dic$dic[1]
dic
result$model.random

save.image("fitting_result.RData")



####4 prediction####
nsample=500
samp<-inla.posterior.sample(n = nsample, result, use.improved.mean = TRUE, seed=1000)#int.stratey
h<-t(sapply(samp, function(x) x$latent))

#get random i
a<-paste(c("i.","0"),collapse="")
b<-paste(c("i.",dim(Ap_bino)[2]-1),collapse="")
beg<-grep(a, rownames(samp[[1]]$latent))
end<-grep(b, rownames(samp[[1]]$latent))
samp.s<-t(h[,beg:end])

#get coefficient
beg1<-grep("Intercept", rownames(samp[[1]]$latent))
end1<-grep("Diag2", rownames(samp[[1]]$latent))
samp.coef<-t(h[,beg1:end1])

#predicted data
data.whole<-data.whole
data.whole$Diag1<-0
data.whole$Diag2<-0

data.whole$tpop<-NA
data.whole[which(data.whole$NAME_0=="Thailand"),"tpop"]<-data.whole[which(data.whole$NAME_0=="Thailand"),"population"]*68658000/69269551
data.whole[which(data.whole$NAME_0=="Laos"),"tpop"]<-data.whole[which(data.whole$NAME_0=="Laos"),"population"]*6664000/6847193
data.whole[which(data.whole$NAME_0=="Vietnam"),"tpop"]<-data.whole[which(data.whole$NAME_0=="Vietnam"),"population"]*93572000/91039071
data.whole[which(data.whole$NAME_0=="Cambodia"),"tpop"]<-data.whole[which(data.whole$NAME_0=="Cambodia"),"population"]*15518000/16004622

#population2018
data.whole$tpop1<-NA
data.whole[which(data.whole$NAME_0=="Thailand"),"tpop1"]<-data.whole[which(data.whole$NAME_0=="Thailand"),"tpop"]*exp(3*0.22/100)
data.whole[which(data.whole$NAME_0=="Laos"),"tpop1"]<-data.whole[which(data.whole$NAME_0=="Laos"),"tpop"]*exp(3*1.45/100)
data.whole[which(data.whole$NAME_0=="Vietnam"),"tpop1"]<-data.whole[which(data.whole$NAME_0=="Vietnam"),"tpop"]*exp(3*1.00/100)
data.whole[which(data.whole$NAME_0=="Cambodia"),"tpop1"]<-data.whole[which(data.whole$NAME_0=="Cambodia"),"tpop"]*exp(3*1.49/100)

coords.pred <- cbind(data.whole$Longitude,data.whole$Latitude)
colnames(coords.pred)<-c("x","y")
plot(mesh)
points(boundaries,lwd=3)
points(as.matrix(coords),pch=20,cex=1,col=2)
points(as.matrix(coords.pred),pch=20,cex=1,col=3)

#Prediction of points
Apred_p = inla.spde.make.A(mesh, loc=matrix(rep(t(as.matrix(coords.pred)),length(1978:2018)),ncol= ncol(as.matrix(coords.pred)),byrow=TRUE), 
                           group.mesh=mesh.t, group = rep(1978:2018,each=nrow(data.whole)))
pred.s1<-Apred_p%*%samp.s

path<-"F:\\result_path\\"#path for result

for (t_start in 0:40){
  print(t_start)
  
  t_end<-t_start+1
  a1=t_start*dim(data.whole)[1]+1
  a2=t_end*dim(data.whole)[1]
  pred.s<-pred.s1[a1:a2,]
  
  env<-data.whole[,c("Intercept","SurveyType","slst_day11","slst_day12","shii","swb_new","selevation","scity","Diag1","Diag2")]
  pred.env<-as.matrix(env)%*%samp.coef
  pred.linear<-pred.env+pred.s
  
  ###transfer to prevalence
  pred.prev<-exp(pred.linear)/(1+exp(pred.linear))
  pred.prev<-as.matrix(pred.prev)
  pred.prev<-data.frame(pred.prev)
  data.whole<-data.frame(data.whole)
  pred.samp<-cbind(data.whole$Longitude,data.whole$Latitude,pred.prev,
                   data.whole$tpop1,data.whole$NAME_0,data.whole$ID_0,
                   data.whole$NAME_1,data.whole$ID_1,data.whole$NAME_2,data.whole$ID_2)
  pred.samp<-data.frame(pred.samp)
  
  ###summary low, median, up, sd
  n<-nrow(data.whole)
  pred.prev<-as.matrix(pred.prev)
  library(matrixStats)
  lower<-rowQuantiles(pred.prev, probs=0.025)
  upper<-rowQuantiles(pred.prev, probs=0.975)
  median<-rowQuantiles(pred.prev, probs=0.5)
  pred.prev<-as.matrix(pred.prev)
  stdev<-rowSds(pred.prev)
  results<-cbind(coords.pred[,1:2],lower,median,upper,stdev,data.whole$tpop1)
  colnames(results)<-c("long","lat","lower","median","upper","stdev","tpop1")
  results<-data.frame(results)
  results<-cbind(results,data.whole$ID_0,data.whole$NAME_0,data.whole$ID_1,data.whole$NAME_1)
  
  write.dbf(results,paste(path,"result_year",t_start+1978,".dbf",sep=""))
}

save.image(paste(path,"prediction_result.RData",sep=""))
