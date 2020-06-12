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
##### data #####
setwd()
ADM1<-read_excel("data_adm1.xlsx")
ADM1<-as.data.frame(ADM1)
ADM2<-read_excel("data_adm2.xlsx")
ADM2<-as.data.frame(ADM2)
ADM3<-read_excel("data_adm3.xlsx")
ADM3<-as.data.frame(ADM3)
point<-read_excel("data_point.xlsx")
point<-as.data.frame(point)
data.whole<-read_excel("data_whole.xlsx")
data.whole<-as.data.frame(data.whole)


##### model set #####
library(INLA)
library(maps)
library(MASS)
library(MBA)
library(fields)
library(combinat)
library(spdep)
library(gtools)
library(matrixStats)

validate<-c("NumPositiv")
exam<-c("NumExamine")
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

# coordinates 
coords <- cbind(data_p$Longitude,data_p$Latitude)
colnames(coords)<-c("x","y")
coords.a1 <- cbind(data_a1$Longitude,data_a1$Latitude)
colnames(coords.a1)<-c("x","y")
coords.a2 <- cbind(data_a2$Longitude,data_a2$Latitude)
colnames(coords.a2)<-c("x","y")
coords.a3 <- cbind(data_a3$Longitude,data_a3$Latitude)
colnames(coords.a3)<-c("x","y")
coords.pred <- cbind(data.whole$Longitude,data.whole$Latitude)
colnames(coords.pred)<-c("x","y")

a<-shapefile("map_adm0.shp")
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

loc<-data_a1[,c("Longitude","Latitude")]#make sure mesh points existed in most district
bnd <- inla.nonconvex.hull(coords.pred,convex=-0.06)#boundary
mesh <- inla.mesh.2d(loc=loc, boundary = bnd, offset=c(0.1, 0.4),
                     max.edge=c(0.35,0.5))
plot(mesh)
points(boundaries,lwd=3)
points(as.matrix(coords),pch=20,cex=1,col=2)

# Matern SPDE model object
nu <- 1 #Matern smoothness parameter
set.seed(1)
dist<-dist(as.matrix(coords.pred), method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
ran.pr <- median(dist)

#between prediction grids
kap.pr <- sqrt(8*nu)/ran.pr
spde <- inla.spde2.matern(mesh=mesh, alpha=2,
                          B.tau=matrix(c(0, 1, 0),nrow=1,ncol=3),
                          B.kappa=matrix(c(0, 0, 1),nrow=1,ncol=3),
                          theta.prior.mean=c(0, log(kap.pr)),
                          theta.prior.prec=c(1, 1)#specifying the joint prior precision for theta
)
hyper=spde$f$hyper.default
s.indices = inla.spde.make.index("s", mesh$n, n.group=length(table(data_p$year)))#SPDE model index vector generation


##### point data set#####
data_p<-rbind(data_p,data_a2[,as.vector(colnames(data_p))],data_a3[,as.vector(colnames(data_p))])
area.poly <- shapefile("map_adm1.shp")
coords <- cbind(data_p$Longitude,data_p$Latitude)
colnames(coords)<-c("x","y")

point.poly.no <- rep(NA, nrow(coords))
sloc<-SpatialPoints(coords,proj4string=CRS(proj4string(area.poly)))
for(i in 1:length(area.poly)){
  t<-(over(sloc,area.poly[i,]))
  x<-as.vector(which(!is.na(t[1])))
  if(length(x)!=0){
    print(i)
    for(j in 1:length(x)){
      point.poly.no[x[j]] <-i
    }
  }
}
data_p$ID=point.poly.no
data_p$SurveyType<-as.numeric(data_p$SurveyType)

Ap = inla.spde.make.A(mesh, loc=as.matrix(coords), group=data_p$year)

stack.p = inla.stack(tag="point",
                     data = list(validate=data_p$NumPositiv,exam=data_p$NumExamine),
                     A = list(Ap, 1),
                     effects = list(c(s.indices,list(Intercept=1)), list(
                       data.frame(slst_day11=data_p$slst_day11, slst_day12=data_p$slst_day12, slst_day=data_p$slst_day,
                                  slst_night11=data_p$slst_night11, slst_night12=data_p$slst_night12, slst_night=data_p$slst_night,
                                  sbio12_11=data_p$sbio12_11, sbio12_12=data_p$sbio12_12, sbio12=data_p$sbio12, 
                                  shii11=data_p$shii11, shii12=data_p$shii12, shii=data_p$shii, 
                                  swb11=data_p$swb11, swb12=data_p$swb12, swb=data_p$swb, swater1=data_p$swater1, swater2=data_p$swater2, urban=data_p$urban, 
                                  landcover1=data_p$landcover1, landcover2=data_p$landcover2, landcover4=data_p$landcover4, landcover5=data_p$landcover5, 
                                  selevation11=data_p$selevation11,selevation12=data_p$selevation12, selevation=data_p$selevation, 
                                  svi11=data_p$svi11, svi12=data_p$svi12, svi=data_p$svi, 
                                  scity11=data_p$scity11, scity12=data_p$scity12, scity=data_p$scity,
                                  year1=data_p$year1, year2=data_p$year2, year3=data_p$year3,
                                  SurveyType=data_p$SurveyType,year=data_p$year))))


##### area data set#####
mesh.coord.in <- NULL
for(i in 1:length(area.poly)){
  mesh.coord <- SpatialPoints(mesh$loc,proj4string=CRS(proj4string(area.poly)))
  x <- as.vector(which(!is.na(over(mesh.coord, area.poly[i,]))))
  x <- x[which(x<=length(mesh.coord))]
  mesh.coord.in <- rbind(mesh.coord.in , mesh$loc[x,])
}

mesh.coord.in <- mesh.coord.in[,1:2]
block <- rep(0, nrow(mesh.coord.in))
print(length(block))
for(i in 1:length(area.poly)){
  locin<-SpatialPoints(mesh.coord.in,proj4string=CRS(proj4string(area.poly)))
  x<-as.vector(which(!is.na((over(locin,area.poly[i,])))))
  if(length(x)!=0){
    for(j in 1:length(x)){
      if(x[j]<=nrow(mesh.coord.in)){
        block[x[j]] <- i
      }
    }
  }
}

data_a1<-ADM1
data_a1<-data_a1[-which(is.na(data_a1$NumExamine)==TRUE),]
data_random_a1t<-data_a1
data_random_a1t$surveyid<-seq(1,nrow(data_random_a1t))
data_random_a1t<-data_random_a1t[with(data_random_a1t,order(data_random_a1t$ID,data_random_a1t$surveyid)),]
data_random_a1t$bs_id<-cumsum(!duplicated(data_random_a1t[,c("ID","surveyid")]))

locint<-cbind(data.frame(locin),block)
colnames(locint)<-c("long_loc","lat_loc","ID")
merg<-merge(locint,data_random_a1t)
merg<-merg[with(merg,order(merg$bs_id)),]

Aa <- inla.spde.make.A(mesh=mesh,loc=as.matrix(data.frame(merg[,c("long_loc","lat_loc")])),group=merg$year,block=merg$bs_id, block.rescale="sum")
stack.a <- inla.stack(tag='areal',
                      data = list(validate=data_random_a1t$NumPositiv,exam=data_random_a1t$NumExamine),
                      A=list(Aa, 1),
                      effects = list(c(s.indices,list(Intercept=1)),
                                     list(data.frame(slst_day11=data_random_a1t$slst_day11, slst_day12=data_random_a1t$slst_day12, slst_day=data_random_a1t$slst_day, 
                                                     slst_night11=data_random_a1t$slst_night11, slst_night12=data_random_a1t$slst_night12, slst_night=data_random_a1t$slst_night,
                                                     sbio12_11=data_random_a1t$sbio12_11, sbio12_12=data_random_a1t$sbio12_12, sbio12=data_random_a1t$sbio12, 
                                                     shii11=data_random_a1t$shii11, shii12=data_random_a1t$shii12, shii=data_random_a1t$shii, 
                                                     swb11=data_random_a1t$swb11, swb11=data_random_a1t$swb12, swb=data_random_a1t$swb, swater1=data_random_a1t$swater1, swater2=data_random_a1t$swater2, 
                                                     landcover1=data_random_a1t$landcover1, landcover2=data_random_a1t$landcover2, 
                                                     landcover4=data_random_a1t$landcover4, landcover5=data_random_a1t$landcover5, urban=data_random_a1t$urban, 
                                                     selevation11=data_random_a1t$selevation11, selevation12=data_random_a1t$selevation12, selevation=data_random_a1t$selevation, 
                                                     svi11=data_random_a1t$svi11, svi12=data_random_a1t$svi12, svi=data_random_a1t$svi, 
                                                     scity11=data_random_a1t$scity11, scity12=data_random_a1t$scity12, scity=data_random_a1t$scity,
                                                     year1=data_random_a1t$year1, year2=data_random_a1t$year2,year3=data_random_a1t$year3,
                                                     SurveyType=data_random_a1t$SurveyType,
                                                     year=data_random_a1t$year))))


##### model fitting process #####
for1<-c("-1","Intercept","f(s, model=spde,group=s.group, control.group=list(model='ar1'))")
para2<-c("SurveyType","scity11+scity12","shii","swb11+swb12","sbio12")

formula <- as.formula(paste("validate ~ ", paste(paste(for1, collapse= "+"),paste(para2, collapse= "+"),sep="+")))
stack.full <- inla.stack(stack.p, stack.a)

link<-rep(NA,length(stack.full$data$data$exam))
link[which(is.na(link))]<-1

result<- inla(formula, data=inla.stack.data(stack.full), family="binomial",
              control.family = list(link="logit"),
              Ntrials=stack.full$data$data$exam,
              control.predictor=list(compute=TRUE, A=inla.stack.A(stack.full), quantiles=c(0.025,0.5,0.975), link=link),
              control.compute = list(dic=TRUE,config=TRUE,cpo=TRUE), 
              control.inla = list(int.strategy = "eb"),
              verbose=TRUE)


res.s <- inla.spde.result(result,"s",spde,do.transform = TRUE)
# eg for the range
inla.qmarginal(c(0.025,0.5,0.975),res.s$marginals.range.nominal[[1]])
# for point spatial variance
inla.qmarginal(c(0.025,0.5,0.975), res.s$marginals.variance.nominal[[1]])
summary(result)
log.score<--mean(log(result$cpo$cpo),na.rm = TRUE)
log.score
dic<-result$dic$dic[1]
dic
result$model.random


##### prediction and mapping ######
set.seed(100)
nsample=500
samp<-inla.posterior.sample(n = nsample, result, use.improved.mean = TRUE, seed=1000)#int.stratey
h<-t(sapply(samp, function(x) x$latent))

###get random s
a<-paste(c("s.","0"),collapse="")
b<-paste(c("s.",dim(Ap)[2]-1),collapse="")
beg<-grep(a, rownames(samp[[1]]$latent))
end<-grep(b, rownames(samp[[1]]$latent))
samp.s<-t(h[,beg:end])

###get coefficient 
beg1<-grep("Intercept", rownames(samp[[1]]$latent))
end1<-grep("sbio12", rownames(samp[[1]]$latent))
samp.coef<-t(h[,beg1:end1])

#read predicted data
data.whole<-data.whole

#adjust population to 2018
data.whole$tpop<-NA
data.whole[which(data.whole$NAME_0=="Thailand"),"tpop"]<-data.whole[which(data.whole$NAME_0=="Thailand"),"population"]*68658000/69269551
data.whole[which(data.whole$NAME_0=="Laos"),"tpop"]<-data.whole[which(data.whole$NAME_0=="Laos"),"population"]*6664000/6847193
data.whole[which(data.whole$NAME_0=="Vietnam"),"tpop"]<-data.whole[which(data.whole$NAME_0=="Vietnam"),"population"]*93572000/91039071
data.whole[which(data.whole$NAME_0=="Cambodia"),"tpop"]<-data.whole[which(data.whole$NAME_0=="Cambodia"),"population"]*15518000/16004622

data.whole$tpop1<-NA
data.whole[which(data.whole$NAME_0=="Thailand"),"tpop1"]<-data.whole[which(data.whole$NAME_0=="Thailand"),"tpop"]*exp(3*0.22/100)
data.whole[which(data.whole$NAME_0=="Laos"),"tpop1"]<-data.whole[which(data.whole$NAME_0=="Laos"),"tpop"]*exp(3*1.45/100)
data.whole[which(data.whole$NAME_0=="Vietnam"),"tpop1"]<-data.whole[which(data.whole$NAME_0=="Vietnam"),"tpop"]*exp(3*1.00/100)

coords.pred <- cbind(data.whole$Longitude,data.whole$Latitude)
colnames(coords.pred)<-c("x","y")
plot(mesh)
points(boundaries,lwd=3)
points(as.matrix(coords),pch=20,cex=1,col=2)
points(as.matrix(coords.pred),pch=20,cex=1,col=3)

# Prediction of points
Apred_p = inla.spde.make.A(mesh, loc=matrix(rep(t(as.matrix(coords.pred)),4),ncol= ncol(as.matrix(coords.pred)),byrow=TRUE), group = rep(1:4,each=nrow(data.whole)))
#select the recent year period,if year=3,a1=2/4*dim(samp.s)[1]+1,a2=3/4*dim(samp.s)[1]
#if year=2,a1=1/4*dim(samp.s)[1]+1,a2=a1=2/4*dim(samp.s)[1]...
pred.s1<-Apred_p%*%samp.s
a1=3*dim(data.whole)[1]+1
a2=4*dim(data.whole)[1]
pred.s<-pred.s1[a1:a2,]

#for predictors
env<-data.whole[,c("Intercept","SurveyType","scity11","scity12","shii","swb11","swb12","sbio12")]
head(env)
pred.env<-as.matrix(env)%*%samp.coef
samp.coef_new<-data.frame(samp.coef)
m<-as.vector(rep(0,8))
for (i in 1:dim(samp.coef)[1]){
  m[i]=0
  for (j in 1:dim(samp.coef)[2]){
    if(samp.coef_new[i,j]>0){m[i]=m[i]+1}
  }
}
m2<-m/500

lower_coef<-rowQuantiles(samp.coef, probs=0.025)
upper_coef<-rowQuantiles(samp.coef, probs=0.975)
median_coef<-rowQuantiles(samp.coef, probs=0.5)
coef_beta<-cbind(median_coef,lower_coef,upper_coef)

pred.linear<-as.matrix(env)%*%samp.coef+pred.s
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
x<-results[,4]*results[,7]
x<-data.frame(x)
sum(x[,1])/sum(data.whole$tpop1)
results<-cbind(results,x)
colnames(results)<-c("long","lat","lower","median","upper","stdev","tpop1","inf")
results<-cbind(results,data.whole$ID_0,data.whole$NAME_0,data.whole$ID_1,data.whole$NAME_1)

write.dbf(results,"result_year4.dbf")