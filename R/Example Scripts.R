#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
# Scripts for Levi Newediuk
# Brett Jesmer, PhD
# brettjesmer@vt.edu
# 08/04/2022
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#

#_#_#_#_#_#_#_#_#_#_#_#_#
# DISTANCE TO RASTER----
#_#_#_#_#_#_#_#_#_#_#_#_#

#step 1, reclassify habitat raster

library(raster)

beginCluster(type="SOCK")

willow<-clusterR(NLCD, reclassify, args=list(rcl=matrix(c(11,NA, 
                                                          12,NA, 
                                                          21,NA, 
                                                          22,NA, 
                                                          23,NA, 
                                                          24,NA, 
                                                          31,NA, 
                                                          41,NA,
                                                          42,NA, 
                                                          43,NA, 
                                                          51,NA, 
                                                          52,NA, 
                                                          71,NA, 
                                                          72,NA, 
                                                          73,NA, 
                                                          74,NA, 
                                                          81,NA, 
                                                          82,NA, 
                                                          90,1, 
                                                          95,NA),
                                                        20,2, byrow=T)),
                 filename = "F:/NLCD/Utah/Utah_willow.tif", format="GTiff", 
                 datatype="INT2U", overwrite=TRUE, progress='text')

#step 2, calculate distance to raster

willow.dist<-distance(willow, filename="F:/NLCD/Utah/Utah_Willow_dist.tif", format="GTiff", dataType="FLT4S", progress='text', overwrite=TRUE)

#step 3, repeat for each habitat class of interest

#_#_#_#_#_#_#_#_#_#_#_#_#
# FOCAL STAT RASTER----
#_#_#_#_#_#_#_#_#_#_#_#_#

#step 1, reclassify habitat raster

wetland.bin<-clusterR(NLCD, reclassify, args=list(rcl=matrix(c(11,0, 
                                                               12,0, 
                                                               21,0, 
                                                               22,0, 
                                                               23,0, 
                                                               24,0, 
                                                               31,0, 
                                                               41,0, 
                                                               42,0,
                                                               43,0, 
                                                               51,0, 
                                                               52,0, 
                                                               71,0, 
                                                               72,0, 
                                                               73,0, 
                                                               74,0, 
                                                               81,0, 
                                                               82,0, 
                                                               90,0, 
                                                               95,1),
                                                             20,2, byrow=T)), 
                      filename="F:/NLCD/Utah/Utah_wetland_bin.tif", format="GTiff", 
                      datatype="INT2U", overwrite=TRUE, progress='text')

#step 2, calculate zonal statistics 

w<-focalWeight(NLCD, d=1000, type="circle")
### can change size of circular window, can also change cirdular window to rectangular or other shape

willow.foc<-focal(willow.bin, w=w, fun=sum, na.rm=TRUE, progress='text',
                  filename="F:/NLCD/Utah/Utah_willow_foc.tif", 
                  format="GTiff", overwrite=TRUE)
### here we calcuate amount of given habitat within 1km circular moving window

#step 3, repeat for each habitat class of interest


#_#_#_#_#_#_#_#_#_#_#_#_#
# HIERARCHICAL CLUSTERING----
#_#_#_#_#_#_#_#_#_#_#_#_#

#find dietary/habitat clusters

library(stats)
library(vegan)
library(cluster)
library(devtools)
install_github("vqv/ggbiplot")

#build dendrogram based on Ward's minimum variance clustering, uses MANOVA to help define groups
stand<-data[,5:19] #select only columns with habitat covariate data
norm<-decostand(stand,"normalize") #standardize/normalize covariate distributions
dist<-vegdist(norm, method="euc") #euclidean dist if negative values present
# see book for other distance options
strat<-hclust(dist, method="ward") #fit hierarchical cluster model
strat$height<-sqrt(strat$height) #standardize values for plotting
plot(strat) #plot

### use silhouette width to determine optimum number of clusters (approach 1)
asw<-numeric(nrow(stand))

for(k in 2:(nrow(stand.wint.F.BC)-1)){
  sil<-silhouette(cutree(strat.wint.F.BC,k=k), dist.wint.F.BC)
  asw[k]<-summary(sil)$avg.width
}

k.best<-which.max(asw)

par(mar=c(7,4,4,2))
plot(1:nrow(stand), asw, type="h", xlab="k (Number of groups)", ylab="Ave Silhouette Width")
points(k.best, max(asw),pch=16, col="red")
axis(1,k.best,paste("optimum",k.best,sep="="),col="red",col.axis="red", font=2, las=2)

#what is the optimum number of clusters?

### use silhouette width, combine dist and binary matrices (approach 2)

#function to compute binary distance matrix from groups
grpdist<-function(x){
  require(cluster)
  gr<-as.data.frame(as.factor(x))
  distgr<-daisy(gr,"gower")
  distgr
}

kt<-data.frame(k=1:nrow(stand), r=0)

for(i in 2:(nrow(stand)-1)){
  gr<-cutree(strat, i)
  distgr<-grpdist(gr)
  mt<-cor(dist, distgr, method="pearson")
  kt[i,2]<-mt
}

(k.best<-which.max(kt$r))

plot(kt$k, kt$r, type="h", xlab="k (number of groups", ylab="Pearson's Correlation")
points(k.best, max(kt$r),pch=16, col="red")
axis(1,k.best,paste("optimum",k.best,sep="="),col="red",col.axis="red", font=2, las=2)

#is there agreement between the two methods????

#plot final dendrogram with n=2 clusters
setwd("C:/Users/brettjesmer/My Drive/Literature/books/Numerical Ecology/NEwR_updated_material_R332/NEwR_updated_material_R332-NEwR_ed1/")
source("hcoplot.R")
library(gclus)

# setwd("/Users/BrettJesmer/Box Sync/Jesmer-Moose/Data/Moose_Data/diet data/2017")
# tiff("wint_diet_clusters_male.tiff",
#      type = "cairo",
#      height = 6,
#      width = 15, units = "in",
#      res= 300,
#      antialias = "subpixel",
#      family = "Arial")
# 
# par(mar=c(6,4,4,2))
hcoplot(strat, dist, k=2)

# dev.off()

#add grouping/cluster variable to winter df
data$group<-cutree(strat, k=2)

#_#_#_#_#_#
# LDA----
#_#_#_#_#_#

library(car)
library(MASS)
library(ggbiplot)
library(dplyr)

data.lda<-data %>% 
  select(group, forest, crop, etc) %>%  #columns containing habitat variables, and grouping variable
  mutate(forest=scale(forest),
         crop=scale(crop),
         etc=scale(etc)) #must scale and center variables for LDA and PCA

lda.mod<-lda(group~., data=data.lda) #fit LDA
lda.mod
plot(lda.mod)

g <- ggbiplot(lda.mod, choices=c(1,2), var.scale=FALSE,
              groups = as.factor(data.lda$group), ellipse = TRUE, 
              ellipse.prob = 0.95, circle = FALSE, var.axes=TRUE)
# g <- g + scale_color_manual(name = '', values=c("black","green","orange","blue","darkturquoise","deeppink"))
# g <- g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                panel.background = element_blank(), axis.line = element_line(colour = "black"),
#                axis.text=element_text(face="bold",size=11, colour="black"), 
#                axis.title = element_text(face="bold", size=14))
# g <- g + xlab("PC1")+ylab("PC2")
g <- g + scale_x_continuous(limits = c(-1.5, 1.5)) + scale_y_continuous(limits = c(-5, 0.5))
# g <- g + geom_vline(xintercept=0, linetype="longdash") + geom_hline(yintercept=0, linetype="longdash")
# g <- g + theme(legend.direction = 'horizontal', legend.position = 'top', 
#                legend.text=element_text(size=10, colour="black",face="bold"),
#                legend.key = element_rect(colour = NA, fill = NA))
print(g)


