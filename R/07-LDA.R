
library(tidyverse)
library(vegan)

# Load data
sample_dat <- readRDS('derived_data/uid_prop_use_samples.rds')

lda_ord <- sample_dat %>% 
  mutate(across(crop_prop:wet_prop, list(sc = scale)))
  select(crop_prop:wet_prop) %>%  #columns containing habitat variables
  mutate(forest=scale(forest),
         crop=scale(crop),
         etc=scale(etc)) #must scale and center variables for LDA and PCA

lda.mod<-lda(group~., data=data.lda)
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




