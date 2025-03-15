# 702-project-
---
title: "Binf702 Project"
date: "2024-04-15"
output:
  pdf_document:
    latex_engine: xelatex
---
```{r}
library(plsgenomics)
data(SRBCT)      # 83 patients/samples and 2308 
#genessum
sum(SRBCT$Y==1)   # 29 cases of Ewing sarcoma (EWS)
sum(SRBCT$Y==2)   # 11 cases of Burkitt lymphoma (BL)
sum(SRBCT$Y==3)   # 18 cases of neuroblastoma (NB)
sum(SRBCT$Y==4)   # 25 cases of rhabdomyosarcoma (RMS)

# check data structure
str(SRBCT)
khan_df <- as.data.frame(SRBCT$X); dim(khan_df)   # 83 2308
khan_fac <- factor(SRBCT$Y,levels=1:4, labels=c("EWS","BL","NB","RMS"))
table(khan_fac)
gnames <- as.data.frame(SRBCT$gene.names)

# check SRBCT$gene.names
length(gnames$Image.Id.)           # 2308
length(unique(gnames$Image.Id.))   # 2277
dup_gnames <- gnames$Image.Id.[duplicated(gnames$Image.Id.)]# 31
dup_gnames
na_Image.Id. <- subset(gnames,is.na(gnames$Image.Id.)==TRUE)   # 5 have no Image.Id. and no Gene.Description
na_Image.Id.
na_Gene.Description <- subset(gnames,is.na(gnames$Gene.Description)==TRUE | trimws(gnames$Gene.Description)=="") # 78
na_Gene.Description
# There are 31 duplicated Image.Id. in which 5 have no Image.Id. and 78 have no Gene.Description
# Create new id's, and will address each gene using this new unique gene id.
gnames$id <- paste0("g",1:nrow(gnames))
colnames(khan_df) <- gnames$id

# order patients by cancer type
t1 <- cbind(khan_df,khan_fac)
t2 <- t1[order(t1$khan_fac),]
khan_fac <- t2$khan_fac           # EWS=1-29; BL=30-40; NB=41-58; RMS=59-83
khan_df <- t2[,1:ncol(t2)-1]
khan_df_t <- as.data.frame(t(khan_df))
colnames(khan_df_t) <- paste0("p",1:ncol(khan_df_t))
rownames(khan_df) <- colnames(khan_df_t)

set.seed(123)
training <- sample(1:83, 65, replace = FALSE)
testing <- setdiff(1:83,training) 
train_khan_df <- khan_df[training,]
test_khan_df <- khan_df[testing,]


```

DATA ASSESSMENT – HYPOTHESIS TESTING:

```{r}
# Check normality:

EWS.shapiro.pval <- apply(khan_df[khan_fac=="EWS",], 2, function(x) shapiro.test(x)$p.value)
BL.shapiro.pval <- apply(khan_df[khan_fac=="BL",], 2, function(x) shapiro.test(x)$p.value)
NB.shapiro.pval <- apply(khan_df[khan_fac=="NB",], 2, function(x) shapiro.test(x)$p.value)
RMS.shapiro.pval <- apply(khan_df[khan_fac=="RMS",], 2, function(x) shapiro.test(x)$p.value)

# proportion of normally distributed
sum(EWS.shapiro.pval>0.05)/length(EWS.shapiro.pval)   # 0.3461872
sum(BL.shapiro.pval>0.05)/length(BL.shapiro.pval)     # 0.7677643
sum(NB.shapiro.pval>0.05)/length(NB.shapiro.pval)     # 0.6078856
sum(RMS.shapiro.pval>0.05)/length(RMS.shapiro.pval)   # 0.4458406

# Both EWS and RMS have less than 50% normality, and EWS is only about 35% normal.  BL is the most normal.

```

```{r}
# Check outliers:
  
library(outliers)
EWS.grubbs.pval <- apply(khan_df[khan_fac=="EWS",], 2, function(x) grubbs.test(x)$p.value)
BL.grubbs.pval <- apply(khan_df[khan_fac=="BL",], 2, function(x) grubbs.test(x)$p.value)
NB.grubbs.pval <- apply(khan_df[khan_fac=="NB",], 2, function(x) grubbs.test(x)$p.value)
RMS.grubbs.pval <- apply(khan_df[khan_fac=="RMS",], 2, function(x) grubbs.test(x)$p.value)

# proportion of positive outliers
sum(EWS.grubbs.pval<0.05)/length(EWS.grubbs.pval)   # 0.4813692
sum(BL.grubbs.pval<0.05)/length(BL.grubbs.pval)     # 0.2318024
sum(NB.grubbs.pval<0.05)/length(NB.grubbs.pval)     # 0.3041594
sum(RMS.grubbs.pval<0.05)/length(RMS.grubbs.pval)   # 0.4441075

# EWS and RMS have more outliers which are consistent with their lower normality rate.  BL has the least outlier.

```

Dendrogram with all dimensions:

```{r}

khan_df_sd <- scale(khan_df)  # standardize the variables to have mean zero and standard deviation one
plot(hclust(dist(khan_df_sd,method="euclidian"),method="complete"), main="Dendrogram", 
     xlab="Clustering of patients by gene expression", sub="",ylab="Distance"   #, hang=-1
     )
khan_hc_out <- hclust(dist(khan_df_sd,method="euclidian"),method="complete")
khan_hc_clusters <- cutree(khan_hc_out,h=79) 
hc_cutree_num_clusters <- max(khan_hc_clusters)
hc_cutree_num_clusters    ## 6 clusters
plot(hclust(dist(khan_df_sd,method="euclidian"),method="complete"), main="", 
     xlab="Clustering of patients by gene expression", sub="",ylab="Distance"
     #, hang=-1
     )
rect.hclust(khan_hc_out,k = length(unique(khan_hc_clusters)))

```

CART  (Classification and Regression Trees) with all dimensions

```{r}
library(rpart); library(rpart.plot)
rpartFit <- rpart(khan_fac~., method="class", data=khan_df)            # linear model
prp(rpartFit,branch.lwd=4,branch.col="blue",extra=101) 

rpartPredict <- predict(rpartFit, type="class")
table(rpartPredict, khan_fac)       # confusion table
3/(3+28+10+18+24)    # error rate = 0.036
```

SAMPLE VARIABILITY AND PRINCIPAL COMPONENT ANALYSIS:

Variances:

```{r}

eigenval <- eigen(cor(khan_df_t))$values[1:10] # first ten components
sum(eigenval/nrow(khan_df)*100)  # The first ten components represent 81% variance.
sum((eigen(cor(khan_df_t))$values[1:9])/nrow(khan_df)*100) # The first nine components represent 80% variance.
sum((eigen(cor(khan_df_t))$values[1:2])/nrow(khan_df)*100)  # The first two components represent 63% variance.
sum((eigen(cor(khan_df_t))$values[1])/nrow(khan_df)*100)    # The first component represents 57% variance.


```
Estimate 95% Confidence Interval using Bootstrap for the eigenvalues:

```{r}

dt <- khan_df_t
nboot<-1000
boot_eigenval <- array(dim=c(nboot, ncol(dt)))  
set.seed(1)
for (i in 1:nboot){
  dat.star <- dt[sample(1:nrow(dt), replace=TRUE),]
  boot_eigenval[i,] <- eigen(cor(dat.star))$values
  }

for (j in 1:ncol(dt)) print(quantile(boot_eigenval[, j], c(0.025,0.975)))
for (j in 1:10) cat(j, as.numeric(quantile(boot_eigenval[,j], c(0.025,0.975))),"\n" )

# The null hypothesis of the eigenvalue being equal to one is rejected for the first nine components but not rejected for the tenth component. 
# The tenth eigenvalue represents less variance than an individual variable.

```

```{r}
# weights of the first two eigenvectors
-eigen(cor(khan_df_t))$vec[, 1:2]

# All weights of the first eigenvector are positive while the weights of the second principal component contain both positive and negative.  From 30 to 40 (BL), almost all their weights are negative.  

```

Print the names of the genes with the largest expression values from the second principal component:

```{r}
pca <- princomp(khan_df_t, cor=TRUE, scores=TRUE)
pcs <- as.data.frame(pca$scores)
pcs9 <- pcs[,1:9]

o2 <- order(pca$scores[,2])                   # second component scores
gnames[o2[1:10],2:3]                          # first ten genes
gnames[o2[(nrow(gnames)-10):nrow(gnames)],2:3]  # last ten genes

# These biomarkers could likely indicate BL.

```

```{r}

biplot(princomp(khan_df_t, cor=TRUE), pc.biplot=TRUE, cex=0.5, expand=1.0, scale=0.5)   # see Biplot.png

```



```{r}
# Use prcomp():

pr_out <- prcomp(khan_df, scale=TRUE)              # scale the variables (genes)

Cols= function (vec){
           cols= rainbow(length(unique(vec)))
return (cols[as.numeric(as.factor(vec))])
}

par(mfrow=c(1,3))
plot(pr_out$x[,1:2],    col=Cols(khan_fac), pch=19, xlab="Z1", ylab="Z2")  # plot first and second principal components
plot(pr_out$x[,c(1,3)], col=Cols(khan_fac), pch=19, xlab="Z1", ylab="Z3")   # plot first and third principal components
plot(pr_out$x[,c(1,4)], col=Cols(khan_fac), pch=19, xlab="Z1", ylab="Z4")   # plot first and fourth principal components

# Z1 and Z4 give better clustering.

```

```{r}

summary(pr_out)

# pr_out <- prcomp(khan_df, scale=TRUE)   # scale the variables (genes)
pve <- 100*pr_out$sdev^2/sum(pr_out$sdev^2)    # proportion of variance explained
par(mfrow=c(1,2))
plot(pve, type="o", ylab="PVE", xlab="Principal Component",col ="blue")     # scree plot
plot(cumsum(pve), type="o", ylab="Cumulative PVE", xlab="Principal Component", col ="brown3")

summary(pr_out)$importance[2,]    # pve
summary(pr_out)$importance[3,]    # cumsum(pve)

# The 12th PC looks like the elbow.  The first twelve principal components explain about 58% of the variance.  

```

CLUSTERING WITH REDUCED DIMENSIONALITY:

```{r}
library(mclust)
pca_df9 <- pr_out$x[,1:9]  # 9 PC's explain about 51% of variance.

BIC <- mclustBIC(pca_df9)
plot(BIC)  
summary(BIC)
# The VEI model best fits the data.

mod1 <- Mclust(pca_df9, x = BIC)
summary(mod1, parameters = TRUE)
plot(mod1, what = "classification")   # see mclust_BIC_plot.png

```

```{r}
mclust_pca9 <- Mclust(pca_df9)              # to evaluate 1–9 clusters and select the optimal number of components based on BIC
clus9 <- mclust_pca9$classification

plot(pca_df9[,1], pca_df9[,2], pch=19, col=clus9, xlab="PC1",ylab="PC2", main="Clusters for 9 PC's")
legend("topleft", legend=unique(clus9), col=unique(clus9), pch=19, title="Cluster")
text(pca_df9[, 1], pca_df9[, 2], labels=rownames(pca_df9), pos=1, cex=0.8)

plot(pca_df9[,1], pca_df9[,3], pch=19, col=clus9, xlab="PC1",ylab="PC3", main="Clusters for 9 PC's")
legend("topleft", legend=unique(clus9), col=unique(clus9), pch=19, title="Cluster")
text(pca_df9[, 1], pca_df9[, 3], labels=rownames(pca_df9), pos=1, cex=0.8)

plot(pca_df9[,1], pca_df9[,4], pch=19, col=clus9, xlab="PC1",ylab="PC4", main="Clusters for 9 PC's")
legend("topleft", legend=unique(clus9), col=unique(clus9), pch=19, title="Cluster")
text(pca_df9[, 1], pca_df9[, 4], labels=rownames(pca_df9), pos=1, cex=0.8)

```


UNSUPERVISED ML:

kmeans - non-parametric:

```{r}

set.seed(6)
cl6 <- kmeans(pca_df9, centers=6, iter.max = 10) 
plot(pca_df9,pch=3,col=cl6$cluster )
text(pca_df9, labels=cl6$cluster, pos=3, cex=0.8, col =cl6$cluster)
points(cl6$centers, col=1:6, pch = 8, cex=2)

```

USE GENES WITH DIFFERENT MEANS:

```{r}
# select the genes with an ANOVA p-value smaller than 0.000001 (genes with different means over the patient groups).
anova.pValue <- apply(khan_df_t, 1, function(x) anova(lm(x ~ khan_fac))$Pr[1])
probeData <- khan_df_t[which(anova.pValue<0.000001),]

if(FALSE){
probeData_1 <- khan_df_t[which(anova.pValue<0.0000001),] # 123x83
probeData <- khan_df_t[which(anova.pValue<0.000001),]    # 186x83  
probeData1 <- khan_df_t[which(anova.pValue<0.00001),]    # 276x83
probeData2 <- khan_df_t[which(anova.pValue<0.0001),]     # 399x83
probeData3 <- khan_df_t[which(anova.pValue<0.001),]      # 594x83
probeData4 <- khan_df_t[which(anova.pValue<0.01),]       # 924x83
}

```

Parallel Coordinates Plot:

```{r}
library(MASS)

parcoord( t(probeData),
         # (t(probeData))[,1:30],      # first 30 genes
         main = "Parallel Coordinates Plot - all 186 significant genes", lwd = 2,
         col=as.integer(khan_fac)+2,
         oma=c(3,3,1,5),
         cex = 1.1 ) 
par(xpd=TRUE)
legend("bottomright", legend=levels(khan_fac),col=unique(as.integer(khan_fac)+2),lwd=2,title="Cancer")

# Genes of high expression values are visible per cancer  type.

```

```{r}
parcoord(# t(probeData),
         (t(probeData))[,1:30],      # first 30 genes
         main = "Parallel Coordinates Plot - first 30 significant genes", lwd = 2,
         col=as.integer(khan_fac)+2,
         oma=c(3,3,1,5),
         cex = 1.1 ) 
par(xpd=TRUE)
legend("bottomright", legend=levels(khan_fac),col=unique(as.integer(khan_fac)+2),lwd=2,title="Cancer")

```
Boxplots for 10 significant genes:

```{r}
sig10 <- as.data.frame(t(probeData)[,1:10])
par(mfrow=c(1,5))
boxplot(sig10$g1 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(sig10$g2 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(sig10$g17 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(sig10$g29 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(sig10$g52 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
par(mfrow=c(1,5))
boxplot(sig10$g54 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(sig10$g74 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(sig10$g85 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(sig10$g119 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(sig10$g123 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
par(mfrow=c(1,1))

# Visual Differentiation:
# EWS -  g29, g52
# BL -  g1, g85, g123
# NB -  g17
# RMS -  g2
(gnames[gnames$id %in% c("g29","g52","g1","g85","g123","g17","g2"),])[c("id","Gene.Description")]

```


```{r}
biplot(princomp(probeData, cor=TRUE), pc.biplot=TRUE, cex=0.5, expand=1.0, scale=0.5)   # see ana_Biplot.png

```

```{r}
ana_pr_out <- prcomp(t(probeData), scale=TRUE)              # scale the variables (genes)

Cols= function (vec){
           cols= rainbow(length(unique(vec)))
return (cols[as.numeric(as.factor(vec))])
}

par(mfrow=c(1,3))
plot(ana_pr_out$x[,1:2],    col=Cols(khan_fac), pch=19, xlab="Z1", ylab="Z2")  # plot first and second principal components

#plot(ana_pr_out$x[,c(1,3)], col=Cols(khan_fac), pch=19, xlab="Z1", ylab="Z3")   # plot first and third principal components
#plot(ana_pr_out$x[,c(1,4)], col=Cols(khan_fac), pch=19, xlab="Z1", ylab="Z4")   # plot first and fourth principal components

```

```{r}

summary(ana_pr_out)

ana_pve <- 100*ana_pr_out$sdev^2/sum(ana_pr_out$sdev^2)    # proportion of variance explained
par(mfrow=c(1,2))
plot(ana_pve, type="o", ylab="PVE", xlab="Principal Component",col ="blue")     # scree plot
plot(cumsum(ana_pve), type="o", ylab="Cumulative PVE", xlab="Principal Component", col ="brown3")

summary(ana_pr_out)$importance[2,]    # pve
summary(ana_pr_out)$importance[3,]    # cumsum(pve)

# The 6th PC looks like the elbow.  The first six principal components explain about 57% of the variance.  

```

```{r}
library(mclust)

ana_pca_df9 <- ana_pr_out$x[,1:9]                   # 9 PC's explain about 64% of variance.

ana_BIC <- mclustBIC(ana_pca_df9)
plot(ana_BIC)  
summary(ana_BIC)

ana_mod1 <- Mclust(ana_pca_df9, x = ana_BIC)
summary(ana_mod1, parameters = TRUE)
plot(ana_mod1, what = "classification")   # see ana_mclust_BIC_plot.png


ana_mclust_pca9 <- Mclust(ana_pca_df9)              # to evaluate 1–9 clusters and select the optimal number of components based on BIC
ana_clus9 <- ana_mclust_pca9$classification

plot(ana_pca_df9[,1], ana_pca_df9[,2], pch=19, col=ana_clus9,xlab="PC1",ylab="PC2", main="Clusters for 9 PC's")
legend("topleft", legend=unique(ana_clus9), col=unique(ana_clus9), pch=19, title="Cluster")
text(ana_pca_df9[, 1], ana_pca_df9[, 2], labels=rownames(ana_pca_df9), pos=1, cex=0.8)

```

```{r}

set.seed(7)
ana_cl7 <- kmeans(ana_pca_df9, centers=7, iter.max = 10)  
plot(ana_pca_df9,pch=3,col=ana_cl7$cluster )
text(ana_pca_df9, labels=ana_cl7$cluster, pos=3, cex=0.8, col = ana_cl7$cluster)
points(ana_cl7$centers, col=1:7, pch = 8, cex=2)

```

SUPERVISED ML:

#random forest - variable importance plot

```{r}
library(randomForest)
X <- t(probeData)
Y <- khan_fac
rf1 <- randomForest(X,Y,ntree=1000, importance=TRUE, proximity=TRUE)
rf1
varImpPlot(rf1,n.var = 15,pch=19,main=NULL,col="red",gcolor="blue",lcolor="darkgreen")

```

```{r}
# VIP top ten genes with description
imp <- as.data.frame(rf1$importance)
imp1 <- imp[order(imp$MeanDecreaseGini,decreasing=TRUE),]

top_ten_vip <- imp1[1:10,]
vip_gnames10 <- gnames[gnames$id %in% rownames(top_ten_vip),]
vip_gnames10

vip_probeData10 <- (t(probeData))[,colnames(t(probeData)) %in% vip_gnames10$id]
parcoord( vip_probeData10,
         main = "Parallel Coordinates Plot - VIP Top 10 Genes", lwd = 2,
         col=as.integer(khan_fac)+2,
         oma=c(3,3,1,5),
         cex = 1.1 ) 
par(xpd=TRUE)
legend("bottomright", legend=levels(khan_fac),col=unique(as.integer(khan_fac)+2),lwd=2,title="Cancer")

## Top ten genes do not show BL having high expression value.
## Several highly expressed BL genes are showed in the plot "Parallel Coordinates Plot - first 30 significant genes".

top_20_vip <- imp1[1:20,]
vip_gnames20 <- gnames[gnames$id %in% rownames(top_20_vip),]
vip_gnames20

vip_probeData20 <- (t(probeData))[,colnames(t(probeData)) %in% vip_gnames20$id]
parcoord( vip_probeData20,
         main = "Parallel Coordinates Plot - VIP Top 20 Genes", lwd = 2,
         col=as.integer(khan_fac)+2,
         oma=c(3,3,1,5),
         cex = 1.1 ) 
par(xpd=TRUE)
legend("bottomright", legend=levels(khan_fac),col=unique(as.integer(khan_fac)+2),lwd=2,title="Cancer")

## Not even top 20 genes reveal BL's high expression genes


```

```{r}
par(mfrow=c(1,5))
boxplot(khan_df$g246 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(khan_df$g509 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(khan_df$g545 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(khan_df$g742 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(khan_df$g842 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
par(mfrow=c(1,5))
boxplot(khan_df$g1003 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(khan_df$g1389 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(khan_df$g1954 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(khan_df$g1955 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
boxplot(khan_df$g2050 ~ khan_fac, col=unique(as.integer(khan_fac)+2))
par(mfrow=c(1,1))

# Visual Differentiation:
# EWS -  g246, g545
# BL - ?? g742 (low expression)
# NB -  g742
# RMS -  g509

(gnames[gnames$id %in% c("g246","g545","g742","g509"),])[c("id","Gene.Description")]

```


```{r}
# Top 20 genes with description
top_ten_vip20 <- imp1[1:20,]
vip_gnames20 <- gnames[gnames$id %in% rownames(top_ten_vip20),]
vip_gnames20

vip_probeData20 <- (t(probeData))[,colnames(t(probeData)) %in% vip_gnames20$id]
parcoord( vip_probeData20,
         main = "Parallel Coordinates Plot - Top 20 Genes", lwd = 2,
         col=as.integer(khan_fac)+2,
         oma=c(3,3,1,5),
         cex = 1.1 ) 
par(xpd=TRUE)
legend("bottomright", legend=levels(khan_fac),col=unique(as.integer(khan_fac)+2),lwd=2,title="Cancer")

```

SVM

```{r}
library(e1071)
set.seed(0)
svm_training <- sample(1:nrow(khan_df), 55, replace = FALSE)
svm_testing <- setdiff(1:nrow(khan_df),training)  

Xt <- khan_df[svm_training,]
Yt <- khan_fac[svm_training]
set.seed(1)
svm_train <- svm(Xt, Yt, type = "C-classification", kernel = "linear", probability = TRUE)
pt <- predict(svm_train, Xt, probability=TRUE)

Xv <- khan_df[svm_testing,]
Yv <- khan_fac[svm_testing]
pv <- predict(svm_train, Xv, probability=TRUE)

table(pt, Yt)
table(pv, Yv)

## Error rate if the test set is 0.023.  
1/(1+16+6+7+14)

```

Multi-layer Perceptron - Non-linear:

```{r}
library(nnet)
set.seed(0)
nn_training <- sample(1:nrow(khan_df), floor(0.67 * nrow(khan_df)), replace = FALSE)
nn_testing <- setdiff(1:nrow(khan_df), nn_training)  
nn_fac <- as.data.frame(khan_fac)
nn_data <- cbind(t(probeData),nn_fac)
nn_data <-nn_data[order(as.integer(nn_data$khan_fac)),]
set.seed(1)
nnest <- nnet(khan_fac ~., data=nn_data, subset=nn_training, size = 5, maxit=500, decay = 0.01)

prednnt <- predict(nnest, nn_data[training,], type = "class")
table(prednnt, Yt=nn_data$khan_fac[training])

prednnv <- predict(nnest, nn_data[testing,], type = "class")
table(prednnv, Yv= nn_data$khan_fac[testing])   

# Both training and testing sets have perfect predictions.

```

Linear Discriminate Analysis (Linear Classifier):

```{r}

library(MASS)

lda_fac <- as.data.frame(khan_fac)
lda_data <- cbind(khan_df,lda_fac)

set.seed(123)
lda_fit <- lda(khan_fac~., data=lda_data )      
lda_pd <- predict(lda_fit, lda_data) 
lda_pd_class <- lda_pd$class
table(pd=lda_pd_class, tru=lda_data$khan_fac)   

# LDA has perfect prediction

```

K-Nearest Neighbor - Non-parametric:

```{r}
library(class)
set.seed(1)
knn_training <- sample(1:nrow(khan_df), floor(0.67 * nrow(khan_df)), replace = FALSE)
knn_testing <- setdiff(1:nrow(khan_df), knn_training)
train <- khan_df[knn_training,]
test <- khan_df[knn_testing,]
train_fac <- khan_fac[knn_training]
test_fac <- khan_fac[knn_testing]

set.seed(2)
knn.pred4 = knn(train, test, cl=train_fac, k=4)
table(knn.pred4, test_fac)

# Error rate for k=4 is 0.07
2/(2+11+3+3+9)


```





