################################################################################
######################### - GWAS - HEADING DATE - ############################## 
################################################################################

### - 0. Setting the worked directory - ###

setwd("C:/Users/guilh/OneDrive/Documents/WorkshopGWAS")

#------------------------------------------------------------------------------#

### - 1. Read phenotypes and check the assumptions of the models - ###

pheno <- read.delim("PhenotypicData.txt", header = TRUE, sep = "\t", 
                    check.names = FALSE)
?read.delim

### - a. Checking the phenotypic data - ###

dim(pheno)
head(pheno)
str(pheno)
summary(pheno) # GID and ENV are char - we need factors here #

### - b. Transforming GID and ENV in factors - ####

pheno$GID <- as.factor(pheno$GID)
pheno$ENV <- as.factor(pheno$ENV)

str(pheno) # Now, GID and ENV are factors #

#------------------------------------------------------------------------------#

### - 2. Read phenotypes and check the assumptions of the models - ###

geno <- read.delim("GenotypicData.txt", header = TRUE, sep = "\t", 
                   check.names = FALSE)
map <- read.delim("GenotypicMap.txt", header = TRUE, sep = "\t")

### - a. Transposing the geno - ###

library(data.table)

genoT <- transpose (geno)
colnames(geno) <- rownames(geno)
rownames(genoT) <- colnames(geno)

### - b. The column name is not the markers, so..lve that! - ###
names(genoT) <- as.matrix(genoT[1, ])
genoT <- genoT[-1, ]
genoT[] <- lapply(genoT, function(x) type.convert(as.character(x))); genoT

### - c. Filtering the genoT - ###

filter.fun <- function(genoT,IM,MM,H){
  #Remove individuals with more than a certain %
  individual.missing <- apply(genoT,1, function(x){
    return(length(which(is.na(x)))/ncol(genoT))
  })
  #Remove markers with certain % missing data
  marker.missing <- apply(genoT,2,function(x)
  {return(length(which(is.na(x)))/nrow(genoT))
  })
  length(which(marker.missing>0.5))
  #Remove individuals with high heterozygous calls.
  heteroz <- apply(genoT,1,function(x){
    return(length(which(x==0))/length(!is.na(x)))
  })
  filter1 <- genoT[which(individual.missing<IM),
                       which(marker.missing<MM)]
  filter2 <- filter1[,(heteroz<H)]
  return(filter2)
}
genoTF <- filter.fun(genoT[,1:15353],0.5,0.50,0.1)
dim(genoTF)

#------------------------------------------------------------------------------#

### - 3. Imputation of the gwasTF, and performing the Kinship matrix.- ###

library(rrBLUP)

Imputation <- A.mat(genoTF,impute.method="EM", 
                    return.imputed=T,min.MAF=0.05)
K.mat <- Imputation$A ; dim(K.mat) # Kinship matrix#
genoTFI <- Imputation$imputed; 
dim(GgenoTFI) #NEW geno data.
genoTFI[1:5,1:5]
K.mat[1:5,1:5]

#------------------------------------------------------------------------------#

### - 4. Look for population structure effects.- ###

geno.scale <- scale(GwasGenoT.filtered.imputation,center=T
                    ,scale=F) # Data needs to be center #
Svdgeno <- svd(geno.scale)
PCA <- geno.scale%*%Svdgeno$v # Principal components #
colnames(PCA) <- paste("PCA",1:ncol(PCA),sep="")
PCA[1:5,1:5]
plot(round((Svdgeno$d)^2/sum((Svdgeno$d)^2),d=7)[1:10],type="o",
     main="Screeplot",xlab="PCAs",ylab="% variance")
PCA1 <- 100*round((Svdgeno$d[1])^2/sum((Svdgeno$d)^2),d=3); PCA1
PCA2 <- 100*round((Svdgeno$d[2])^2/sum((Svdgeno$d)^2),d=3); PCA2
Eucl <- dist(GwasGenoT.filtered.imputation) # Euclinean distance
Fit <- hclust(Eucl,method="ward") # makes clusters with same size #
groups2 <- cutree(Fit,k=2) # Selecting two clusters.
table(groups2) #Number of individuals per cluster.
plot(PCA[,1],PCA[ ,2],xlab=paste("Pcomp:",PCA1,"%",sep=""),
     ylab=paste("Pcomp:",PCA2,"%",sep=""),pch=20,cex=0.7,col=groups2)

#------------------------------------------------------------------------------#

pheno=pheno[pheno$GID%in%rownames(genoTFI),]
pheno$GID<-factor(as.character(pheno$GID),levels=rownames(genoTFI))
X <- model.matrix(~-1+ENV, data=pheno)
phenoGWAS <- data.frame(GID=pheno$GID,X,HD=pheno$HD)
head(phenoGWAS)

genoTFI <- genoTFI[rownames(genoTFI)%in%
                            phenoGWAS$GID,]
phenoGWAS <- phenoGWAS[phenoGWAS$GID%in%rownames(genoTFI),]
genoTFI <- genoTFI[rownames(genoTFI)%in%
                             rownames(K.mat),]
K.mat <- K.mat[rownames(K.mat)%in%rownames(genoTFI),colnames(K.mat)%in%
                 rownames(genoTFI)]
phenoGWAS <- phenoGWAS[phenoGWAS$GID%in%
                               rownames(K.mat),]

genoTFI <-genoTFI[,match(map$marker,colnames(genoTFI))] 
head(map)

genoTFI <- genoTFI[,colnames(genoTFI)%in%map$marker]
map <- map[map$marker%in%colnames(genoTFI),]
genoTFI2 <- data.frame(mark=colnames(genoTFI),chr=map$chr,pos=map$pos,
                       t(genoTFI))
dim(genoTFI2)

colnames(genoTFI2)[4:ncol(genoTFI2)] <-
  rownames(genoTFI)

#------------------------------------------------------------------------------#

### - 6. Perform GWAS - rrBLUP (population structure and Kinship matrix- ###


HDGWASRESULTS<-GWAS(phenoGWAS,genoTFI2,fixed=colnames(phenoGWAS)[2:4],
                    K=K.mat, plot=T,n.PC = 3)
str(HDGWASRESULTS)


write.csv(HDGWASRESULTS, file = "HDGWASRESULTS")
