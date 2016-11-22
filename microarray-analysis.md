# microarray analysis
Rich Chien  
October 23, 2016  



### Setup

Load packages


```r
library(affy)   
library(limma)
library(GEOquery)
library(inSilicoMerging)
```


Exploring GEO datasets  
The GEOmetadb library allows a searching through GEO data via SQlite database.

```r
library(GEOmetadb)

# Download SQlite file with all meta data (large file)
if (!file.exists("GEOmetadb.sqlite")) {
    # Download database only if it's not done already
    getSQLiteFile()
}

geo <- dbConnect(SQLite(), "GEOmetadb.sqlite")

dbListTables(geo)
dbListFields(geo)

# Find a dataset (GSE12276) via query
dbGetQuery(geo_con, "SELECT gse.ID, gse.title, gse.gse FROM gse WHERE gse.gse='GSE12276';")
```


Full set of GEO databases used in this study  
\* Duplicate samples were removed from TAM, UNT, UPP


|Set.Name   |Source               |Survival  |Sample.. |
|:----------|:--------------------|:---------|:--------|
|BC         |GEO:GSE21653         |RFS       |266      |
|DFHCC      |GEO: GSE19615        |DMFS      |115      |
|EMC2       |GEO: GSE12276        |DMFS      |204      |
|E-TABM-158 |AE: E-TABM-158       |DMFS, RFS |118      |
|HATZIS     |GEO:GSE25055         |DMFS      |310      |
|HATZIS2    |GEO:GSE25065         |DMFS      |198      |
|MAINZ      |GEO: GSE11121        |DMFS      |200      |
|MDA5       |GEO: GSE17705        |DMFS      |298      |
|MSK        |GEO: GSE2603         |DMFS      |99       |
|TAM        |GEO: GSE6532/GSE9195 |DMFS, RFS |251*     |
|TRANSBIG   |GEO: GSE7390         |DMFS, RFS |198      |
|UNT        |GEO: GSE2990         |DMFS, RFS |83*      |
|UPP        |GEO: GSE3494         |RFS       |190*     |
|VDX        |GEO: GSE2034/GSE5327 |DMFS      |344      |
|VDX3       |GEO: GSE12093        |DMFS      |136      |


Download GEO Raw files and unzip  
Example using the EMC dataset GSE12276


```r
if(length(dir("Data/GEO/",pattern="GSE12276"))==0)
{
  getGEOSuppFiles("GSE12276", makeDirectory = TRUE, baseDir = "./Data/GEO/")
  untar("./Data/GEO/GSE12276/GSE12276_RAW.tar", exdir="./Data/GEO/GSE12276/", tar=Sys.getenv("TAR"))
}

# (optional) gunzip
cels <- list.files("./Data/GEO/", pattern = "[gz]")
sapply(paste("./Data/GEO/", cels, sep="/"), gunzip)
```

Download phenotype data from source


RMA normalization


```r
celfiles <- system.file("extdata", package="arrays")
eset <- justRMA(phenoData=phenoData, celfile.path=celfiles)
```

Build combined database with COMBAT batch effect removal 


```r
esetlist = list(list_of_eset)
memory.limit(size=6000)
combinedset = merge(esetlist, method="COMBAT")
```

When multiple probes match same gene
* take mean of all probes matching to the gene. (http://www.inside-r.org/packages/cran/WGCNA/docs/collapseRows)
* take the probe with the largest variance. (nsfilter from genefilter)

filter off dups and affy probes

```r
fcombinedset <- nsFilter(combinedset, require.entrez=FALSE, require.GOBP=FALSE, require.GOCC=FALSE, require.GOMF=FALSE, require.CytoBand=FALSE, remove.dupEntrez=TRUE, var.func=IQR, var.cutoff=0.5, var.filter=FALSE, filterByQuantile=FALSE, feature.exclude="^AFFX")
```


### Quality Metrics  

array Quality Metrics can be used to generate a quality report for microarray data


```r
library(arrayQualityMetrics)
arrayQualityMetrics(fcombined,force = FALSE,do.logtransform = FALSE,spatial = FALSE)
```


### Breast cancer subtype generation  

PAM50 clustering provided by the genefu package can be applied to predict breast cancer subtypes.


```r
library(genefu)

# expression data
expdata <- t(exprs(fcombinedset))

# meta data
pdata <- pData(fcombinedset)

# probe data
fdata <- fData(fcombinedset)
# rename column
names(fdata)[names(fdata)=="ENTREZID"] <- "EntrezGene.ID"
fdata$probe <- rownames(fdata)

# make prediction
pred.pam50 <- intrinsic.cluster.predict(
  sbt.model = pam50.robust, data=expdata, annot=fdata, 
  do.mapping = TRUE, do.prediction.strength = FALSE, 
  verbose = FALSE)

table(pred.pam50$subtype)

pdata$pam50 <- pred.pam50$subtype
```


### Differential gene expression  

Explore differentially expressed genes in Basal type breast cancer.


```r
# Set up design matrix
type <- as.factor(pdata$pam50)
design.mat <- model.matrix(~0 + type)

# rename columns to easier names
colnames(design.mat) <- gsub("type", "", colnames(design.mat))

# Set up contrast matrix
contrast.mat <- makeContrasts(
  Basal_LumA = Basal - LumA, 
  Basal_LumB = Basal - LumB, 
  levels=design)

# fit limma model to expression data
library(limma)

fit <- lmFit(exprs(fcombinedset, design.mat))
cfit <- contrasts.fit(fit, contrast.mat)
efit <- eBayes(cfit)
# return the top 10 results for contrasts
# coef=1 is first contrast, coef=2 is 2nd contrast
res <- topTable(efit, number=10, coef=1)

# annotate with gene names
biocLite("hgu133plus2.db")
library(hgu133plus2.db)
library(annotate)

gene <- getSYMBOL(res$ID, "hgu133plus2")
res <- cbind(res, gene)
```




