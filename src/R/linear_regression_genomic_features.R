
## 26/06/2017
###############

## Read input data
###################
rt = read.table("/Users/brodriguez/Research/References/Annotations/H.sapiens/hg19/GenomicFeatures/retrotransposition_genetic_features.tsv", header=1, sep="\t", stringsAsFactors=F)
head(rt)

hist(rt$nbL1,100) # L1 follows a poisson distribution. So use poisson linear regression

## Make poisson linear regression
##################################
## All the features
model = glm(nbL1 ~ medianRT + log(pmax(medianExpr,10)) + geneDensity, family = poisson, data = rt)
summary(model)
1 - (model$deviance/model$null.deviance)
# [1] 0.4620262 ** 46.20%

## Only replication time
model = glm(nbL1 ~ medianRT, family = poisson, data = rt)
summary(model)
1 - (model$deviance/model$null.deviance)
# [1] 0.4461445 ** 44.61%

## Only expression
model = glm(nbL1 ~ log(pmax(medianExpr,10)), family = poisson, data = rt)
summary(model)
1 - (model$deviance/model$null.deviance)
# [1] 0.05840325 ** 5.84%

## Only gene density
model = glm(nbL1 ~ geneDensity, family = poisson, data = rt)
summary(model)
1 - (model$deviance/model$null.deviance)
# [1] 0.05166693 ** 5.16%
