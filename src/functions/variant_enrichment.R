##############################################
### Variant enrichment analysis functions ###
#############################################




#Aggregation function for data dcasting for enrichments.
# Returns 0 or 1 based on whether the annotation exists for the individual or not 
# (i.e. whether the individual has a rare variant for the gene)
cat.agg <- function(x){
  return((length(x) > 0) + 0)
}

# Function for enrichment analysis: 
#Log odds ratios and 95% Wald confidence intervals from logistic regression models of outlier class 
#as a function of feachures (variant class)

enrich <- function(data, outlier, features, scale = T){
  estims = c()
  data <- as.data.frame(data)
  for(i in 1:length(features)){
    if(scale){
      mod = summary(glm(data[, outlier]  ~ scale(data[, features[i]]), data = data, family = binomial))$coefficients
      estims.temp = data.frame(Cat = features[i], Estim = mod[2, 1], Std = mod[2, 2], Pval = mod[2, 4], stringsAsFactors = F)
      estims = rbind(estims, estims.temp)
    }else{
      mod = summary(glm(data[, outlier]   ~ data[, features[i]], data = data, family = binomial))$coefficients
      estims.temp = data.frame(Cat = features[i], Estim = mod[2, 1], Std = mod[2, 2], Pval = mod[2, 4], stringsAsFactors = F)
      estims = rbind(estims, estims.temp)
    }
  }
  estims = estims[order(estims$Estim, -estims$Pval), ]
  cat.names = as.character(estims$Cat)
  estims$Cat = factor(as.character(estims$Cat), levels = as.character(estims$Cat), labels = cat.names)
  return(estims)
}


