#library(MSstatsTMT)
library(limma)
library(qvalue)
library(data.table)
library(tidyr)
library(dplyr)
library(edgeR)

groupComparison <- function(data,
                            method,
                            contrast.matrix = "pairwise",
                            moderated = TRUE, 
                            adj.method = "BH"){
  
  ## check the option for method
  method.list <- c("msstatstmt", "limma", "anova", "ttest", "edgeR",
                   "limma+twoway", 
                   "anova+twoway",
                   "wilcox",
                   "limma+timecourse",
                   "limma+timecourse+twoway")
  
  if (sum(method == method.list) != 1) {
    stop(" 'method' must be one of the following : 'msstatstmt', 'limma',
             'anova', 'ttest', 'edgeR, 'limma+twoway', 'anova+twoway', 'wilcox', 'limma+timecourse', 'limma+timecourse+twoway'. ")
  }
  
  if(method == "msstatstmt"){
    data.res <- groupComparisonTMT.v2(data = data, 
                                      contrast.matrix = contrast.matrix,
                                      moderated = moderated, # do moderated t test
                                      adj.method = adj.method) # multiple comparison adjustment
  } 
  
  if(method == "limma"){
    data.res <- .limma.model(data = data, 
                             model = "oneway",
                             adj.method = adj.method)
  }
   
  if(method == "limma+twoway"){
    data.res <- .limma.model(data = data, 
                             model = "twoway",
                             adj.method = adj.method)
  }
  
  if(method == "anova"){
    data.res <- .anova.model(data = data, 
                             model = "oneway",
                             contrast.matrix = contrast.matrix,
                             adj.method = adj.method)
  } 
  
  if(method == "anova+twoway"){
    data.res <- .anova.model(data = data, 
                             model = "twoway",
                             contrast.matrix = contrast.matrix,
                             adj.method = adj.method)
  }
  
  if(method == "ttest"){
    data.res <- .ttest.model(data = data, 
                             adj.method = adj.method)
  }
  
  if(method == "edgeR"){
    data.res <- .edgeR.model(data = data,
                             adj.method = adj.method)
    
  }
  
  if(method == "wilcox"){
    data.res <- .wilcox.model(data = data,
                             adj.method = adj.method)
    
  }
  
  if(method == "limma+timecourse"){
    data.res <- .limma.timecourse.model(data = data, 
                                        model = "oneway",
                                        contrast.matrix = contrast.matrix,
                                        adj.method = adj.method)
  }
  
  if(method == "limma+timecourse+twoway"){
    data.res <- .limma.timecourse.model(data = data, 
                                        model = "twoway",
                                        contrast.matrix = contrast.matrix,
                                        adj.method = adj.method)
  }
  
  ### check column name in order to use groupComparisonPlot from MSstats
  colnames(data.res)[colnames(data.res) == 'Comparison'] <- 'Label'
  colnames(data.res)[colnames(data.res) == 'adjusted.pvalue'] <- 'adj.pvalue'
  
  return(data.res)
}

.limma.model <- function(data,
                         model = "oneway",
                         adj.method = "BH") {
  
  data$Subject <- paste(data$Run, data$Channel, sep = "_")
  annotation <- as.data.frame(unique(data[,c("Run",  "Condition", "Subject")]))
  rownames(annotation) <- annotation$Subject
  
  wide_data <- data %>% 
    dplyr::select(Subject, Abundance, Protein) %>% 
    spread(Subject, Abundance)
  rns <- wide_data$Protein
  wide_data <- wide_data %>% dplyr::select(-Protein)
  wide_data <- as.matrix(wide_data)
  rownames(wide_data) <- rns
  
  # limma, factorial design; 
  tr <- annotation[colnames(wide_data), "Condition"]
  ex <- annotation[colnames(wide_data), "Run"]
  
  if(model == "oneway"){
    design <- model.matrix(~ 0 + tr)
  } else{
    design <- model.matrix(~ 0 + tr + ex)
  }

  n <- dim(wide_data)[1]
  fit <- lmFit(wide_data, design)
  
  groups <- paste("tr", unique(tr), sep = "")
  # pairwise comparison
  ## store inference results
  resList <- list()
  for(j in 1:(length(groups)-1)){
    for(k in (j+1):length(groups)){
      g1 <- groups[j]
      g2 <- groups[k]
      comp <- paste(g1, g2, sep="-")
      contrast.matrix <- makeContrasts(contrasts=comp, levels=design)
      fit2 <- contrasts.fit(fit, contrast.matrix)
      fit2 <- eBayes(fit2)
      temp <- as.data.frame(fit2)
      
      proteins <- rownames(fit2$coefficients)
      pvalue <- as.vector(fit2$p.value)
      log2FC <- as.vector(fit2$coefficients)
      DF <- as.vector(fit2$df.total)
      SE <- as.vector(sqrt(fit2$s2.post) * fit2$stdev.unscaled)
      adjusted.pvalue <- p.adjust(pvalue, adj.method)
      resList[[paste(groups[j], groups[k],sep="-")]] <- data.frame("Protein" = proteins,
                                                                   "log2FC" = log2FC,
                                                                   "pvalue" = pvalue,
                                                                   "SE" = SE,
                                                                   "DF" = DF,
                                                                   "adjusted.pvalue" = adjusted.pvalue)
    }
  }
  ## Finalize the inference results
  res <- rbindlist(resList, use.names=TRUE, idcol = "Comparison")
  #res$adjusted.pvalue <- p.adjust(res$pvalue, adj.method)
  res$Comparison <- gsub("tr", "", res$Comparison)
  data.res<-res[, c("Protein",
                    "Comparison",
                    "log2FC",
                    "SE",
                    "DF",
                    "pvalue",
                    "adjusted.pvalue")]
  
  ### check column name in order to use groupComparisonPlot from MSstats
  colnames(data.res)[colnames(data.res) == 'Comparison'] <- 'Label'
  colnames(data.res)[colnames(data.res) == 'adjusted.pvalue'] <- 'adj.pvalue'
  
  return(data.res)
}


.limma.timecourse.model <- function(data,
                                    model = "oneway",
                                    contrast.matrix = "pairwise",
                                    adj.method = "BH") {
  
  data$Replicate <- paste(data$Run, data$Channel, sep = "_")
  annotation <- as.data.frame(unique(data[,c("Run",  "Condition", "BioReplicate", "Replicate")]))
  rownames(annotation) <- annotation$Replicate
  
  wide_data <- data %>% 
    dplyr::select(Replicate, Abundance, Protein) %>% 
    tidyr::spread(Replicate, Abundance)
  rns <- wide_data$Protein
  wide_data <- wide_data %>% dplyr::select(-Protein)
  wide_data <- as.matrix(wide_data)
  rownames(wide_data) <- rns
  
  # limma, factorial design; 
  tr <- as.factor(as.character(annotation[colnames(wide_data), "Condition"]))
  ex <- as.factor(as.character(annotation[colnames(wide_data), "Run"]))
  sub <- as.factor(as.character(annotation[colnames(wide_data), "BioReplicate"]))
  
  if(model == "oneway"){
    design <- model.matrix(~ 0 + tr)
    if(is.matrix(contrast.matrix)){
      colnames(contrast.matrix) <- paste("tr", colnames(contrast.matrix), sep = "")
    }
  } else{
    design <- model.matrix(~ 0  + tr + ex )
    if(is.matrix(contrast.matrix)){
      colnames(contrast.matrix) <- paste("tr", colnames(contrast.matrix), sep = "")
      contrast.matrix.add.cols <- colnames(design)[grepl("ex", colnames(design))]
      contrast.matrix.add <- matrix(0, nrow = nrow(contrast.matrix), ncol=length(contrast.matrix.add.cols))
      colnames(contrast.matrix.add) <- contrast.matrix.add.cols
      contrast.matrix <- cbind(contrast.matrix, contrast.matrix.add)
    }
  }
  dupcor <- duplicateCorrelation(wide_data,design,block=sub)
  fit <- lmFit(wide_data,design,block=sub,correlation=dupcor$consensus)
  n <- dim(wide_data)[1]
  
  ## store inference results
  resList <- list()
  if(is.matrix(contrast.matrix)){
    for(j in seq_len(nrow(contrast.matrix))){
  
      contrast.matrix.single <- as.matrix(contrast.matrix[j,])
      colnames(contrast.matrix.single) <- rownames(contrast.matrix)[j]
      fit2 <- contrasts.fit(fit, contrast.matrix.single)
      fit2 <- eBayes(fit2)
      
      proteins <- rownames(fit2$coefficients)
      pvalue <- as.vector(fit2$p.value)
      log2FC <- as.vector(fit2$coefficients)
      DF <- as.vector(fit2$df.total)
      SE <- as.vector(sqrt(fit2$s2.post) * fit2$stdev.unscaled)
      adjusted.pvalue <- p.adjust(pvalue, adj.method)
      
      resList[[rownames(contrast.matrix)[j]]] <- data.frame("Protein" = proteins,
                                                            "log2FC" = log2FC,
                                                            "pvalue" = pvalue,
                                                            "SE" = SE,
                                                            "DF" = DF,
                                                            "adjusted.pvalue" = adjusted.pvalue)
    }
  } else{
    groups <- paste("tr", unique(tr), sep = "")
    # pairwise comparison
    for(j in 1:(length(groups)-1)){
      for(k in (j+1):length(groups)){
        g1 <- groups[j]
        g2 <- groups[k]
        comp <- paste(g1, g2, sep="-")
        contrast.matrix <- makeContrasts(contrasts=comp, levels=design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        
        proteins <- rownames(fit2$coefficients)
        pvalue <- as.vector(fit2$p.value)
        log2FC <- as.vector(fit2$coefficients)
        DF <- as.vector(fit2$df.total)
        SE <- as.vector(sqrt(fit2$s2.post) * fit2$stdev.unscaled)
        adjusted.pvalue <- p.adjust(pvalue, adj.method)
        resList[[paste(substring(groups[j], 3), substring(groups[k], 3),sep="-")]] <- data.frame("Protein" = proteins,
                                                                                                 "log2FC" = log2FC,
                                                                                                 "pvalue" = pvalue,
                                                                                                 "SE" = SE,
                                                                                                 "DF" = DF,
                                                                                                 "adjusted.pvalue" = adjusted.pvalue)
      }
    }
  }
  
  ## Finalize the inference results
  res <- rbindlist(resList, use.names=TRUE, idcol = "Comparison")
  data.res<-res[, c("Protein",
                    "Comparison",
                    "log2FC",
                    "SE",
                    "DF",
                    "pvalue",
                    "adjusted.pvalue")]
  
  ### check column name in order to use groupComparisonPlot from MSstats
  colnames(data.res)[colnames(data.res) == 'Comparison'] <- 'Label'
  colnames(data.res)[colnames(data.res) == 'adjusted.pvalue'] <- 'adj.pvalue'
  data.res$consensus.correlation <- dupcor$consensus.correlation
  
  return(data.res)
}


.anova.model <- function(data,
                         model = "oneway",
                         contrast.matrix = "pairwise",
                         adj.method = "BH") {
  
  ## remove the rows with NA intensities
  data <- data[!is.na(data$Abundance),]
  colnames(data)[colnames(data) == 'BioReplicate'] <- 'Subject'
  colnames(data)[colnames(data) == 'Condition'] <- 'Group'
  
  Abundance <- Group <- Protein <- NULL
  
  groups <- as.character(unique(data$Group)) # groups
  if(length(groups) < 2){
    stop("Please check the Condition column in annotation file. There must be at least two conditions!")
  }
  
  ## contrast matrix can be matrix or character vector.
  if(is.matrix(contrast.matrix)){
    # comparison come from contrast.matrix
    if (!all(colnames(contrast.matrix) %in% groups)) {
      stop("Please check the contrast.matrix. Column names of contrast.matrix must be matched with conditions!")
    }
  } else{  
    # create constrast matrix for pairwise comparison
    contrast.matrix <- .makeContrast(groups)
  }
  ncomp <- nrow(contrast.matrix)
  
  proteins <- unique(as.character(data$Protein))
  num.protein <- length(proteins)
  res <- as.data.frame(matrix(rep(NA, 6  * num.protein * ncomp), ncol = 6)) ## store the inference results
  colnames(res) <- c("Protein", "Comparison", "log2FC", "pvalue", "SE", "DF")
  data$Group <- as.factor(data$Group) # make sure group is factor
  data$Run <- as.factor(data$Run)
  nrun <- length(unique(data$Run)) # check the number of MS runs in the data
  count <- 0
  for(i in seq_along(proteins)){
    message(paste("Testing for Protein :", proteins[i] , "(", i, " of ", num.protein, ")"))
    
    ## get the data for protein i
    sub_data <- data %>% dplyr::filter(Protein == proteins[i]) ## data for protein i
    ## record the contrast matrix for each protein
    sub.contrast.matrix <- contrast.matrix
    
    sub_groups <- as.character(unique(sub_data$Group)) # groups in the sub data
    sub_groups <- sort(sub_groups) # sort the groups based on alphabetic order
    
    if(model == "oneway"){
      fit <- suppressMessages(try(lm(Abundance ~ 1 + Group, data = sub_data), TRUE))
    } else{
      fit <- suppressMessages(try(lm(Abundance ~ 1 + Group + Run, data = sub_data), TRUE))
    }
    
    if(!inherits(fit, "try-error")){
      ## Estimate the group variance from fixed model
      av <- anova(fit)
      # use error variance for testing
      df.post <- av["Residuals", "Df"]
      
      ## check the model is fittable 
      # single run case 
      ## Get estimated fold change from mixed model
      coeff <- coef(fit)
      coeff <- coeff[!is.na(coeff)]
      
      ## Compare one specific contrast
      # perform testing for required contrasts
      for(j in seq_len(nrow(sub.contrast.matrix))){
        count <- count + 1
        res[count, "Protein"] <- proteins[i] ## protein names
        res[count, "Comparison"] <- row.names(sub.contrast.matrix)[j] ## comparison
        
        # groups with positive coefficients
        positive.groups <- colnames(sub.contrast.matrix)[sub.contrast.matrix[j,]>0]
        # groups with negative coefficients
        negative.groups <- colnames(sub.contrast.matrix)[sub.contrast.matrix[j,]<0]
        # make sure at least one group from each side of the contrast exist
        if(any(positive.groups %in% sub_groups) & 
           any(negative.groups %in% sub_groups)){
          
          contrast.matrix.single <- as.vector(sub.contrast.matrix[j,])
          names(contrast.matrix.single) <- colnames(sub.contrast.matrix)
          
          cm <- .make.contrast.single(fit, contrast.matrix.single, sub_data)
          
          ## logFC
          FC <- (cm%*%coeff)[,1]
          
          ## variance and df
          temp1 <- summary.lm(fit)
          variance <- diag(t(cm) %*% temp1$cov.unscaled %*% cm)*(temp1$sigma ^ 2)
          
          ## calculate the t statistic
          t <- FC/sqrt(variance)
          
          ## calculate p-value
          p <- 2*pt(-abs(t), df = df.post)
          res[count, "pvalue"] <- p
          
          ## save testing results
          res[count, "log2FC"] <- FC
          res[count, "SE"] <- sqrt(variance)
          res[count, "DF"] <- df.post
          
        } else{
          # at least one condition is missing
          res[count, "log2FC"] <- NA
          res[count, "pvalue"] <- NA
          res[count, "SE"] <- NA
          res[count, "DF"] <- NA
          
        }
      } # for constrast matrix
    } else {
      # very few measurements so that the model is unfittable
      for (j in 1:nrow(sub.contrast.matrix)) {
        count <- count + 1
        res[count, "Protein"] <- proteins[i] ## protein names
        res[count, "Comparison"] <- row.names(sub.contrast.matrix)[j] ## comparison
        res[count, "log2FC"] <- NA
        res[count, "pvalue"] <- NA
        res[count, "SE"] <- NA
        res[count, "DF"] <- NA
        
      } # end loop for comparison    
    } # if the linear model is fittable
  } # for each protein
  
  res <- as.data.frame(res[seq_len(count),])
  res$Protein <- as.factor(res$Protein)
  res$log2FC <- as.numeric(as.character(res$log2FC))
  res$pvalue <- as.numeric(as.character(res$pvalue))
  res$adjusted.pvalue <- NA
  comps <- unique(res$Comparison)
  
  ## Adjust multiple tests for each comparison
  for(i in seq_along(comps)){
    res[res$Comparison == comps[i], "adjusted.pvalue"] <- p.adjust(res[res$Comparison == comps[i], "pvalue"], adj.method)
  }
  
  res <- res[, c("Protein",
                 "Comparison",
                 "log2FC",
                 "SE",
                 "DF",
                 "pvalue",
                 "adjusted.pvalue")]
  return(res)
}

.ttest.model <- function(data,
                         adj.method = "BH") {
  
  ## remove the rows with NA intensities
  data <- data[!is.na(data$Abundance),]
  colnames(data)[colnames(data) == 'BioReplicate'] <- 'Subject'
  colnames(data)[colnames(data) == 'Condition'] <- 'Group'
  
  Abundance <- Group <- Protein <- NULL
  
  groups <- as.character(unique(data$Group)) # groups
  if(length(groups) < 2){
    stop("Please check the Condition column in annotation file. There must be at least two conditions!")
  }
  ncomp <- length(groups)*(length(groups)-1)/2
  
  proteins <- unique(as.character(data$Protein))
  num.protein <- length(proteins)
  res <- as.data.frame(matrix(rep(NA, 7 * num.protein * ncomp), ncol = 7)) ## store the inference results
  colnames(res) <- c("Protein", "Comparison", "log2FC", "pvalue", "SE", "DF", "issue")
  data$Group <- as.factor(data$Group) # make sure group is factor
  data$Run <- as.factor(data$Run)
  nrun <- length(unique(data$Run)) # check the number of MS runs in the data
  count <- 0
  for(i in seq_along(proteins)){
    message(paste("Testing for Protein :", proteins[i] , "(", i, " of ", num.protein, ")"))
    
    ## get the data for protein i
    sub_data <- data %>% dplyr::filter(Protein == proteins[i]) ## data for protein i
     
    for(j in 1:(length(groups)-1)){
      for(k in (j+1):length(groups)){
        
        count <- count + 1
        res[count, "Protein"] <- proteins[i] ## protein names
        res[count, "Comparison"] <- paste(groups[j], groups[k], sep="-") ## comparison
        
        sub_data_group1 <- sub_data %>% dplyr::filter(Group == groups[j]) 
        sub_data_group2 <- sub_data %>% dplyr::filter(Group == groups[k])
        
        if(nrow(sub_data_group1)< 2 & nrow(sub_data_group2)< 2){
          res[count, "log2FC"] <- NA
          res[count, "pvalue"] <- NA
          res[count, "SE"] <- NA
          res[count, "DF"] <- NA
          res[count, "issue"] <- "completeMissing"
          
        } else{
          if(nrow(sub_data_group1)< 2){
            res[count, "log2FC"] <- (-Inf)
            res[count, "pvalue"] <- NA
            res[count, "SE"] <- NA
            res[count, "DF"] <- NA
            res[count, "issue"] <- "oneConditionMissing"
            
          } else{
            if(nrow(sub_data_group2)< 2){
              res[count, "log2FC"] <- (Inf)
              res[count, "pvalue"] <- NA
              res[count, "SE"] <- NA
              res[count, "DF"] <- NA
              res[count, "issue"] <- "oneConditionMissing"
            } else{
              ttest <- t.test(sub_data_group1$Abundance, sub_data_group2$Abundance)
              
              res[count, "log2FC"] <- ttest$estimate[1]-ttest$estimate[2]
              res[count, "pvalue"] <- ttest$p.value
              res[count, "SE"] <- (ttest$estimate[1]-ttest$estimate[2])/ttest$statistic
              res[count, "DF"] <- ttest$parameter
              res[count, "issue"] <- NA
            }
          }
        }
      }
    }
  } # for each protein
  
  res <- as.data.frame(res[seq_len(count),])
  res$Protein <- as.factor(as.character(res$Protein))
  res$log2FC <- as.numeric(as.character(res$log2FC))
  res$pvalue <- as.numeric(as.character(res$pvalue))
  res$adjusted.pvalue <- NA
  comps <- unique(res$Comparison)
  
  ## Adjust multiple tests for each comparison
  for(i in seq_along(comps)){
    res[res$Comparison == comps[i], "adjusted.pvalue"] <- p.adjust(res[res$Comparison == comps[i], "pvalue"], adj.method)
  }
  
  res <- res[, c("Protein",
                 "Comparison",
                 "log2FC",
                 "SE",
                 "DF",
                 "pvalue",
                 "adjusted.pvalue",
                 "issue")]
  return(res)
}

.wilcox.model <- function(data,
                         adj.method = "BH") {
  
  ## remove the rows with NA intensities
  data <- data[!is.na(data$Abundance),]
  colnames(data)[colnames(data) == 'BioReplicate'] <- 'Subject'
  colnames(data)[colnames(data) == 'Condition'] <- 'Group'
  
  Abundance <- Group <- Protein <- NULL
  
  groups <- as.character(unique(data$Group)) # groups
  if(length(groups) < 2){
    stop("Please check the Condition column in annotation file. There must be at least two conditions!")
  }
  ncomp <- length(groups)*(length(groups)-1)/2
  
  proteins <- unique(as.character(data$Protein))
  num.protein <- length(proteins)
  res <- as.data.frame(matrix(rep(NA, 7 * num.protein * ncomp), ncol = 7)) ## store the inference results
  colnames(res) <- c("Protein", "Comparison", "log2FC", "pvalue", "SE", "DF", "issue")
  data$Group <- as.factor(data$Group) # make sure group is factor
  data$Run <- as.factor(data$Run)
  nrun <- length(unique(data$Run)) # check the number of MS runs in the data
  count <- 0
  for(i in seq_along(proteins)){
    message(paste("Testing for Protein :", proteins[i] , "(", i, " of ", num.protein, ")"))
    
    ## get the data for protein i
    sub_data <- data %>% dplyr::filter(Protein == proteins[i]) ## data for protein i
    
    for(j in 1:(length(groups)-1)){
      for(k in (j+1):length(groups)){
        
        count <- count + 1
        res[count, "Protein"] <- proteins[i] ## protein names
        res[count, "Comparison"] <- paste(groups[j], groups[k], sep="-") ## comparison
        
        sub_data_group1 <- sub_data %>% dplyr::filter(Group == groups[j]) 
        sub_data_group2 <- sub_data %>% dplyr::filter(Group == groups[k])
        
        if(nrow(sub_data_group1)< 2 & nrow(sub_data_group2)< 2){
          res[count, "log2FC"] <- NA
          res[count, "pvalue"] <- NA
          res[count, "SE"] <- NA
          res[count, "DF"] <- NA
          res[count, "issue"] <- "completeMissing"
          
        } else{
          if(nrow(sub_data_group1)< 2){
            res[count, "log2FC"] <- (-Inf)
            res[count, "pvalue"] <- NA
            res[count, "SE"] <- NA
            res[count, "DF"] <- NA
            res[count, "issue"] <- "oneConditionMissing"
            
          } else{
            if(nrow(sub_data_group2)< 2){
              res[count, "log2FC"] <- (Inf)
              res[count, "pvalue"] <- NA
              res[count, "SE"] <- NA
              res[count, "DF"] <- NA
              res[count, "issue"] <- "oneConditionMissing"
            } else{
              wilcoxtest <- wilcox.test(sub_data_group1$Abundance, sub_data_group2$Abundance, conf.int = TRUE)
              
              res[count, "log2FC"] <- wilcoxtest$estimate
              res[count, "pvalue"] <- wilcoxtest$p.value
              res[count, "SE"] <- NA
              res[count, "DF"] <- wilcoxtest$statistic
              res[count, "issue"] <- NA
            }
          }
        }
      }
    }
  } # for each protein
  
  res <- as.data.frame(res[seq_len(count),])
  res$Protein <- as.factor(as.character(res$Protein))
  res$log2FC <- as.numeric(as.character(res$log2FC))
  res$pvalue <- as.numeric(as.character(res$pvalue))
  res$adjusted.pvalue <- NA
  comps <- unique(res$Comparison)
  
  ## Adjust multiple tests for each comparison
  for(i in seq_along(comps)){
    res[res$Comparison == comps[i], "adjusted.pvalue"] <- p.adjust(res[res$Comparison == comps[i], "pvalue"], adj.method)
  }
  
  res <- res[, c("Protein",
                 "Comparison",
                 "log2FC",
                 "SE",
                 "DF",
                 "pvalue",
                 "adjusted.pvalue",
                 "issue")]
  return(res)
}

## make contrast matrix for pairwise comparisons
#' @keywords internal
.makeContrast <- function(groups) {
  
  ncomp <- length(groups) * (length(groups) - 1) / 2 # Number of comparison
  contrast.matrix <- matrix(rep(0, length(groups) * ncomp), ncol = length(groups))
  colnames(contrast.matrix) <- groups
  
  count <- 0
  contrast.matrix.rownames <- NULL
  for(j in seq_len(length(groups)-1)){
    for(k in (j+1):length(groups)){
      
      count <- count + 1
      # save row name
      contrast.matrix.rownames <- c(contrast.matrix.rownames, paste(groups[j], groups[k], sep = "-"))
      # set constrast value
      contrast.matrix[count, groups[j]] <- 1
      contrast.matrix[count, groups[k]] <- -1
    }
  }
  rownames(contrast.matrix) <- contrast.matrix.rownames
  
  return(contrast.matrix)
}

.eb.fit <- function(dat, design, tr,
                    adj.method = "BH"){
  
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  
  # count the sample size per group per protein
  samples <- paste("tr", tr, sep = "")
  unique_groups <- unique(samples)
  # samples_group_table <- t(apply(dat, 1, function(x) ifelse(is.na(x), NA, samples)))
  # groups_size <- t(apply(samples_group_table, 1, table))
  # groups_size <- groups_size[rownames(fit.eb$coef), ]
  
  # calculate the group mean
  coeff <- fit.eb$coef
  coeff[, colnames(coeff) != "(Intercept)"] <- 
    coeff[, colnames(coeff) != "(Intercept)"] + coeff[, colnames(coeff) == "(Intercept)"]
  # Find the group name for baseline
  colnames(coeff)[colnames(coeff) == "(Intercept)"] <- setdiff(as.character(unique_groups), colnames(coeff))
  coeff <- coeff[rownames(fit.eb$coef), ]
  
  # pairwise comparison
  ## store inference results
  resList <- list()
  for(j in 1:(length(unique_groups)-1)){
    for(k in (j+1):length(unique_groups)){
      g1 <- unique_groups[j]
      g2 <- unique_groups[k]
      comp <- paste(g1, g2, sep="-")
      
      logFC <- coeff[, g1] - coeff[, g2]
      stdev.all <- fit.eb$stdev.unscaled
      if(g1 %in% colnames(stdev.all) & g2 %in% colnames(stdev.all)){
        stdev.pairwise <-  apply(stdev.all[,c(g1,g2)], 1, max)
      } else{
        if(g1 %in% colnames(stdev.all) & !g2 %in% colnames(stdev.all)){
          stdev.pairwise <- stdev.all[,g1]
        } else{
          stdev.pairwise <- stdev.all[,g2]
        }
      }
      
      df.0 <- rep(fit.eb$df.prior, n)
      df.r <- fit.eb$df.residual
      s2.0 <- rep(fit.eb$s2.prior, n)
      s2 <- (fit.eb$sigma)^2
      t.ord <- logFC/fit.eb$sigma/stdev.pairwise
      p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
      
      # s2.post <- fit.eb$s2.post
      # t.mod <- fit.eb$t[, "tr2"]
      # p.mod <- fit.eb$p.value[, "tr2"]
      
      s2.post <- (df.r*s2 +df.0*s2.0)/(df.r+df.0)
      t.mod <- logFC/sqrt(s2.post)/stdev.pairwise
      p.mod <- 2*pt(-abs(t.mod), df.r+df.0) 
      
      proteins <- rownames(fit.eb$coef)
      pvalue <- p.mod
      log2FC <- as.vector(logFC)
      DF <- as.vector(df.r+df.0)
      SE <- as.vector(sqrt(s2.post)*stdev.pairwise)
      #adjusted.pvalue <- qvalue(pvalue)$q
      adjusted.pvalue <- p.adjust(pvalue, adj.method)
      resList[[comp]] <- data.frame("Protein" = proteins,
                                    "log2FC" = log2FC,
                                    "pvalue" = pvalue,
                                    "SE" = SE,
                                    "DF" = DF,
                                    "adjusted.pvalue" = adjusted.pvalue)
    }
  }
  ## Finalize the inference results
  res <- rbindlist(resList, use.names=TRUE, idcol = "Comparison")
  #res$adjusted.pvalue <- p.adjust(res$pvalue, adj.method)
  res$Comparison <- gsub("tr", "", res$Comparison)
  data.res<-res[, c("Protein",
                    "Comparison",
                    "log2FC",
                    "SE",
                    "DF",
                    "pvalue",
                    "adjusted.pvalue")]
  
  ### check column name in order to use groupComparisonPlot from MSstats
  colnames(data.res)[colnames(data.res) == 'Comparison'] <- 'Label'
  colnames(data.res)[colnames(data.res) == 'adjusted.pvalue'] <- 'adj.pvalue'
  
  return(data.res)
}


## make constrast
#	MSstats
#' @importFrom stats coef
#' @importFrom lme4 fixef
#' @keywords internal
.make.contrast.single <- function(fit, contrast, sub_data) {
  
  ## when there are some groups which are all missing
  sub_groups <- as.character(unique(sub_data[, c("Group")]))
  
  # groups with positive coefficients
  positive.groups <- names(contrast)[contrast>0]
  # groups with negative coefficients
  negative.groups <- names(contrast)[contrast<0]
  
  # if some groups not exist in the protein data
  if(!(all(positive.groups %in% sub_groups) & 
       all(negative.groups %in% sub_groups))){
    
    contrast.single <- contrast[sub_groups]
    
    ## tune the coefficients of positive groups so that their summation is 1
    temp <- contrast.single[contrast.single > 0]
    temp <- temp*(1/sum(temp, na.rm = TRUE))
    contrast.single[contrast.single > 0] <- temp
    
    ## tune the coefficients of positive groups so that their summation is 1
    temp2 <- contrast.single[contrast.single < 0]
    temp2 <- temp2*abs(1/sum(temp2, na.rm = TRUE))
    contrast.single[contrast.single < 0] <- temp2
    
    ## set the coefficients of non-existing groups to zero
    contrast[] <- 0
    contrast[sub_groups] <- contrast.single
  }
  
  if (inherits(fit, "lm")) {
    coef_name <- names(stats::coef(fit))
  } else {
    coef_name <- names(lme4::fixef(fit))
  }
  
  ## intercept
  temp <- coef_name[grep("Intercept", coef_name)]
  intercept_c <- rep(0, length(temp))
  names(intercept_c) <- temp
  if (length(temp) == 0) {
    intercept_c <- NULL
  }
  
  ## group
  temp <- coef_name[grep("Group", coef_name)]
  tempcontrast <- contrast[sub_groups]
  group_c <- tempcontrast[gsub("Group", "", temp)] 
  names(group_c) <- temp
  if (length(temp) == 0) {
    group_c<-NULL
  }
  
  ## run
  temp <- coef_name[grep("Run", coef_name)]
  if (length(temp) > 0) {
    run_c <- rep(0, length(temp))
    names(run_c) <- temp
  } else {
    run_c <- NULL
  }
  
  ## combine all
  newcontrast <- c(intercept_c, group_c, run_c)
  if(inherits(fit, "lm")) {
    contrast1 <- newcontrast[!is.na(stats::coef(fit))]
  } else {
    contrast1 <- newcontrast[!is.na(lme4::fixef(fit))]
  }
  
  return(contrast1)
}

.edgeR.model <- function(data,
                         adj.method = "BH") {
  
  data$Subject <- paste(data$Run, data$Channel, sep = "_")
  annotation <- unique(data[,c("Run",  "Condition", "Subject")])
  rownames(annotation) <- annotation$Subject
  
  wide_data <- data %>% 
    dplyr::select(Subject, Abundance, Protein) %>% 
    spread(Subject, Abundance)
  wide_data <- na.omit(wide_data)
  rns <- wide_data$Protein
  wide_data <- wide_data %>% dplyr::select(-Protein)
  wide_data <- as.matrix(wide_data)
  wide_data <- 2^(wide_data)
  rownames(wide_data) <- rns
  
  # set up the sample mapping
  group <- annotation[colnames(wide_data),"Condition"]
  wide_data <- wide_data[,group!="Norm"]
  group <- group[group!="Norm"]
  group <- as.factor(as.character(group))
  # create a DGEList object with the IRS data
  y_irs <- DGEList(counts = wide_data, group = group)
  # we will skip running TMM (using the calcNormFactors function)
  # we need to estimate the dispersion terms (global and local)
  y_irs <- estimateDisp(y_irs)
  
  # the exact test object has columns like fold-change, CPM, and p-values
  gs <- unique(group)
  resList <- list()
  for(j in 1:(length(gs)-1)){
    for(k in (j+1):length(gs)){
      g1 <- gs[j]
      g2 <- gs[k]
      et_irs <- exactTest(y_irs, pair = c(g1, g2))
      
      # the topTags function adds the BH FDR values to an exactTest data frame. Make sure we do not change the row order!
      tt_irs <- topTags(et_irs, n = Inf, sort.by = "none", adjust.method=adj.method)
      tt_irs <- tt_irs$table # tt_sl is a list. We just need the data frame table
      # add the protein name
      tt_irs$Protein <- rownames(tt_irs)
      resList[[paste(g1, g2,sep="-")]] <- data.frame("Protein" = tt_irs$Protein,
                                                     "log2FC" = tt_irs$logFC,
                                                     "pvalue" = tt_irs$PValue,
                                                     "adj.pvalue" = tt_irs$FDR)
    }
  }
  ## Finalize the inference results
  res <- rbindlist(resList, use.names=TRUE, idcol = "Label")
  data.res<-res[, c("Protein",
                    "Label",
                    "log2FC",
                    "pvalue",
                    "adj.pvalue")]
  
  return(data.res)
}