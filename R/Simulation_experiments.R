library(dplyr)
library(tidyr)
library(invgamma)

setwd("~/Dropbox/Ting-thesis/Complex_design/results/Simulation/github")
source('R/modelEvaluation.R')
source("R/groupComparisonTMT.R")
source("R/utils_group_comparison.R")

# select the mixtures for simulation
selected_mixtures <- c("161117_SILAC_HeLa_UPS1_TMT10_Mixture1_01raw",
                       "161117_SILAC_HeLa_UPS1_TMT10_Mixture2_01raw",
                       "161117_SILAC_HeLa_UPS1_TMT10_Mixture3_01raw",
                       "161117_SILAC_HeLa_UPS1_TMT10_Mixture4_01raw",
                       "161117_SILAC_HeLa_UPS1_TMT10_Mixture5_01raw")

load("protein-quant/TMT-Controlledmix-MS3-PD-quant.rda")
data <- processed.quant$ProteinLevelData %>% filter(Run %in% selected_mixtures) 
data <- data %>% filter(!Condition %in% c("Norm", "Empty"))
data$Condition <- as.character(data$Condition)
data$BioReplicate <- paste(data$Mixture, data$Channel, sep = "_")

# assume different spiked-in concentrations are different time points
# assume each mixture and condition has twp biological replicates
data <- data %>% 
  select(Mixture, Channel, BioReplicate, Condition) %>% 
  unique() %>% 
  group_by(Mixture, Condition) %>% 
  mutate(Subject_within_mixture = row_number()) %>%
  right_join(data) %>% 
  mutate(BioReplicate = paste(Mixture, Subject_within_mixture, sep="_")) %>% 
  select(-Subject_within_mixture)

########################################################################
########################################################################
# simulate experiments
subjects <- as.character(unique(data$BioReplicate))
proteins <- as.character(unique(data$Protein))
sim_sub_sd <- c(0.05, 0.1, 0.15, 0.2, 0.4) #simulated biological variance
count <- 0
for(i in seq_along(sim_sub_sd)){
  subject.effect <- list()
  for(k in seq_along(proteins)){
    count <- count + 1
    set.seed(20220628+count)
    
    sub_subject.effect <- data.frame(BioReplicate = subjects,
                                     Sim_noise = rnorm(length(subjects), mean = 0, sd = sqrt(rinvgamma(1, shape = 1/(sim_sub_sd[i]^2)+1))))
    subject.effect[[proteins[k]]] <- sub_subject.effect
  }
  subject.effect <- rbindlist(subject.effect, use.names = TRUE, idcol = "Protein")
  
  # mean(subject.effect$Sim_noise)
  # sd(subject.effect$Sim_noise)
  
  data.sim <- left_join(data, subject.effect)
  data.sim$Abundance <- data.sim$Abundance + data.sim$Sim_noise
  data.sim <- data.sim %>% dplyr::select(-Sim_noise)
  data.sim <- as.data.frame(data.sim)
  
  save(data.sim, file = paste0("simulation-experiments/controlledmixtures-sub_sd=", sim_sub_sd[i], ".abun.rda"))
  
  rm(subject.effect)
  rm(data.sim)
}


########################################################################
########################################################################
# run statistical modeling with simulated sub var
sim_sub_sd <- c(0.05, 0.1, 0.15, 0.2, 0.4)

fx1 <- formula("Abundance ~ 1 + Group + (1|Mixture) + (1|Mixture:Subject)")
fx2 <- formula("Abundance ~ 1 + Group + (1|Subject)")
fx3 <- formula("Abundance ~ 1 + Group + (1|Mixture)")
fx4 <- formula("Abundance ~ 1 + Group")

methods <- c("CMS", "CS", "CM", "C")

model_list <- c(fx1, fx2, fx3, fx4)
for(i in seq_along(sim_sub_sd)){
  load(paste0("simulation-experiments/controlledmixtures-", "sub_sd=", sim_sub_sd[i], ".abun.rda"))
  for(k in seq_along(model_list)){
    
    # without moderation
    data.res <- groupComparisonTMT.v2(data = data.sim,
                                      formula = model_list[[k]],
                                      contrast.matrix = "pairwise",
                                      moderated = FALSE)
    filename1 <- paste0('simulation-experiments/controlledmixtures-sub_sd=', sim_sub_sd[i], '-by-moderated', methods[k], ".testing.rda")
    save(data.res, file=filename1)
    rm(data.res)
    rm(filename1)
    
    # with moderation
    data.res <- groupComparisonTMT.v2(data = data.sim,
                                      formula = model_list[[k]],
                                      contrast.matrix = "pairwise",
                                      moderated = TRUE)

    filename2 <- paste0('simulation-experiments/controlledmixtures-sub_sd=', sim_sub_sd[i], '-by-moderated', methods[k], ".testing-EB.rda")
    save(data.res, file=filename2)
    rm(data.res)
    rm(filename2)
  }
  
  data.res <- groupComparison(data = data.sim,
                              method = "limma+timecourse",
                              contrast.matrix = "pairwise",
                              moderated = TRUE, 
                              adj.method = "BH")
  filename3 <- paste0('simulation-experiments/controlledmixtures-sub_sd=', sim_sub_sd[i], '-by-limmaCS.testing.rda')
  save(data.res, file=filename3)
  rm(data.res)
  rm(filename3)
  
  # run two-way limma
  data.res <- groupComparison(data = data.sim,
                              method = "limma+timecourse+twoway",
                              contrast.matrix = "pairwise",
                              moderated = TRUE, 
                              adj.method = "BH")
  filename4 <- paste0('simulation-experiments/controlledmixtures-sub_sd=', sim_sub_sd[i], '-by-limmaCMS.testing.rda')
  save(data.res, file=filename4)
  rm(data.res)
  rm(filename4)
}


########################################################################
########################################################################
# estimate variance components for CS model
sim_sub_sd <- c(0.05, 0.1, 0.15, 0.2, 0.4)

for(i in seq_along(sim_sub_sd)){
  load(paste0("simulation-experiments/controlledmixtures-", "sub_sd=", sim_sub_sd[i], ".abun.rda"))
  
  proteins <- as.character(unique(data.sim$Protein))
  varscomp <- matrix(rep(0, 2*length(proteins)), ncol = 2)
  terms <- c("BioReplicate", "Residual")
  for(j in 1:length(proteins)) {
    message("Protein: ", j)
    sub_data <- data.sim[data.sim$Protein == proteins[j],]
    ## linear mixed model
    fit.mixed <- try(lmer(Abundance ~ 1 + Condition + (1|BioReplicate), data=sub_data), TRUE)
    if(!inherits(fit.mixed, "try-error")){
      vc <- as.data.frame(VarCorr(fit.mixed, comp="Variance"))
      rownames(vc) <- vc$grp
      varscomp[j,] <- vc[terms, "vcov"]
    } else{
      varscomp[j,] <- NA
    }
  }
  colnames(varscomp) <- terms
  rownames(varscomp) <- proteins
  save(varscomp, file = paste0('simulation-experiments/controlledmixtures-sub_sd=', sim_sub_sd[i], '-CS-variance-comp.rda'))
  
  data.res <- groupComparison(data = data.sim,
                              method = "limma+timecourse",
                              contrast.matrix = "pairwise",
                              moderated = TRUE, 
                              adj.method = "BH")
  filename3 <- paste0('simulation-experiments/controlledmixtures-sub_sd=', sim_sub_sd[i], '-by-limmaCS.testing.rda')
  save(data.res, file=filename3)
  rm(data.res)
  rm(filename3)
  
  # run two-way limma
  data.res <- groupComparison(data = data.sim,
                              method = "limma+timecourse+twoway",
                              contrast.matrix = "pairwise",
                              moderated = TRUE, 
                              adj.method = "BH")
  filename4 <- paste0('simulation-experiments/controlledmixtures-sub_sd=', sim_sub_sd[i], '-by-limmaCMS.testing.rda')
  save(data.res, file=filename4)
  rm(data.res)
  rm(filename4)
}