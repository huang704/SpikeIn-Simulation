library(dplyr)
library(tidyr)
library(invgamma)
library(lme4)

setwd("~/Dropbox/Ting-thesis/Complex_design/results/Simulation/github")
source('R/modelEvaluation.R')
source("R/groupComparisonTMT.R")
source("R/utils_group_comparison.R")

fx1 <- formula("Abundance ~ 1 + Group + (1|Mixture) + (1|Mixture:Subject)")
fx2 <- formula("Abundance ~ 1 + Group + (1|Subject)")
fx3 <- formula("Abundance ~ 1 + Group + (1|Mixture)")
fx4 <- formula("Abundance ~ 1 + Group")

########################################################################
########################################################################
selected_mixtures <- c("161117_SILAC_HeLa_UPS1_TMT10_Mixture1_01raw",
                       "161117_SILAC_HeLa_UPS1_TMT10_Mixture2_01raw",
                       "161117_SILAC_HeLa_UPS1_TMT10_Mixture3_01raw",
                       "161117_SILAC_HeLa_UPS1_TMT10_Mixture4_01raw",
                       "161117_SILAC_HeLa_UPS1_TMT10_Mixture5_01raw")

load("protein-quant/TMT-Controlledmix-MS3-PD-quant.rda")
data <- processed.quant$ProteinLevelData %>% filter(Run %in% selected_mixtures) 
data <- data %>% filter(!Condition %in% c("Norm", "Empty"))
data$Condition <- as.character(data$Condition)

# assume different spiked-in concentrations are different time points
# assume there are two conditions (treatments) which do not have protein abundance changes
# assume each mixture and condition corresponds to a biological replicate
data <- data %>% 
  dplyr::select(Mixture, Channel, BioReplicate, Condition) %>% 
  unique() %>% 
  dplyr::group_by(Mixture, Condition) %>% 
  dplyr::mutate(Treatment = paste0("G", row_number())) %>%
  dplyr::right_join(data) %>% 
  dplyr::mutate(BioReplicate = paste(Mixture, Treatment, sep="_"), Time=Condition, Condition = paste(Treatment, Time, sep="_"))

# perform pairwise comparison
levels(as.factor(data$Condition))
# 'Norm' will be removed during tesing and should be not considered in the contrast
comparison1 <-matrix(c(-1/4,-1/4,-1/4,-1/4,
                       1/4,1/4,1/4,1/4),nrow=1)

comparison2 <-matrix(c(-0.5,0.5,0,0,
                       -0.5,0.5,0,0),nrow=1)

comparison3 <-matrix(c(-0.5,0,0.5,0,
                       -0.5,0,0.5,0),nrow=1)

comparison4 <-matrix(c(-0.5,0,0,0.5,
                       -0.5,0,0,0.5),nrow=1)

comparison5 <-matrix(c(0,-0.5,0.5,0,
                       0,-0.5,0.5,0),nrow=1)

comparison6 <-matrix(c(0,-0.5,0,0.5,
                       0,-0.5,0,0.5),nrow=1)

comparison7 <-matrix(c(0,0,-0.5,0.5,
                       0,0,-0.5,0.5),nrow=1)


comparison <- rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6, comparison7)
# Set the column names
colnames(comparison)<- levels(as.factor(data$Condition))
# Set the names of each row
row.names(comparison)<-c("G2-G1", 
                         "0.5-0.125",
                         "0.667-0.125",
                         "1-0.125",
                         "0.667-0.5",
                         "1-0.5",
                         "1-0.667")

comparison

########################################################################
########################################################################
# simulate experiments
subjects <- as.character(unique(data$BioReplicate))
proteins <- as.character(unique(data$Protein))
sim_sub_sd <- c(0.05, 0.1, 0.15, 0.2, 0.4)
count <- 0
for(i in seq_along(sim_sub_sd)){
  subject.effect <- list()
  for(k in seq_along(proteins)){
    count <- count + 1
    set.seed(23210628+count)
    sub_subject.effect <- data.frame(BioReplicate = subjects,
                                     Sim_noise = rnorm(length(subjects), mean = 0, sd = sqrt(rinvgamma(1, shape = 1/(sim_sub_sd[i]^2)+1))))
    subject.effect[[proteins[k]]] <- sub_subject.effect
  }
  subject.effect <- rbindlist(subject.effect, use.names = TRUE, idcol = "Protein")
  
  mean(subject.effect$Sim_noise)
  sd(subject.effect$Sim_noise)
  
  data.sim <- left_join(data, subject.effect)
  data.sim$Abundance <- data.sim$Abundance + data.sim$Sim_noise
  data.sim <- data.sim %>% dplyr::select(-Sim_noise)
  
  save(data.sim, file = paste0("simulation-experiments-hybrid/sub_sd=", sim_sub_sd[i], ".abun.rda"))
  
  rm(subject.effect)
  rm(data.sim)
}

########################################################################
########################################################################
# run statistical modeling with simulated sub var
model_list <- c(fx1, fx2, fx3, fx4)
sim_sub_sd <- c(0.05, 0.1, 0.15, 0.2, 0.4)

methods <- c("CMS", "CS", "CM", "C")

for(i in seq_along(sim_sub_sd)){
  load(paste0("simulation-experiments-hybrid/sub_sd=", sim_sub_sd[i], ".abun.rda"))
  for(k in seq_along(model_list)){
    # # without moderation
    # data.res <- groupComparisonTMT.v2(data = as.data.frame(data.sim),
    #                                   formula = model_list[[k]],
    #                                   contrast.matrix = comparison,
    #                                   moderated = FALSE)
    # 
    # filename1 <- paste0('simulation-experiments-hybrid/controlledmixtures-sub_sd=', sim_sub_sd[i], '-by-moderated', methods[k], ".testing.rda")
    # save(data.res, file=filename1)
    # rm(data.res)
    # rm(filename1)

    # with moderation
    data.res <- groupComparisonTMT.v2(data = as.data.frame(data.sim),
                                      formula = model_list[[k]],
                                      contrast.matrix = comparison,
                                      moderated = TRUE)

    filename2 <- paste0('simulation-experiments-hybrid/controlledmixtures-sub_sd=', sim_sub_sd[i], '-by-moderated', methods[k], ".testing-EB.rda")
    save(data.res, file=filename2)
    rm(data.res)
    rm(filename2)
    }

  data.res <- groupComparison(data = as.data.frame(data.sim),
                              method = "limma+timecourse",
                              contrast.matrix = comparison,
                              moderated = TRUE,
                              adj.method = "BH")
  filename4 <- paste0('simulation-experiments-hybrid/controlledmixtures-sub_sd=', sim_sub_sd[i], '-by-limmaCS.testing.rda')
  save(data.res, file=filename4)
  rm(data.res)
  rm(filename4)

  data.res <- groupComparison(data = as.data.frame(data.sim),
                              method = "limma+timecourse+twoway",
                              contrast.matrix = comparison,
                              moderated = TRUE,
                              adj.method = "BH")
  filename5 <- paste0('simulation-experiments-hybrid/controlledmixtures-sub_sd=', sim_sub_sd[i], '-by-limmaCMS.testing.rda')
  save(data.res, file=filename5)
  rm(data.res)
  rm(filename5)

  rm(data.sim)
}

########################################################################
########################################################################
# estimate variance components for CS model
sim_sub_sd <- c(0.05, 0.1, 0.15, 0.2, 0.4)

for(i in seq_along(sim_sub_sd)){
  load(paste0("simulation-experiments-hybrid/sub_sd=", sim_sub_sd[i], ".abun.rda"))
  
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
  save(varscomp, file = paste0('simulation-experiments-hybrid/controlledmixtures-sub_sd=', sim_sub_sd[i], '-CS-variance-comp.rda'))
}