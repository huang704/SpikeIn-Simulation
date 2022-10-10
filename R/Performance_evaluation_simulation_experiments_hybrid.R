library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pROC)
library(RColorBrewer)

sim_sub_sd <- c(0.05, 0.1, 0.15, 0.2, 0.4)
methods <- c("moderatedCMS", "moderatedCS", "moderatedCM", "moderatedC", "limmaCS", "limmaCMS")
new_methods <- c("moderatedTMS", "moderatedTS", "moderatedTM", "moderatedT", "limmaTS", "limmaTMS")
total <- 0.006659+0.00133+0.0131
sim_sub_perc <- round(sim_sub_sd^2/(sim_sub_sd^2+total),4)


all_data <- list()
all_data_summary <- list()
count = 0
for(j in seq_along(sim_sub_sd)){
  
  message("sim_sub_sd: ", sim_sub_sd[j])
  
  UPS.Proteins <- NULL
  Proteins <- list()
  for(i in seq_along(methods)){
    if(i <= 4){
      load(paste0('simulation-experiments-hybrid/controlledmixtures-sub_sd=', sim_sub_sd[j], '-by-', methods[i], ".testing-EB.rda"))
      data.res <-data.res$ComparisonResult
    }
    else{
      load(paste0('simulation-experiments-hybrid/controlledmixtures-sub_sd=', sim_sub_sd[j], '-by-', methods[i], ".testing.rda"))
    }
    data.res <- as.data.frame(na.omit(data.res[, c("Protein", "Label", "log2FC", "adj.pvalue")]))
    #data.res <- data.res %>% filter(!Protein %in% rm_proteins)
    pro <- as.character(unique(data.res$Protein))
    Proteins[[i]] <- pro
    
    ups.pro <- pro[grepl("ups",pro)]
    UPS.Proteins <- c(UPS.Proteins, ups.pro)
  }
  UPS.Proteins <- unique(UPS.Proteins) 
  shared.prots <- Reduce(intersect,  Proteins) 
  
  ############################### evaluate the results #######################
  res <- list()
  AUC <- NULL
  qvalueTPs <- NULL
  qvalueFPs <- NULL
  met <- NULL
  ntests <- NULL
  qvalue_thres <- 0.05
  for(i in 1:length(methods)){
    count <- count + 1
    if(i <= 4){
      load(paste0('simulation-experiments-hybrid/controlledmixtures-sub_sd=', sim_sub_sd[j], '-by-', methods[i], ".testing-EB.rda"))
      data.res <-data.res$ComparisonResult
    }
    else{
      load(paste0('simulation-experiments-hybrid/controlledmixtures-sub_sd=', sim_sub_sd[j], '-by-', methods[i], ".testing.rda"))
    }
    data.res <- data.res[, c("Protein", "Label", "log2FC", "SE", "DF", "pvalue",  "adj.pvalue")]
    data.res <- data.res %>% filter(!is.na(adj.pvalue) & Label != "G2-G1")
    data.res <- as.data.frame(data.res)
    data.res$Label <- as.character(data.res$Label)
    data.res$Protein <- as.character(data.res$Protein)
    data.res$Method <- new_methods[i]
    
    data.res$ref <- ifelse(data.res$Protein %in% UPS.Proteins, 1, 0)
    res[[count]] <- data.res
    
    met <- c(met, new_methods[i])
    ntests <- c(ntests, nrow(data.res))
    AUC <- c(AUC, roc(as.factor(data.res$ref), data.res$adj.pvalue, quiet = T)$auc)
    
    output3 <- data.res %>% filter(adj.pvalue <= qvalue_thres)
    if(nrow(output3) == 0){
      qvalueTPs <- c(qvalueTPs, 0)
      qvalueFPs <- c(qvalueFPs, 0) 
    } else{
      qvalueTPs <- c(qvalueTPs, sum(output3$ref))
      qvalueFPs <- c(qvalueFPs, (nrow(output3) - sum(output3$ref)))
    }
    rm(data.res)
  }
  
  res_summary <- data.table(Method = met, 
                            qvalueTPs = qvalueTPs, 
                            qvalueFPs= qvalueFPs, 
                            auc= AUC, 
                            ntests = ntests)
  
  res1 <- data.table::rbindlist(res)
  res1$log2FC <- as.numeric(as.character(res1$log2FC))
  res1$adj.pvalue <- as.numeric(as.character(res1$adj.pvalue))
  res1 <- res1 %>% filter(!is.na(adj.pvalue))
  res1$Sub_sd <- sim_sub_perc[j]
  all_data[[j]] <- res1
  
  totalFPs <- res1 %>% filter(ref == 0 ) %>% group_by(Method) %>% dplyr::summarise(totalfps = n(), totalbg = n_distinct(Protein))
  totalTPs <- res1 %>% filter(ref == 1 ) %>% group_by(Method) %>% dplyr::summarise(totaltps = n(), totalups = n_distinct(Protein))
  res_summary <- left_join(res_summary, totalFPs)
  res_summary <- left_join(res_summary, totalTPs)
  res_summary$sensitivity <- round(res_summary$qvalueTPs/res_summary$totaltps, 4)
  res_summary$specificity <- round(1-res_summary$qvalueFPs/res_summary$totalfps, 4)
  res_summary$ppv <- round(res_summary$qvalueTPs/(res_summary$qvalueFPs+res_summary$qvalueTPs), 4)
  res_summary$npro <- res_summary$totalbg + res_summary$totalups
  res_summary$np <- res_summary$qvalueTPs + res_summary$qvalueFPs
  res_summary$auc <- round(res_summary$auc, 4)
  res_summary
  
  res_summary$fdr <- 1-res_summary$ppv
  output <- res_summary[,.(Method, npro, qvalueTPs, qvalueFPs, ppv, fdr, sensitivity, specificity, auc)]
  
  output$Sub_sd <- sim_sub_perc[j]
  
  print(output)
  
  all_data_summary[[count]] <- output
}

all_data_summary <- do.call(rbind, all_data_summary)
all_data_summary <- within(all_data_summary, Method <- factor(Method,  
                                                              levels = new_methods))

write.csv(all_data_summary, file =  paste0("simulation-experiments-hybrid/res.csv")) 

all_data_summary <- all_data_summary %>% filter(Method != "moderatedCMS")

blues_colors <- brewer.pal(n = 9, name = "Blues")
purple_colors <- brewer.pal(n = 9, name = "PuRd")
colors <- c(purple_colors[9], blues_colors[6], blues_colors[8], purple_colors[7], purple_colors[4])

pdf("simulation-experiments-hybrid/Compare_different_sim_vars_TPs.pdf",width=3.5,height=3)
g4 <- all_data_summary %>% ggplot(aes(x=Sub_sd, y= qvalueTPs, group=Method, color = Method, shape = Method, size = Method))+
  geom_line(size=1)+
  geom_point() +
  scale_size_manual(values=c(3, 2, 2, 2, 2))+
  scale_shape_manual(values=c(15, 16, 16, 16, 16))+
  ylim(0,100)+
  scale_color_manual(values=colors) + 
  labs(x=("% of bio variance in total variance"), y = ("# of TPs"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'none')

g4
dev.off()

pdf("simulation-experiments-hybrid/Compare_different_sim_vars_FPs.pdf",width=3.5,height=3)
g4 <- all_data_summary %>% ggplot(aes(x=Sub_sd, y= qvalueFPs, group=Method, color = Method, shape = Method, size = Method))+
  geom_line(size=1)+
  geom_point() +
  scale_size_manual(values=c(3, 2, 2, 2, 2))+
  scale_shape_manual(values=c(15, 16, 16, 16, 16))+
  ylim(0,50)+
  scale_color_manual(values=colors) + 
  labs(x=("% of bio variance in total variance"), y = ("# of FPs"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'none')

g4
dev.off()


pdf("simulation-experiments-hybrid/Compare_different_sim_vars_FDR.pdf",width=3.5,height=3)
g4 <- all_data_summary %>% ggplot(aes(x=Sub_sd, y= fdr, group=Method, color = Method, shape = Method, size = Method))+
  geom_line(size=1)+
  geom_point() +
  scale_size_manual(values=c(3, 2, 2, 2, 2))+
  scale_shape_manual(values=c(15, 16, 16, 16, 16))+
  ylim(0,1)+
  scale_color_manual(values=colors) + 
  labs(x=("% of bio variance in total variance"), y = ("Empirical false discovery rate"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'none')

g4
dev.off()

pdf("simulation-experiments-hybrid/Compare_different_sim_vars_sensitivity.pdf",width=3.5,height=3)
g4 <- all_data_summary %>% ggplot(aes(x=Sub_sd, y= sensitivity, group=Method, color = Method, shape = Method, size = Method))+
  geom_line(size=1)+
  geom_point() +
  scale_size_manual(values=c(3, 2, 2, 2, 2))+
  scale_shape_manual(values=c(15, 16, 16, 16, 16))+
  ylim(0.5,1)+
  scale_color_manual(values=colors) + 
  labs(x=("% of bio variance in total variance"), y = ("Sensitivity"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'none')

g4
dev.off()

pdf("simulation-experiments-hybrid/Compare_different_sim_vars_specificity.pdf",width=3.5,height=3)
g4 <- all_data_summary %>% ggplot(aes(x=Sub_sd, y= specificity, group=Method, color = Method, shape = Method, size = Method))+
  geom_line(size=1)+
  geom_point() +
  scale_size_manual(values=c(3, 2, 2, 2, 2))+
  scale_shape_manual(values=c(15, 16, 16, 16, 16))+
  ylim(0.5,1)+
  scale_color_manual(values=colors) + 
  labs(x=("% of bio variance in total variance"), y = ("Specificity"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'none')

g4
dev.off()

pdf("simulation-experiments-hybrid/Compare_different_sim_vars_auc.pdf",width=3.5,height=3)
g4 <- all_data_summary %>% ggplot(aes(x=Sub_sd, y= auc, group=Method, color = Method, shape = Method, size = Method))+
  geom_line(size=1)+
  geom_point() +
  scale_size_manual(values=c(3, 2, 2, 2, 2))+
  scale_shape_manual(values=c(15, 16, 16, 16, 16))+
  ylim(0.5,1)+
  scale_color_manual(values=colors) + 
  labs(x=("% of bio variance in total variance"), y = ("AUC"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'none')

g4
dev.off()

# plot the estimated fold changes, SE, and DF
all_data <- do.call(rbind, all_data)
all_data <- within(all_data, Method <- factor(Method,  
                                              levels = new_methods))
all_data$pred <- ifelse(all_data$adj.pvalue <= 0.05, 1, 0)

res1 <- all_data
save(res1, file = "simulation-experiments-hybrid/res1.rda")

res1 <- res1 %>% filter(Method != "moderatedTMS" & Label == "0.5-0.125")
res1[res1$Method == "moderatedTS", "Method"] <- "moderatedTS (proposed)"
res1 <- within(res1, Method <- factor(Method, levels = c("moderatedTS (proposed)", "moderatedTM", "moderatedT", "limmaTS", "limmaTMS")))

# Distribution of log2FC
pdf("simulation-experiments-hybrid/Compare_different_sim_vars_se.pdf",width=11,height=4)
g4 <- res1 %>% ggplot(aes(x=as.factor(Sub_sd), y= SE, color= Method))+
  #geom_jitter(position=position_jitter(0.1))+ # v1 : width=0.5, v2 n v3 : width=1
  geom_boxplot(lwd=0.6, outlier.alpha = 0.02, notch = TRUE)+ 
  scale_color_manual(values=colors) + 
  ylim(0,0.5)+
  labs(x=("% of bio variance in total variance"), y = ("Standard error"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'bottom')

g4
dev.off()


pdf("simulation-experiments-hybrid/Compare_different_sim_vars_logfc.pdf",width=11,height=4)
g4 <- res1 %>% ggplot(aes(x=as.factor(Sub_sd), y= log2FC, fill= Method))+
  #geom_jitter(position=position_jitter(0.1))+ # v1 : width=0.5, v2 n v3 : width=1
  geom_boxplot(lwd=0.3, outlier.alpha = 0.02)+ 
  scale_fill_manual(values=colors) + 
  ylim(-0.3,0.3)+
  labs(x=("% of bio variance in total variance"), y = ("Estimated log2-fold change"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'bottom')

g4
dev.off()



pdf("simulation-experiments-hybrid/Compare_different_sim_vars_DF.pdf",width=11,height=4)
g4 <- res1 %>% ggplot(aes(x=as.factor(Sub_sd), y= DF, fill= Method))+
  #geom_jitter(position=position_jitter(0.1))+ # v1 : width=0.5, v2 n v3 : width=1
  geom_violin(lwd=0.3)+ 
  scale_fill_manual(values=colors) +
  #geom_hline(yintercept=39, linetype="dashed", color = "grey") +
  #geom_hline(yintercept=27, linetype="dashed", color = "grey") +
  #geom_hline(yintercept=32, linetype="dashed", color = "grey") +
  labs(x=("% of bio variance in total variance"), y = ("Degree of freedoms"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'bottom',
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"))

g4
dev.off()


########################################################################
########################################################################
sim_sub_sd <- c(0.05, 0.1, 0.15, 0.2, 0.4)
methods <- c("moderatedCMS", "moderatedCS", "moderatedCM", "moderatedC", "limmaCS", "limmaCMS")
new_methods <- c("moderatedTMS", "moderatedTS", "moderatedTM", "moderatedT", "limmaTS", "limmaTMS")
#rm_proteins <- c("P01112")
total <- 0.006659+0.00133+0.0131
sim_sub_perc <- round(sim_sub_sd^2/(sim_sub_sd^2+total),4)

res <- list()
count = 0

for (i in seq_along(sim_sub_sd)){
  for(k in seq_along(methods)){
    count <- count + 1
    if(k <= 4){
      load(paste0('simulation-experiments-hybrid/controlledmixtures-sub_sd=', sim_sub_sd[i], '-by-', methods[k], ".testing-EB.rda"))
      data.res <-data.res$ComparisonResult
    }
    else{
      load(paste0('simulation-experiments-hybrid/controlledmixtures-sub_sd=', sim_sub_sd[i], '-by-', methods[k], ".testing.rda"))
    }
    data.res <- data.res[, c("Protein", "Label", "log2FC", "SE", "DF", "pvalue",  "adj.pvalue")]
    data.res <- data.res %>% filter(Label == "G2-G1"|Label == "G1-G2")
    data.res <- as.data.frame(data.res)
    data.res$Label <- as.character(data.res$Label)
    data.res$Protein <- as.character(data.res$Protein)
    data.res$Method <- new_methods[k]
    data.res$Sub_sd <- sim_sub_perc[i]
    res[[count]] <- data.res
  }
}

res2 <- data.table::rbindlist(res)
res2$log2FC <- as.numeric(as.character(res2$log2FC))
res2$adj.pvalue <- as.numeric(as.character(res2$adj.pvalue))

res2 <- within(res2, Method <- factor(Method, levels = new_methods))

res2$pred <- ifelse(res2$adj.pvalue <= 0.05, 1, 0)
unique(res2$Label)

res2 <- res2 %>% filter(Method != "moderatedTMS")
res2[res2$Method == "moderatedTS", "Method"] <- "moderatedTS (proposed)"
res2 <- within(res2, Method <- factor(Method, levels = c("moderatedTS (proposed)", "moderatedTM", "moderatedT", "limmaTS", "limmaTMS")))


num_testable_prots <- res2 %>% group_by(Sub_sd, Method) %>%
  dplyr::summarise(n= n_distinct(Protein)) %>% 
  spread(Sub_sd, n)
num_testable_prots

sig <- res2 %>% filter(adj.pvalue <= 0.05)
perf <- sig %>% group_by(Sub_sd, Method) %>%
  dplyr::summarise(n= n_distinct(Protein)) %>% 
  spread(Sub_sd, n)
perf

write.csv(perf, file =  paste0("simulation-experiments-hybrid/res_compare_conditions.csv")) 


# plot the estimated fold changes, SE, and DF

# Distribution of log2FC
pdf("simulation-experiments-hybrid/Compare_different_sim_vars_null_se.pdf",width=11,height=4)
g4 <- res2 %>% ggplot(aes(x=as.factor(Sub_sd), y= SE, color= Method))+
  #geom_jitter(position=position_jitter(0.1))+ # v1 : width=0.5, v2 n v3 : width=1
  geom_boxplot(lwd=0.6, outlier.alpha = 0.02, notch = TRUE)+ 
  scale_color_manual(values=colors) + 
  ylim(0,0.5)+
  labs(x=("% of bio variance in total variance"), y = ("Standard error"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'bottom')

g4
dev.off()


pdf("simulation-experiments-hybrid/Compare_different_sim_vars_null_logfc.pdf",width=11,height=4)
g4 <- res2 %>% ggplot(aes(x=as.factor(Sub_sd), y= log2FC, fill= Method))+
  #geom_jitter(position=position_jitter(0.1))+ # v1 : width=0.5, v2 n v3 : width=1
  geom_boxplot(lwd=0.3, outlier.alpha = 0.02)+ 
  scale_fill_manual(values=colors) + 
  ylim(-0.3,0.3)+
  labs(x=("% of bio variance in total variance"), y = ("Estimated log2-fold change"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'bottom')

g4
dev.off()



pdf("simulation-experiments-hybrid/Compare_different_sim_vars_null_DF.pdf",width=11,height=4)
g4 <- res2 %>% ggplot(aes(x=as.factor(Sub_sd), y= DF, fill= Method))+
  #geom_jitter(position=position_jitter(0.1))+ # v1 : width=0.5, v2 n v3 : width=1
  geom_violin(lwd=0.3)+ 
  scale_fill_manual(values=colors) +
  #geom_hline(yintercept=39, linetype="dashed", color = "grey") +
  #geom_hline(yintercept=27, linetype="dashed", color = "grey") +
  #geom_hline(yintercept=32, linetype="dashed", color = "grey") +
  labs(x=("% of bio variance in total variance"), y = ("Degree of freedoms"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour=NA, fill=NA),
        legend.position = 'bottom',
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"))

g4
dev.off()



res2$ref <- 0
res <- rbind(res1, res2)

pdf("simulation-experiments-hybrid/Compare_moderatedTS_se.pdf",width=15,height=3.5)
g4 <- res %>% filter(Method %in% c("moderatedTS (proposed)")) %>% 
  select(Protein, Label, SE, Method, Sub_sd) %>% 
  pivot_wider( names_from = Label, values_from = SE) %>% 
  ggplot(aes(x=`0.5-0.125`, y= `G2-G1`))+
  geom_point() + 
  facet_grid(~ Sub_sd) +
  xlim(0, 1)+
  ylim(0, 1)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(angle=0, size=12),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12))

g4
dev.off()


pdf("simulation-experiments-hybrid/Compare_limmaTS_se.pdf",width=15,height=3.5)
g4 <- res %>% filter(Method %in% c("limmaTS")) %>% 
  select(Protein, Label, SE, Method, Sub_sd) %>% 
  pivot_wider( names_from = Label, values_from = SE) %>% 
  ggplot(aes(x=`0.5-0.125`, y= `G2-G1`))+
  geom_point() + 
  facet_grid(~ Sub_sd) +
  xlim(0, 1)+
  ylim(0, 1)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(angle=0, size=12),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12))

g4
dev.off()


pdf("simulation-experiments-hybrid/Compare_moderatedTS_df.pdf",width=15,height=3.5)
g4 <- res %>% filter(Method %in% c("moderatedTS (proposed)")) %>% 
  select(Protein, Label, DF, Method, Sub_sd) %>% 
  pivot_wider( names_from = Label, values_from = DF) %>% 
  ggplot(aes(x=`0.5-0.125`, y= `G2-G1`))+
  geom_point() + 
  facet_grid(~ Sub_sd) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(angle=0, size=12),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12))

g4
dev.off()


pdf("simulation-experiments-hybrid/Compare_limmaTS_df.pdf",width=15,height=3.5)
g4 <- res %>% filter(Method %in% c("limmaTS")) %>% 
  select(Protein, Label, DF, Method, Sub_sd) %>% 
  pivot_wider( names_from = Label, values_from = DF) %>% 
  ggplot(aes(x=`0.5-0.125`, y= `G2-G1`))+
  geom_point() + 
  facet_grid(~ Sub_sd) +
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(angle=0, size=12),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12))

g4
dev.off()


pdf("simulation-experiments-hybrid/Compare_moderatedTS_se_0.3217.pdf",width=4,height=3.5)
g4 <- res %>% filter(Method %in% c("moderatedTS (proposed)") & Sub_sd == "0.3217") %>% 
  select(Protein, Label, SE, Method, Sub_sd) %>% 
  pivot_wider( names_from = Label, values_from = SE) %>% 
  ggplot(aes(x=`0.5-0.125`, y= `G2-G1`))+
  geom_point() + 
  geom_abline(intercept = 0, linetype = "dashed") +
  xlim(0, 1)+
  ylim(0, 1)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(angle=0, size=12),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12))

g4
dev.off()


pdf("simulation-experiments-hybrid/Compare_limmaTS_se_0.3217.pdf",width=4,height=3.5)
g4 <- res %>% filter(Method %in% c("limmaTS") & Sub_sd == "0.106") %>% 
  select(Protein, Label, SE, Method, Sub_sd) %>% 
  pivot_wider( names_from = Label, values_from = SE) %>% 
  ggplot(aes(x=`0.5-0.125`, y= `G2-G1`))+
  geom_point() + 
  geom_abline(intercept = 0, linetype = "dashed") +
  xlim(0, 1)+
  ylim(0, 1)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(angle=0, size=12),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12))

g4
dev.off()


pdf("simulation-experiments-hybrid/Compare_moderatedTS_df_0.3217.pdf",width=4,height=3.5)
g4 <- res %>% filter(Method %in% c("moderatedTS (proposed)") & Sub_sd == "0.3217") %>% 
  select(Protein, Label, DF, Method, Sub_sd) %>% 
  pivot_wider( names_from = Label, values_from = DF) %>% 
  ggplot(aes(x=`0.5-0.125`, y= `G2-G1`))+
  geom_point() + 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(angle=0, size=12),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12))

g4
dev.off()


pdf("simulation-experiments-hybrid/Compare_limmaTS_df_0.3217.pdf",width=4,height=3.5)
g4 <- res %>% filter(Method %in% c("limmaTS") & Sub_sd == "0.3217") %>% 
  select(Protein, Label, DF, Method, Sub_sd) %>% 
  pivot_wider( names_from = Label, values_from = DF) %>% 
  ggplot(aes(x=`0.5-0.125`, y= `G2-G1`))+
  geom_point() + 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(angle=0, size=12),
        strip.text.x = element_text(size=12), # angle=75
        strip.text.y = element_text(size=12))

g4
dev.off()

pdf("simulation-experiments-hybrid/Compare_se_0.8835.pdf",width=16,height=3.5)
g4 <- res %>% filter(Sub_sd == "0.8835") %>% 
  select(Protein, Label, SE, Method, Sub_sd) %>% 
  pivot_wider( names_from = Label, values_from = SE) %>% 
  ggplot(aes(x=`0.5-0.125`, y= `G2-G1`))+
  geom_point(alpha = 1/20) + 
  geom_abline(intercept = 0, linetype = "dashed") +
  xlim(0, 1)+
  ylim(0, 1)+
  facet_grid(~ Method) +
  labs(x=("Standard error of comparing T2-T1"), y = ("Standard error of comparing C2-C1"), size=12)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(angle=0, size=12, face = "bold"),
        strip.text.x = element_text(size=14, face = "bold"), # angle=75
        strip.text.y = element_text(size=12))

g4
dev.off()

pdf("simulation-experiments-hybrid/Compare_df_0.8835.pdf",width=16,height=3.5)
g4 <- res %>% filter(Sub_sd == "0.8835") %>% 
  select(Protein, Label, DF, Method, Sub_sd) %>% 
  pivot_wider( names_from = Label, values_from = DF) %>% 
  ggplot(aes(x=`0.5-0.125`, y= `G2-G1`))+
  geom_point(alpha = 1/20) + 
  theme_bw()+
  facet_grid(~ Method) +
  labs(x=("Degrees of freedom of comparing T2-T1"), y = ("Degrees of freedom of comparing C2-C1"), size=12)+ 
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.y=element_text(angle=0, size=12), #, vjust = 0.43
        axis.text.x=element_text(angle=0, size=12),
        axis.title=element_text(angle=0, size=12, face = "bold"),
        strip.text.x = element_text(size=14, face = "bold"), # angle=75
        strip.text.y = element_text(size=12))

g4
dev.off()