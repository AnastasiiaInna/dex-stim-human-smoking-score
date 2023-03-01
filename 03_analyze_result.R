# Analyze the result from differently normalized DS

library(dplyr)
library(ggplot2)
library(reshape2) 

# Set up variables

src.data.pre  <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/"

score.qn.fn      <- "smoking_score_DexStim_EPIC_2020_QN.csv"
score.bmiq.fn    <- "smoking_score_DexStim_EPIC_2020_BMIQ.csv"
score.qn.bmiq.fn <- "smoking_score_DexStim_EPIC_2020_QN_BMIQ.csv"

# Load data

score.qn.df      <- read.csv(paste0(src.data.pre, score.qn.fn), sep = ";", header = T, dec = ",")
score.bmiq.df    <- read.csv(paste0(src.data.pre, score.bmiq.fn), sep = ";", header = T, dec = ",")
score.qn.bmiq.df <- read.csv(paste0(src.data.pre, score.qn.bmiq.fn), sep = ";", header = T, dec = ",")

score.illig.df <- data.frame(cbind(Individual = score.bmiq.df$Individual, 
                                   Score_BMIQ = score.bmiq.df$smokingScoreIllig, 
                                   Score_QN = score.qn.df$smokingScoreIllig, 
                                   Score_QN_BMIQ = score.qn.bmiq.df$smokingScoreIllig)) %>% 
                  distinct()

score.illig.df[, 2:4] <- apply(score.illig.df[, 2:4], 2, as.numeric)
summary(score.illig.df)

score.df <- score.illig.df %>% melt()

ggplot(score.df, aes(x = value)) + 
  geom_density(aes(group = variable, colour = variable))

ggplot(score.df, aes(x = Individual , y = value,group = variable)) + 
  geom_line(aes(colour = variable)) + 
  theme(legend.position = "bottom", axis.text.x=element_blank())

# Analyze predicted status

smoking.status.df <- data.frame(cbind(Individual = score.bmiq.df$Individual, 
                              Score_BMIQ = score.bmiq.df$PredictedSmokingStatus, 
                              Score_QN = score.qn.df$PredictedSmokingStatus, 
                              Score_QN_BMIQ = score.qn.bmiq.df$PredictedSmokingStatus)) %>% 
                     distinct()

apply(smoking.status.df[, 2:4], 2, table)
summary(smoking.status.df)
