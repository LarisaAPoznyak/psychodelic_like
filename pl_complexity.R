library(reshape2)
library(data.table)
library(ggplot2)
library(lme4)
# library("ggpubr")
library(emmeans)
library(lmerTest)
library(stringi)
library(stringr)
library(dplyr)
library(purrr)
library(tidyverse)
library(LMERConvenienceFunctions)
library(rstatix)
library(plotrix)
library(ggpubr)
#############################veg


df <- read.csv( "D:/hse/psychodelic_like_experience/data_processing/complexity/c_baselined_alpha.csv")

df <- subset(df, subject != "YT50" ) #exclude sleepy participant

df <- df %>%
  mutate(hyper = if_else(cond %in% c("Fractal", "HoneyComb"), 'hyperbolic', 'nonhyperbolic'))
#data2 <- data2[, !colnames(data2) %in% c('X')]

data3 <- subset(df, cond %in% c("CubesControl", "HoneyComb"))
data4 <- subset(df, cond %in% c("Fractal", "kaleidoscope"))

#

sensors<- unique(df$sensor)

p_vals <- data.table()
#cols <- c("psd")
############### for green heads (main_effects) ##############
for (i in 1:length(sensors)) {
  temp <-subset(data3, sensor == sensors[i])
  m1 <- lmer(comp_bl ~ hyper + (1|subject), data = temp)
  d<-romr.fnc(m1, temp, trim = 3) ##### remove outliers
  data<- d$data 
  m <- lmer(comp_bl ~ hyper + (1|subject), data = data)
  an <- anova(m)
  an <- data.table(an,keep.rownames = TRUE)
  an_cols <- c('rn','Pr(>F)') 
  an <- an[, ..an_cols]
  an$p_value <- format(an$`Pr(>F)`, digits = 3)
  #an$interval <- j
  #an$interval <- gsub('beta power','',an$interval)
  #an <- dcast(an,formula = interval~rn,value.var = 'Pr(>F)')
  an$sensor <- sensors[i] 
  #an$sensor_name <- files[sensor==i]$Name
  p_vals <- rbind(p_vals,an)
  
}
setwd('D:/hse/psychodelic_like_experience/data_processing/stats/')
write.csv(p_vals, "aov_complexity_HCCC_alpha.csv")




p_vals <- data.table()
#cols <- c("psd")
############### for green heads (main_effects) ##############
for (i in 1:length(sensors)) {
  temp <-subset(data4, sensor == sensors[i])
  m1 <- lmer(comp_bl ~ hyper + (1|subject), data = temp)
  d<-romr.fnc(m1, temp, trim = 3) ##### remove outliers
  data<- d$data 
  m <- lmer(comp_bl ~ hyper + (1|subject), data = data)
  an <- anova(m)
  an <- data.table(an,keep.rownames = TRUE)
  an_cols <- c('rn','Pr(>F)') 
  an <- an[, ..an_cols]
  an$p_value <- format(an$`Pr(>F)`, digits = 3)
  #an$interval <- j
  #an$interval <- gsub('beta power','',an$interval)
  #an <- dcast(an,formula = interval~rn,value.var = 'Pr(>F)')
  an$sensor <- sensors[i] 
  #an$sensor_name <- files[sensor==i]$Name
  p_vals <- rbind(p_vals,an)
  
}
setwd('D:/hse/psychodelic_like_experience/data_processing/stats/')
write.csv(p_vals, "aov_complexity_fk_alpha.csv")


pvalsensor <- subset(p_vals, rn == 'hyper' & p_value <0.05)

temp <-subset(data3, sensor  %in% c( 'CP1'))
m <- lmer(comp_bl ~ hyper + (1|subject), data = temp)
marginal_em <- emmeans(m, ~ as.factor(hyper),level = 0.95,lmer.df = "satterthwaite")
marginal_em<-as.data.frame(marginal_em)
Tuk<-data.table(summary(emmeans(m, pairwise ~ hyper, adjust = 'tukey',lmer.df = "satterthwaite"))$contrasts)
Tuk <- Tuk[, group1:=gsub(' -.*', '', contrast)][, group2:=gsub('.*- ', '', contrast)]
Tuk <- Tuk[p.value<0.1, p_significant:=format(p.value, digits = 3)]

n <- Tuk[!is.na(p_significant), .N]

Tuk[p.value<0.001, stars:='***']
Tuk[p.value<0.01 & p.value>0.001 , stars:='**']
Tuk[p.value<0.05 & p.value>0.01 , stars:='*']
Tuk[p.value>0.05 & p.value<0.1 , stars:='#']


plot <- ggplot(marginal_em, aes(x = hyper, y = emmean, group = hyper, color = hyper, fill = hyper)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.5, position = position_dodge(width = 0.2), size = 2) +
  #geom_bar(binwidth = 5,aes(group = hyper, fill = hyper), stat = "identity")+
  #geom_line(size = 2, aes(group = hyper), position = position_dodge(width = 0.2)) +
  geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.2), aes(fill = hyper)) +
  
  scale_color_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                     values=c("magenta2","olivedrab3")) +
  scale_fill_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                    values=c("white","white")) +
  
  # Facet the plot
  #facet_wrap(~set, ncol = 2) +
  
  stat_pvalue_manual(Tuk, label = 'stars', size = 7, bracket.size = 1.5, tip.length = 0.01,
                     y.position =c(max(marginal_em$emmean)+max(marginal_em$SE) + 0.002),inherit.aes = FALSE) +
  # Customize the appearance of the plot
  theme_bw() +
  labs(
    title = "Normalized complexity - 'Fz', 'F8', 'FC2'",
    x = "",
    y = "Marginal means",
    caption = "CubesControl vs HoneyComb"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(angle = 0, size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    legend.position = c(.98, .98),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.box.background = element_rect(color="black", size=1),
    
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 15),
    plot.caption = element_text(hjust = 0, size = 12)
  )
plot <- ggpar(plot,
              ylim = c(-0.002, 0.015)
)
# Print the plot
print(plot)

temp <-subset(data3, sensor  %in% c( 'P7', 'P4', 'O1', 'Oz'))
m <- lmer(comp_bl ~ hyper + (1|subject), data = temp)
marginal_em <- emmeans(m, ~ as.factor(hyper),level = 0.95,lmer.df = "satterthwaite")
marginal_em<-as.data.frame(marginal_em)
Tuk<-data.table(summary(emmeans(m, pairwise ~ hyper, adjust = 'tukey',lmer.df = "satterthwaite"))$contrasts)
Tuk <- Tuk[, group1:=gsub(' -.*', '', contrast)][, group2:=gsub('.*- ', '', contrast)]
Tuk <- Tuk[p.value<0.1, p_significant:=format(p.value, digits = 3)]

n <- Tuk[!is.na(p_significant), .N]

Tuk[p.value<0.001, stars:='***']
Tuk[p.value<0.01 & p.value>0.001 , stars:='**']
Tuk[p.value<0.05 & p.value>0.01 , stars:='*']
Tuk[p.value>0.05 & p.value<0.1 , stars:='#']


plot <- ggplot(marginal_em, aes(x = hyper, y = emmean, group = hyper, color = hyper, fill = hyper)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.5, position = position_dodge(width = 0.2), size = 2) +
  #geom_bar(binwidth = 5,aes(group = hyper, fill = hyper), stat = "identity")+
  #geom_line(size = 2, aes(group = hyper), position = position_dodge(width = 0.2)) +
  geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.2), aes(fill = hyper)) +
  
  scale_color_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                     values=c("magenta2","olivedrab3")) +
  scale_fill_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                    values=c("white","white")) +
  
  # Facet the plot
  #facet_wrap(~set, ncol = 2) +
  
  stat_pvalue_manual(Tuk, label = 'stars', size = 7, bracket.size = 1.5, tip.length = 0.01,
                     y.position =c(max(marginal_em$emmean)+max(marginal_em$SE) + 0.002),inherit.aes = FALSE) +
  # Customize the appearance of the plot
  theme_bw() +
  labs(
    title = "Normalized complexity - 'P7', 'P4', 'O1', 'Oz'",
    x = "",
    y = "Marginal means",
    caption = "CubesControl vs HoneyComb"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(angle = 0, size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    legend.position = c(.98, .98),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.box.background = element_rect(color="black", size=1),
    
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 15),
    plot.caption = element_text(hjust = 0, size = 12)
  )
plot <- ggpar(plot,
              ylim = c(0, 0.025)
)
# Print the plot
print(plot)



temp <-subset(data4, sensor  %in% c( 'Fp1'))
m <- lmer(comp_bl ~ hyper + (1|subject), data = temp)
marginal_em <- emmeans(m, ~ as.factor(hyper),level = 0.95,lmer.df = "satterthwaite")
marginal_em<-as.data.frame(marginal_em)
Tuk<-data.table(summary(emmeans(m, pairwise ~ hyper, adjust = 'tukey',lmer.df = "satterthwaite"))$contrasts)
Tuk <- Tuk[, group1:=gsub(' -.*', '', contrast)][, group2:=gsub('.*- ', '', contrast)]
Tuk <- Tuk[p.value<0.1, p_significant:=format(p.value, digits = 3)]

n <- Tuk[!is.na(p_significant), .N]

Tuk[p.value<0.001, stars:='***']
Tuk[p.value<0.01 & p.value>0.001 , stars:='**']
Tuk[p.value<0.05 & p.value>0.01 , stars:='*']
Tuk[p.value>0.05 & p.value<0.1 , stars:='#']


plot <- ggplot(marginal_em, aes(x = hyper, y = emmean, group = hyper, color = hyper, fill = hyper)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.5, position = position_dodge(width = 0.2), size = 2) +
  #geom_bar(binwidth = 5,aes(group = hyper, fill = hyper), stat = "identity")+
  #geom_line(size = 2, aes(group = hyper), position = position_dodge(width = 0.2)) +
  geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.2), aes(fill = hyper)) +
  
  scale_color_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                     values=c("magenta2","olivedrab3")) +
  scale_fill_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                    values=c("white","white")) +
  
  # Facet the plot
  #facet_wrap(~set, ncol = 2) +
  
  stat_pvalue_manual(Tuk, label = 'stars', size = 7, bracket.size = 1.5, tip.length = 0.01,
                     y.position =c(max(marginal_em$emmean)+max(marginal_em$SE) + 0.002),inherit.aes = FALSE) +
  # Customize the appearance of the plot
  theme_bw() +
  labs(
    title = "Normalized complexity - 'Fp1'",
    x = "",
    y = "Marginal means",
    caption = "Fractal vs kaleidoscope"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(angle = 0, size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    legend.position = c(.98, .98),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.box.background = element_rect(color="black", size=1),
    
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 15),
    plot.caption = element_text(hjust = 0, size = 12)
  )
#plot <- ggpar(plot,
 #             ylim = c(-0.005, 0.015)
#)
# Print the plot
print(plot)

