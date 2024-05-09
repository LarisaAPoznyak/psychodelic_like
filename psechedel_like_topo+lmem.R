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




path<- "D:/hse/psychodelic_like_experience/data_processing/psd_df/"


read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}

df <-list.files(path,pattern = "*.csv", full.names = T) %>% 
  map_df(~read_plus(.))
#write.csv(df, "D:/hse/psychodelic_like_experience/data_processing/df_full_0905.csv")
#df <- read.csv( "D:/hse/psychodelic_like_experience/data_processing/df_full.csv", colClasses = c( "factor", "factor", "factor", "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric", "factor" ))
head(df, 3)
df <- df %>%
  filter(freq <= 30) %>%
  mutate(freq = case_when(
    freq >= 0 & freq < 4 ~ "delta",
    freq >= 4 & freq < 8 ~ "theta",
    freq >= 8 & freq < 12 ~ "alpha",
    freq >= 12 & freq <= 30 ~ "beta",
    TRUE ~ as.character(freq) # Keep unchanged if none of the conditions match
  ))

df$subject <- substr(df$filename, start = 59, stop = 62)
colnames(df)

df1 <- df[, !colnames(df) %in% c( '...1', 'filename')]
melted_data <- melt(df1, id=c("condition", "epoch", "freq",  "subject"))
colnames(melted_data)
head(melted_data, 3)

colnames(melted_data)[5] <- "sensor"
colnames(melted_data)[6] <- "psd"

table <- aggregate(melted_data$psd, by=list(melted_data$subject, melted_data$freq, melted_data$epoch, melted_data$sensor, melted_data$condition), FUN=function(x) c(mean=mean(x))) # spectral density

head(table, 3) #data check

colnames(table) <- c('subject', 'freq',  'epoch', 'sensor', 'condition', 'psd')

table <- subset(table, subject != "YT50" ) #exclude sleepy participant

#table$psd_db <- 10 * log10(table$psd) #convert to dB

z_scores <- table %>%
  group_by(subject,freq) %>% ###########???
  mutate(value = scale(psd))

z_scores <- separate(z_scores, condition, into = c("type", "color"), sep = "/") #separate type and color
colnames(df)



z_sc <- z_scores[, !colnames(z_scores) %in% c( 'psd')]
result2 <- spread(z_sc, key = sensor, value = value)
##############################################################################
####################################
######################
numeric_cols <- sapply(result2, is.numeric)

result3 <- aggregate(result2[, numeric_cols], 
                     by = list(result2$freq, result2$type), 
                     FUN=function(x) c(mean=mean(x)))
colnames(result3)[2]<-'type'
colnames(result3)[1]<-'freq'

#result2$type<- as.factor(resul2$type)
#levels(result1$type)
CC <- subset(result3, type == "CubesControl" )
CC <- CC[, !colnames(CC) %in% c( 'epoch', 'type')]
CC  <- as.data.frame(t(CC))
colnames(CC) <- CC[1, ]
CC <- CC[-1, ]
CC <- CC[, c("delta", "theta", "alpha", "beta")]


HC<- subset(result3, type == "HoneyComb" )
HC <- HC[, !colnames(HC) %in% c( 'epoch', 'type')]
HC  <- as.data.frame(t(HC))
colnames(HC) <- HC[1, ]
HC <- HC[-1, ]
HC <- HC[, c("delta", "theta", "alpha", "beta")]

Fr<- subset(result3, type == "Fractal" )
Fr <- Fr[, !colnames(Fr) %in% c( 'epoch', 'type')]
Fr  <- as.data.frame(t(Fr))
colnames(Fr) <- Fr[1, ]
Fr <- Fr[-1, ]
Fr <- Fr[, c("delta", "theta", "alpha", "beta")]

k<- subset(result3, type == "kaleidoscope" )
k <- k[, !colnames(k) %in% c( 'epoch', 'type')]
k  <- as.data.frame(t(k))
colnames(k) <- k[1, ]
k <- k[-1, ]
k <- k[, c("delta", "theta", "alpha", "beta")]


write.csv(CC, "D:/hse/psychodelic_like_experience/data_processing/CC_topo_0905_.csv")
write.csv(HC, "D:/hse/psychodelic_like_experience/data_processing/HC_topo_0905_.csv")
write.csv(Fr, "D:/hse/psychodelic_like_experience/data_processing/Fr_topo_0905_.csv")
write.csv(k, "D:/hse/psychodelic_like_experience/data_processing/k_topo_0905_.csv")



#colors decoding
data2 <- z_scores %>%
  mutate(color = case_when(
    color == 'GreenPurple'  ~ "PurpleGreen",
    color == 'BlueRed'  ~ "RedBlue",
    color == 'PinkYellow'  ~ "OrangePink",
    color == 'PinkBlue'  ~ "BluePink",
    TRUE ~ as.character(color) # Keep unchanged if none of the conditions match
  ))
colnames(data2)[8]<-'z_score'

#epoch decoding
data2 <- data2 %>%
  mutate(epoch = case_when(
    epoch >= 0 & epoch < 4 ~ "b",
    epoch >= 4 & epoch < 15 ~ "m",
    epoch >= 15 & epoch <= 19 ~ "e",
    TRUE ~ as.character(epoch) # Keep unchanged if none of the conditions match
  ))
write.csv(data2, "D:/hse/psychodelic_like_experience/data_processing/new_data_0905__.csv")
data2 <- read.csv("D:/hse/psychodelic_like_experience/data_processing/new_data_0905__.csv")
data2 <- data2[, !colnames(data2) %in% c( 'X')]
#colnames(data2)[]<-'z_score'

#data2 <-melted_data
data2 <- data2 %>%
  mutate(hyper = if_else(type %in% c("Fractal", "HoneyComb"), 'hyperbolic', 'nonhyperbolic'))
#data2 <- data2[, !colnames(data2) %in% c('X')]

data3 <- subset(data2, type %in% c("CubesControl", "HoneyComb"))
data4 <- subset(data2, type %in% c("Fractal", "kaleidoscope"))
#theta$hyper <- as.factor(theta$hyper)
#theta$color <- as.factor(theta$color)
#theta$subject <- as.factor(theta$subject)

#levels(theta$subject)
delta <- subset(data3, freq == "delta" )
theta <- subset(data3, freq == "theta" )
alpha <- subset(data3, freq == "alpha" )
beta <- subset(data3, freq == "beta" )

delta <- subset(data4, freq == "delta" )
theta <- subset(data4, freq == "theta" )
alpha <- subset(data4, freq == "alpha" )
beta <- subset(data4, freq == "beta" )
levels()


sensors<- unique(data2$sensor)

p_vals <- data.table()
#cols <- c("psd")
############### for green heads (main_effects) ##############
for (i in 1:length(sensors)) {
  temp <-subset(delta, sensor == sensors[i])
  m1 <- lmer(z_score ~ hyper*color + (1|subject), data = temp)
  d<-romr.fnc(m1, temp, trim = 3) ##### remove outliers
  data<- d$data 
  m <- lmer(z_score ~ hyper*color + (1|subject), data = data)
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
write.csv(p_vals, "aov_results_delta_fk_z_sc_0905_.csv")



p_vals <- data.table()
#cols <- c("psd")
############### for green heads (main_effects) ##############
for (i in 1:length(sensors)) {
  temp <-subset(theta, sensor == sensors[i])
  m1 <- lmer(z_score ~ hyper*color + (1|subject), data = temp)
  d<-romr.fnc(m1, temp, trim = 3) ##### remove outliers
  data<- d$data 
  m <- lmer(z_score ~ hyper*color + (1|subject), data = data)
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
write.csv(p_vals, "aov_results_theta_fk_z_sc_0905_.csv")

p_vals <- data.table()
#cols <- c("psd")
############### for green heads (main_effects) ##############
for (i in 1:length(sensors)) {
  temp <-subset(alpha, sensor == sensors[i])
  m1 <- lmer(z_score ~ hyper*color + (1|subject), data = temp)
  d<-romr.fnc(m1, temp, trim = 3) ##### remove outliers
  data<- d$data 
  m <- lmer(z_score ~ hyper*color + (1|subject), data = data)
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
write.csv(p_vals, "aov_results_alpha_fk_z_sc_0905_.csv")


p_vals <- data.table()
#cols <- c("psd")
############### for green heads (main_effects) ##############
for (i in 1:length(sensors)) {
  temp <-subset(beta, sensor == sensors[i])
  m1 <- lmer(z_score ~ hyper*color + (1|subject), data = temp)
  d<-romr.fnc(m1, temp, trim = 3) ##### remove outliers
  data<- d$data 
  m <- lmer(z_score ~ hyper*color + (1|subject), data = data)
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
write.csv(p_vals, "aov_results_beta_fk_z_sc_0905_.csv")

Tuk<-data.table(summary(emmeans(m, pairwise ~ US_type|CS_type, adjust = 'tukey',lmer.df = "satterthwaite"))$contrasts)
Tuk <- Tuk[, group1:=gsub(' -.*', '', contrast)][, group2:=gsub('.*- ', '', contrast)]
Tuk <- Tuk[p.value<0.1, p_significant:=format(p.value, digits = 3)]

emm_options(lmerTest.limit = 10000)
emm_options(pbkrtest.limit = 9000)

alpha <- subset(data3, freq == "alpha" )

temp <-subset(alpha, sensor == 'Fp2')
m <- lmer(z_score ~ hyper*color + (1|subject), data = temp)
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
  geom_line(size = 2, aes(group = hyper), position = position_dodge(width = 0.2)) +
  geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.2), aes(fill = hyper)) +
  
  scale_color_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                     values=c("magenta2","olivedrab3")) +
  scale_fill_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                    values=c("white","white")) +
  
  # Facet the plot
  #facet_wrap(~set, ncol = 2) +
  
  stat_pvalue_manual(Tuk, label = 'stars', size = 7, bracket.size = 1.5, tip.length = 0.01,
                     y.position =c(max(marginal_em$emmean)+max(marginal_em$SE) + 0.1),inherit.aes = FALSE) +
  # Customize the appearance of the plot
  theme_bw() +
  labs(
    title = "Normalized power - Fp2",
    x = "",
    y = "Marginal means",
    caption = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(angle = 0, size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    legend.position = c(.98, .02),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.box.background = element_rect(color="black", size=1),
    
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 15),
    plot.caption = element_text(hjust = 0, size = 10)
  )
plot <- ggpar(plot,
            ylim = c(-1, 1)
            )
# Print the plot
print(plot)
print(sensors)


temp <-subset(alpha, sensor  %in% c( 'F8', 'FC6', 'T8'))
m <- lmer(z_score ~ hyper*color + (1|subject), data = temp)
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
  geom_line(size = 2, aes(group = hyper), position = position_dodge(width = 0.2)) +
  geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.2), aes(fill = hyper)) +
  
  scale_color_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                     values=c("magenta2","olivedrab3")) +
  scale_fill_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                    values=c("white","white")) +
  
  # Facet the plot
  #facet_wrap(~set, ncol = 2) +
  
  stat_pvalue_manual(Tuk, label = 'stars', size = 7, bracket.size = 1.5, tip.length = 0.01,
                     y.position =c(max(marginal_em$emmean)+max(marginal_em$SE) + 0.1),inherit.aes = FALSE) +
  # Customize the appearance of the plot
  theme_bw() +
  labs(
    title = "Normalized power - Temporal group",
    x = "",
    y = "Marginal means",
    caption = ""
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
    plot.caption = element_text(hjust = 0, size = 10)
  )
plot <- ggpar(plot,
              ylim = c(-1, 1)
)
# Print the plot
print(plot)



temp <-subset(alpha, sensor  %in% c( 'CP2', 'CP2', 'Pz', 'Cz', 'FC2'))
m <- lmer(z_score ~ hyper*color + (1|subject), data = temp)
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
  geom_line(size = 2, aes(group = hyper), position = position_dodge(width = 0.2)) +
  geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.2), aes(fill = hyper)) +
  
  scale_color_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                     values=c("magenta2","olivedrab3")) +
  scale_fill_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                    values=c("white","white")) +
  
  # Facet the plot
  #facet_wrap(~set, ncol = 2) +
  
  stat_pvalue_manual(Tuk, label = 'stars', size = 7, bracket.size = 1.5, tip.length = 0.01,
                     y.position =c(max(marginal_em$emmean)+max(marginal_em$SE) + 0.1),inherit.aes = FALSE) +
  # Customize the appearance of the plot
  theme_bw() +
  labs(
    title = "Normalized power - Central group",
    x = "",
    y = "Marginal means",
    caption = ""
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
    plot.caption = element_text(hjust = 0, size = 10)
  )
plot <- ggpar(plot,
              ylim = c(-1, 1)
)
# Print the plot
print(plot)


delta <- subset(data4, freq == "delta" )
pvalsensor <- subset(p_vals, rn == 'hyper' & p_value <0.05)
pvalsensor<- unique(pvalsensor$sensor)
print(pvalsensor)
# Fp1 Fz  F7  FC1 C3  CP1 Pz  P3  P7  P4  P8  CP6 CP2 Cz  C4  T8  FC2 F8  Fp2
temp <-subset(delta, sensor  %in% c('Fz', 'Fp1', 'Fp2', 'F7', 'F8', 'FC1', 'FC2'))
m <- lmer(z_score ~ hyper*color + (1|subject), data = temp)
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
  geom_line(size = 2, aes(group = hyper), position = position_dodge(width = 0.2)) +
  geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.2), aes(fill = hyper)) +
  
  scale_color_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                     values=c("magenta2","olivedrab3")) +
  scale_fill_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                    values=c("white","white")) +
  
  # Facet the plot
  #facet_wrap(~set, ncol = 2) +
  
  stat_pvalue_manual(Tuk, label = 'stars', size = 7, bracket.size = 1.5, tip.length = 0.01,
                     y.position =c(max(marginal_em$emmean)+max(marginal_em$SE) + 0.1),inherit.aes = FALSE) +
  # Customize the appearance of the plot
  theme_bw() +
  labs(
    title = "Normalized power - Frontal group",
    x = "",
    y = "Marginal means",
    caption = ""
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
    plot.caption = element_text(hjust = 0, size = 10)
  )
plot <- ggpar(plot,
              ylim = c(-1, 1)
)
# Print the plot
print(plot)


# Fp1 Fz  F7  FC1 C3  CP1 Pz  P3  P7  P4  P8  CP6 CP2 Cz  C4  T8  FC2 F8  Fp2
temp <-subset(delta, sensor  %in% c('С3','CP1', 'CP6', 'CP2', 'Cz', 'C4'))
m <- lmer(z_score ~ hyper*color + (1|subject), data = temp)
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
  geom_line(size = 2, aes(group = hyper), position = position_dodge(width = 0.2)) +
  geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.2), aes(fill = hyper)) +
  
  scale_color_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                     values=c("magenta2","olivedrab3")) +
  scale_fill_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                    values=c("white","white")) +
  
  # Facet the plot
  #facet_wrap(~set, ncol = 2) +
  
  stat_pvalue_manual(Tuk, label = 'stars', size = 7, bracket.size = 1.5, tip.length = 0.01,
                     y.position =c(max(marginal_em$emmean)+max(marginal_em$SE) + 0.1),inherit.aes = FALSE) +
  # Customize the appearance of the plot
  theme_bw() +
  labs(
    title = "Normalized power - Central group",
    x = "",
    y = "Marginal means",
    caption = ""
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
    plot.caption = element_text(hjust = 0, size = 10)
  )
plot <- ggpar(plot,
              ylim = c(-1, 1)
)
# Print the plot
print(plot)


# Fp1 Fz  F7  FC1 C3  CP1 Pz  P3  P7  P4  P8  CP6 CP2 Cz  C4  T8  FC2 F8  Fp2
temp <-subset(delta, sensor  %in% c('Pz',  'P3',  'P7',  'P4',  'P8'))
m <- lmer(z_score ~ hyper*color + (1|subject), data = temp)
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
  geom_line(size = 2, aes(group = hyper), position = position_dodge(width = 0.2)) +
  geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.2), aes(fill = hyper)) +
  
  scale_color_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                     values=c("magenta2","olivedrab3")) +
  scale_fill_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                    values=c("white","white")) +
  
  # Facet the plot
  #facet_wrap(~set, ncol = 2) +
  
  stat_pvalue_manual(Tuk, label = 'stars', size = 7, bracket.size = 1.5, tip.length = 0.01,
                     y.position =c(max(marginal_em$emmean)+max(marginal_em$SE) + 0.1),inherit.aes = FALSE) +
  # Customize the appearance of the plot
  theme_bw() +
  labs(
    title = "Normalized power - Parietal group",
    x = "",
    y = "Marginal means",
    caption = ""
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
    plot.caption = element_text(hjust = 0, size = 10)
  )
plot <- ggpar(plot,
              ylim = c(-1, 1)
)
# Print the plot
print(plot)






# Fp1 Fz  F7  FC1 C3  CP1 Pz  P3  P7  P4  P8  CP6 CP2 Cz  C4  T8  FC2 F8  Fp2
temp <-subset(delta, sensor  %in% c('T8'))
m <- lmer(z_score ~ hyper*color + (1|subject), data = temp)
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
  geom_line(size = 2, aes(group = hyper), position = position_dodge(width = 0.2)) +
  geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.2), aes(fill = hyper)) +
  
  scale_color_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                     values=c("magenta2","olivedrab3")) +
  scale_fill_manual(breaks = c('hyperbolic', 'nonhyperbolic'),
                    values=c("white","white")) +
  
  # Facet the plot
  #facet_wrap(~set, ncol = 2) +
  
  stat_pvalue_manual(Tuk, label = 'stars', size = 7, bracket.size = 1.5, tip.length = 0.01,
                     y.position =c(max(marginal_em$emmean)+max(marginal_em$SE) + 0.1),inherit.aes = FALSE) +
  # Customize the appearance of the plot
  theme_bw() +
  labs(
    title = "Normalized power - Temporal group",
    x = "",
    y = "Marginal means",
    caption = ""
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
    plot.caption = element_text(hjust = 0, size = 10)
  )
plot <- ggpar(plot,
              ylim = c(-1, 1)
)
# Print the plot
print(plot)
