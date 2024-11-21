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
data<- read.csv("D:/hse/psychodelic_like_experience/data_processing/psd/ap_analysis.csv")
data1 <- subset(data, er < 0.11)

data2 <- separate(data1, condition, into = c("type", "color"), sep = "/")
table1 <- data2[, !colnames(data2) %in% c('X1', 'X')]
#numeric_cols <- sapply(table1, is.numeric)
table1 <- subset(table1, a2 != 'NaN')

result3 <- aggregate(table1$z_score, 
                     by = list(table1$type), 
                     FUN=function(x) c(mean=mean(x), std=std.error(x)))
colnames(result3) = c('type', 'x')
result3 <- result3 %>%
  mutate(hyper = if_else(type %in% c("Fractal", "HoneyComb"), 'hyperbolic', 'nonhyperbolic'))%>%
  mutate(group = if_else(type %in% c("Fractal", "kaleidoscope"), 'Fractal/kaleidoscope', 'HoneyComb/CubesControl'))

result3$type <- factor(result3$type, levels = c("CubesControl",    "HoneyComb", "Fractal" ,     "kaleidoscope")) 

plot <- ggplot(result3, aes(x = group, y = x[,"mean"], color = hyper)) +
  geom_point(size = 3, position = position_dodge(width = 0.2), shape = 23) +  
  scale_color_manual(values = c("hyperbolic" = "magenta2", "nonhyperbolic" = "olivedrab3")) +
  geom_errorbar(aes(ymin = x[,"mean"] - x[,"std"], ymax = x[,"mean"] + x[,"std"], color = factor(hyper, levels = c("hyperbolic", "nonhyperbolic"))), width = 0.2, position = position_dodge(width = 0.2), size = 2) +
  #facet_grid( ~ hyper) +  # Facet by group
  labs(title = "Periodic component of alpha-band", x = "Type", y = "Scale mean of peak power ", caption = 'Errorbars represent SE') +  # Labels
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),  # Adjust x-axis title size
    axis.title.y = element_text(size = 16),  # Adjust y-axis title size
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    #panel.grid.major.x = element_blank(),
    legend.position = "bottom"  # Adjust as needed
  )

# Show the plot
print(plot)

table1 <- table1 %>%
  group_by(subj) %>% ###########???
  mutate(z_score = scale(a2))

table1 <- table1 %>%
  mutate(hyper = if_else(type %in% c("Fractal", "HoneyComb"), 'hyperbolic', 'nonhyperbolic'))

data3 <- subset(table1, type %in% c("CubesControl", "HoneyComb"))
data4 <- subset(table1, type %in% c("Fractal", "kaleidoscope"))

m1 <- lmer(z_score ~ hyper + (1|subj), data = data3)
d<-romr.fnc(m1, data3, trim = 3) ##### remove outliers
data<- d$data 
m <- lmer(z_score ~ hyper + (1|subj), data = data3)
an <- anova(m)
an <- data.table(an,keep.rownames = TRUE)
an_cols <- c('rn','Pr(>F)') 
an <- an[, ..an_cols]
an$p_value <- format(an$`Pr(>F)`, digits = 3)


data<- read.csv("D:/hse/psychodelic_like_experience/data_processing/psd/ap_analysis.csv")
data1 <- subset(data, er < 0.11)

data2 <- separate(data1, condition, into = c("type", "color"), sep = "/")
table1 <- data2[, !colnames(data2) %in% c('X1', 'X')]

table1 <- subset(table1, ap3 != 'NaN') #exponent

result3 <- aggregate(table1$ap3, 
                     by = list(table1$type), 
                     FUN=function(x) c(mean=mean(x), std=std.error(x)))
colnames(result3) = c('type', 'x')
result3 <- result3 %>%
  mutate(hyper = if_else(type %in% c("Fractal", "HoneyComb"), 'hyperbolic', 'nonhyperbolic'))%>%
  mutate(group = if_else(type %in% c("Fractal", "kaleidoscope"), 'Fractal/kaleidoscope', 'HoneyComb/CubesControl'))

result3$type <- factor(result3$type, levels = c("CubesControl",    "HoneyComb", "Fractal" ,     "kaleidoscope")) 

plot <- ggplot(result3, aes(x = group, y = x[,"mean"], color = hyper)) +
  geom_point(size = 3, position = position_dodge(width = 0.2), shape = 23) +  
  scale_color_manual(values = c("hyperbolic" = "magenta2", "nonhyperbolic" = "olivedrab3")) +
  geom_errorbar(aes(ymin = x[,"mean"] - x[,"std"], ymax = x[,"mean"] + x[,"std"], color = factor(hyper, levels = c("hyperbolic", "nonhyperbolic"))), width = 0.2, position = position_dodge(width = 0.2), size = 2) +
  #facet_grid( ~ hyper) +  # Facet by group
  labs(title = "Aperiodic component (exponent)", x = "Type", y = "Mean", caption = 'Errorbars represent SE') +  # Labels
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),  # Adjust x-axis title size
    axis.title.y = element_text(size = 16),  # Adjust y-axis title size
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    #panel.grid.major.x = element_blank(),
    legend.position = "bottom"  # Adjust as needed
  )

# Show the plot
print(plot)




#################AP_on_sensors
data<- read.csv("D:/hse/psychodelic_like_experience/data_processing/psd/ap_analysis_sensors.csv")
data1 <- subset(data, er < 0.15)

data2 <- separate(data1, condition, into = c("type", "color"), sep = "/")
table1 <- data2[, !colnames(data2) %in% c('X1', 'X')]
#numeric_cols <- sapply(table1, is.numeric)
table1 <- subset(table1, a2 != 'NaN')

result3 <- aggregate(table1$a2, 
                     by = list(table1$type, table1$sensor ), 
                     FUN=function(x) c(mean=mean(x)))
colnames(result3) = c('type','sensor', 'x')


CC <- subset(result3, type == "CubesControl" )
CC <- CC[, !colnames(CC) %in% c('type')]
CC  <- t(CC)
colnames(CC) <- CC[1, ]
CC <- CC[-1, ]


HC<- subset(result3, type == "HoneyComb" )
HC <- HC[, !colnames(HC) %in% c( 'type')]
HC  <- t(HC)
colnames(HC) <- HC[1, ]
HC <- HC[-1, ]


Fr<- subset(result3, type == "Fractal" )
Fr <- Fr[, !colnames(Fr) %in% c( 'type')]
Fr  <- t(Fr)
colnames(Fr) <- Fr[1, ]
Fr <- Fr[-1, ]


k<- subset(result3, type == "kaleidoscope" )
k <- k[, !colnames(k) %in% c('type')]
k  <- t(k)
colnames(k) <- k[1, ]
k <- k[-1, ]


numeric_cols <- sapply(HC, is.numeric)

#HCCC <- HC[, numeric_cols] - CC[, numeric_cols]
#Fk <- Fr[, numeric_cols] - k[, numeric_cols]

#HCCC$freq <- CC$freq
#Fk$freq <- k$freq
write.csv(CC, "D:/hse/psychodelic_like_experience/data_processing/CC_alpha.csv")
write.csv(HC, "D:/hse/psychodelic_like_experience/data_processing/HC_alpha.csv")
write.csv(Fr, "D:/hse/psychodelic_like_experience/data_processing/Fr_alpha.csv")
write.csv(k, "D:/hse/psychodelic_like_experience/data_processing/k_alpha.csv")








result3 <- result3 %>%
  mutate(hyper = if_else(type %in% c("Fractal", "HoneyComb"), 'hyperbolic', 'nonhyperbolic'))%>%
  mutate(group = if_else(type %in% c("Fractal", "kaleidoscope"), 'Fractal/kaleidoscope', 'HoneyComb/CubesControl'))

result3$type <- factor(result3$type, levels = c("CubesControl",    "HoneyComb", "Fractal" ,     "kaleidoscope")) 

plot <- ggplot(result3, aes(x = group, y = x[,"mean"], color = hyper)) +
  geom_point(size = 3, position = position_dodge(width = 0.2), shape = 23) +  
  scale_color_manual(values = c("hyperbolic" = "magenta2", "nonhyperbolic" = "olivedrab3")) +
  geom_errorbar(aes(ymin = x[,"mean"] - x[,"std"], ymax = x[,"mean"] + x[,"std"], color = factor(hyper, levels = c("hyperbolic", "nonhyperbolic"))), width = 0.2, position = position_dodge(width = 0.2), size = 2) +
  #facet_grid( ~ hyper) +  # Facet by group
  labs(title = "Periodic component of alpha-band", x = "Type", y = "Mean of peak power", caption = 'Errorbars represent SE') +  # Labels
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),  # Adjust x-axis title size
    axis.title.y = element_text(size = 16),  # Adjust y-axis title size
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    #panel.grid.major.x = element_blank(),
    legend.position = "bottom"  # Adjust as needed
  )

# Show the plot
print(plot)

table1 <- table1 %>%
  mutate(color = case_when(
    color == 'GreenPurple'  ~ "PurpleGreen",
    color == 'BlueRed'  ~ "RedBlue",
    color == 'PinkYellow'  ~ "OrangePink",
    color == 'PinkBlue'  ~ "BluePink",
    TRUE ~ as.character(color) # Keep unchanged if none of the conditions match
  ))
sensors<- unique(data2$sensor)
table1 <- table1 %>%
  group_by(subj) %>% ###########???
  mutate(z_score = scale(a2))

table1 <- table1 %>%
  mutate(hyper = if_else(type %in% c("Fractal", "HoneyComb"), 'hyperbolic', 'nonhyperbolic'))

data3 <- subset(table1, type %in% c("CubesControl", "HoneyComb"))
data4 <- subset(table1, type %in% c("Fractal", "kaleidoscope"))


p_vals <- data.table()
#cols <- c("psd")
############### for green heads (main_effects) ##############
for (i in 1:length(sensors)) {
  temp <-subset(data3, sensor == sensors[i])
  m1 <- lmer(z_score ~ hyper + (1|subj), data = temp)
  d<-romr.fnc(m1, temp, trim = 3) ##### remove outliers
  data<- d$data 
  m <- lmer(z_score ~ hyper + (1|subj), data = data)
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

data<- read.csv("D:/hse/psychodelic_like_experience/data_processing/psd/ap_analysis_sensors.csv")
data1 <- subset(data, er < 0.11)

data2 <- separate(data1, condition, into = c("type", "color"), sep = "/")
table1 <- data2[, !colnames(data2) %in% c('X1', 'X')]

table1 <- subset(table1, ap3 != 'NaN') #exponent

result3 <- aggregate(table1$ap3, 
                     by = list(table1$type), 
                     FUN=function(x) c(mean=mean(x), std=std.error(x)))
colnames(result3) = c('type', 'x')
result3 <- result3 %>%
  mutate(hyper = if_else(type %in% c("Fractal", "HoneyComb"), 'hyperbolic', 'nonhyperbolic'))%>%
  mutate(group = if_else(type %in% c("Fractal", "kaleidoscope"), 'Fractal/kaleidoscope', 'HoneyComb/CubesControl'))

result3$type <- factor(result3$type, levels = c("CubesControl",    "HoneyComb", "Fractal" ,     "kaleidoscope")) 

plot <- ggplot(result3, aes(x = group, y = x[,"mean"], color = hyper)) +
  geom_point(size = 3, position = position_dodge(width = 0.2), shape = 23) +  
  scale_color_manual(values = c("hyperbolic" = "magenta2", "nonhyperbolic" = "olivedrab3")) +
  geom_errorbar(aes(ymin = x[,"mean"] - x[,"std"], ymax = x[,"mean"] + x[,"std"], color = factor(hyper, levels = c("hyperbolic", "nonhyperbolic"))), width = 0.2, position = position_dodge(width = 0.2), size = 2) +
  #facet_grid( ~ hyper) +  # Facet by group
  labs(title = "Aperiodic component (exponent)", x = "Type", y = "Mean", caption = 'Errorbars represent SE') +  # Labels
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),  # Adjust x-axis title size
    axis.title.y = element_text(size = 16),  # Adjust y-axis title size
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    #panel.grid.major.x = element_blank(),
    legend.position = "bottom"  # Adjust as needed
  )

# Show the plot
print(plot)

table1 <- table1 %>%
  mutate(color = case_when(
    color == 'GreenPurple'  ~ "PurpleGreen",
    color == 'BlueRed'  ~ "RedBlue",
    color == 'PinkYellow'  ~ "OrangePink",
    color == 'PinkBlue'  ~ "BluePink",
    TRUE ~ as.character(color) # Keep unchanged if none of the conditions match
  ))
sensors<- unique(table1$sensor)
table1 <- table1 %>%
  group_by(subj) %>% ###########???
  mutate(z_score = scale(ap3))

table1 <- table1 %>%
  mutate(hyper = if_else(type %in% c("Fractal", "HoneyComb"), 'hyperbolic', 'nonhyperbolic'))

data3 <- subset(table1, type %in% c("CubesControl", "HoneyComb"))
data4 <- subset(table1, type %in% c("Fractal", "kaleidoscope"))


p_vals <- data.table()
#cols <- c("psd")
############### for green heads (main_effects) ##############
for (i in 1:length(sensors)) {
  temp <-subset(data4, sensor == sensors[i])
  m1 <- lmer(z_score ~ hyper + (1|subj), data = temp)
  d<-romr.fnc(m1, temp, trim = 3) ##### remove outliers
  data<- d$data 
  m <- lmer(z_score ~ hyper + (1|subj), data = data)
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

