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
#############################veg
path<- "D:/hse/psychodelic_like_experience/data_processing/pletism_base/"
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
  mutate(filename = flnm)
}
df_kgr <-list.files(path,pattern = "*.csv", full.names = T) %>% 
  map_df(~read_plus(.))

replace_values <- function(x) {
  x <- ifelse(substr(x, 1, 1) == "1", "Fractal", 
              ifelse(substr(x, 1, 1) == "2", "kaleidoscope",
                     ifelse(substr(x, 1, 1) == "3", "CubesControl",
                            ifelse(substr(x, 1, 1) == "4", "HoneyComb", x))))
  return(x)
}

# Apply the function to the values
df_kgr$Condition <- replace_values(df_kgr$Condition)
df_kgr <- df_kgr %>%
  mutate(color = case_when(
    Label == 101 ~ "Blue",
    Label == 102 ~ "Purple",
    Label == 103 ~ "Green",
    Label == 104 ~ "Red",
    Label == 105 ~ "Yellow",
    Label == 201 ~ "Blue",
    Label == 202 ~ "Purple",
    Label == 203 ~ "Green",
    Label == 204 ~ "Red",
    Label == 205 ~ "Yellow",
    Label == 301 ~ "BluePink",
    Label == 302 ~ "GreenPurple",
    Label == 303 ~ "RedBlue",
    Label == 304 ~ "OrangePink",
    Label == 305 ~ "GreenOrange",
    Label == 401 ~ "BluePink",
    Label == 402 ~ "GreenPurple",
    Label == 403 ~ "RedBlue",
    Label == 404 ~ "OrangePink",
    Label == 405 ~ "GreenOrange",
    TRUE ~ as.character(Label) # Handle other cases if needed
  ))

df_kgr <- df_kgr %>%
  mutate(subj = substr(filename, nchar(filename) - 7, nchar(filename) - 4))


#df_kgr <- separate(df_kgr, stim, into = c("type", "color"), sep = "_")
#df_kgr$kgr <- df_kgr$mean - df_kgr$base
table <- aggregate(df_kgr$PPG_Rate_Mean, by=list(df_kgr$Condition, df_kgr$color), FUN=function(x) c(mean=mean(x), std=std.error(x))) # spectral density
colnames(table) <- c( 'type', 'color', 'x')

levels(table$color)
table2 <- table %>%
  mutate(hyper = if_else(type %in% c("Fractal", "HoneyComb"), 'psechedelic_like', 'control'))%>% 
  mutate(set = if_else(type %in% c("Fractal", "kaleidoscope"), 'Mystic', 'Hyperbolic'))
table2$type <- factor(table2$type, levels = c("Fractal", "HoneyComb", "kaleidoscope", "CubesControl"))
table2$color <- factor(table2$color, levels = c("Blue","BluePink","Purple","GreenPurple","Green","GreenOrange", "Red","RedBlue","Yellow","OrangePink"))


plot <- ggplot(table2, aes(x = type, y = x[,"mean"], group = color, color = color, fill = hyper)) +
  #geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.2), aes(fill = color)) +
  geom_errorbar(aes(ymin = x[,"mean"] - x[,"std"], ymax = x[,"mean"] + x[,"std"]), 
                 width = 1, position = position_dodge(width = 0.5), size = 2) +
  geom_line(size = 2, aes(group = color), position = position_dodge(width = 0.5)) +
  geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.5), aes(fill = hyper)) +
  
  scale_color_manual(breaks = c("Blue","BluePink","Purple","GreenPurple","Green","GreenOrange", "Red","RedBlue","Yellow","OrangePink"),
                    values=c("blue","blue","purple","purple","forestgreen","forestgreen","firebrick1","firebrick1","gold1","gold1")) +
  #scale_fill_manual(breaks = c("Blue","Purple","Green","Red","Yellow","BluePink","GreenPurple","RedBlue","OrangePink","GreenOrange"),
   #                  values=c("blue","purple","green","red","yellow","blue","green","red","magenta","orange")) +
  scale_fill_manual(breaks = c('psechedelic_like', 'control'),
                                      values=c("black","white")) +
                    
  # Facet the plot
  facet_wrap(~set, ncol = 2) +
  

# Customize the appearance of the plot
theme_bw() +
  labs(
    title = "PPG mean values - Baseline (-10 - 0)",
    x = "Stimulus",
    y = "Mean value",
    caption = "Baseline (-10 - 0), SE between subjects"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    plot.caption = element_text(hjust = 0, size = 10)
  )

# Print the plot
print(plot)


table <- aggregate(df_kgr$PPG_Rate_Mean, by=list(df_kgr$Condition), FUN=function(x) c(mean=mean(x), std=std.error(x))) # spectral density
colnames(table) <- c( 'type',  'x')

levels(table$color)
table2 <- table %>%
  mutate(hyper = if_else(type %in% c("Fractal", "HoneyComb"), 'psechedelic_like', 'control'))%>% 
  mutate(set = if_else(type %in% c("Fractal", "kaleidoscope"), 'Mystic', 'Hyperbolic'))
table2$type <- factor(table2$type, levels = c("Fractal", "HoneyComb", "kaleidoscope", "CubesControl"))
table2$color <- factor(table2$color, levels = c("Blue","BluePink","Purple","GreenPurple","Green","GreenOrange", "Red","RedBlue","Yellow","OrangePink"))


plot <- ggplot(table2, aes(x = type, y = x[,"mean"], group = hyper, color = hyper, fill = hyper)) +
  geom_errorbar(aes(ymin = x[,"mean"] - x[,"std"], ymax = x[,"mean"] + x[,"std"]), 
                width = 0.5, position = position_dodge(width = 0.2), size = 2) +
  geom_line(size = 2, aes(group = hyper), position = position_dodge(width = 0.2)) +
  geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.2), aes(fill = hyper)) +
  
  scale_color_manual(breaks = c('psechedelic_like', 'control'),
                    values=c("magenta2","olivedrab3")) +
    scale_fill_manual(breaks = c('psechedelic_like', 'control'),
                    values=c("white","white")) +
  
  # Facet the plot
  facet_wrap(~set, ncol = 2) +
  
  
  # Customize the appearance of the plot
  theme_bw() +
  labs(
    title = "PPG mean values - Baseline (-10 - 0)",
    x = "Stimulus",
    y = "Mean value",
    caption = "Baseline (-10 - 0), SE between trials"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    plot.caption = element_text(hjust = 0, size = 10)
  )

# Print the plot
print(plot)

df_kgr <- df_kgr %>%
  mutate(hyper = if_else(Condition %in% c("Fractal", "HoneyComb"), 'psechedelic_like', 'control'))%>% 
  mutate(set = if_else(Condition %in% c("Fractal", "kaleidoscope"), 'Mystic', 'Hyperbolic'))
sets <- unique(df_kgr$set)
p_vals <- data.table()
#cols <- c("psd")
############### for green heads (main_effects) ##############
for (i in 1:length(sets)) {
  temp <-subset(df_kgr, set == sets[i])
  m1 <- lmer(PPG_Rate_Mean ~ hyper*color + (1|subj), data = temp)
  d<-romr.fnc(m1, temp, trim = 3) ##### remove outliers
  data<- d$data 
  m <- lmer(PPG_Rate_Mean ~ hyper*color + (1|subj), data = data)
  an <- anova(m)
  an <- data.table(an,keep.rownames = TRUE)
  an_cols <- c('rn','Pr(>F)') 
  an <- an[, ..an_cols]
  an$p_value <- format(an$`Pr(>F)`, digits = 3)
  #an$interval <- j
  #an$interval <- gsub('beta power','',an$interval)
  #an <- dcast(an,formula = interval~rn,value.var = 'Pr(>F)')
  an$set <- sets[i] 
  #an$sensor_name <- files[sensor==i]$Name
  p_vals <- rbind(p_vals,an)
  
}
setwd('D:/hse/psychodelic_like_experience/data_processing/stats/')
write.csv(p_vals, "ppg_models.csv")
emm_options(lmerTest.limit = 6000)

Tuk1<- NULL
Tuk1<-data.table(summary(emmeans(m, pairwise ~ color|hyper, adjust = 'Tukey',lmer.df = "satterthwaite"))$contrasts)
Tuk1 <- Tuk1[, group1:=gsub(' -.*', '', contrast)][, group2:=gsub('.*- ', '', contrast)]
Tuk1 <- Tuk1[p.value<0.1, p_significant:=format(p.value, digits = 3)]

n <- Tuk1[!is.na(p_significant), .N]

Tuk1[p.value<0.001, stars:='***']
Tuk1[p.value<0.01 & p.value>0.001 , stars:='**']
Tuk1[p.value<0.05 & p.value>0.01 , stars:='*']
Tuk1[p.value>0.05 & p.value<0.1 , stars:='#']

#if (n>1){
#  Tuk <- Tuk[!is.na(p_significant), y.position := seq((thr1+0.01), (thr1+0.3), 0.29/(n-1))]
#} else {
#  Tuk <- Tuk[!is.na(p_significant), y.position := thr1+0.1]
#}
#y.position<-Tuk$y.position
temp <-subset(df_kgr, set == sets[1])
m1 <- lmer(PPG_Rate_Mean ~ hyper*color + (1|subj), data = temp)
d<-romr.fnc(m1, temp, trim = 3) ##### remove outliers
data<- d$data 
m <- lmer(PPG_Rate_Mean ~ hyper*color + (1|subj), data = data)

#Tuk$emmean<-y.position
Tuk2<- NULL
#thr1 <- max(data_test[, mean(data_test) + sterr(data_test), by=c('stimuli', 're')]$V1) 
#thr1 <- thr1+0.02 #for RT

#thr1_min <- min(means[!is.na(mean_beta), mean(mean_beta) - sterr(mean_beta), by=c('trial_type')]$V1) 

Tuk2<-data.table(summary(emmeans(m, pairwise ~ hyper|color, adjust = 'Tukey',lmer.df = "satterthwaite"))$contrasts)
Tuk2 <- Tuk2[, group1:=gsub(' -.*', '', contrast)][, group2:=gsub('.*- ', '', contrast)]
Tuk2 <- Tuk2[p.value<0.1, p_significant:=format(p.value, digits = 3)]

n <- Tuk2[!is.na(p_significant), .N]

Tuk2[p.value<0.001, stars:='***']
Tuk2[p.value<0.01 & p.value>0.001 , stars:='**']
Tuk2[p.value<0.05 & p.value>0.01 , stars:='*']
Tuk2[p.value>0.05 & p.value<0.1 , stars:='#']





z_scores <- table %>%
  group_by(subject) %>%
  mutate(z_score = scale(kgr))
table1 <- aggregate(z_scores$z_score, by=list(z_scores$type, z_scores$color), FUN=function(x) c(mean=mean(x))) # spectral density
colnames(table1) <- c('type',  'color', 'z_score')

table2 <- table1 %>%
  mutate(hyper = if_else(type %in% c("Fractal", "HoneyComb"), 'hyperbolic', 'nonhyperbolic'))%>% 
  mutate(set = if_else(type %in% c("Fractal", "kaleidoscope"), 'set1', 'set2'))
plot <- ggplot(table1, aes(x = type, y = z_score, group = color, color = color)) +
  geom_point(size = 3, position = position_dodge(width = 0.2)) +
  # geom_errorbar(aes(ymin = x[,"mean"] - x[,"std"], ymax = x[,"mean"] + x[,"std"]), 
  #               width = 0.2, position = position_dodge(width = 0.2)) +
  #geom_line(aes(group = color), position = position_dodge(width = 0.2)) +
  #scale_color_manual(breaks = c("compl_plus", "pictu_plus",  "sound_plus", "compl_minus",  "pictu_minus",  "sound_minus"),
  #                  values=c("darkred", "darkgreen", "darkblue", "red", "green", "blue")) +
  #scale_color_manual(breaks = c("non_reinforced", "reinforced"),
  #                   values=c("green", "red")) +
  
  # Facet the plot
  #facet_wrap(~set, ncol = 1) +
  

# Customize the appearance of the plot
theme_minimal() +
  labs(
    title = "GSR mean values - Baseline (-5 - 0)",
    x = "Type",
    y = "Mean",
    caption = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )

# Print the plot
print(plot)


#########################################
path<- "D:/hse/psychodelic_like_experience/data_processing/gsr/"
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}
df_kgr <-list.files(path,pattern = "*.csv", full.names = T) %>% 
  map_df(~read_plus(.))

replace_values <- function(x) {
  x <- ifelse(substr(x, 1, 1) == "1", "Fractal", 
              ifelse(substr(x, 1, 1) == "2", "kaleidoscope",
                     ifelse(substr(x, 1, 1) == "3", "CubesControl",
                            ifelse(substr(x, 1, 1) == "4", "HoneyComb", x))))
  return(x)
}

# Apply the function to the values
df_kgr$Condition <- replace_values(df_kgr$Condition)
df_kgr <- df_kgr %>%
  mutate(color = case_when(
    Label == 101 ~ "Blue",
    Label == 102 ~ "Purple",
    Label == 103 ~ "Green",
    Label == 104 ~ "Red",
    Label == 105 ~ "Yellow",
    Label == 201 ~ "Blue",
    Label == 202 ~ "Purple",
    Label == 203 ~ "Green",
    Label == 204 ~ "Red",
    Label == 205 ~ "Yellow",
    Label == 301 ~ "BluePink",
    Label == 302 ~ "GreenPurple",
    Label == 303 ~ "RedBlue",
    Label == 304 ~ "OrangePink",
    Label == 305 ~ "GreenOrange",
    Label == 401 ~ "BluePink",
    Label == 402 ~ "GreenPurple",
    Label == 403 ~ "RedBlue",
    Label == 404 ~ "OrangePink",
    Label == 405 ~ "GreenOrange",
    TRUE ~ as.character(Label) # Handle other cases if needed
  ))

df_kgr <- df_kgr %>%
  mutate(subj = substr(filename, nchar(filename) - 7, nchar(filename) - 4))


table <- aggregate(df_kgr$EDA_SCR, by=list(df_kgr$Condition, df_kgr$color), FUN=function(x) c(sum=sum(x))) # spectral density
colnames(table) <- c( 'type', 'color', 'x')

table2 <- table %>%
  mutate(hyper = if_else(type %in% c("Fractal", "HoneyComb"), 'psechedelic_like', 'control'))%>% 
  mutate(set = if_else(type %in% c("Fractal", "kaleidoscope"), 'Mystic', 'Hyperbolic'))
table2$type <- factor(table2$type, levels = c("Fractal", "HoneyComb", "kaleidoscope", "CubesControl"))
#table2$color <- factor(table2$color, levels = c("Blue","BluePink","Purple","GreenPurple","Green","GreenOrange", "Red","RedBlue","Yellow","OrangePink"))


plot <- ggplot(table2, aes(x = color, y = x, group = hyper, color = hyper, fill = hyper)) +
  #geom_errorbar(aes(ymin = x[,"mean"] - x[,"std"], ymax = x[,"mean"] + x[,"std"]), 
  #              width = 0.5, position = position_dodge(width = 0.2), size = 2) +
  geom_bar(binwidth = 5,aes(group = hyper, fill = hyper), stat = "identity")+
  geom_line(size = 2, aes(group = hyper), position = position_dodge(width = 0.2)) +
  #geom_point(size = 2.5, shape = 23, position = position_dodge(width = 0.2), aes(fill = hyper)) +
  
  scale_color_manual(breaks = c('psechedelic_like', 'control'),
                     values=c("magenta2","olivedrab3")) +
  scale_fill_manual(breaks = c('psechedelic_like', 'control'),
                    values=c("white","white")) +
  
  # Facet the plot
  #facet_wrap(~set, ncol = 2) +
  
  
  # Customize the appearance of the plot
  theme_bw() +
  labs(
    title = "Number of SCR",
    x = "Color",
    y = "Number",
    caption = "EDA, SE between trials"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    legend.position = "top",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    plot.caption = element_text(hjust = 0, size = 10)
  )

# Print the plot
print(plot)




data <- data.frame(Factor1, Factor2, Outcome)

# Проведение теста на равенство долей для каждой комбинации факторов
test_result <- prop.test(table(data$Outcome, data$Factor1, data$Factor2), by = c("Factor1", "Factor2"))
print(test_result)



################опроснички
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm,
           number = row_number())
}
path<- "D:/hse/psychodelic_like_experience/data_processing/otvety/"

otvety <-list.files(path,pattern = "*.csv", full.names = T) %>% 
  map_df(~read_plus(.))
otvety <- otvety %>%
  mutate(subj = substr(filename, nchar(filename) - 7, nchar(filename) - 4))
colnames(otvety) <- c( "Subject",   "Condition", "Otvet",    "filename",  "number",   "subj")

otvety <- otvety[, c( "Subject",   "Condition", "Otvet",    "number")]



replacements <- list(
  "HoneyComb" = "HoneyComb/",
  "kaleidoscope" = "kaleidoscope/",
  "Fractal" = "Fractal/",
  "CubesControl" = "CubesControl/"
)
library(stringr)


for (pattern in names(replacements)) {
  replacement <- replacements[[pattern]]
  
  otvety$Condition <- gsub(pattern, replacement, otvety$Condition)
}
otvety <- separate(otvety, Condition, into = c("type", "color"), sep = "/")
write.csv(otvety, "D:/hse/psychodelic_like_experience/data_processing/otvety.csv")











