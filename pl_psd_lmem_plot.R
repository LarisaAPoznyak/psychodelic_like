library(reshape2)
library(data.table)
library(ggplot2)
library(lme4)
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




path<- "D:/hse/psychodelic_like_experience/data_processing/psd_df/"


read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}

df <-list.files(path,pattern = "*.csv", full.names = T) %>% 
  map_df(~read_plus(.))

head(df, 3) #data check

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


z_scores <- table %>%
  group_by(subject,freq) %>% #normalize data to reduce intersubject differences and balance powers 
  mutate(value = scale(psd))

z_scores <- separate(z_scores, condition, into = c("condition", "color"), sep = "/") #separate type and color
colnames(df)
colnames(z_scores)[8]<-'z_score'


#fix messy color codes
data2 <- z_scores %>%
  mutate(color = case_when(
    color == 'GreenPurple'  ~ "PurpleGreen",
    color == 'BlueRed'  ~ "RedBlue",
    color == 'PinkYellow'  ~ "OrangePink",
    color == 'PinkBlue'  ~ "BluePink",
    TRUE ~ as.character(color)
  ))


data2 <- data2 %>%
  mutate(hyper = if_else(condition %in% c("Fractal", "HoneyComb"), 'hyperbolic', 'nonhyperbolic'))


data3 <- subset(data2, condition %in% c("CubesControl", "HoneyComb"))
data4 <- subset(data2, condition %in% c("Fractal", "kaleidoscope"))

# Chose set
#№1 
delta <- subset(data3, freq == "delta" )
theta <- subset(data3, freq == "theta" )
alpha <- subset(data3, freq == "alpha" )
beta <- subset(data3, freq == "beta" )

#№2
delta <- subset(data4, freq == "delta" )
theta <- subset(data4, freq == "theta" )
alpha <- subset(data4, freq == "alpha" )
beta <- subset(data4, freq == "beta" )


sensors<- unique(data2$sensor)

p_vals <- data.table()
############### for main_effects ##############
for (i in 1:length(sensors)) {
  temp <-subset(alpha, sensor == sensors[i])
  m1 <- lmer(z_score ~ hyper*color + (1|subject), data = temp)
  d<-romr.fnc(m1, temp, trim = 3) ##### remove outliers
  data<- d$data 
  m <- lmer(z_score ~ hyper*color + (1|subject), data = data)
  an <- anova(m)
  an <- data.table(an,keep.rownames = TRUE)
  an$p_value <- format(an$`Pr(>F)`, digits = 3)
  an$sensor <- sensors[i] 
  p_vals <- rbind(p_vals,an)
  
}
setwd('D:/hse/psychodelic_like_experience/data_processing/stats/')
write.csv(p_vals, "aov_results_alpha_honey_0108.csv")




emm_options(lmerTest.limit = 100000)
emm_options(pbkrtest.limit = 90000)

sens <- subset(p_vals, rn == "color")
sens <- subset(sens, p_value < 0.05)
ses <- c(sens$sensor)
alpha_color <- subset(alpha, sensor %in% ses)
alpha_color


m1 <- lmer(z_score ~ hyper*color + (1|subject), data = alpha_color)

#means within the model
marginal_em <- emmeans(m1, ~ as.factor(color|hyper),level = 0.95,lmer.df = "satterthwaite")
marginal_em<-as.data.frame(marginal_em)


#post-hoc comparisons
m1 <- lmer(z_score ~ hyper*color + (1|subject), data = alpha_color)
Tuk<-data.table(summary(emmeans(m1, pairwise ~ color|1, adjust = 'tukey',lmer.df = "satterthwaite"))$contrasts)
Tuk <- Tuk[, group1:=gsub(' -.*', '', contrast)][, group2:=gsub('.*- ', '', contrast)]
Tuk <- Tuk[p.value<0.1, p_significant:=format(p.value, digits = 3)]
n <- Tuk[!is.na(p_significant), .N]
Tuk[p.value<0.001, stars:='***']
Tuk[p.value<0.01 & p.value>0.001 , stars:='**']
Tuk[p.value<0.05 & p.value>0.01 , stars:='*']
Tuk[p.value>0.05 & p.value<0.1 , stars:='#']

#plot
plot <- ggplot(marginal_em, aes(x = color, y = emmean, group = hyper, color = hyper, fill = hyper)) +
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
  
  #stat_pvalue_manual(Tuk, label = 'stars', size = 7, bracket.size = 1.5, tip.length = 0.01,
  #                   y.position =c(max(marginal_em$emmean)+max(marginal_em$SE) + 0.1),inherit.aes = FALSE) +
  # Customize the appearance of the plot
  theme_bw() +
  labs(
    title = "",
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

# Print the plot
print(plot)


