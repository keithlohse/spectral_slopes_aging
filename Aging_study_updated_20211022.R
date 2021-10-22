# Aging study analysis
## Opening libraries -----------------------------------------------------------
library("lme4"); library("lmerTest"); library("car"); 
library("psych"); library("mediation"); library("tidyverse");
library("ggdark"); library("patchwork")

## Set working Directory ------------------------------------------------
getwd()

setwd("~/GitHub/spectral_slopes_aging")

# let's see what is in the data folder
list.files()


## Read csv file ---------------------------------------------------------
MASTER <- read.csv("./MASTER_EO_and_EC_EEG_KRL.csv", 
                   header = TRUE,na.strings=c("","NA","na"),
                   stringsAsFactors = TRUE)
head(MASTER)
summary(as.factor(MASTER$Hz))
summary(MASTER$condition)
summary(MASTER$subID)


MASTER$Frontal<-(MASTER$F7+MASTER$F3+MASTER$Fz+MASTER$F4+MASTER$F8)/5
MASTER$Central<-(MASTER$C3+MASTER$Cz+MASTER$C4)/3
MASTER$Parietal<-(MASTER$P3+MASTER$P7+MASTER$Pz+MASTER$P4+MASTER$P8)/5
MASTER$Occipital<-(MASTER$O1+MASTER$Oz+MASTER$O2)/3


DATA<-subset(MASTER, condition=="ec")
DATA<-subset(DATA, group!="mc")

head(DATA)

write.csv(DATA, "./MASTER_FOOOF.csv")

xtabs(~Hz+subID+condition, data=MASTER)

DATA$Labels <- fct_recode(DATA$group, 
                              "Older Adults" = "oa", 
                              "Younger Adults" = "ya")

DATA$Labels <- fct_relevel(DATA$Labels, 
                               "Younger Adults", 
                               "Older Adults")

DATA <- subset(DATA, Hz >=2)
DATA <- subset(DATA, Hz <=25)
summary(DATA$Hz)

# Figure 2A ----
ggplot(data = DATA, 
       mapping = aes(x = log(Hz), y = log(Frontal))) +
  geom_line(aes(col=subID, lty=Labels), lwd=0.5)+
  stat_summary(aes(lty=Labels), fun="mean", geom="line", lwd=1.5)+
  facet_grid(~Labels)+
  scale_x_continuous(name = "Frequency (ln(Hz))")+
  scale_y_continuous(name = "Frontal Power (ln(uV^2))")+
  theme_bw() + 
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0.5,
                                lineheight=1.2),  # title
        panel.spacing = unit(2, "lines"),
        plot.caption=element_text(size=16, face = "bold"),  # caption
        strip.text = element_text(size=16),
        axis.title.x=element_text(size=16, face = "bold"),  # X axis title
        axis.title.y=element_text(size=16, face = "bold"),  # Y axis title
        axis.text.x=element_text(size=16, color="black"),  # X axis text
        axis.text.y=element_text(size=16, color="black"),
        legend.title=element_text(size=16, face = "bold"),
        legend.text=element_text(size=16),  # X axis text
        legend.position = "none")


### FOOOF outputs ----------------------------------------------------

### Calculate exponents and offsets using FOOOF from the 
### MASTER_FOOOF csv file for Frontal, Central, Parietal and Occipital averaged channels
### save the FOOOF output variables in a separate csv titled Aging_study_FOOOF_var


DAT1 <- read.csv("./Aging_study_FOOOF_var.csv", 
                   header = TRUE,na.strings=c("","NA","na"),
                   stringsAsFactors = TRUE)
head(DAT1)
tail(DAT1)

DAT1$group<-factor(DAT1$group)
summary(DAT1$group)

summary(DAT1$R2)



### Are slopes different between the two age groups? -----------------------------------
## Are there age related differences in different regions

### Age related differences in Aperiodic slopes
M1 <- lmer(Exponent ~ group*Area+(1|SubID), data = DAT1, REML=FALSE)
Anova(M1, Type="III") 
anova(M1)
summary(M1)

#### Since we have a significant interaction, we will test each channel
Frontal <- subset(DAT1, Area =="Frontal")
Central <- subset(DAT1, Area =="Central")
Parietal <- subset(DAT1, Area =="Parietal")
Occipital  <- subset(DAT1, Area =="Occipital")


### Frontal channel 
mod1<- lm(Exponent ~ 1+group, data = Frontal)
anova(mod1)


### Central Channel
mod2<- lm(Exponent ~ 1+group, data = Central)
anova(mod2)


### Parietal Channel
mod3<- lm(Exponent ~ 1+group, data = Parietal)
anova(mod3)


### Occipital Channel
mod4<- lm(Exponent ~ 1+ group, data = Occipital)
anova(mod4)

# Age-related differences in Intercepts
M2 <- lmer(Offset ~ group*Area+(1|SubID), data = DAT1, REML=FALSE)
Anova(M2, Type="III")  ### the difference in slopes between both hemispheres is non-sig
anova(M2)
summary(M2)




#### Exponent and offsets in frontal channels ------------------
DAT1$Labels <- fct_recode(DAT1$group, 
                          "Older Adults" = "OA", 
                          "Younger Adults" = "YA")

DAT1$Labels <- fct_relevel(DAT1$Labels, 
                           "Younger Adults", 
                           "Older Adults")

DAT1$group <- fct_relevel(DAT1$group, 
                           "YA", 
                           "OA")

DAT1$Area <- fct_relevel(DAT1$Area, 
                          "Frontal", "Central", 
                          "Parietal", "Occipital")

# Figure 2B ----
ggplot(data = DAT1, 
             mapping = aes(x = group, y = Exponent)) +
  geom_point(aes(fill=Labels), pch=21, size=2, position=position_jitter(width=0.1)) +
  geom_boxplot(aes(fill=Labels), col="black", alpha=0.2, outlier.shape = NA) +
  facet_wrap(~Area, ncol=2)+
  scale_fill_manual(values=c("white", "grey20"))+
  scale_x_discrete(name = "Age Group")+
  scale_y_continuous(name = "Exponent")+
  theme_bw() + 
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0.5,
                                lineheight=1.2),  # title
        panel.spacing = unit(2, "lines"),
        plot.caption=element_text(size=16, face = "bold"),  # caption
        strip.text = element_text(size=16),
        axis.title.x=element_text(size=16, face = "bold"),  # X axis title
        axis.title.y=element_text(size=16, face = "bold"),  # Y axis title
        axis.text.x=element_text(size=16, color="black"),  # X axis text
        axis.text.y=element_text(size=16, color="black"),
        legend.title=element_text(size=16, face = "bold"),
        legend.text=element_text(size=16),  # X axis text
        legend.position = "none")



head(DAT1)
summary(DAT1$Exponent)
summary(DAT1$Offset)

# Analysis #2: Age Related Differences in the RBANS ----------------------------
list.files()

BEHAVIOR <- read.csv("./RBANS_aging_study_08262021.csv", 
                     header = TRUE,na.strings=c("","NA","na"),
                     stringsAsFactors = TRUE)

summary(BEHAVIOR$subID)

# Drop subjects who are missing EEG data:
BEHAVIOR <- subset(BEHAVIOR, subID != "ya20")
BEHAVIOR <- subset(BEHAVIOR, subID != "oa22")


# Reliability across RBANS subtest items:
head(BEHAVIOR[,c(7:18)])
psych::alpha(BEHAVIOR[,c(7:18)])
# There is a reasonably high reliability across items, justifying the 
# normalized average scores

summary(BEHAVIOR$group)


BEHAVIOR$im_mem <- ((BEHAVIOR$list_learn/40+BEHAVIOR$story_mem/24)/2)*100
head(BEHAVIOR)


BEHAVIOR$vis <- ((BEHAVIOR$fig_copy/20 + BEHAVIOR$line_orn/20)/2)*100 
head(BEHAVIOR)

BEHAVIOR$lang <- ((BEHAVIOR$pic_name/10 + BEHAVIOR$sem_flu/40)/2)*100
head(BEHAVIOR)

BEHAVIOR$att <- ((BEHAVIOR$dig_spn/16 + BEHAVIOR$coding/89)/2)*100
head(BEHAVIOR)


BEHAVIOR$dm <- ((BEHAVIOR$list_rec/10 + BEHAVIOR$list_recog/20+ BEHAVIOR$story_rec/12+ BEHAVIOR$fig_rec/20)/4)*100
head(BEHAVIOR)



# Reformatting the data for graphing:
df<- BEHAVIOR %>% 
  select(subID, group, im_mem, vis, lang, att, dm, raw_ave)%>%
  mutate(raw_ave = raw_ave*100)
head(df)
df<-gather(df, scale, score, im_mem:raw_ave, factor_key=TRUE)
head(df)
summary(df$scale)

SUMMARY<-df %>% group_by(group, scale) %>%
  summarize(mean = mean(score, na.rm=TRUE),
            std = sd(score, na.rm=TRUE),
            count = n(),
            ll = mean(score, na.rm=TRUE)-(sd(score, na.rm=TRUE)/sqrt(n()))*qt(0.975, n()-1, lower.tail = TRUE),
            ul = mean(score, na.rm=TRUE)+(sd(score, na.rm=TRUE)/sqrt(n()))*qt(0.975, n()-1, lower.tail = TRUE))

write.csv(SUMMARY, "./descriptive_group_stats.csv")

# Plot of RBANS scores by group and test
scale_names <- c(`im_mem` = "Immediate Memory", `vis` = "Visuospatial/Constructional", 
                 `lang` = "Language", `att` = "Attention",
                 `dm` = "Delayed Memory", `raw_ave` = "Average")

df$Labels <- fct_recode(df$group, 
                        "OA" = "OA", 
                        "YA" = "YA",
                        "MCI" = "MC")

df$Labels <- fct_relevel(df$Labels, "YA", "OA", "MCI")

summary(df$group)
df_MCI <- subset(df, group == "MC")
df <- subset(df, group != "MC")

# Figure 3: RBANS Subtests
head(df)

summary(df$scale)
ggplot(df, aes(x = Labels, y = score)) +
  geom_point(aes(fill=Labels), pch=21, size=1, stroke=1, col="black",
             position = position_jitterdodge()) +
  geom_boxplot(aes(fill=Labels), col="black", alpha=0.2, outlier.shape = NA)+
  facet_wrap(~scale, ncol=3, scales="free", labeller = as_labeller((scale_names)))+
  scale_fill_manual(values=c("white", "grey20"))+
  scale_x_discrete(name = "Age Group") +
  scale_y_continuous(name = "Domain Score (%)")+
  theme_bw() + 
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0.5,
                                lineheight=1.2),  # title
        #plot.background = element_rect(fill="grey10"),
        #panel.background = element_rect(fill="grey20"),
        plot.caption=element_text(size=16, face = "bold"),  # caption
        strip.text = element_text(size=16),
        axis.title.x=element_text(size=16, face = "bold"),  # X axis title
        axis.title.y=element_text(size=16, face = "bold"),  # Y axis title
        axis.text.x=element_text(size=16, color="black"),  # X axis text
        axis.text.y=element_text(size=16, color="black"),
        legend.title=element_text(size=16, face = "bold"),
        legend.text=element_text(size=16),  # X axis text
        legend.position = "none")


# T-tests comparing Older and Younger Adults ----
head(BEHAVIOR)

t.test(im_mem ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less") # Non-sig
t.test(vis ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less") # Non-Sig
t.test(lang ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less") # Non-sig
t.test(att ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less") # Sig
t.test(dm ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less") # Sig
t.test(raw_ave ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less") # Sig


# Analysis #3: Does the spectral slope mediate age-related differences? --------

# Showing the relationship between RBANS scores and Frontal Slopes:
BEHAVIOR$Labels <- fct_recode(BEHAVIOR$group, 
                              "OA" = "OA", 
                              "YA" = "YA",
                              "MCI" = "MC")

BEHAVIOR$Labels <- fct_relevel(BEHAVIOR$Labels, "YA", "OA", "MCI")

BEHAVIOR <- subset(BEHAVIOR, group!="MC")

head(BEHAVIOR)
# Figure 3A: Raw Average
BEHAVIOR$raw_ave <- BEHAVIOR$raw_ave*100
head(BEHAVIOR)
summary(BEHAVIOR$Labels)

fig3A <- ggplot(data = BEHAVIOR, aes(x = Frontal_Exponent, y = raw_ave)) +
  geom_point(aes(fill=Labels), col="black", pch = 21, size = 2, alpha=0.80) + 
  stat_smooth(aes(col=Labels, lty=Labels), method="lm", se=FALSE, lwd=1) +
  scale_x_continuous(name = "Exponent (Frontal)") +
  scale_y_continuous(name = "Average (%)")+
  scale_color_grey(start=0.0, end=0.4)+
  scale_fill_manual(values=c("white", "grey20"))+
  theme_bw() + 
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0.5,
                                lineheight=1.2),  # title
        #plot.background = element_rect(fill="grey10"),
        #panel.background = element_rect(fill="grey20"),
        plot.caption=element_text(size=14, face = "bold"),  # caption
        strip.text = element_text(size=14),
        axis.title.x=element_text(size=14, face = "bold"),  # X axis title
        axis.title.y=element_text(size=14, face = "bold"),  # Y axis title
        axis.text.x=element_text(size=12, color="black"),  # X axis text
        axis.text.y=element_text(size=12, color="black"),
        legend.title=element_text(size=14, face = "bold"),
        legend.text=element_text(size=14),  # X axis text
        legend.position = "none")

head(BEHAVIOR)
# Figure 3B: Attention
fig3B <- ggplot(data = BEHAVIOR, aes(x = Frontal_Exponent, y = att)) +
  geom_point(aes(fill=Labels), col="black", pch = 21, size = 2, alpha=0.80) + 
  stat_smooth(aes(col=Labels, lty=Labels), method="lm", se=FALSE, lwd=1) +
  scale_color_grey(start=0.0, end=0.4)+
  scale_fill_manual(values=c("white", "grey20"))+
  scale_x_continuous(name = "Exponent (Frontal)") +
  scale_y_continuous(name = "Attention (%)")+
  theme_bw() + 
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0.5,
                                lineheight=1.2),  # title
        #plot.background = element_rect(fill="grey10"),
        #panel.background = element_rect(fill="grey20"),
        plot.caption=element_text(size=14, face = "bold"),  # caption
        strip.text = element_text(size=14),
        axis.title.x=element_text(size=14, face = "bold"),  # X axis title
        axis.title.y=element_text(size=14, face = "bold"),  # Y axis title
        axis.text.x=element_text(size=12, color="black"),  # X axis text
        axis.text.y=element_text(size=12, color="black"),
        legend.title=element_text(size=14, face = "bold"),
        legend.text=element_text(size=14),  # X axis text
        legend.position = "bottom")

# Figure 3C: Delayed Memory
fig3C <- ggplot(data = BEHAVIOR, aes(x = Frontal_Exponent, y = dm)) +
  geom_point(aes(fill=Labels), col="black", pch = 21, size = 2, alpha=0.80) + 
  stat_smooth(aes(col=Labels, lty=Labels), method="lm", se=FALSE, lwd=1) +
  scale_x_continuous(name = "Exponent (Frontal)") +
  scale_y_continuous(name = "Delayed Memory (%)")+
  scale_color_grey(start=0.0, end=0.4)+
  scale_fill_manual(values=c("white", "grey20"))+
  theme_bw() + 
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0.5,
                                lineheight=1.2),  # title
        #plot.background = element_rect(fill="grey10"),
        #panel.background = element_rect(fill="grey20"),
        plot.caption=element_text(size=14, face = "bold"),  # caption
        strip.text = element_text(size=14),
        axis.title.x=element_text(size=14, face = "bold"),  # X axis title
        axis.title.y=element_text(size=14, face = "bold"),  # Y axis title
        axis.text.x=element_text(size=12, color="black"),  # X axis text
        axis.text.y=element_text(size=12, color="black"),
        legend.title=element_text(size=14, face = "bold"),
        legend.text=element_text(size=14),  # X axis text
        legend.position = "none")



(fig3A|fig3B|fig3C)


# Mediation Analysis: RBANS Average --------------------------------------------
# Y = raw_ave
# Mediator = Exponent (Frontal)
# X = Age Group

# Y ~ X
mod1 <- lm(raw_ave~group,
           data=BEHAVIOR[BEHAVIOR$raw_ave!="NA",])
summary(mod1)

# M ~ X
mod2 <- lm(Frontal_Exponent~group,
           data=BEHAVIOR[BEHAVIOR$raw_ave!="NA",])
summary(mod2)

# Y ~ X + M
mod3 <- lm(raw_ave~group+Frontal_Exponent,
           data=BEHAVIOR[BEHAVIOR$raw_ave!="NA",])
summary(mod3)

# Bootstrapping Mediation Analysis:
set.seed(1)
raw_ave.med.result<- mediate(model.m=mod2, model.y=mod3, 
                             treat='group', mediator='Frontal_Exponent',
                             boot=TRUE, sims=5000)
summary(raw_ave.med.result)

plot(raw_ave.med.result)


# Mediation Analysis: Attention ---------------------------------------------------
# Y = Attention
# Mediator = Exponent (Frontal)
# X = Age Group

# Y ~ X
mod1 <- lm(att~group,
           data=BEHAVIOR[BEHAVIOR$att!="NA",])
summary(mod1)

# M ~ X
mod2 <- lm(Frontal_Exponent~group,
           data=BEHAVIOR[BEHAVIOR$att!="NA",])
summary(mod2)

# Y ~ X + M
mod3 <- lm(att~group+Frontal_Exponent,
           data=BEHAVIOR[BEHAVIOR$att!="NA",])
summary(mod3)

# Bootstrapping Mediation Analysis:
set.seed(1)
att.med.result<- mediate(model.m=mod2, model.y=mod3, 
                            treat='group', mediator='Frontal_Exponent',
                            boot=TRUE, sims=5000)
summary(att.med.result)
plot(att.med.result)


# Mediation Analysis: Delayed Memory ----------------------------------------------
# Y = Delayed Memory
# Mediator = Spectral Slope (Frontal)
# X = Age Group

# Y ~ X
mod1 <- lm(dm~group,
           data=BEHAVIOR[BEHAVIOR$dm!="NA",])
summary(mod1)

# M ~ X
mod2 <- lm(Frontal_Exponent~group,
           data=BEHAVIOR[BEHAVIOR$dm!="NA",])
summary(mod2)

# Y ~ X + M
mod3 <- lm(dm~group+Frontal_Exponent,
           data=BEHAVIOR[BEHAVIOR$dm!="NA",])
summary(mod3)

# Bootstrapping Mediation Analysis:
set.seed(1)
dm.med.result<- mediate(model.m=mod2, model.y=mod3, 
                                treat='group', mediator='Frontal_Exponent',
                                boot=TRUE, sims=5000)
summary(dm.med.result)
plot(dm.med.result)




## Follow-Up Mediation looking at Coding and Digit Span Alone
# Mediation Analysis: Coding ---------------------------------------------------
# Y = Coding
# Mediator = Exponent (Frontal)
# X = Age Group

BEHAVIOR$coding.p <- BEHAVIOR$coding/89*100

# Y ~ X
mod1 <- lm(coding.p~group,
           data=BEHAVIOR[BEHAVIOR$att!="NA",])
summary(mod1)

# M ~ X
mod2 <- lm(Frontal_Exponent~group,
           data=BEHAVIOR[BEHAVIOR$att!="NA",])
summary(mod2)

# Y ~ X + M
mod3 <- lm(coding.p~group+Frontal_Exponent,
           data=BEHAVIOR[BEHAVIOR$att!="NA",])
summary(mod3)

# Bootstrapping Mediation Analysis:
set.seed(1)
coding.med.result<- mediate(model.m=mod2, model.y=mod3, 
                         treat='group', mediator='Frontal_Exponent',
                         boot=TRUE, sims=5000)
summary(coding.med.result)
plot(coding.med.result)



ggplot(data = BEHAVIOR, aes(x = Frontal_Exponent, y = (coding/89)*100)) +
  geom_point(aes(fill=Labels), col="black", pch = 21, size = 2, alpha=0.80) + 
  stat_smooth(aes(col=Labels, lty=Labels), method="lm", se=FALSE, lwd=1) +
  scale_x_continuous(name = "Exponent (Frontal)") +
  scale_y_continuous(name = "Coding (%)")+
  scale_color_grey(start=0.0, end=0.4)+
  scale_fill_manual(values=c("white", "grey20"))+
  theme_bw() + 
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0.5,
                                lineheight=1.2),  # title
        #plot.background = element_rect(fill="grey10"),
        #panel.background = element_rect(fill="grey20"),
        plot.caption=element_text(size=14, face = "bold"),  # caption
        strip.text = element_text(size=14),
        axis.title.x=element_text(size=14, face = "bold"),  # X axis title
        axis.title.y=element_text(size=14, face = "bold"),  # Y axis title
        axis.text.x=element_text(size=12, color="black"),  # X axis text
        axis.text.y=element_text(size=12, color="black"),
        legend.title=element_text(size=14, face = "bold"),
        legend.text=element_text(size=14),  # X axis text
        legend.position = "none")









# Mediation Analysis: Digit Span -----------------------------------------------
# Y = Digit Span
# Mediator = Exponent (Frontal)
# X = Age Group

BEHAVIOR$dig_spn.p <- BEHAVIOR$dig_spn/16*100

# Y ~ X
mod1 <- lm(dig_spn.p~group,
           data=BEHAVIOR[BEHAVIOR$att!="NA",])
summary(mod1)

# M ~ X
mod2 <- lm(Frontal_Exponent~group,
           data=BEHAVIOR[BEHAVIOR$att!="NA",])
summary(mod2)

# Y ~ X + M
mod3 <- lm(dig_spn.p~group+Frontal_Exponent,
           data=BEHAVIOR[BEHAVIOR$att!="NA",])
summary(mod3)

# Bootstrapping Mediation Analysis:
set.seed(1)
ds.med.result<- mediate(model.m=mod2, model.y=mod3, 
                         treat='group', mediator='Frontal_Exponent',
                         boot=TRUE, sims=5000)
summary(ds.med.result)
plot(ds.med.result)




ggplot(data = BEHAVIOR, aes(x = Frontal_Exponent, y = (dig_spn/16)*100)) +
  geom_point(aes(fill=Labels), col="black", pch = 21, size = 2, alpha=0.80) + 
  stat_smooth(aes(col=Labels, lty=Labels), method="lm", se=FALSE, lwd=1) +
  scale_x_continuous(name = "Exponent (Frontal)") +
  scale_y_continuous(name = "Digit Span (%)")+
  scale_color_grey(start=0.0, end=0.4)+
  scale_fill_manual(values=c("white", "grey20"))+
  theme_bw() + 
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0.5,
                                lineheight=1.2),  # title
        #plot.background = element_rect(fill="grey10"),
        #panel.background = element_rect(fill="grey20"),
        plot.caption=element_text(size=14, face = "bold"),  # caption
        strip.text = element_text(size=14),
        axis.title.x=element_text(size=14, face = "bold"),  # X axis title
        axis.title.y=element_text(size=14, face = "bold"),  # Y axis title
        axis.text.x=element_text(size=12, color="black"),  # X axis text
        axis.text.y=element_text(size=12, color="black"),
        legend.title=element_text(size=14, face = "bold"),
        legend.text=element_text(size=14),  # X axis text
        legend.position = "none")

