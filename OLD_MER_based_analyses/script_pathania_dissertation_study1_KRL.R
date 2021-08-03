# Aging study analysis
## Opening libraries -----------------------------------------------------------
library("lme4"); library("lmerTest"); library("car"); 
library("psych"); library("mediation"); library("tidyverse");
library("ggdark"); library("patchwork")

## Set working Directory ------------------------------------------------
getwd()

setwd("C:/Users/kelop/Box/Trainees/Anupriya Pathania/Dissertation")

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


DATA<-subset(MASTER, condition=="ec")
# DATA<-subset(DATA, group!="mc")

DATA %>% group_by(subID, group) %>% summarize(count=length(Hz))


DATA$group<-factor(DATA$group)
summary(DATA$group)

summary(DATA$subID)

# Average out power across all frontal electrodes -------------------------
DATA$Frontal<-(DATA$F7+DATA$F3+DATA$Fz+DATA$F4+DATA$F8)/5
DATA$Central<-(DATA$C3+DATA$Cz+DATA$C4)/3
DATA$Parietal<-(DATA$P3+DATA$P7+DATA$Pz+DATA$P4+DATA$P8)/5
DATA$Occipital<-(DATA$O1+DATA$Oz+DATA$O2)/3

summary(DATA$Frontal)
summary(DATA$Central)
summary(DATA$Parietal)
summary(DATA$Occipital)

head(DATA)
DF<-DATA[,c("subID", "group","Hz", "Frontal", "Central", "Parietal", "Occipital")]
head(DF)

DF_LONG <- gather(DF, Region, Power, Frontal:Occipital, factor_key=TRUE)
head(DF_LONG)

DF_LONG$lgHz <- log(DF_LONG$Hz)
DF_LONG$lgPower <- log(DF_LONG$Power)

# Filter data between 2 and 30 Hz -------------------------------
FILTERED<-subset(DF_LONG, Hz>=2) # Lose everything below 2
FILTERED<-subset(FILTERED, Hz<=32) # Lose everything above 32 (drop gamma)


head(FILTERED)
FILTERED$Labels <- fct_recode(FILTERED$group, 
                              "Older Adults" = "oa", 
                              "Younger Adults" = "ya",
                              "MCI" = "mc")

FILTERED$Labels <- fct_relevel(FILTERED$Labels, 
                               "Younger Adults", 
                               "Older Adults",
                               "MCI")

DF_LONG_MEAN <- FILTERED %>%
  group_by(group, Labels, Region, Hz, lgHz) %>%
  summarise(Power=mean(Power, na.rm=TRUE),
            lgPower=mean(lgPower, na.rm=TRUE))

head(DF_LONG_MEAN)

NO_MCI <- subset(DF_LONG_MEAN, Labels != "MCI")

summary(FILTERED$Group)

FRONTAL <- subset(FILTERED, Region=="Frontal")
CENTRAL <- subset(FILTERED, Region=="Central")
PARIETAL <- subset(FILTERED, Region=="Parietal")
OCCIPITAL <- subset(FILTERED, Region=="Occipital")


# Analysis #1: Age Related Differences in the Spectral Slope
# # FIGURE 00: PSD Showing the Raw Data across positions
# ggplot(data = FRONTAL, 
#        mapping = aes(x = Hz, y = Power)) +
#   geom_line(aes(group=subID, lty=Labels, col=subID), lwd=1) +
#   #geom_point(aes(fill=subID), pch = 21, size = 2) + 
#   facet_wrap(~Labels)+
#   scale_x_continuous(name = "Frequency (Hz)") +
#   scale_y_continuous(name = "Frontal Power (uV^2)")+
#   theme_bw() + 
#   theme(plot.title=element_text(size=15, 
#                                 face="bold", 
#                                 hjust=0.5,
#                                 lineheight=1.2),  # title
#         #plot.background = element_rect(fill="grey10"),
#         #panel.background = element_rect(fill="grey20"),
#         plot.caption=element_text(size=16, face = "bold"),  # caption
#         strip.text = element_text(size=16),
#         axis.title.x=element_text(size=16, face = "bold"),  # X axis title
#         axis.title.y=element_text(size=16, face = "bold"),  # Y axis title
#         axis.text.x=element_text(size=16),  # X axis text
#         axis.text.y=element_text(size=16),
#         legend.position = "none") +
#   geom_line(data=DF_LONG_MEAN[DF_LONG_MEAN$Region=="Frontal",], aes(lty=Labels), col="black", lwd=1.25)
#   


# FIGURE 1A: PSD Showing the Raw Data across positions----
summary(DF_LONG_MEAN$group)

Fig2A<-ggplot(data = FRONTAL[FRONTAL$Labels!="MCI",], 
       mapping = aes(x = lgHz, y = lgPower)) +
  geom_line(aes(group=subID, lty=Labels, col=subID), lwd=0.5) + 
  #geom_point(aes(fill=subID), pch = 21, size = 2) + 
  facet_wrap(~Labels)+
  scale_x_continuous(name = "Frequency (ln(Hz))", ) +
  scale_y_continuous(name = "Frontal Power (ln(uV^2))")+
  theme_bw() + 
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0.5,
                                lineheight=1.2),  # title
        panel.spacing = unit(2, "lines"),
        #plot.background = element_rect(fill="grey10"),
        #panel.background = element_rect(fill="grey20"),
        plot.caption=element_text(size=16, face = "bold"),  # caption
        strip.text = element_text(size=16),
        axis.title.x=element_text(size=16, face = "bold"),  # X axis title
        axis.title.y=element_text(size=16, face = "bold"),  # Y axis title
        axis.text.x=element_text(size=16, color="black"),  # X axis text
        axis.text.y=element_text(size=16, color="black"),
        legend.position = "none") +
  geom_line(data=NO_MCI[NO_MCI$Region=="Frontal",], aes(lty=Labels), col="black", lwd=1.5)
Fig2A

# FIGURE 1B: Spectra by Region in YA and OA----
head(NO_MCI)

Fig2B<-ggplot(data = NO_MCI, 
       mapping = aes(x = lgHz, y = lgPower)) +
  annotate("rect", xmin=log(8), xmax=log(13), ymin=-Inf, ymax=Inf, alpha=0.6, fill="grey")+
  geom_line(aes(group=Region, lty=Region, col=Region), lwd=1) +
  #geom_point(aes(fill=subID), pch = 21, size = 2) + 
  #stat_smooth(aes(lty=Group), method="gam", se=FALSE, col="black")+
  facet_wrap(~Labels)+
  scale_x_continuous(name = "Frequency (ln(Hz))") +
  scale_y_continuous(name = "Power (ln(uV^2))", limits=c(-1.5,2))+
  scale_color_grey(start=0.0, end=0.5)+
  #guides(lty=TRUE)+
  theme_bw() + 
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0.5,
                                lineheight=1.2),  # title
        #plot.background = element_rect(fill="grey10"),
        #panel.background = element_rect(fill="grey20"),
        panel.spacing = unit(2, "lines"),
        plot.caption=element_text(size=16, face = "bold"),  # caption
        strip.text = element_text(size=16),
        axis.title.x=element_text(size=16, face = "bold"),  # X axis title
        axis.title.y=element_text(size=16, face = "bold"),  # Y axis title
        axis.text.x=element_text(size=16, color="black"),  # X axis text
        axis.text.y=element_text(size=16, color="black"),
        legend.title=element_text(size=16, face = "bold"),
        legend.text=element_text(size=16),  # X axis text
        legend.position = "bottom")
Fig2B

# FIGURE 1C: Spectra by Region in YA and OA----
head(FILTERED)


# Excluding the Alpha Band
# Exclude alpha
NO_ALPHA<-subset(FILTERED, Hz <= 7 | Hz >= 13)
summary(as.factor(NO_ALPHA$Hz))


# Excluding the MCI group for statistical comparisons.
NO_ALPHA_NO_MCI <- subset(NO_ALPHA, group!="mc")
NO_ALPHA_NO_MCI$group <- factor(NO_ALPHA_NO_MCI$group)
NO_ALPHA_NO_MCI$subID <- factor(NO_ALPHA_NO_MCI$subID)

DAT1<-NO_ALPHA_NO_MCI %>% group_by(subID, group) %>%
  summarize(count=length(Hz))
summary(DAT1$group)


NO_ALPHA_NO_MCI$Group.c <- NO_ALPHA_NO_MCI$group
contrasts(NO_ALPHA_NO_MCI$Group.c) <- contr.poly(2)

NO_ALPHA_NO_MCI$Region.c <- NO_ALPHA_NO_MCI$Region
contrasts(NO_ALPHA_NO_MCI$Region.c) <- contr.poly(4)

# Extracting slopes for the frontal region
rand_mod<-lmer(lgPower~
                 # Fixed-effects
                 1+lgHz+
                 # Random-effects
                 (1+lgHz|subID), data=NO_ALPHA_NO_MCI[NO_ALPHA_NO_MCI$Region=="Frontal",], REML=TRUE,
               control=lmerControl(optimizer="Nelder_Mead",
                                   optCtrl=list(maxfun=2e5)))

FRONTAL_SLOPES<-data.frame(coef(rand_mod)$subID)
head(FRONTAL_SLOPES)
FRONTAL_SLOPES <- cbind(subID = rownames(FRONTAL_SLOPES), FRONTAL_SLOPES)
rownames(FRONTAL_SLOPES) <- 1:nrow(FRONTAL_SLOPES)
FRONTAL_SLOPES$Group <-FRONTAL_SLOPES$subID <- factor(substr(FRONTAL_SLOPES$subID, 
                                                             start = 1, 
                                                             stop = 2))
FRONTAL_SLOPES$Group <- fct_recode(FRONTAL_SLOPES$Group, 
                              "Older Adults" = "oa", 
                              "Younger Adults" = "ya")

FRONTAL_SLOPES$Group <- fct_relevel(FRONTAL_SLOPES$Group, 
                               "Younger Adults", 
                               "Older Adults")

FRONTAL_SLOPES <- rename(FRONTAL_SLOPES, Intercepts=X.Intercept., Slopes=lgHz)
head(FRONTAL_SLOPES)

FRONTAL_SLOPES <- FRONTAL_SLOPES %>% pivot_longer(!c("subID", "Group"), names_to = "Variable", values_to="Values")
head(FRONTAL_SLOPES)


Fig2C<-ggplot(data = FRONTAL_SLOPES, 
       mapping = aes(x = Variable, y = Values)) +
  geom_point(fill="black", pch=21, size=2, position=position_jitter(width=0.1)) +
  geom_boxplot(col="black", fill="grey", alpha=0.2, outlier.shape = NA) +
  facet_wrap(~Group)+
  scale_x_discrete(name = NULL)+
  scale_y_continuous(name = NULL)+
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

Fig2C

## Analysis Comparing Younger to Older Adults to calculate individual intercepts and slopes ---------
m1<-lmer(lgPower~
           # Fixed-effects
           1+lgHz*Region.c*Group.c+
           # Random-effects
           (1+lgHz|subID), data=NO_ALPHA_NO_MCI, REML=TRUE,
         control=lmerControl(optimizer="Nelder_Mead",
                             optCtrl=list(maxfun=2e5)))

anova(m1)
summary(m1)

# Post-Hoc Tests within Age Groups:
OA1<-lmer(lgPower~
           # Fixed-effects
           1+lgHz*Region.c+
           # Random-effects
           (1+lgHz|subID), data=NO_ALPHA_NO_MCI[NO_ALPHA_NO_MCI$Group=="oa",], REML=TRUE,
          control=lmerControl(optimizer="Nelder_Mead",
                              optCtrl=list(maxfun=2e5)))

anova(OA1)
summary(OA1)


YA1<-lmer(lgPower~
            # Fixed-effects
            1+lgHz*Region.c+
            # Random-effects
            (1+lgHz|subID), data=NO_ALPHA_NO_MCI[NO_ALPHA_NO_MCI$Group=="ya",], REML=TRUE,
          control=lmerControl(optimizer="Nelder_Mead",
                              optCtrl=list(maxfun=2e5)))

anova(YA1)
summary(YA1)



# Post-Hoc Tests within Regions:
Frontal1<-lmer(lgPower~
            # Fixed-effects
            1+lgHz*Labels+
            # Random-effects
            (1+lgHz|subID), data=NO_ALPHA_NO_MCI[NO_ALPHA_NO_MCI$Region=="Frontal",], REML=TRUE,
          control=lmerControl(optimizer="Nelder_Mead",
                              optCtrl=list(maxfun=2e5)))
anova(Frontal1)
summary(Frontal1)

Central1<-lmer(lgPower~
                 # Fixed-effects
                 1+lgHz*Labels+
                 # Random-effects
                 (1+lgHz|subID), data=NO_ALPHA_NO_MCI[NO_ALPHA_NO_MCI$Region=="Central",], REML=TRUE,
               control=lmerControl(optimizer="Nelder_Mead",
                                   optCtrl=list(maxfun=2e5)))
anova(Central1)
summary(Central1)



Parietal1<-lmer(lgPower~
                 # Fixed-effects
                 1+lgHz*Labels+
                 # Random-effects
                 (1+lgHz|subID), data=NO_ALPHA_NO_MCI[NO_ALPHA_NO_MCI$Region=="Parietal",], REML=TRUE,
               control=lmerControl(optimizer="Nelder_Mead",
                                   optCtrl=list(maxfun=2e5)))
anova(Parietal1)
summary(Parietal1)

Occipital1<-lmer(lgPower~
                  # Fixed-effects
                  1+lgHz*Labels+
                  # Random-effects
                  (1+lgHz|subID), data=NO_ALPHA_NO_MCI[NO_ALPHA_NO_MCI$Region=="Occipital",], REML=TRUE,
                control=lmerControl(optimizer="Nelder_Mead",
                                    optCtrl=list(maxfun=2e5)))
anova(Occipital1)
summary(Occipital1)









# Analysis #2: Age Related Differences in the RBANS ----------------------------
list.files()

BEHAVIOR <- read.csv("./RBANS_aging_study_06242020.csv", 
                   header = TRUE,na.strings=c("","NA","na"),
                   stringsAsFactors = TRUE)

# Drop subjects who are missing EEG data:
BEHAVIOR <- subset(BEHAVIOR, subID != "ya20")
BEHAVIOR <- subset(BEHAVIOR, subID != "oa22")


# Reliability across RBANS subtest items:
head(BEHAVIOR[,c(7:18)])
psych::alpha(BEHAVIOR[,c(7:18)])
# There is a reasonably high reliability across items, justifying the 
# normalized average scores

summary(BEHAVIOR$group)
# Reformatting the data for graphing:
df<- BEHAVIOR %>% 
  select(subID, group, list_learn, story_mem, fig_copy, line_orn, pic_name,
         sem_flu, dig_spn, coding, list_rec, list_recog, story_rec, fig_rec, raw_ave)%>%
  mutate(raw_ave = raw_ave*100)
head(df)
df<-gather(df, scale, score, list_learn:raw_ave, factor_key=TRUE)
head(df)
summary(df$scale)

SUMMARY<-df %>% group_by(group, scale) %>%
  summarize(mean = mean(score, na.rm=TRUE),
            std = sd(score, na.rm=TRUE),
            count = n(),
            ll = mean(score, na.rm=TRUE)-(sd(score, na.rm=TRUE)/sqrt(n()))*qt(0.975, n()-1, lower.tail = TRUE),
            ul = mean(score, na.rm=TRUE)+(sd(score, na.rm=TRUE)/sqrt(n()))*qt(0.975, n()-1, lower.tail = TRUE))

write.csv(SUMMARY, "./descriptive group stats.csv")

# Plot of RBANS scores by group and test
scale_names <- c(`list_learn` = "LL", `story_mem` = "SM", `fig_copy` = "FC", `line_orn` = "LO",
                 `pic_name` = "PN", `sem_flu` = "SF", `dig_spn` = "DS", `coding` = "Code",
                 `list_rec` = "L Recall", `list_recog` = "L Recog",`story_rec` = "S Recall", `fig_rec` = "F Recall",
                 `raw_ave` = "Average")

df$Labels <- fct_recode(df$group, 
                              "OA" = "OA", 
                              "YA" = "YA",
                              "MCI" = "MC")

df$Labels <- fct_relevel(df$Labels, "YA", "OA", "MCI")

summary(df$group)
df_MCI <- subset(df, group == "MC")
df <- subset(df, group != "MC")

# Figure 2: RBANS Subtests
ggplot(df[df$scale != "raw_ave",], aes(x = Labels, y = score)) +
  geom_point(aes(fill=Labels), pch=21, size=1, stroke=1, col="black",
             position = position_jitterdodge()) +
  geom_boxplot(aes(fill=Labels), col="black", alpha=0.2, outlier.shape = NA)+
  facet_wrap(~scale, nrow=3, scales="free", labeller = as_labeller((scale_names)))+
  scale_x_discrete(name = "Group") +
  scale_y_continuous(name = "Subtest Scores")+
  scale_fill_grey()+
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

t.test(raw_ave   ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less")
t.test(list_learn~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less")
t.test(story_mem ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less")
t.test(fig_copy  ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less")
t.test(line_orn  ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less")
t.test(pic_name  ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less")
t.test(sem_flu  ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less")
t.test(dig_spn  ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less")
t.test(coding  ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less")
t.test(list_rec ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less")
t.test(list_recog  ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less")
t.test(story_rec  ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less")
t.test(fig_rec   ~group, BEHAVIOR[BEHAVIOR$group!="MC",], paired=FALSE, var.equal=FALSE, alternative = "less")




# Analysis #3: Does the spectral slope mediate age-related differences? --------

# Extracting Frontal Spectral Slopes for all Individuals
F1<-lmer(lgPower~
           # Fixed-effects
           1+lgHz+
           # Random-effects
           (1+lgHz|subID), data=NO_ALPHA[NO_ALPHA$Region=="Frontal",], REML=TRUE,
         control=lmerControl(optimizer="Nelder_Mead",
                             optCtrl=list(maxfun=2e5)))
coef(F1)$subID[,1]

SLOPES <- NO_ALPHA %>% group_by(subID) %>%
  summarize()

SLOPES$Frontal_Intercepts<-coef(F1)$subID[,1]
SLOPES$Frontal_Slopes <- coef(F1)$subID[,2]
SLOPES

head(BEHAVIOR)
BEHAVIOR <-merge(BEHAVIOR, SLOPES, by="subID", all=FALSE)
head(BEHAVIOR)

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

fig3A <- ggplot(data = BEHAVIOR, aes(x = Frontal_Slopes, y = raw_ave)) +
  geom_point(aes(fill=Labels), col="black", pch = 21, size = 2) + 
  stat_smooth(aes(col=Labels, lty=Labels), method="lm", se=FALSE, lwd=1) +
  scale_x_continuous(name = "Spectral Slope (Frontal)") +
  scale_y_continuous(name = "Normalized Average (%)")+
  scale_color_grey(start=0.0, end=0.4)+scale_fill_grey()+
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

# Figure 3B: Coding
fig3B <- ggplot(data = BEHAVIOR, aes(x = Frontal_Slopes, y = coding)) +
  geom_point(aes(fill=Labels), col="black", pch = 21, size = 2) + 
  stat_smooth(aes(col=Labels, lty=Labels), method="lm", se=FALSE, lwd=1) +
  scale_color_grey(start=0.0, end=0.4)+scale_fill_grey()+
  scale_x_continuous(name = "Spectral Slope (Frontal)") +
  scale_y_continuous(name = "Coding Score")+
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

# Figure 3C: List Recall
fig3C <- ggplot(data = BEHAVIOR, aes(x = Frontal_Slopes, y = list_rec)) +
  geom_point(aes(fill=Labels), col="black", pch = 21, size = 2) + 
  stat_smooth(aes(col=Labels, lty=Labels), method="lm", se=FALSE, lwd=1) +
  scale_x_continuous(name = "Spectral Slope (Frontal)") +
  scale_y_continuous(name = "List Recall Score")+
  scale_color_grey(start=0.0, end=0.4)+scale_fill_grey()+
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


# Figure 3D: List Recognition
fig3D <- ggplot(data = BEHAVIOR, aes(x = Frontal_Slopes, y = list_recog)) +
  geom_point(aes(fill=Labels), col="black", pch = 21, size = 2) + 
  stat_smooth(aes(col=Labels, lty=Labels), method="lm", se=FALSE, lwd=1) +
  scale_x_continuous(name = "Spectral Slope (Frontal)") +
  scale_y_continuous(name = "List Recognition Score")+
  scale_color_grey(start=0.0, end=0.4)+scale_fill_grey()+
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


# Figure 3E: Figure Recall
fig3E <- ggplot(data = BEHAVIOR, aes(x = Frontal_Slopes, y = fig_rec)) +
  geom_point(aes(fill=Labels), col="black", pch = 21, size = 2) + 
  stat_smooth(aes(col=Labels, lty=Labels), method="lm", se=FALSE, lwd=1) +
  scale_x_continuous(name = "Spectral Slope (Frontal)") +
  scale_y_continuous(name = "List Recognition Score")+
  scale_color_grey(start=0.0, end=0.4)+scale_fill_grey()+
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



(fig3A|fig3B|fig3C)
(fig3D|fig3E)





# Removing the MCI subjects before statistical analysis:
NO_MCI <- subset(BEHAVIOR, group!="MC")

summary(NO_MCI$group)

# Mediation Analysis: RBANS Average --------------------------------------------
# Y = raw_ave
# Mediator = Spectral Slope (Frontal)
# X = Age Group

# Y ~ X
mod1 <- lm(raw_ave~group,
           data=NO_MCI[NO_MCI$raw_ave!="NA",])
summary(mod1)

# M ~ X
mod2 <- lm(Frontal_Slopes~group,
           data=NO_MCI[NO_MCI$raw_ave!="NA",])
summary(mod2)

# Y ~ X + M
mod3 <- lm(raw_ave~group+Frontal_Slopes,
           data=NO_MCI[NO_MCI$raw_ave!="NA",])
summary(mod3)

# Bootstrapping Mediation Analysis:
set.seed(1)
raw_ave.med.result<- mediate(model.m=mod2, model.y=mod3, 
                             treat='group', mediator='Frontal_Slopes',
                             boot=TRUE, sims=5000)
summary(raw_ave.med.result)





# Mediation Analysis: Coding ---------------------------------------------------
# Y = Coding
# Mediator = Spectral Slope (Frontal)
# X = Age Group

# Y ~ X
mod1 <- lm(coding~group,
           data=NO_MCI[NO_MCI$coding!="NA",])
summary(mod1)

# M ~ X
mod2 <- lm(Frontal_Slopes~group,
           data=NO_MCI[NO_MCI$coding!="NA",])
summary(mod2)

# Y ~ X + M
mod3 <- lm(coding~group+Frontal_Slopes,
           data=NO_MCI[NO_MCI$coding!="NA",])
summary(mod3)

# Bootstrapping Mediation Analysis:
set.seed(1)
coding.med.result<- mediate(model.m=mod2, model.y=mod3, 
                     treat='group', mediator='Frontal_Slopes',
                     boot=TRUE, sims=5000)
summary(coding.med.result)
plot(coding.med.result)


# Mediation Analysis: List Recall ----------------------------------------------
# Y = List Recall
# Mediator = Spectral Slope (Frontal)
# X = Age Group

# Y ~ X
mod1 <- lm(list_rec~group,
           data=NO_MCI[NO_MCI$list_rec!="NA",])
summary(mod1)

# M ~ X
mod2 <- lm(Frontal_Slopes~group,
           data=NO_MCI[NO_MCI$list_rec!="NA",])
summary(mod2)

# Y ~ X + M
mod3 <- lm(list_rec~group+Frontal_Slopes,
           data=NO_MCI[NO_MCI$list_rec!="NA",])
summary(mod3)

# Bootstrapping Mediation Analysis:
set.seed(1)
list_recal.med.result<- mediate(model.m=mod2, model.y=mod3, 
                            treat='group', mediator='Frontal_Slopes',
                            boot=TRUE, sims=5000)
summary(list_recal.med.result)
plot(list_recal.med.result)




# Mediation Analysis: List Recognition -----------------------------------------
# Y = List Recognition
# Mediator = Spectral Slope (Frontal)
# X = Age Group

# Y ~ X
mod1 <- lm(list_recog~group,
           data=NO_MCI[NO_MCI$list_recog!="NA",])
summary(mod1)

# M ~ X
mod2 <- lm(Frontal_Slopes~group,
           data=NO_MCI[NO_MCI$list_recog!="NA",])
summary(mod2)

# Y ~ X + M
mod3 <- lm(list_recog~group+Frontal_Slopes,
           data=NO_MCI[NO_MCI$list_recog!="NA",])
summary(mod3)

# Bootstrapping Mediation Analysis:
set.seed(1)
list_recog.med.result<- mediate(model.m=mod2, model.y=mod3, 
                                treat='group', mediator='Frontal_Slopes',
                                boot=TRUE, sims=5000)
summary(list_recog.med.result)
plot(list_recog.med.result)




# Mediation Analysis: Figure Recall --------------------------------------------
# Y = raw_ave
# Mediator = Spectral Slope (Frontal)
# X = Age Group

# Y ~ X
mod1 <- lm(fig_rec~group,
           data=NO_MCI[NO_MCI$fig_rec!="NA",])
summary(mod1)

# M ~ X
mod2 <- lm(Frontal_Slopes~group,
           data=NO_MCI[NO_MCI$fig_rec!="NA",])
summary(mod2)

# Y ~ X + M
mod3 <- lm(fig_rec~group+Frontal_Slopes,
           data=NO_MCI[NO_MCI$fig_rec!="NA",])
summary(mod3)

# Bootstrapping Mediation Analysis:
set.seed(1)
raw_ave.med.result<- mediate(model.m=mod2, model.y=mod3, 
                                treat='group', mediator='Frontal_Slopes',
                                boot=TRUE, sims=5000)
summary(raw_ave.med.result)
