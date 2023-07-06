##### Opening libraries ----------------------------------------------
library(tidyverse); library(car); library(lme4); library(lmerTest); 
library(patchwork); library(ggpubr);library(reshape2);
library(emmeans); library(ez); library(reticulate) 

#Set working directory 
# E.g.,:
setwd("~/GitHub/spectral_slopes_aging/clark_alpha_modulation/")


# 1.0 Visualizing the Data ----
# Before we plot the FOOOF'd data, can we include some plots of the original power spectra?

# 1.1. Plotting the raw power spectra ----
SPECTRA <- read.csv("https://raw.githubusercontent.com/keithlohse/spectral_slopes_aging/main/MASTER_EO_and_EC_EEG_KRL.csv",
                    stringsAsFactors = TRUE, na.strings=c("NA","NaN"," ",""))

# Filter the data to include only those participants who had r^2>0.80 from FOOOF in this project:
summary(SPECTRA$subID)

# Filter the power spectra to 2-30 Hz and just O1, Oz, and O2
head(SPECTRA)
SPECTRA <- SPECTRA %>% select(subID, group, condition,  band, Hz, O1, Oz, O2) %>%
  filter(Hz>=1) %>% # Truncate the low end
  filter(Hz<=31) %>% #truncate the high end
  rowwise() %>% 
  mutate(Occipital=mean(c(O1, Oz, O2), na.rm=TRUE),
         log_Occipital=log(Occipital, base=10))


SPECTRA$condition <- factor(ifelse(SPECTRA$condition == "ec", "Eyes Closed", "Eyes Open"))


# Extracting FOOOF parameters from Pathania et al., 2022 data ------------------
head(SPECTRA)

# if you don't already have these installed
#reticulate::py_install("fooof")

# Import python packages
fooof <- import("fooof")
np <- import("numpy")
plt <- import("matplotlib.pyplot")


# same as python syntax, but use $ instead of . to initialize/access methods and attributes of a class
# see how R types are converted to Python here: https://rstudio.github.io/reticulate/
fm <- fooof$FOOOF(
  peak_width_limits = c(2, 8), 
  max_n_peaks = 8, 
  min_peak_height = 0.1, 
  peak_threshold = 2,
  aperiodic_mode = "fixed",
  verbose = T
)


# need to use np$ravel to get data into format required for fooof module
summary(SPECTRA$subID)
summary(SPECTRA$condition)

# OA04 ----
OA04_EC <- SPECTRA %>% filter(subID == "oa04" && condition =="Eyes Closed")
head(OA04_EC)
OA04_EO <- SPECTRA %>% filter(subID == "oa04" && condition =="Eyes Open")
head(OA04_EO)

fm$fit(freqs = np$ravel(OA04_EC$Hz), 
       power_spectrum = np$ravel(OA04_EC$Occipital), 
       freq_range = c(1.5, 30.5))

results <- fm$get_results()
results

fm$plot()
plt$show()

OA04_EC_results <- tibble::tibble(
  a_params = list(results$aperiodic_params),
  p_params = list(results$peak_params),
  g_params = list(results$gaussian_params),
  r_squared = results$r_squared,
  error = results$error,
  setting_am = fm$aperiodic_mode,
  setting_mnp = fm$max_n_peaks,
  setting_mph = fm$min_peak_height,
  setting_pt = fm$peak_threshold,
  setting_pwl = list(fm$peak_width_limits)
)

fm$report()


OA04_EC_results$a_params
OA04_EC_results$p_params



# OA05 ----
OA05_EC <- SPECTRA %>% filter(subID == "oa05" && condition =="Eyes Closed")
head(OA05_EC)
OA05_EO <- SPECTRA %>% filter(subID == "oa05" && condition =="Eyes Open")
head(OA05_EO)

fm$fit(freqs = np$ravel(OA05_EC$Hz), 
       power_spectrum = np$ravel(OA05_EC$Occipital), 
       freq_range = c(1.5, 30.5))

results <- fm$get_results()
results

fm$plot()
plt$show()

OA05_EC_results <- tibble::tibble(
  a_params = list(results$aperiodic_params),
  p_params = list(results$peak_params),
  g_params = list(results$gaussian_params),
  r_squared = results$r_squared,
  error = results$error,
  setting_am = fm$aperiodic_mode,
  setting_mnp = fm$max_n_peaks,
  setting_mph = fm$min_peak_height,
  setting_pt = fm$peak_threshold,
  setting_pwl = list(fm$peak_width_limits)
)

fm$report()






# Setting up the for-loop ---- 
head(SPECTRA)

EO <- SPECTRA %>% filter(condition=="Eyes Open")
EC <- SPECTRA %>% filter(condition=="Eyes Closed")



EO_RESULTS <- vector(mode="list", length=0)
EC_RESULTS <- vector(mode="list", length=0)

sub_list <- c(unique(EO$subID))
sub_list

for(s in sub_list) {
  print(s)
  
  SUB <- EO %>% filter(subID == s)
  
  fm <- fooof$FOOOF(
      peak_width_limits = c(2, 8), 
      max_n_peaks = 8, 
      min_peak_height = 0.1, 
      peak_threshold = 2,
      aperiodic_mode = "fixed",
      verbose = T
    )
    
    fm$fit(freqs = np$ravel(SUB$Hz), 
           power_spectrum = np$ravel(SUB$Occipital), 
           freq_range = c(2, 30))
    
    results <- fm$get_results()
    
    # nested tibble of model information
    EO_RESULTS[[s]] <- tibble::tibble(
      ps = "EO",
      a_params = list(results$aperiodic_params),
      p_params = list(results$peak_params),
      g_params = list(results$gaussian_params),
      r_squared = results$r_squared,
      error = results$error,
      setting_am = fm$aperiodic_mode,
      setting_mnp = fm$max_n_peaks,
      setting_mph = fm$min_peak_height,
      setting_pt = fm$peak_threshold,
      setting_pwl = list(fm$peak_width_limits)
    )
}

# Subject Level
RESULTS_DF <- bind_rows(EO_RESULTS, .id = "column_label")

RESULTS_DF <-RESULTS_DF %>% unnest_wider(a_params, names_sep = "_") %>%
    rename_at(vars(contains("a_params")), ~ c("exponent", "offset"))

# max(lengths(RESULTS_DF$p_params))
# max(lengths(RESULTS_DF$p_params))/3
# 
# a<-c("CF_", "PW_", "BW_")
# b<-c(sort(rep(seq(from=1, 
#           to = max(lengths(RESULTS_DF$p_params))/3),
#          3)))
# paste(a, b, sep="")

RESULTS_DF <-RESULTS_DF %>% unnest_wider(p_params, names_sep = "_") %>%
  rename_at(vars(contains("p_params")), ~ c(paste("Peak_",
                                                  seq(from=1, 
                                                      to = max(lengths(RESULTS_DF$p_params))/3),
                                                  sep=""))
  )

  
RESULTS_DF <-RESULTS_DF %>% unnest_wider(setting_pwl, names_sep = "_") %>%
    rename_at(vars(contains("setting_pwl")), ~ c("peak_width_LL", "peak_width_UL"))
  

write.csv(RESULTS_DF %>% select(-g_params), "./EO_results.csv")

head(SPECTRA)
ggplot(data=SPECTRA[SPECTRA$subID=="oa25",], aes(x=Hz, y=log_Occipital)) +
  geom_rect(aes(xmin=8, xmax=12, ymin=-Inf, ymax=Inf), fill="grey90")+
  geom_line(aes(group=subID), col="black", lwd=0.5) +
  facet_wrap(~condition) +
  scale_x_continuous(name = "Frequency (Hz)") +
  scale_y_continuous(name = "Power log_10(uV^2)") +
  #scale_color_manual(values=c("black", "red"))+
  theme_bw()+
  theme(axis.text=element_text(size=10, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "bottom")




