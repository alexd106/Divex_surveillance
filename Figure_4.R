# Plot predicted copy number for samples collected monthly 2017 - 2019 at 
# four different sites.

# created : # Fri Jun 19 13:05:52 2020 ------------------------------
# last modified: # Tue Jun 22 13:24:58 2021 ------------------------------
# created by: Alex Douglas

# code for Figure 4

# required packages
library(ednar)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(lubridate)
library(RColorBrewer)

# data import 
# data for calibration curve
calib_dat <- read.table('data/dvex_eDNA_calib.txt', header = TRUE, sep = '\t',
                        stringsAsFactors = TRUE)
# field data 
field_dat <- read.table('data/dvex_eDNA_field.txt', header = TRUE, sep = '\t',
                        stringsAsFactors = TRUE)

# data from master calibration curve for LOD, LOD3 and LOQ estimates
pcr <- read.table('data/positive_cont.txt', header = TRUE, sep = "\t", 
                  stringsAsFactors = TRUE)

# LOD, LOD3 and LOQ estimates for each site
divex_lod <- calib_lod(calib_dat)

# warning: Potential outliers detected: 
# Target   SQ rep.number     Cq
# 340 90720203 5548          1 36.751

# There was an error in data transcription. The Ct value should be 26.751 not 
# 36.751. Confirmed with IM Fri Jun 19 2020
calib_dat$Cq[340] <- 26.751

# re-estimate LOD, LOD3 and LOQ
divex_lod <- calib_lod(calib_dat)

# LOD, LOD3 and LOQ based on master calibration curve
divex_lod_master <- calib_lod(pcr)
LOD <- divex_lod_master$assaySum[1, 6]
LOD3 <- divex_lod_master$assaySum[1, 9]
LOQ <- divex_lod_master$assaySum[1, 7]

# predict CN from Ct for field samples
preds <- calib_predict(calib_df = calib_dat, ct_df = field_dat)

# All locations plotted with log10 mean CN axis
filt_dat <- preds %>%
  mutate(labs = paste(Location, Year, Month, Sampling.site, Filter.number, sep = '_')) %>%
  filter(Ct.value < 41) %>%
  add_count(labs) %>%
  mutate(rate = n/3) %>%
  filter(rate > 0.5)  # remove when detections < 2

filt_ave <- filt_dat %>%
  group_by(Location, Year, Month, Sampling.site, Filter.number) %>%
  summarise(mean_CN = mean(CN.back, na.rm = TRUE)) %>%
  mutate(Time = as.Date(ymd(paste0(Year,Month,"01")), format ="%m.%Y")) %>%
  mutate(Month_n = month(Time))

colourCount <- length(unique(preds$Sampling.site))
getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))

p_all <- ggplot(filt_ave, aes(x = Month_n, y = mean_CN)) +
  geom_jitter(aes(colour = Sampling.site), width = 0.2, size = 1) +
  facet_grid(Location ~ Year) + 
  geom_hline(yintercept = LOD3, color = "darkgrey", linetype = 2, size = 0.7) +
  geom_hline(yintercept = LOQ, color = "darkgrey", linetype = 3, size = 0.7) +
  scale_colour_manual(values = getPalette(colourCount)) +
  xlab("Month") +
  ylab(expression(mean~copy~number~(copies~reaction^{-1}))) +
  scale_x_continuous(breaks = 2:11, labels=c("Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")) +
  scale_y_continuous(trans = "log10") +
  # coord_trans(y = "log2") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.spacing.x = unit(0.2, "lines"), panel.spacing.y = unit(0.2, "lines")) +
  theme(strip.background = element_rect(colour = "darkgrey", fill = "white")) +
  theme(strip.text = element_text(colour = "black", size = 13)) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.title = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.position = "none")
p_all
ggsave('output/figures_pdf/all_locations_CN_log10_filtered.pdf', device = 'pdf', dpi = 600)
ggsave('output/figures_tiff/all_locations_CN_log10_filtered.tiff', device = 'tiff', dpi = 100)

# individual plots for each sampling location 
# Creran data 
filt1_creran_dat <- preds %>%
  mutate(labs = paste(Location, Year, Month, Sampling.site, Filter.number, sep = '_')) %>%
  filter(Ct.value < 41) %>%
  add_count(labs) %>%
  mutate(rate = n/3) %>%
  filter(rate > 0.5 & Location == "Creran")

filt1_creran_ave <- filt1_creran_dat %>%
  group_by(Location, Year, Month, Sampling.site, Filter.number) %>%
  summarise(mean_CN = mean(CN.back, na.rm = TRUE)) %>%
  mutate(Time = as.Date(ymd(paste0(Year,Month,"01")), format ="%m.%Y")) %>%
  mutate(Month_n = month(Time))

# plot for Creran
p <- ggplot(filt1_creran_ave, aes(x = Month_n, y = mean_CN, colour = Sampling.site)) +
  geom_jitter(width = 0.2, size = 1) +
  facet_grid(~ Year) + 
  geom_hline(yintercept = LOD3, color = "darkgrey", linetype = 2, size = 1) +
  scale_colour_brewer(palette = "Dark2") +
  xlab("Month") +
  ylab(expression(mean~copy~number~(copies~reaction^{-1}))) +
  scale_x_continuous(breaks = c(5, 6, 7, 8, 9, 10, 11), labels=c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")) +
  scale_y_continuous(trans = "log10") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.spacing.x=unit(0.2, "lines")) +
  theme(strip.background =element_rect(colour="darkgrey", fill="lightgrey")) +
  theme(strip.text = element_text(colour = "black", size = 13)) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.title = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.position = "none")
p
ggsave('output/figures_pdf/creran_CN_log10_filtered.pdf', device = 'pdf', dpi = 600)
ggsave('output/figures_tiff/creran_CN_log10_filtered.tiff', device = 'tiff', dpi = 100)

# Fairlie data
filt1_Fairlie_dat <- preds %>%
  mutate(labs = paste(Location, Year, Month, Sampling.site, Filter.number, sep = '_')) %>%
  filter(Ct.value < 41) %>%
  add_count(labs) %>%
  mutate(rate = n/3) %>%
  filter(rate > 0.5 & Location == "Fairlie")

filt1_Fairlie_ave <- filt1_Fairlie_dat %>%
  group_by(Location, Year, Month, Sampling.site, Filter.number) %>%
  summarise(mean_CN = mean(CN.back, na.rm = TRUE)) %>%
  mutate(Time = as.Date(ymd(paste0(Year,Month,"01")), format ="%m.%Y")) %>%
  mutate(Month_n = month(Time))

p_f <- ggplot(filt1_Fairlie_ave, aes(x = Month_n, y = mean_CN, colour = Sampling.site)) +
  geom_jitter(width = 0.2, size = 3) +
  facet_grid(~ Year) + 
  geom_hline(yintercept = LOD3, color = "darkgrey", linetype = 2, size = 1) +
  scale_colour_brewer(palette = "Dark2") +
  xlab("Month") +
  ylab(expression(mean~copy~number~(copies~reaction^{-1}))) +
  scale_x_continuous(breaks = 2:10, labels=c("Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")) +
  scale_y_continuous(trans = "log10") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.spacing.x=unit(0.2, "lines")) +
  theme(strip.background =element_rect(colour="darkgrey", fill="lightgrey")) +
  theme(strip.text = element_text(colour = "black", size = 13)) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.title = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.position = "none")
p_f
ggsave('output/figures_pdf/fairlie_CN_log10_filtered.pdf', device = 'pdf', dpi = 600)
ggsave('output/figures_tiff/fairlie_CN_log10_filtered.tiff', device = 'tiff', dpi = 100)

# Portavadie data
filt1_Portavadie_dat <- preds %>%
  mutate(labs = paste(Location, Year, Month, Sampling.site, Filter.number, sep = '_')) %>%
  filter(Ct.value < 41) %>%
  add_count(labs) %>%
  mutate(rate = n/3) %>%
  filter(rate > 0.5 & Location == "Portavadie")

filt1_Portavadie_ave <- filt1_Portavadie_dat %>%
  group_by(Location, Year, Month, Sampling.site, Filter.number) %>%
  summarise(mean_CN = mean(CN.back, na.rm = TRUE)) %>%
  mutate(Time = as.Date(ymd(paste0(Year,Month,"01")), format ="%m.%Y")) %>%
  mutate(Month_n = month(Time))

p_p <- ggplot(filt1_Portavadie_ave, aes(x = Month_n, y = mean_CN, colour = Sampling.site)) +
  geom_jitter(width = 0.2, size = 3) +
  facet_grid(~ Year) + 
  geom_hline(yintercept = LOD3, color = "darkgrey", linetype = 2, size = 1) +
  scale_colour_brewer(palette = "Dark2") +
  xlab("Month") +
  ylab(expression(mean~copy~number~(copies~reaction^{-1}))) +
  scale_x_continuous(breaks = 2:11, labels=c("Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")) +
  scale_y_continuous(trans = "log10") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.spacing.x=unit(0.2, "lines")) +
  theme(strip.background =element_rect(colour="darkgrey", fill="lightgrey")) +
  theme(strip.text = element_text(colour = "black", size = 13)) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.title = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.position = "none")
p_p
ggsave('output/figures_pdf/portavadie_CN_log10_filtered.pdf', device = 'pdf', dpi = 600)
ggsave('output/figures_tiff/portavadie_CN_log10_filtered.tiff', device = 'tiff', dpi = 100)

# Largs data
filt1_Largs_dat <- preds %>%
  mutate(labs = paste(Location, Year, Month, Sampling.site, Filter.number, sep = '_')) %>%
  filter(Ct.value < 41) %>%
  add_count(labs) %>%
  mutate(rate = n/3) %>%
  filter(rate > 0.5 & Location == "Largs")

filt1_Largs_ave <- filt1_Largs_dat %>%
  group_by(Location, Year, Month, Sampling.site, Filter.number) %>%
  summarise(mean_CN = mean(CN.back, na.rm = TRUE)) %>%
  mutate(Time = as.Date(ymd(paste0(Year,Month,"01")), format ="%m.%Y")) %>%
  mutate(Month_n = month(Time))

p_l <- ggplot(filt1_Largs_ave, aes(x = Month_n, y = mean_CN, colour = Sampling.site)) +
  geom_jitter(width = 0.2, size = 3) +
  facet_grid(~ Year) + 
  geom_hline(yintercept = LOD3, color = "darkgrey", linetype = 2, size = 1) +
  scale_colour_brewer(palette = "Dark2") +
  xlab("Month") +
  ylab(expression(mean~copy~number~(copies~reaction^{-1}))) +
  scale_x_continuous(breaks = 2:10, labels=c("Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")) +
  scale_y_continuous(trans = "log10") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(panel.spacing.x=unit(0.2, "lines")) +
  theme(strip.background =element_rect(colour="darkgrey", fill="lightgrey")) +
  theme(strip.text = element_text(colour = "black", size = 13)) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  theme(axis.title = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.position = "none")
p_l
ggsave('output/figures_pdf/largs_CN_log10_filtered.pdf', device = 'pdf', dpi = 600)
ggsave('output/figures_tiff/largs_CN_log10_filtered.tiff', device = 'tiff', dpi = 100)

# output data files
# all predicted copy number data where CT.value greater than 41 and less than 
# two technical replicates excluded 
filt1_all_dat <- preds %>%
  mutate(labs = paste(Location, Year, Month, Sampling.site, Filter.number, sep = '_')) %>%
  filter(Ct.value < 41) %>%
  add_count(labs) %>%
  mutate(rate = n/3) %>%
  filter(rate > 0.5) # add to remove n < 2 replicates

write.csv(filt1_all_dat, 'output/data_summaries/filtered_allData.csv')

# average predicted copy number (averaged over technical replicates) where 
# CT.value greater than 41 and less than two technical replicates excluded 
filt1_all_ave <- filt1_all_dat %>%
  group_by(Location, Year, Month, Sampling.site, Filter.number) %>%
  summarise(mean_CN = mean(CN.back, na.rm = TRUE)) %>%
  mutate(Time = as.Date(ymd(paste0(Year,Month,"01")), format ="%m.%Y")) %>%
  mutate(Month_n = month(Time))

write.csv(filt1_all_ave, 'output/data_summaries/filtered_meanData.csv')

# lod calcs to output summary
write.csv(divex_lod$assaySum, 'output/data_summaries/field_data_assay_summary.csv')

# all standard curve plots for each sampling occasion
tar_name <- unique(calib_dat$Target)[1:8]
divex_1_8 <- calib_plot_all(calib_dat[calib_dat$Target %in% tar_name,])
ggsave('output/figures_pdf/calibration_curves_1_8.pdf', divex_1_8, device = 'pdf')
ggsave('output/figures_tiff/calibration_curves_1_8.tiff', divex_1_8, device = 'tiff', dpi = 100)

tar_name <- unique(calib_dat$Target)[9:16]
divex_9_16 <- calib_plot_all(calib_dat[calib_dat$Target %in% tar_name,])
ggsave('output/figures_pdf/calibration_curves_9_16.pdf', divex_9_16, device = 'pdf')
ggsave('output/figures_tiff/calibration_curves_9_16.tiff', divex_9_16, device = 'tiff', dpi = 100)

tar_name <- unique(calib_dat$Target)[17:21]
divex_17_19 <- calib_plot_all(calib_dat[calib_dat$Target %in% tar_name,])
mod_17_19 <- divex_17_19 + plot_spacer() + plot_spacer() + plot_spacer() 
ggsave('output/figures_pdf/calibration_curves_17_21.pdf', mod_17_19, device = 'pdf')
ggsave('output/figures_tiff/calibration_curves_17_21.tiff', mod_17_19, device = 'tiff', dpi = 100)

