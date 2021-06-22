# calibration curve for divex water experiment (file: calibration_curve.txt)
# construct calibration curve and fit linear model

# Predict copy number from C_t values from divex water 
# samples (file:divex_water_samples.txt) using inverse prediction
# with the inverse.predict() function from the chemCal package

# plot the estimated C_t values from water sample against dilution
# and add lines for LOD and LOQ

# LOD and LOQ were estimated from the calibration curve data 
# using the LoD-calculator.R script (changed DAT to calib_t) 

# code for Figure 3

# created: # Fri Jun 19 15:53:17 2020 ------------------------------
# last modified: # # Tue Jun 23 11:03:45 2020 ------------------------------

# required packages
library(ednar)
library(ggplot2)
library(dplyr)
library(chemCal)

# data import
# data for calibration curve
calib <- read.table('data/calibration_curve.txt', header =TRUE, sep = "\t",
                    stringsAsFactors = TRUE)
# data from qPCR water tanks
divex <- read.table('data/divex_water_samples.txt',header = TRUE, sep = "\t", 
                    stringsAsFactors = TRUE)

# data from master calibration curve for LOD, LOD3 and LOQ estimates
pcr <- read.table('data/positive_cont.txt', header = TRUE, sep = "\t", 
                  stringsAsFactors = TRUE)

# data wrangling for calibration curve construction
names(calib)[c(1, 3, 5)] <- c("SQ", "Cq", "Target") 
names(divex)[c(5, 6)] <- c("Ct.value", "calib.curve")

# plot calibration curve
calib_plot(calib, target ='01052020IM1')

# Alternative plot from scratch 
ggplot(calib, aes(x = SQ, y = Cq)) +
  geom_point(alpha = 0.65, size = 2) +
  scale_x_continuous(trans = 'log10', labels = scales::comma) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  xlab(expression(standard~concentrations~(copies~reaction^{-1}))) +
  ylab("Cq - value") +
  theme_light() +
  theme(axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        plot.margin = unit(c(5.5, 15, 5.5, 5.5), "points"))

# estimate predicted copy number from Cq values using the calib_predict
# function from ednar package
preds <- calib_predict(calib_df = calib, ct_df = divex)

# generate LOD, LOD3 and LOQ values using the calib_lod function from the 
# ednar package
lod_estimates <- calib_lod(pcr)
LOD <- lod_estimates$assaySum[1, 6]
LOD3 <- lod_estimates$assaySum[1, 9]
LOQ <- lod_estimates$assaySum[1, 7]

# plot estimated copy number against dilution for the divex
# water data
ggplot(preds, aes(x = sample_dilution, y = CN.back)) +
  geom_point(alpha = 0.65, size = 2) +
  scale_x_log10(labels = c("0.00001", "0.0001", 0.001, 0.01, 0.1, 1), breaks = c(0.00001,0.0001, 0.001, 0.01, 0.1, 1)) +
  scale_y_log10(labels = c(1, 10, 100, 1000, 10000), breaks = c(1, 10, 100, 1000, 10000)) +
  xlab("mesocosm DNA dilution") +
  ylab(expression(estimated~copy~number~(copies~reaction^{-1}))) +
  geom_hline(yintercept = LOD, color = "black", size = 0.7) +
  geom_hline(yintercept = LOQ, color = "red", size = 0.7, linetype = "dashed") +
  geom_hline(yintercept = LOD3, color = "black", size = 0.7) +
  annotate(geom = "text", x = 0.9, y = 13, label = "LOD", colour = "black",  size = 4) +
  annotate(geom = "text", x = 0.9, y = 27, label = "LOQ", colour = "red",  size = 4) +
  annotate(geom = "text", x = 0.9, y = 1.7, label = "LOD3", colour = "black",  size = 4) +
  theme_light() +
  theme(axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        plot.margin = unit(c(5.5, 15, 10, 5.5), "points"))

ggsave("output/est_CN_vs_water_dilution.pdf", device = "pdf")