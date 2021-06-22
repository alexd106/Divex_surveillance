# create standard curve plots from qPCR data 
# Add LOD, LOD3 and LOQ values onto each plot.
# code for Figure 2

# LOD, LOD3 and LOQ values estimated using ednar package
# https://alexd106.github.io/ednar/

# created : # Fri Jun 19 13:05:52 2020 ------------------------------
# last modified: # Tue Jun 22 11:35:03 2021 ------------------------------
# created by: Alex Douglas

# required packages
library(ggplot2)
library(ednar)

# data import
pcr <- read.table('data/positive_cont.txt', header = TRUE, sep = "\t", 
                  stringsAsFactors = TRUE)

# generate LOD, LOD3 and LOQ values for the plot annotation using the calib_lod
# function from the ednar package
lod_est <- calib_lod(pcr)

# write this file as output
write.csv(lod_est$assaySum, 'output/data_summaries/assay_summary.csv')

# extract LOD, LOD3 and LOQ values
LOD <- lod_est$assaySum[1, 6]
LOD3 <- lod_est$assaySum[1, 9]
LOQ <- lod_est$assaySum[1, 7]

##  Generate Figure 2
SC.plot <- ggplot(pcr, aes(x = SQ, y = Cq)) +
  geom_point(alpha = 0.65, size = 2) +
  scale_x_continuous(trans = 'log10', labels = scales::comma) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  geom_vline(xintercept = LOD, color = "black", size = 0.7) +
  geom_vline(xintercept = LOD3, color = "black", size = 0.7) +
  geom_vline(xintercept = LOQ, color = "red", linetype = "dashed", size = 0.7) +
  annotate(geom = "text", x = 11.5, y = 43, label = "LOD", colour = "black", angle = 270, size = 5) +
  annotate(geom = "text", x = 1.5, y = 43, label = "LOD3", colour = "black", angle = 270, size = 5) +
  annotate(geom = "text", x = 24.5, y = 43, label = "LOQ", colour = "red", angle = 270, size = 5) +
  xlab(expression(standard~concentrations~(copies~reaction^{-1}))) +
  ylab("Cq - value") +
  theme_light() +
  theme(axis.title = element_text(size = 13),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11))

ggsave("output/figures_pdf/standard_curve_LOD_LOD3_LOQ.pdf", device = "pdf")
ggsave("output/figures_tiff/standard_curve_LOD_LOD3_LOQ.tiff", device = "tiff", dpi = 150)
