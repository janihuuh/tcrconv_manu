
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(RColorBrewer)

me=system("whoami", intern = T)
setwd(paste0("/Users/", me, "/Dropbox/covid19/"))

theme_set(theme_classic(base_size = 12))

add_guide    <- guides(colour = guide_legend(override.aes = list(size=5)))
facets_nice3 <- theme(strip.text.x = element_text(size = 8, angle = 0, hjust = 0), strip.background = element_rect(colour="white", fill="white"))

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))

## Run all fun_* codes
for(code in list.files("src/", "fun", full.names = T, recursive = T)){
  
  message(code)
  source(code)
  
}
