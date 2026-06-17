# Graphics, tables and summary statistics for publication

# Packages ----------------------------------------------------------------
packages <- c("dplyr","tidyr","magrittr", "terra", "tidyterra","viridis", "ggplot2", 
              "stringr", "cowplot", "patchwork", "ggExtra", "purrr")
sapply(packages, FUN = library, character.only = T)


# Font --------------------------------------------------------------------
library(showtext)
font_add_google("Inter", "inter")
showtext_auto()

# data --------------------------------------------------------------------
dt.f <- rast("data/ForestArea_germany_100m.tif")
bay.f <- rast("data/ForestArea_bavaria_100m.tif")

germany <- vect("data/DE_Verwaltungsgebiete5000/vg5000_ebenen_1231/VG5000_KRS.shp") %>% 
  aggregate() %>% 
  project(., dt.f)
bavaria <- vect("data/DE_Verwaltungsgebiete5000/vg5000_ebenen_1231/VG5000_LAN.shp") %>% 
  filter(GEN == "Bayern") %>% 
  project(.,dt.f)
ofr <- vect("data/DE_Verwaltungsgebiete5000/vg5000_ebenen_1231/VG5000_KRS.shp") %>% 
  filter(GEN %in% c("Coburg", "Kronach", "Hof", "Lichtenfels", "Kulmbach", "Wunsiedel i.Fichtelgebirge", "Bamberg", "Bayreuth", "Forchheim")) %>% 
  project(.,dt.f)


# Species
source("script/Species_select.R")
species.final <- readRDS("output/Fits/Sapling/h50d7_Germany/Sapling_model_final.rds")$species
species.tab <- read.csv("data/DE_BWI3_regeneration_species.csv")

# Colors
sunset = colorRampPalette(c("#FFEC9DFF", "#F2AF4AFF", "#EB7F54FF", "#C36377FF", "#61599DFF", "#1D457F", "#191F40FF", "black"))
cult.col = colorRampPalette(c("grey90", "#FCFD8F","#F3CE65","#EB9F3C","#9A3F07"))
div = colorRampPalette(c("grey90","#FFEC9DFF", "#F2AF4AFF", "#EB7F54FF", "#9A3F07"))


# GERMANY no hist ---------------------------------------------------------

# SPECIES RICHNESS ---------------------------------------------------------------
# species richness rule BaySF Waldbauhandbuch Baumartenwahl 2020
# 4 species per hectare for species with minimum abundances of >=5 %
# 3 species Lindner et al 2020 (policy brief)
sprich.df <- readRDS("output/Graphs/Speciesrichness.rds")
sprich.rast <- rast("output/Graphs/Speciesrichness.tif")


## Plot --------------------------------------------------------------------

### Class -------------------------------------------------------------------

p.div <- 
  ggplot() +
  theme_void() +
  geom_spatraster(data = sprich.rast$class) +
  scale_fill_manual(      
    "Species richness",
    na.translate = FALSE,
    labels = c("1-2", "3-4", "\u22655", ""),
    na.value = "transparent",
    values = c("#9A3F07","#F2AF4AFF","grey90"),
    guide = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) +
  geom_spatvector(data = germany, fill = "transparent", colour = "black", linewidth = 0.1)+ 
  theme(text = element_text(family = "inter"))

p.div

ggsave(filename = "output/Graphs/Figure4_nohist.pdf",
       height = 8, width = 10, units = "cm", dpi = 900,
       bg = "white",
       device="pdf")


# OBERFRANKEN -------------------------------------------------------------

# REGENERATION DISTRIBUTION ------------------------------------------
# Visualising predicted species distribution Maps

## Distribution Species ----------------------------------------------------------------
sunset = colorRampPalette(c("#FFEC9DFF", "#F2AF4AFF", "#EB7F54FF", "#C36377FF", "#61599DFF", "#1D457F", "#191F40FF", "black"))

### Plot --------------------------------------------------------------------
sapling.map <- function(species.vect=species.vect, source, scale, scale.plot, max.count, text.size){
  for(species in species.vect){
    print(species)
    # load
    pred <- rast(paste0("output/Predictions/Regeneration_",species,".tif")) %>% 
      rename(count_ha = paste0(species)) %>% 
      select(count_ha)
    # Sapling
    if (scale != "Germany") {
      pred %<>% terra::crop(., scale.plot, mask = T)}
    # plot
    p <-
      ggplot()+
      theme_void()+
      geom_spatraster(data = pred) +
      scale_fill_gradientn(
        "Regeneration\ndensity [ha\u207B\u00B9]",
        na.value = "transparent",
        colors = sunset(7),
        space = "Lab",
        trans = "log1p",
        breaks = c(0,10,100,1000,10000, 100000, 1e6,1e7),#Sapling
        labels = scales::comma(c(0,10,100,1000,10000, 100000,1e6,1e7)),#Sapling
        limits = c(0, max.count))+
      guides(fill = guide_colourbar(barwidth = 1))+
      theme(legend.title = element_text(size = 8),
            legend.text = element_text(size = 7))+
      geom_spatvector(data = scale.plot, fill = "transparent", colour = "black", linewidth = 0.1)+
      annotate("text", x = -Inf, y = Inf, hjust=-0.1, vjust = 1, size = text.size, label = paste(species.tab[species.tab$name.id==species,]$name.scient), fontface = 'bold.italic')
    
    assign(paste0("p.",species), p , envir = .GlobalEnv)
    
    # save
    # ggsave(plot = p,
    #        filename = paste0("output/Graphs/Regeneration_",scale,"_",species,".png"),
    #        height = 8, width = 8, units = "cm", dpi = 900,
    #        bg = "white",
    #        device=grDevices::png)
  }
}



## Selected species --------------------------------------------------------
species.sel = c("Picea.abies","Abies.alba", "Fagus.sylvatica")

# Regeneration stack
regstack <- list()
for (species in species.final){
  regstack[[species]] <- rast(paste0("output/Predictions/Regeneration_",species,".tif")) %>% 
    select(all_of(species))
}
regstack <- rast(regstack) # Turn list into multilayer



# Regstack of selected species
regstack.sel <- regstack %>%
  select(all_of(species.sel))

# Oberfranken
sapling.map(species.vect = c("Picea.abies","Abies.alba", "Fagus.sylvatica"),
            text.size = 3, scale = "OFR", scale.plot = ofr, max.count = 56075) # keep max of germany

fag <- p.Fagus.sylvatica + theme(text = element_text(family = "inter"))
pic <- p.Picea.abies     + theme(text = element_text(family = "inter"))
ab <- p.Abies.alba      + theme(text = element_text(family = "inter"))

final.ofr = fag + pic + ab +  plot_layout(guides = 'collect')

ggsave(final.ofr, filename = "output/Graphs/Figure2_ofr.pdf",
       height = 7, width = 25, units = "cm", dpi = 900,
       bg = "white",
       device = "pdf")


# FUTURE SUITABILITY ---------------------------------------------------------
## Species -----------------------------------------------------------------
## Load -----------------------------------------------------------------
species.tab <- read.csv("data/DE_BWI3_regeneration_species.csv")

species.final <- readRDS("output/Fits/Sapling/h50d7_Germany/Sapling_model_final.rds")$species
# species.final <- "Abies.alba" # for test reasons
# species = "Abies.alba"

cult.col = colorRampPalette(c("grey90", "#FCFD8F","#F3CE65","#EB9F3C","#9A3F07"))

sap.risk <- rast("output/Suitability/Regeneration_suitability.tif") %>% 
  select(higher) %>% 
  terra::crop(., ofr, mask = T)
  
crs(sap.risk)==crs(ofr)


## Stats-------------------------------------------------

sap.risk.df <- 
  as.data.frame(sap.risk.ofr)


sap.risk.df.nona <- sap.risk.df %>% 
  drop_na()


## Plot --------------------------------------------------------------
# Map
b.risk <-
  ggplot()+
  theme_void()+
  geom_spatvector(data = ofr, fill = "white", colour = "black", linewidth = 0.1) +
  geom_spatraster(data = sap.risk) +
  scale_fill_gradientn("Proportion of\nregeneration of\nlow future suitability\n[%]",
                       colours = cult.col(50),
                       na.value = NA,
                       limits = c(0, 100)) +
  geom_spatvector(data = ofr, fill = "transparent", colour = "black", linewidth = 0.1) +
  theme(text = element_text(family = "inter")) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA))

ggsave(b.risk,
       filename = "output/Graphs/Risk_ofr.pdf",
       height = 8, width = 14, units = "cm", dpi = 900,
       bg = "transparent",
       device="pdf")



# BOX BAVARIA: TOTAL, DIVERSITY and CULT RISK -------------------------------------------------
## Cult risk ---------------------------------------------------------------

# Map
b.risk <-
  ggplot()+
  theme_void()+
  geom_spatraster(data = sap.risk) +
  scale_fill_gradientn("Proportion of\nregeneration of\nlow future suitability\n[%]",
                       colours = cult.col(50),
                       na.value = "transparent",
                       limits = c(0, 100)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) +
  geom_spatvector(data = ofr, fill = "transparent", colour = "black", linewidth = 0.1) +
  theme(text = element_text(family = "inter"))

# Histogram
b.risk.h <-
  ggplot(sap.risk.df.nona, aes(x = higher)) +
  theme_classic() +
  geom_histogram(binwidth = 5,
                 boundary = 0,#-0.5
                 fill=cult.col(20), color="black") +
  xlab("Proportion of regeneration of\nlow future suitability [%]") +
  ylab("Area [10\u00B3 ha]")+
  scale_y_continuous(limits=c(0,7e4),
                     breaks = c(0,2e4,4e4,6e4),
                     labels = c(0,20,40,60)) +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)) +
  annotate("segment", x = round(median(sap.risk.df.nona$higher, na.rm =T), digits=0), xend = 0, y = 65e3, yend = 65e3,
           arrow = arrow(ends = "both", angle = 90, length = unit(.15,"cm")),
           colour = "#1D457F") +
  annotate("text", x = round(median(sap.risk.df.nona$higher, na.rm =T))+ 2, y = 65e3, 
           hjust = 0, vjust = 0.5,
           label = paste(strwrap(paste0("50% of the forest area has <", round(median(sap.risk.df.nona$higher, na.rm =T), digits=1),"% of low future suitability"), 26), collapse = "\n"),
           colour = "#1D457F",
           size = 2) +
  annotate("segment", x = 100, xend = 75, y = 2e4, yend = 2e4,
           arrow = arrow(ends = "both", angle = 90, length = unit(.15,"cm")),
           colour = "#1D457F")+
  annotate("text", x = 100, y = 2.5e4, 
           hjust = 1, vjust = -0.1,
           label = paste(strwrap(paste0(round(table(sap.risk.df.nona$higher >= 75)[2]/dim(sap.risk.df.nona)[1]*100, digits=1),"% forest area has a high proportion of regeneration of low suitability (\u226575%)"), 30), collapse = "\n"),
           colour = "#1D457F",
           size = 2)+
  theme(text = element_text(family = "inter"))

## Total Density ---------------------------------------------------------------------
regtot.bay <- rast("output/Graphs/Density.tif") %>% 
  project(.,sap.risk) %>% 
  crop(., sap.risk, mask = T)

# Categorize the raster values into classes
regtot.bay$class <- classify(regtot.bay, c(0, 1000, 2000, Inf), right= F)

b.tot <- 
  ggplot() + 
  theme_void() +
  geom_spatraster(data = regtot.bay$class)+
  scale_fill_manual(
    "Regeneration\ndensity [ha\u207B\u00B9]",
    na.translate = FALSE,
    labels = c("0-1000", "1000-2000", "\u22652000"),
    values = c("#9A3F07","#F2AF4AFF","grey90"),
    guide = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) +
  geom_spatvector(data = ofr, fill = "transparent", colour = "black", linewidth = 0.1) +
  theme(text = element_text(family = "inter"))

regtot.bay.df <- as.data.frame(regtot.bay)
length(regtot.bay.df[regtot.bay.df$count_tot_ha <= 1000,]$count_tot_ha)/length(regtot.bay.df$count_tot_ha)*100


## Sp rich --------------------------------------------------------
sprich.bay <- rast("output/Graphs/Speciesrichness.tif") %>% 
  project(.,sap.risk) %>% 
  crop(., sap.risk, mask = T)
sprich.bay$class <- classify(sprich.bay$sprich, c(0, 2, 4, Inf))

b.sprich <- 
  ggplot() +
  theme_void() +
  geom_spatraster(data = sprich.bay$class) +
  scale_fill_manual(
    "Species richness",
    na.translate = FALSE,
    labels = c("1-2", "3-4", "\u22655"),
    values = c("#9A3F07","#F2AF4AFF","grey90"),
    guide = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) +
  geom_spatvector(data = ofr, fill = "transparent", colour = "black", linewidth = 0.1) +
  theme(text = element_text(family = "inter"))

# Stats
sprich.bay.df <- as.data.frame(sprich.bay)
table(sprich.bay.df$sprich <= 2)[2]/length(sprich.bay.df$sprich)*100

## All plots ---------------------------------------------------------------
b.tot + b.sprich + b.risk + b.risk.h +
  plot_layout(ncol = 2)+
  plot_annotation(tag_levels = "A")

ggsave(filename = paste0("output/Graphs/Figure5_ofr.pdf"),
       height = 13.5, width = 21, units = "cm", dpi = 900,
       bg = "white",
       device="pdf")
