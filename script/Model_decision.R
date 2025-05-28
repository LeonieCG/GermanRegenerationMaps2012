# Model decision analysis
packages <- c("tidyr", "dplyr","magrittr", "ggplot2", "tidyverse")
sapply(packages, FUN = library, character.only = T)


# Model assumptions
# DHARMA, there is A problem when:
# disp.p < 0.05: significant ratio > 1 indicates overdispersion, a significant ratio < 1 underdispersion (?testDispersion)
# zeroinfl.p < 0.05: a value < 1 means that the observed data has less zeros than expected,
#                    a value > 1 means that it has more zeros than expected (aka zero-inflation) (?testZeroInflation)
# spatautocorr.p < 0.05
#-> check them also visually 


# Model performance check ----------------------------------------------------------------
ms_reg <- read.csv2("output/Fits/Sapling/h50d7_Germany/Sapling_model_summary.csv") 

ms_reg %<>% mutate(mae.rule = ifelse(mae.relative.median <= 2, "passed", "failed"),
                   rsq.rule = ifelse(rsq.test.median >= 0.1, "passed", "failed"))

ms_reg_final <- ms_reg %>%
  filter(mae.rule == "passed" & rsq.rule == "passed")

ms_reg_out <- ms_reg %>%
  filter(mae.rule == "failed" | rsq.rule == "failed" | is.na(cv.method))


# DHARMa checks p-value and visual ------------------------------------------
# ms_reg_final %>%
#   filter(disp.p <= 0.05 | zeroinfl.p <= 0.05 | spatautocorr.p <= 0.05)

# Visual interpretation: when red line leaves distribution = not so good, but it depends on mangitude of deviation from the mean!

# Dispersion:
ms_reg_final %>% 
  filter(disp.p <= 0.05) %>% 
  select(species)

# Tilia within bars
# Quercus.rubra within bars
# Pinus.sylvestris 0.5
# Picea.abies 0.5
# Abies.alba 0.5
# Acer.campestre within bars
# Carpinus.betulus 0.5
# Fraxinus.excelsior 0.75
# Alnus.glutinosa within bars
# Fagus.sylvatica 0.4
# Prunus.serotina within bars
# Prunus.avium within bars
# -> no major deviation!

# Zero inflation:
ms_reg_final %>% 
  filter(zeroinfl.p <= 0.05)
# -> none!

# Spatial autocorrelation
ms_reg_final %>%
  filter(spatautocorr.p <= 0.05)

ms_reg %>% # before model selection
  filter(spatautocorr.p <= 0.05)


# Save ----------------------------------------------------
saveRDS(ms_reg_final,"output/Fits/Sapling/h50d7_Germany/Sapling_model_final.rds")


# Plot --------------------------------------------------------------------
 ggplot()+
  theme_bw()+
  scale_x_continuous(expand=c(0, 0), limits=c(-1,1)) + #limits=c(-1.1, 1)) +
  scale_y_continuous(expand=c(0, 0), limits=c(-1,1)) + #limits=c(-1.1, 1)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.1) +
  xlab(expression("pseudo-R"["test"]^2)) +
  ylab(expression("pseudo-R"["train"]^2)) +
  geom_polygon(data = data.frame(x =c(-Inf, -Inf, Inf), y= c(-Inf, Inf, Inf)),
               aes(x = x, y = y), fill = "#B49629", alpha = 0.3) + 
  geom_point(data= ms_reg, aes(x = rsq.test.mean, y = rsq.train.mean, colour = mae_rule )) +
  scale_colour_manual(expression(Delta*"mae"), 
                      values = c("passed" = "#D36D39", "failed" = "black"))

  ggplot()+
    theme_bw()+
    scale_x_continuous(expand=c(0, 0), limits=c(-1,1)) + #limits=c(-1.1, 1)) +
    scale_y_continuous(expand=c(0, 0), limits=c(-1,1)) + #limits=c(-1.1, 1)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = 0.1) +
    xlab(expression("pseudo-R"["test"]^2)) +
    ylab(expression("pseudo-R"["train"]^2)) +
    geom_polygon(data = data.frame(x =c(-Inf, -Inf, Inf), y= c(-Inf, Inf, Inf)),
                 aes(x = x, y = y), fill = "#B49629", alpha = 0.3) + 
    geom_point(data= ms_reg, aes(x = rsq.test.median, y = rsq.train.median, colour = mae_rule )) +
    scale_colour_manual(expression(Delta*"mae"), 
                        values = c("passed" = "#D36D39", "failed" = "black"))


