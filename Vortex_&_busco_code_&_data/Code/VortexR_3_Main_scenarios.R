library(vortexR)
library(ggplot2)

setwd("./Vortex_output_manuscript") # set working directory
Pink_pigeon <- collate_dat('Final_vortex_simulations', 1000, save2disk = FALSE) # takes all data files produced by the project Final_vortex_simulations and puts in to alist object called Pink_pigeon

###Visualising the data - qualitative analysis ###
#All plots produced below were done so at four time points 5, 25, 50 100 examine differences over time###
Meta <- subset(Pink_pigeon, Pink_pigeon$pop.name=="Metapopulation")

PW <-pairwise(Pink_pigeon, "Final_vortex_simulations", scenario="Scenario_1",params = c("Nall", "Het"), yrs = c(100), save2disk = FALSE)
#pairwise creates a list of elements so you can access each set of results via an index
#produced SSMD which are better than selection coeff apparently as not influenced by sample size
pval <- PW[[3]]

pval$SSMD_Nall100 <- round(pval$SSMD_Nall100, 4)
pval$SSMD_Het100 <- round(pval$SSMD_Het100, 4)

pval_meta <- subset(pval, pval$Population=="Metapopulation")
pval_meta

###
#Scenario     Population SSMD_Nall100 SSMD_Het100
#13 Scenario_2 Metapopulation       0.3096      0.1595
#14 Scenario_3 Metapopulation       0.0218      0.0000
###

SSMD <- PW[[2]]
SSMD$SSMD_Nall100 <- round(SSMD$SSMD_Nall100, 4)
SSMD$SSMD_Het100 <- round(SSMD$SSMD_Het100, 4)
SSMD_meta <- subset(SSMD, SSMD$Population=="Metapopulation")
SSMD_meta

###
#Scenario     Population SSMD_Nall100 SSMD_Het100
#13 Scenario_2 Metapopulation       0.4970      0.9966
#14 Scenario_3 Metapopulation       2.0185     10.8500
###

###Plot abundance###

####Making nice plots with ggplot###
line_plot_year(data=Meta, project="Final_vortex_simulations", scenario="Scenario_1", plotpops="Metapopulation",params="Nall", save2disk=T) # dot plot showing the change in abundance at 4 time points and so the difference between the three scenarios.
#produces a basic plot but saves data which can be used to build better plot below
#open Final_vortex_simulations_Scenario_1_Nall_plot.rda manually will bei nthe newly created 'Plots' directory

data <- Final_vortex_simulations_Scenario_1_Nall_plot[[1]]
#colour blind friendly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7")
scale_colour_manual(values = cbbPalette)
  
#greyscale as used in manuscript
greyscle <- c("#000000","#999999","#dedede")

ggplot(data, aes(Year, Nall)) +
  geom_line(aes(colour = scen.name), size = 0.5) +
  geom_errorbar(aes(ymin=Nall-SE.Nall., ymax=Nall+SE.Nall., colour = scen.name), width = 1) +
  scale_colour_manual(values = greyscle) +
  scale_y_continuous(breaks=seq(0,600,100), limits = c(0, NA), expand = c(0,0)) +
  scale_x_continuous(limits = c(0, NA), expand = c(0,0)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title = element_text(size = 14),
        axis.line = element_line(size = 0.5),
        legend.position = "none") +
  ylab("Abundance") +
  xlab("Number of years from start of simulation")
