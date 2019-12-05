library(readr)
library(tidyverse)
library(data.table)
library(lme4)
library(ggeffects)
library(sjPlot)
library(stargazer)
library(stats)
library(maps)
library(stringr)
library(viridis)
library(ggrepel)
library(MuMIn)

insect_abundance_subset <- fread("insecta_full.csv", sep="\t", header=TRUE)
# insect_adundance_subset <- read_delim('~/Documents/insects_CR_first200000.csv', delim='\t')

#Kept only Inescta from full dataset and removed any samples 
#that were not classified to species level and rows where elevation was NA
species_na_removed <- insect_abundance_subset %>%
  filter(class == "Insecta", !is.na(species), !is.na(elevation), stateProvince!="Indeterminado")

#Created a dataframe that includes the species richness 
C <- species_na_removed %>%
  group_by(locality, stateProvince, elevation, year) %>%
  summarize(number_species = length(unique(species)))

### HOW TO MAKE THE COSTA RICA MAP ####
cr <- map_data("world", region = "Costa Rica")

labs <- tibble(
  lat = c(10.0159, 9.8638,10.6267,9.9981, 9.9913,9.9778, 9.9281),
  long = c(-84.2142, -83.9162,-85.4437,-84.1198, -83.0415, -84.8294,-84.0907),
  names = c("Alajuela", "Cartago","Guanacastle","Heredia","Limón","Puntarenas","San José"))  

gg1 <- ggplot() + 
  geom_polygon(data = cr, aes(x = long, y = lat, group = group), 
               fill = "violet", color = "purple") + 
  coord_quickmap() 

insectLoc <- unique(insecta_only_abd$locality)

#filtered out outlying points outside of costa rica (possible recording error)
localities <- insecta_only_abd %>%
  filter(!is.na(locality), !is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  filter(decimalLongitude >= -87) %>% 
  dplyr::select(locality, decimalLatitude, decimalLongitude, elevation) %>%
  distinct()

(map <- gg1 + 
    geom_point(data =localities , aes(x = decimalLongitude, 
                                      y = decimalLatitude, colour = elevation), 
               shape = 21, fill = "yellow", size = 0.01) +
    geom_label_repel(data = labs, aes(x = long, y = lat, label = names), 
                     size = 2, hjust = 0.01) + 
    xlab("Longitude") + 
    ylab("Latitude") +
    ggtitle("Data Collection Instances in Costa Rica"))

##How species richness is effected by elevation and time 

#Produces Figure 1
ggplot(species_richness_graph_data)+
  geom_point(aes(x = elevation, y = number_species, colour = year), alpha = 0.5, position = "jitter") +
  facet_wrap(~stateProvince, scales='free') + 
  labs(x='Elevation',y='Species Richness') +  
  scale_colour_viridis() +
  theme_linedraw() + 
  theme(legend.background = element_rect(colour='white', fill='white', linetype='solid'), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        legend.title=element_text(colour='black', size=13, face='bold.italic'), 
        legend.text=element_text(colour='black'), 
        strip.text.x = element_text(size = 9, color = "black", face = "bold.italic"),
        axis.title=element_text(size=13, face='bold.italic'),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank())

#Below are models with species richness as a predictor variable and response variables of time and year, 
#but with differing random effects  

#model (1) with random effect of only Province 
species_richness_model_one <- lmer(log(number_species)~elevation + year +(1|stateProvince), 
                                   data=species_richness_graph_data)
summary(species_richness_model_one)
anova(species_richness_model_one)

pred_model_one <- ggpredict(species_richness_model_one, terms=c('elevation'))

##Check assumptions of model one ##

# Check assumption for homogeneity of variance
par(mfrow=c(2,2))
(plot.model <- plot(species_richness_model_one))

#The residuals of the model are normally distributed.
require("lattice")
qqmath(species_richness_model_one, id=0.05)

#Fitted linear model (1) to data from Figure 1
ggplot(pred_model_one) +
  geom_point(data= species_richness_graph_data, 
             aes(x = elevation, y = number_species, colour = year), alpha = 0.6) +
  geom_line(aes(x=x, y=predicted)) +
  facet_wrap(~stateProvince, scales='free') +
  theme_linedraw() +
  labs(x='Elevation',y='Species Richness') + 
  scale_colour_viridis()+
  theme_linedraw() + 
  theme(legend.background = element_rect(colour='white', fill='white', linetype='solid'), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        legend.title=element_text(colour='black', size=13, face='bold.italic'), 
        legend.text=element_text(colour='black'), 
        strip.text.x = element_text(size = 9, color = "black", face = "bold.italic"),
        axis.title=element_text(size=13, face='bold.italic'),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank())

#Plot of Linear model (1) 
ggpredict(species_richness_model_one, terms=c('elevation','stateProvince'), type='re') %>% 
  plot() 

#Plot of difference of general slope value found in summary model and the estimate for random effect of Province 
re_model_one<- plot_model(species_richness_model_one, type='re', show.values=TRUE)
re_model_one                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

#model (2) with random effect of province and locality 
species_richness_model_two <- lmer(number_species~elevation + year +(1|stateProvince) +(1|locality), data=species_richness_graph_data)
summary(species_richness_model_two)      
anova(species_richness_model_two)

### Check assumptions of model two ###
# Check assumption for homogeneity of variance
par(mfrow=c(2,2))
(plot.model <- plot(species_richness_model_one))

#The residuals of the model are normally distributed.
require("lattice")
qqmath(species_richness_model_two, id=0.05)

pred_model_two <- ggpredict(species_richness_model_two, terms=c('elevation'))

#Fitted linear model (2) to data from Figure 1
ggplot(pred_model_two) +
  geom_point(data= species_richness_graph_data, aes(x = elevation, y = number_species, colour = year)) +
  geom_line(aes(x=x, y=predicted)) + facet_wrap(~stateProvince, scales='free') +
  labs(x='Elevation',y='Species Richness') +
  theme_linedraw() +
  theme(legend.background = element_rect(colour='white', fill='white', linetype='solid'), strip.background = element_rect(colour = "white", fill = "white"), legend.title=element_text(colour='black', size=13, face='bold.italic'), legend.text=element_text(colour='black'), strip.text.x = element_text(size = 9, color = "black", face = "bold.italic"), axis.title=element_text(size=13, face='bold.italic'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())

#Plot of linear model (2)        
ggpredict(species_richness_model_two, terms=c('elevation','stateProvince'), type='re') %>% 
plot()        
   
#Plot of difference of general slope value found in summary model and the estimate for random effect of Province
re_model_two<- plot_model(species_richness_model_two, type='re', show.values=TRUE)
re_model_two

#Must do this to add linear models to regression tables produced by stargazer package
class(species_richness_model_two) <- 'lmerMod'
class(species_richness_model_one) <- 'lmerMod'

#Regression table of model (1) and (2) produced by stargazer package
stargazer(species_richness_model_one, species_richness_model_two, type='text', style='qje', digits=3, star.cutoffs=c(0.05, 0.01, 0.001), dep.var.labels = 'species richness')

#How to calculate R^2 marginal and conditional 
r.squaredGLMM(species_richness_model_one)
r.squaredGLMM(species_richness_model_two)

#How max species richness is effected by elevation and time 

max_species_rich <- species_na_removed %>%
  group_by(stateProvince, elevation, year) %>%
  tally()

max_rich <- max_species_rich %>%
  filter(!is.na(elevation)) %>%
  group_by(year, stateProvince) %>%
  summarise(n = max(n))

d <- semi_join(max_species_rich, max_rich, by = c("stateProvince", "year", "n"))

se <- function(x) {sqrt(var(x)/length(x))}

avg_elevation_df <- d %>%
  filter(!is.na(elevation)) %>%
  group_by(stateProvince) %>%
  summarise(avg_elevation = mean(elevation), se = se(elevation))

#Produces plot of maximum species richness against year and coloured by elevation 
ggplot(d) + 
  geom_point(aes(x = year, y =n, colour =elevation)) +
  facet_wrap(~stateProvince, scales='free') + 
  geom_smooth(aes(x = year, y = n, colour =elevation), colour='black')  labs(x='Year',y='Maximum Species Richness') +
  theme_linedraw() +
  theme(legend.background = element_rect(colour='white', fill='white', linetype='solid'),
        strip.background = element_rect(colour = "white", fill = "white"),
        legend.title=element_text(colour='black', size=13, face='bold.italic'),
        legend.text=element_text(colour='black'), 
        strip.text.x = element_text(size = 9, color = "black", face = "bold.italic"), 
        axis.title=element_text(size=13, face='bold.italic'),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 

#Linear model (1) - predictor variable is maximum species richness, response variables are elevation and year, random effect of Province
max_species_richness_model <- lmer(n~elevation+year +(1|stateProvince) , data=d)
summary(max_species_richness_model)
anova(max_species_richness_model)

### Check assumption for max_species_richness_model ###
# Check assumption for homogeneity of variance
par(mfrow=c(2,2))
(plot.model <- plot(max_species_richness_model))

#The residuals of the model are normally distributed.
require("lattice")
qqmath(max_species_richness_model, id=0.05)

pred_model_max <- ggpredict(max_species_richness_model, terms=c('year'))

##Fitted linear model (1) to data of maximum species irchness against year and coloured by elevation
ggplot(pred_model_max) +
  geom_point(data= d,aes(x = year, y = n, colour = elevation))+ geom_line(aes(x=x, y=predicted)) +
  facet_wrap(~stateProvince, scales='free') +
  labs(x = "Year", y = "Maximum Species Richness") +
  theme_linedraw() +
  theme(legend.background = element_rect(colour='white', fill='white', linetype='solid'), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        legend.title=element_text(colour='black', size=13, face='bold.italic'), 
        legend.text=element_text(colour='black'), 
        strip.text.x = element_text(size = 9, color = "black", face = "bold.italic"), 
        axis.title=element_text(size=13, face='bold.italic'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) 

#Must do this to add linear models to regression tables produced by stargazer package
class(max_species_richness_model_two) <- 'lmerMod'

#Regression table of model (1) produced by stargazer package
stargazer(max_species_richness_model_two, type='text',digits=3, star.cutoffs=c(0.05, 0.01, 0.001), dep.var.labels = 'maximum species richness', no.space = TRUE)

#Produce plot of the average elevation where maximum species richness occurs for each Province 
ggplot(avg_elevation_df) +
  geom_bar(aes(x = stateProvince, y = avg_elevation), stat = "identity", fill = "deepskyblue4") +
  geom_errorbar(aes(x =stateProvince, y = avg_elevation, ymin = avg_elevation-se, 
                    ymax = avg_elevation+se), width = 0.4, size = 1.3, alpha = 0.9)+
  theme_linedraw()+ 
  theme(legend.background = element_rect(colour='white', fill='white', linetype='solid'), 
        strip.background = element_rect(colour = "white", fill = "white"), 
        legend.title=element_text(colour='black', size=10, face='bold.italic'), 
        legend.text=element_text(colour='black'), 
        strip.text.x = element_text(size = 8, color = "black", face = "bold.italic"), 
        axis.title=element_text(size=10, face='bold.italic'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) + labs(x='Province', y='Average Elevation of Maximum Species Richness')            

### VENN DIAGRAM (seen in presentation only) ####
#Below produces venn diagram of the orders shared/unique to each Province
#created a dataframe with two columns: stateProvince and order
orderProv <- insecta_only_abd %>%
  select(stateProvince, order) %>%
  filter(!is.na(order), !is.na(stateProvince))

#created prov, which is a vector containing all the province names (excluding "indetermined)
prov <- unique(orderProv$stateProvince)[1:7]

#created a dataframe with 7 columns (one for each province)
#the number of rows is based on Alajuela, as Alajuela has the most rows in orderProv 
df <- data.frame(matrix(ncol = 7, nrow = length(Alajuela)))
colnames(df) <- prov

#filled data frame with all the orders observed in each province
Alajuela <-  unlist(unname(filter(orderProv, stateProvince == prov[1])[,2]))
df$Alajuela[1:length(Alajuela)] <- Alajuela
Puntarenas <-  unlist(unname(filter(orderProv, stateProvince == prov[2])[,2])) 
df$Puntarenas[1:length(Puntarenas)] <- Puntarenas
Limón <-  unlist(unname(filter(orderProv, stateProvince == prov[3])[,2]))
df$Limón[1:length(Limón)] <- Limón
Guanacaste <-  unlist(unname(filter(orderProv, stateProvince == prov[4])[,2])) 
df$Guanacaste[1:length(Guanacaste)] <- Guanacaste
San.José <-  unlist(unname(filter(orderProv, stateProvince == prov[5])[,2])) 
df$"San José"[1:length(San.José)] <- San.José #unable to use for loop bc/ of mismatch b/w San Jose and San.Jose
Cartago <-  unlist(unname(filter(orderProv, stateProvince == prov[6])[,2])) 
df$Cartago[1:length(Cartago)] <- Cartago
Heredia <-  unlist(unname(filter(orderProv, stateProvince == prov[7])[,2]))
df$Heredia[1:length(Heredia)] <- Heredia

#replaced NAs in the dataframe with ""
df[is.na(df)] <- ""
#assigned myV to be the venn diagram
myV <- plotVenn(as.list(df), nCycles = 140000, systemShow = TRUE)


###### FITTED QUADRATIC MODEL ######
#created linear model to compare to quadratic model
linear.model <- lm(max_rich$n~max_rich$year)
summary(linear.model)

#created quadratic model
quadratic.model <-lm(n ~ year + I(year^2), data = max_rich)
summary(quadratic.model)

#quadratic model had the lower AIC value (), and explained for 10.66% more variance
#in the data as adjusted R-squared .1261 >>> 0.01952
AIC(linear.model, quadratic.model)

#fitted quadratic model to the raw data
quad<-lmer(n~year+I(year^2)+(1|stateProvince), data=max_rich)
summary(quad)

ggplot(max_rich, aes(x = year, y = n)) +  
  geom_point(colour = "blue", alpha = 0.5, position = "jitter") +
  stat_smooth(method = lm, formula = y ~ x + I(x^2), colour = "red", alpha = 0.7) + 
  ylab("Maximum Species Richness") +
  xlab("Year")

# Find number of kingdom, phylum, class, order, family, genus, and species in dataset
length(unique(species_na_removed$kingdom))
length(unique(species_na_removed$phylum))
length(unique(species_na_removed$class))
length(unique(species_na_removed$order))
length(unique(species_na_removed$family))
length(unique(species_na_removed$genus))
length(unique(species_na_removed$species))

