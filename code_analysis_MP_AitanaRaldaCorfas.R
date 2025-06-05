### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## Project: Master Thesis: Environmental challenges and biodiversity loss. 
#                         What makes some communities thrive while others struggle?
## Content: Analysis and Figures 
## Author: Aitana Ralda Corfas
## Date : 02/06/2025
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#Install libraries and set the working directory ####
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(GGally)
library(viridis)
library(dplyr)
library(tidyr)
library(stringr)
library(vegan)
library("beeswarm")
library(geomtextpath)
library("systemfonts")
library(gridExtra)
library(AICcmodavg)
library(car)

setwd("~/")

#Part 1: Coexistence analysis:----
data<-read.csv("rawdata_1_coexistence.csv", sep=";") #import the data
    #calculate Shannon disveristy and total biomass
shannon_diversity <- function(abundances) {
  if (is.character(abundances)) {
    abundances <- as.numeric(strsplit(abundances, ",")[[1]])  # Convert to numeric
  }
  
  abundances <- abundances[!is.na(abundances)]  # Remove NAs
  
  if (length(abundances) == 0 || sum(abundances) == 0) {
    return(NA)  # Avoid log(0) errors
  }
  
  return(vegan::diversity(abundances, index = "shannon"))
}

biomass <- function(abundances) {
  if (is.character(abundances)) {
    abundances <- as.numeric(strsplit(abundances, ",")[[1]])  # Convert to numeric
  }
  
  abundances <- abundances[!is.na(abundances)]  # Remove NAs
  
  if (length(abundances) == 0 || sum(abundances) == 0) {
    return(0)  # Avoid log(0) errors
  }
  
  return(sum(abundances))
}

data <- data %>%
  mutate(
    biomass = sapply(final_abundances, biomass),
    ShannonI = sapply(final_abundances, shannon_diversity))

  ##Plot Figure 3: Richness of communities with different community structure----
    ###Figure 3a- Saturation curves----
##copy the datset to implement the modifciation in the labels for the figure
dataFigure<-data
  # Replace values in specialization column
dataFigure$specialization[dataFigure$specialization==0.5]<-"Random"
dataFigure$specialization[dataFigure$specialization == "0.15"] <- "Generalist"
dataFigure$specialization[dataFigure$specialization == "0.85"] <- "Specialist"
 # Replace values in nestedness column
dataFigure$nestedness[dataFigure$nestedness==-1]<-"No CF"

##create array for inidcating the number of resources given for each combination
a_mean <- dataFigure %>% 
  group_by(Resources.given, specialization, nestedness) %>% 
  summarize(Rg = mean(as.numeric(as.character(Resources.given))))

  # Double-check for updating the names of the specialization column:
a_mean$specialization[a_mean$specialization == "0.15"] <- "Generalist"
a_mean$specialization[a_mean$specialization == "0.85"] <- "Specialist"
  # Add label to a_mean for the dashed line
a_mean$label <- "# Resources given"

## Label mappings for the plot
speci_labels <- c(
  "Generalist" = "Generalist",
  "Random" = "Random",
  "Specialist" = "Specialist"
)

res_labels <- c(
  `8` = "Resources given: 8",
  `16` = "Resources given: 16",
  `32` = "Resources given: 32",
  `64` = "Resources given: 64"
)

# Line for total resources
total_res_line <- data.frame(y = 80, label = "# Total resources")


ggplot(dataFigure, aes(S, nb_coexisting_species, color = as.factor(nestedness))) +
  geom_point(size = 0.2, alpha = 0.2) +
  #geom_smooth() + 
  stat_summary(aes(y = nb_coexisting_species ), fun=mean, linewidth = 1, geom="line")+
  ylim(0, 85) +
  xlab("Sampled species") +
  ylab("Surviving species") +
  
  # Two hlines with legend
  geom_hline(data = total_res_line, aes(yintercept = y, linetype = label), color = "black") +
  geom_hline(data = a_mean, aes(yintercept = Rg, linetype = label), color = "black") +
  
  # Correct facet order and label mapping
  facet_grid( factor(specialization, c("Generalist", "Random", "Specialist")) ~ Resources.given,
              labeller = labeller(Resources.given = res_labels,
                                  specialization = speci_labels)) +
  
  #ggtitle("Coexisting species with different nestedness") +
  labs(color = "Nestedness", linetype = NULL) +
  theme_bw()+
  scale_linetype_manual(values = c(
    "# Total resources" = "solid",
    "# Resources given" = "dashed"
  ))

    ###Figure 3b- Boxplot at S=200----

ggplot(dataFigure[dataFigure$S==200,], aes(nestedness,nb_coexisting_species, fill=nestedness ))+
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  ylab("Surviving species at S= 200")+
  xlab("Nestedness")+
  theme_bw()+
  #ggtitle("Surviving species by community structure") +
  labs(colour = "Nestedness")+
  facet_grid(specialization~Resources.given, labeller =  labeller(Resources.given = res_labels,
                                                                  specialization = speci_labels))

    ###ANOVA ----
  #copy the datset to implement the modifications needed in the analysis
dataTest<-data
dataTest$crossfeeding <- ifelse(dataTest$nestedness == -1, 0, 1) # create a column with presence,abscence crossfeeding 
dataTest$nestedness[dataTest$crossfeeding == 0] <- NA #set as NA for nestedness in non-crossfeeding communities  

      ####a) ANOVA Crossfeeding (all the communities)----
model_C<-aov(nb_coexisting_species ~ specialization + crossfeeding+Resources.given, data=dataTest[dataTest$S==200,])
(summary(model_C))
(coef(model_C))
 
par(mfrow=c(2,2))
plot(model_C) #plot the residuals

      ####b) ANOVA Nestedness (just crossfeeding communities)----
model_N<-aov(nb_coexisting_species ~ specialization + nestedness + Resources.given , data=dataTest[dataTest$S==200,])
(summary(model_N))
coef(model_N)

plot(model_N)#plot the residuals
par(mfrow=c(1,1))







  ##Plot Figure 4: Changes of specialization score from the initial community to the survival one----
    ###Figure 4a- Difference initial- final specialization score----

final05<-read.csv("rawdata_2_05.csv", sep= ";") #load the data of communities with specilaization score of 0.5
dataw05<-rbind(dataFigure[,1:16], final05) #join with the other dataset
dataw05$nestedness[dataw05$nestedness==-1]<-"No CF" #update the No CF label
dataw05$specialization[dataw05$specialization == "0.15"] <- "Generalist" #update specilization label
dataw05$specialization[dataw05$specialization == "0.85"] <- "Specialist" #update specilization label
dataw05$specialization[dataw05$specialization == "0.5"] <- "Mid (0.5)"#update specilization label


#plot the data
    dataw05$A<- dataw05$Speafter-dataw05$Spebefore #Calculate the difference Initial-Final specilization score
    
    ggplot(dataw05, aes((specialization),A, fill=specialization))+
      geom_boxplot()+
      geom_hline(yintercept =0)+
      xlab("")+
      ylab("Difference intial-final specialization score")+
      theme_bw()
      ####T-test Final-initial specialization score ---- 
    t.test(dataw05$Spebefore[dataw05$specialization== "Generalist"], dataw05$Speafter[dataw05$specialization== "Generalist"], alternative= "less")
    #t = -48.508, df = 11921, p-value < 2.2e-16
    t.test(dataw05$Spebefore[dataw05$specialization== "Specialist "], dataw05$Speafter[dataw05$specialization== "Specialist "], alternative= "less")
    # t = 3.7296, df = 13966, p-value = 0.9999
    t.test(dataw05$Spebefore[dataw05$specialization== "Mid (0.5)"], dataw05$Speafter[dataw05$specialization== "Mid (0.5)"], alternative= "less")
    #t = -5.2053, df = 3142.1, p-value = 1.031e-07
    t.test(dataw05$Spebefore[dataw05$specialization== "Random"], dataw05$Speafter[dataw05$specialization== "Random"], alternative= "less")
    #t = -42.071, df = 12479, p-value < 2.2e-16
    
    ###Figure 4b- Initial vs final specialization score global comparison ----
    
    #modify the data to obtain mean +sd 
    BA05 <- dataw05 %>%
      group_by( specialization) %>%
      summarise(
        mean_before = mean(Spebefore, na.rm = TRUE),
        sd_before = sd(Spebefore, na.rm = TRUE),
        mean_after = mean(Speafter, na.rm = TRUE),
        sd_after = sd(Speafter, na.rm = TRUE),
        n = n()
      ) %>%
      ungroup()
    
    BA05_long <- BA05 %>% #transform in long format the mean
      select(specialization, mean_before, mean_after) %>%
      pivot_longer(cols = starts_with("mean"), 
                   names_to = "time", 
                   values_to = "mean") %>%
      mutate(time = ifelse(time == "mean_before", 1, 2))
    
    BA05_sd <- BA05 %>% #transform in long format the errors
      select(specialization, sd_before, sd_after) %>%
      pivot_longer(cols = starts_with("sd"), 
                   names_to = "time", 
                   values_to = "sd") %>%
      mutate(time = ifelse(time == "sd_before", 1, 2))
    
    # Join the two long datasets
    BA05_plot <- left_join(BA05_long, BA05_sd, by = c("specialization", "time"))
    #update the labels
    BA05_plot$specialization[BA05_plot$specialization == "0.15"] <- "Generalist (0.15)"
    BA05_plot$specialization[BA05_plot$specialization == "0.85"] <- "Specialist (0.75)"
    BA05_plot$specialization[BA05_plot$specialization == "0.5"] <- "Mid (0.5)"
    
    # Plot
    ggplot(BA05_plot, aes(x = time, y = mean, group = specialization, color = specialization)) +
      geom_line() +
      ylim(c(0,1))+
      geom_point(size = 3) +
      labs(fill="NA")+
      #scale_fill_manual('Specialization')+
      geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1) +
      scale_x_continuous(breaks = c(1, 2), labels = c("Initial", "Final")) +
      ylab("Specialization score") +
      xlab("") +
      theme_bw()
      
     




  ##Plot Figure 5: Changes in richness for resource production regimes due to resource identity----
    
    K200<- read.csv("rawdata_3_firstRG", sep=";") #load the data 
    K200$specialization[K200$specialization==0.5]<-"Random" #modfiy specialization label
    K200$nestedness[K200$nestedness==-1]<-"No CF" #modify nestedness label

    plot1<-ggplot(K200[K200$Resources.given==16,], aes(nestedness, nb_coexisting_species, fill=nestedness ))+
      geom_boxplot()+
      geom_jitter(alpha=0.5)+
      #geom_smooth(data=data[data$specialization==0.85 & data$Resources.given==8,], linewidth=0.5, alpha= 0.4, linetype=2) + 
      #geom_point(data=data[data$specialization==0.85 & data$Resources.given==8,],  size=0.2, alpha=0.2)+
      ylim(0, 60)+
      xlab("Nestedness")+
      ylab("")+
      ggtitle("First 16 Resources given")+
      theme_bw()+
      theme(legend.position="none")
    
   
    plot2<-ggplot(dataFigure[dataFigure$Resources.given==16 & dataFigure$S==200 & dataFigure$specialization =="Specialist",], aes(nestedness, nb_coexisting_species, fill=nestedness ))+
      geom_boxplot()+
      geom_jitter(alpha=0.5)+
      #geom_smooth(data=data[data$specialization==0.85 & data$Resources.given==8,], linewidth=0.5, alpha= 0.4, linetype=2) + 
      #geom_point(data=data[data$specialization==0.85 & data$Resources.given==8,],  size=0.2, alpha=0.2)+
      ylim(0, 60)+
      xlab("Nestedness")+
      ylab("Surviving species") +
      ggtitle("Random 16 Resources given")+
      theme_bw()+
      theme(legend.position="none")
    
    grid.arrange(plot2, plot1, ncol=2)
 

      ####T-test Random resources given or not    ----
    t.test(data$nb_coexisting_species[data$Resources.given== 16 & data$specialization == 0.85 & data$nestedness==-1 &data$S == 200], K200$nb_coexisting_species[K200$Resources.given== 16 & K200$nestedness== "No CF" ])
    #NO CF_ t = 0.31526, df = 22.082, p-value = 0.7555
    t.test(data$nb_coexisting_species[data$Resources.given== 16 & data$specialization == 0.85 & data$nestedness==0 &data$S == 200], K200$nb_coexisting_species[K200$Resources.given== 16 & K200$nestedness=="0" ])
    #0 -> t = 1.8168, df = 24.812, p-value = 0.08134
    t.test(data$nb_coexisting_species[data$Resources.given== 16 & data$specialization == 0.85 & data$nestedness==0.5 &data$S == 200], K200$nb_coexisting_species[K200$Resources.given== 16 & K200$nestedness=="0.5" ])
    #0.5 -> t = -2.2603, df = 14.457, p-value = 0.03971
    t.test(data$nb_coexisting_species[data$Resources.given== 16 & data$specialization == 0.85 & data$nestedness==1 &data$S == 200], K200$nb_coexisting_species[K200$Resources.given== 16 & K200$nestedness=="1" ])
    #t = -4.5591, df = 11.085, p-value = 0.000802
    hist(K200$nb_coexisting_species[K200$Resources.given== 16 & K200$nestedness=="No CF" ])
    
    
    
#Part 2: Disturbance analysis:----
    disturbance<-read.csv("rawdata_4_disturbance.csv", sep=";") #import the data
    
    ##add initial resourcs given for the disturbances (=K at coex)
    disturbance <- disturbance %>%
      group_by(seed) %>%
      mutate(initial_resources_given = first(Resources.given[type == "coex"])) %>%
      ungroup() %>%
      # Add resources_removed as the difference
      mutate(resources_removed = as.numeric(initial_resources_given) - as.numeric(Resources.given))
    
        #calculate  biomass and shannon diversity of the communities at intial steady state 
    coex_rows <- which(disturbance$type == "coex") #select coex rows
    
    disturbance <- disturbance %>%
          mutate( biomass = sapply(final_abundances_end, biomass),
            ShannonI = sapply(final_abundances_end, shannon_diversity))
        
        
    disturbance$ShannonI[coex_rows] <- sapply( ##apply diversity function
      disturbance$final_abundances[coex_rows],
      shannon_diversity)
    
    disturbance$biomass[coex_rows] <- sapply( #apply biomass function
      disturbance$final_abundances[coex_rows],
      biomass)
    

 #Fitting R50 comparing with gupta
    funfit <- function(x,c50, P) { Ymax<-22 
    return(Ymax/(1 + (x/c50)^P))}
    
    # Initialize a data frame to store the results
    decay_parameter <- data.frame(
      Seed = character(),
      Disturbance = character(),
      Resources.given = numeric(),
      Specialization = numeric(),
      Nestedness = numeric(),
      Estimate = numeric(),
      RSS= numeric(),
      Min.value= numeric(),
      c50=numeric(),
      Rc50=numeric(),
      stringsAsFactors = FALSE
    )
    
    # Get unique seeds and types
    seeds <- unique(disturbance$seed)
    types <- unique(disturbance$type)
    
    
    # Pre-define your output
    decay_parameter <- data.frame()
    
    # Loop over each seed and disturbance type
    for (seed_t in seeds) {
      for (type_t in types) {
        
        # Subset for current seed and disturbance type
        sub_data <- disturbance[disturbance$seed == seed_t & disturbance$type == type_t, ]
        
        if (nrow(sub_data) > 1) {
          
          
          index_min<- sub_data %>%
            filter(nb_coexisting_species_end <= 1) %>%
            summarize(Value = min(Value))
          
          
          Rc50 <- approx(x = sub_data$nb_coexisting_species_end, 
                         y = sub_data$Value, 
                         xout = 11, 
                         ties = mean)$y
          
          # Adjust model with Gupta decay
          model_gupta <- nls(nb_coexisting_species_end ~ funfit(Value, c50, P),
                             data = sub_data,
                             start = list(P = 1, c50= 0.15))
          
          s <- summary(model_gupta)
          estimate <- s$parameters["P", "Estimate"]
          Cestimate <- s$parameters["c50", "Estimate"]
          RSS <- deviance(model_gupta)
          
          # add results to dataset
          decay_parameter <- rbind(decay_parameter, data.frame(
            Seed = seed_t,
            Disturbance = type_t,
            Resources.given = sub_data$initial_resources_given[1],
            Specialization = sub_data$specialization[1],
            Nestedness = sub_data$nestedness[1],
            Estimate = estimate,
            RSS = RSS,
            Min.value = index_min$Value,
            c50=Cestimate,
            Rc50=Rc50
          ))}
        
        
      }
    }
    
    #summary of the results to include mean +sd
    summary_decay <- decay_parameter %>%
      group_by(Resources.given, Nestedness, Disturbance) %>%
      summarise(
        mean_Rc50 = mean(Rc50, na.rm = TRUE), #Rc50 is the interpolated value
        sd_Rc50 = sd(Rc50, na.rm = TRUE),
        n = n()
      ) %>%
      ungroup()
    
  ##Figure 6- Robustness of communities to the three disturbances----
    #Define custom labels for the facets
    facet_labels <- c(
      Ant = "Antibiotic",
      dil = "Dilution",
      R = "Resource removal"
    )
    
    facet_spe <- c(
      "-1" = "Random",
      "0.85" = "Specialist"
    )
    
    ggplot(decay_parameter[decay_parameter$Disturbance != "S",], 
           aes(x = as.factor(Nestedness), y = Rc50, color = as.factor(Resources.given))) +
      geom_point(alpha = 0.9, size = 0.6, position = position_dodge(width = 0.4)) +
      
      geom_errorbar(data = summary_decay[summary_decay$Disturbance != "S",],
                    aes(x = as.factor(Nestedness), 
                        ymin = mean_Rc50 - sd_Rc50, 
                        ymax = mean_Rc50 + sd_Rc50, 
                        color = as.factor(Resources.given)),
                    inherit.aes = FALSE,
                    width = 0.2,
                    position = position_dodge(width = 0.4)) +
      
      geom_point(data = summary_decay[summary_decay$Disturbance != "S",],
                 aes(x = as.factor(Nestedness), y = mean_Rc50, color = as.factor(Resources.given)),
                 inherit.aes = FALSE,
                 position = position_dodge(width = 0.4)) +
      
      xlab("Presence/Absence of crossfeeding") +
      scale_x_discrete(labels=c("No CF", "0", "0.5", "1"))+
      ylab("Disturbance intensity at 50% species loss (R50)")+
      facet_wrap(Specialization~Disturbance, scales = "free", labeller = labeller(Disturbance = facet_labels, Specialization = facet_spe)) +
      labs(color = "Resources given") +
      theme_bw()
      #ggtitle("R50 for different disturbances and community structures")
    
    ###ANOVA R50 ----
    #Modification of the dataset frot he analysis
    decay_analysis<-decay_parameter
    decay_analysis$Resources.given<-as.numeric(as.character(decay_analysis$Resources.given))
    decay_analysis$Specialization<-as.numeric(as.character(decay_analysis$Specialization))
    decay_analysis$Nestedness<-as.numeric(as.character(decay_analysis$Nestedness))
    
    #resetting crossfedign and nestedness labels
    decay_analysis$Crossfeeding <- ifelse(decay_analysis$Nestedness == -1, 0, 1)
    decay_analysis$Nestedness[decay_analysis$Crossfeeding == 0] <- NA
    
    #ANOVA with Nestedness (only crossfeeding communities)
    vad_n<-aov(Rc50 ~ Resources.given+Specialization+Nestedness, decay_analysis[decay_analysis$Disturbance=="dil",])
    van_n<-aov(Rc50 ~ Resources.given+Specialization+Nestedness, decay_analysis[decay_analysis$Disturbance=="Ant",])
    var_n<-aov(Rc50 ~ Resources.given+Specialization+Nestedness, decay_analysis[decay_analysis$Disturbance=="R",])
    
    par(mfrow=c(2,2)) #plot residuals
    plot(vad_n)
    plot(van_n)
    plot(var_n)
    
    #ANOVA with crossfeeding (all the communities)
    vad_c<-aov(Rc50 ~ Resources.given+Specialization+Crossfeeding, decay_analysis[decay_analysis$Disturbance=="dil",])
    van_c<-aov(Rc50 ~ Resources.given+Specialization+Crossfeeding, decay_analysis[decay_analysis$Disturbance=="Ant",])
    var_c<-aov(Rc50 ~ Resources.given+Specialization+Crossfeeding, decay_analysis[decay_analysis$Disturbance=="R",])

    plot(vad_c) #plot residuals 
    plot(van_c)
    plot(var_c)
    par(mfrow=c(1,1))
    
    
    
 
    
  ##Figure 7- Correlation of robustness among disturbances: Antibiotics, Dilution and Resource removal----
    mean_rc50 <- decay_parameter %>%
      filter(Disturbance %in% c("Ant", "R", "dil")) %>%
      group_by(Seed, Disturbance) %>%
      summarise(mean_Rc50 = mean(Rc50), .groups = "drop") %>%
      pivot_wider(names_from = Disturbance, values_from = mean_Rc50, names_prefix = "meanRc50_")
    
    # Plot: x = mean Rc50 for Ant, y = mean Rc50 for R
    p1<-ggplot(mean_rc50, aes(x = meanRc50_Ant, y = meanRc50_R)) +
      geom_point() +
      labs(x = "R50 (Antibiotics)", y = "R50 (Resource removal)") +
      geom_smooth(method= "lm")+
      theme_minimal()
    
    p2<-ggplot(mean_rc50, aes(x = meanRc50_Ant, y = meanRc50_dil)) +
      geom_point() +
      labs(x = " R50 (Antibiotics)", y = " R50 (Dilution)") +
      geom_smooth(method= "lm")+
      theme_minimal()
    
    p3<-ggplot(mean_rc50, aes(x = meanRc50_dil, y = meanRc50_R)) +
      geom_point() +
      labs(x = " R50 (Dilution)", y = " R50 (Resource removal)") +
      geom_smooth(method= "lm")+
      theme_minimal()
    
    grid.arrange(p1,p2,p3, ncol=3)
    
    ##Pearson correlation test
    cor.test(mean_rc50$meanRc50_Ant, mean_rc50$meanRc50_R, method="pearson")
    cor.test(mean_rc50$meanRc50_Ant, mean_rc50$meanRc50_dil, method="pearson")
    cor.test(mean_rc50$meanRc50_dil, mean_rc50$meanRc50_R, method="pearson")

    
       
    
    

#Part 3: Comparison with realistic matrix----
    #load the data
    dalbello<-read.csv("rawdata_5_dbcomparison.csv", sep=",")
    
    dataFigure$nestedness[dataFigure$nestedness==-1]<-"No CF"
    
    
    
    #obtain long dataset with mean and sd
    dalbello_long<- dalbello %>% group_by(Sparcity, Nestedness_pb )  %>%
      summarize(mean_sp= mean(coexisitng.species, na.rm=TRUE),
                sd_sp= sd(coexisitng.species, na.rm=TRUE) )
    
    ##Figure 8- Comparison realistic matrix----
    db_label<-c(
      "0"= "Sparcity: 0.5",
      "0.95"="Sparcity: 0.97" )
    
      ggplot(dalbello, aes(x = as.factor(Nestedness_pb), y = coexisitng.species, color = as.factor(Nestedness_pb))) +
      geom_point(position = position_jitter(width = 0.1, height = 0)) +
      ylim(0,52)+
      geom_point(data = dalbello_long,
                 aes(x = as.factor(Nestedness_pb), y = mean_sp),
                 inherit.aes = FALSE,
                 color = "black", size = 3) +
      
      geom_errorbar(data = dalbello_long,
                    aes(x = as.factor(Nestedness_pb),
                        ymin = mean_sp - sd_sp,
                        ymax = mean_sp + sd_sp),
                    inherit.aes = FALSE,
                    width = 0.1) +
      labs(color = "Nestedness") +
      xlab("Nestedness")+
        ylab("Richness")+
      scale_x_discrete(labels=c( "No CF", "0", "0.5", "1"))+
      theme_bw()+
      theme(legend.position="none")+
     #ggtitle("Richness of communities with different sparcity and nestedness ")+
      facet_grid(~Sparcity, labeller=labeller(Sparcity=db_label))
      
    
    ###ANOVA comparison Dal Bello ----
    dalbello_analysis<-dalbello
    dalbello_analysis$Crossfeeding <- ifelse(dalbello_analysis$Nestedness == -1, 0, 1)
    dalbello_analysis$Nestedness[dalbello_analysis$Crossfeeding == 0] <- NA
    
    #ANOVA Nestedness
    DB_N<-aov(coexisitng.species~Nestedness_pb, dalbello_analysis[dalbello_analysis$Sparcity==0.95,])
    
    #ANOVA Crossfeeding
    DB_C<-aov(coexisitng.species~Crossfeeding, dalbello_analysis[dalbello_analysis$Sparcity==0.95,])


#Suppelementary Figures----

##Figure S3: comparison mid (0.5 specialization score) and random communities
    a_mean <- dataw05 %>% 
      group_by(Resources.given, specialization, nestedness) %>% 
      summarize(Rg = mean(as.numeric(as.character(Resources.given))))
    
    ggplot(dataw05[dataw05$specialization %in% c("Mid (0.5)", "Random"),], aes(S, nb_coexisting_species, color=specialization ))+
      
    geom_point(size = 0.2, alpha = 0.2) +
      #geom_smooth() + 
      stat_summary(aes(y = nb_coexisting_species ), fun=mean, linewidth = 1, geom="line")+
      ylim(0, 85) +
      xlab("Sampled species") +
      ylab("Surviving species") +
    
      # Correct facet order and label mapping
      facet_grid(  nestedness~Resources.given, labeller = label_both)+
      
      #ggtitle("Coexisting species with different nestedness") +
      labs(color = "Specialization", linetype = NULL) +
      theme_bw()
    
##Figure S7: Saturation curve with diversity
ggplot(dataFigure, aes(S, ShannonI, color=as.factor(nestedness)))+
  geom_point(size=0.2, alpha=0.2)+
  geom_smooth() + 
  xlab("Sampled species")+
  ylab("Diversity (Shannon I)")+
  #ggtitle("Diversity with different nestedness") +
  labs(colour = "Nestedness")+
  theme_bw()+
  facet_grid( specialization~Resources.given, labeller = label_both)

##Figure S8: Ratio intial/final specialization score
    dataw05$ratio<- dataw05$Speafter/dataw05$Spebefore
    ggplot(dataw05, aes((specialization),ratio, fill=specialization))+
      geom_boxplot()+
      geom_hline(yintercept =1)+
      xlab("")+
      ylab("Ratio final/initial specialization score")+
      theme_bw

 ##Figure S9- Initial Diversity ----
    speci_labels <- c(
      `-1` = "Random",
     `0.85` = "Specialist"
    )
    
    res_labels <- c(
      `8` = "Resources given: 8",
      `16` = "Resources given: 16",
      `32` = "Resources given: 32",
      `64` = "Resources given: 64"
    )
    
    ggplot(disturbance[disturbance$type =="coex",], aes(as.factor(nestedness) ,ShannonI, fill=as.factor(nestedness) ))+
      geom_boxplot()+
      geom_jitter(color="black", size=0.4, alpha=0.7)+
      ylab("Diversity (Shannon Index)")+
      xlab("Nestedness")+
      #ggtitle("Diversity of communities used for the disturbance analysis") +
      labs(colour = "Nestedness")+
      scale_x_discrete(labels=c( "0", "0.5", "1", "No CF"))+
      theme_bw()+
      theme(legend.position="none")+
      facet_grid(specialization~Resources.given, labeller =  labeller(Resources.given = res_labels,
                                                                      specialization = speci_labels))
    ##Wilcox test on diversity
    ###specialization
    wilcox.test(disturbance$ShannonI[disturbance$type =="coex" & disturbance$specialization == "0.85"], disturbance$ShannonI[disturbance$type =="coex" &disturbance$specialization == "-1"])
    ##crossfeeding
    wilcox.test(disturbance$ShannonI[disturbance$type =="coex" & disturbance$nestedness == "-1"], disturbance$ShannonI[disturbance$type =="coex" &disturbance$nestedness != "-1"])
    #resources given
    wilcox.test(disturbance$ShannonI[disturbance$type =="coex" & disturbance$initial_resources_given == 64], disturbance$ShannonI[disturbance$type =="coex" &disturbance$initial_resources_given == 32])
    #nestedness
    pairwise.wilcox.test(disturbance$ShannonI[disturbance$type =="coex"], disturbance$nestedness[disturbance$type =="coex"])
    
  ##Figure S10- Initial Biomass---- 
    ggplot(disturbance[disturbance$type =="coex",], aes(as.factor(nestedness) ,biomass, fill=as.factor(nestedness) ))+
      geom_boxplot()+
      geom_jitter(color="black", size=0.4, alpha=0.7)+
      ylab("Total abundance")+
      xlab("Nestedness")+
      #ggtitle("Biomass of communities used for the disturbance analysis") +
      labs(colour = "Nestedness")+
      scale_x_discrete(labels=c( "0", "0.5", "1", "No CF"))+
      theme_bw()+
      theme(legend.position="none")+
      facet_grid(specialization~Resources.given, labeller =  labeller(Resources.given = res_labels,
                                                                      specialization = speci_labels))
    ##Wilcox test on initial biomass
    ###specialization
    wilcox.test(disturbance$biomass[disturbance$type =="coex" & disturbance$specialization == "0.85"], disturbance$biomass[disturbance$type =="coex" &disturbance$specialization == "-1"])
    ##crossfeeding
    wilcox.test(disturbance$biomass[disturbance$type =="coex" & disturbance$nestedness == "-1"], disturbance$biomass[disturbance$type =="coex" &disturbance$nestedness != "-1"])
    #resources given
    wilcox.test(disturbance$biomass[disturbance$type =="coex" & disturbance$initial_resources_given == 64], disturbance$biomass[disturbance$type =="coex" &disturbance$initial_resources_given == 32])
    #nestedness
    pairwise.wilcox.test(disturbance$biomass[disturbance$type =="coex"], disturbance$nestedness[disturbance$type =="coex"])
    
    
##Figure S11: Comparison Gupta-Interpolated R50
    ggplot(decay_parameter, aes(c50, Rc50))+
      geom_point()+
      ylab("Interpolated R50")+
      xlab("Van Genuchten-Gupta R50")+
      geom_abline(color="red", )+
      ggtitle("Comparison interpolated and optimized R50")+
      theme_bw()
    
    cor.test(decay_parameter$Rc50, decay_parameter$c50, method="pearson")


