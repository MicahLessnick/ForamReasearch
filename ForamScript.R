####Loading Data and Libraries####
library(vegan)
library(tidyverse)
library(readxl)

setwd("C:/Users/Micah Lessnick/Desktop/Foram research")
#change to project file as needed

Foram_Master_List <- read_excel("Foram Master List.xlsx",
sheet = "PDT - 2")

Foram_Site_Data <- read_excel("ROME_2021_summary.xlsx",
sheet = "SubsetTest")

#NMDS, PCA, Canico

####Prep Datasets####


####Prep for NMDS####

#metaMDS needs community data, grouping data, relative abundance of community data, distance matrix
#tutorial link: https://rpubs.com/CPEL/NMDS

#group data
Foram_group <- Foram_Master_List$LOCATION
#Foram_group <- Foram_Master_List$Group <- used when PDT - 1 is loaded
#Group data is set such that 1 is Rome 1, 2 is Rome 2, etc...

#community data
Foram_comm <- subset(Foram_Master_List, select = -c(LOCATION, DATE, SPECIES))
#Foram_comm <- Foram_Master_List[, 2:length(Foram_Master_List)-1] <- used when PDT - 1 is loaded

#add rownames <- used when PDT - 1 is loaded
#Foram_comm <- as.data.frame(Foram_comm)
#Foram_comm <- Foram_comm %>% remove_rownames() %>% column_to_rownames(var = "Sites")

#relative abundance [using methods from tutorial]
Foram_comm_rel <- decostand(Foram_comm, method = "total")

#calculating distance matrix [using method from tutorial]
Foram_comm_distMat <- vegdist(Foram_comm_rel, method = "bray")

#setting the distance matrix as matrix (and saving it to csv)
Foram_comm_distMat <- as.matrix(Foram_comm_distMat, labels = T)
write.csv(Foram_comm_distMat, file = "Foram_comm_distMat.csv")

####Computing NDMS####
#run NMDS analysis
Foram_comm_NDMS <- metaMDS(Foram_comm_distMat, distance = "bray", k = 3,
                           maxit = 999, trymax = 500, wascores = TRUE)

#save group labels into array
sites=c(Foram_group)

#plot NMDS points
plot(Foram_comm_NDMS, "sites")

#color-coordinate sites according to group and save color selections into separate vector
siteCol <- as.data.frame(sites)%>%
  mutate(Colors = case_when( 
         endsWith(sites, "1")~"red",
         endsWith(sites, "2")~"orange",
         endsWith(sites, "3")~"chartreuse3",
         endsWith(sites, "4")~"blue",
         endsWith(sites, "5")~"purple"))

color<-c(siteCol[,2])

#crease col vector based on Foram_group
orditorp(Foram_comm_NDMS, labels = sites, "sites", col = color)
#draw hulls around related sites
ordihull(Foram_comm_NDMS, groups = sites, draw="polygon",col="grey90",label=F)
#add environmental variable vectors to plot
envData <- select(Foram_Site_Data, -1:-5)
enVec <- envfit(Foram_comm_NDMS,envData, display = "sites", na.rm = TRUE)
plot(enVec)

####PCA####
#remove location, date, species columns
rda_temp <- Foram_Master_List[,-1:-3]
#rda with no parameters performs PCA
Foram_rda <- rda(rda_temp)

#plot PCA
biplot(Foram_rda,
       display = c("sites", 
                   "species"),
       type = c("text",
                "points"))
#generates errors due to 0-length vectors; code skips 0-length vectors
#current PCA is done with species, may need to transpose back to meaningful dataset