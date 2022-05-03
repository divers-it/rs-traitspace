##Script to recode PROTEUS data, removing polymorphisms, etc.

#input discrete data
discdata=read.csv("discrete.csv",colClasses="character",na.strings="")
discdata=discdata[ , colSums(is.na(discdata)) < nrow(discdata)]
#input continuous data
contdata=read.csv("continuous_c.csv")
contdata=contdata[ , colSums(is.na(contdata)) < nrow(contdata)]
#merge
alldata=merge(contdata,discdata,by="X")

#get basic stats
missing=colSums(is.na(alldata))
nspecies=length(alldata$X)
filled=nspecies-missing
filled=as.data.frame(filled)
filled=data.frame(trait=rownames(filled),n=filled[,1])
filled=filled[filled$trait!="X",]
#indicate whether they are DiveRS or eFLOWER traits
filled$origin="PROTEUS"
filled$retain="no"
filled[filled$trait=="Pollinationsyndrome",]$origin="DiveRS"
filled[filled$trait=="Lifehistory",]$origin="DiveRS"
filled[filled$trait=="Leafphenology",]$origin="DiveRS"
filled[filled$trait=="Flower.leaftemporalpattern",]$origin="DiveRS"
filled[filled$trait=="Dispersalsyndrome",]$origin="DiveRS"
filled[filled$trait=="Self.incompatibilitysystem.genetic.",]$origin="DiveRS"
filled[filled$trait=="Autonomousselfing",]$origin="DiveRS"
filled[filled$trait=="Phenotypicmatingsystem",]$origin="DiveRS"
filled[filled$trait=="Inflorescenceattractivecolor",]$origin="DiveRS"
filled[filled$trait=="Inflorescenceattractiveorgan",]$origin="DiveRS"
filled[filled$trait=="Pseudanthiumtype",]$origin="DiveRS"
filled[filled$trait=="Pseudanthiumshape",]$origin="DiveRS"
filled[filled$trait=="Pseudanthia",]$origin="DiveRS"
filled[filled$trait=="Heteranthery",]$origin="DiveRS"
filled[filled$trait=="Fruitfleshiness",]$origin="DiveRS"
filled[filled$trait=="Fruitdehiscence",]$origin="DiveRS"
filled[filled$trait=="Dichogamy",]$origin="DiveRS"
filled[filled$trait=="Herkogamy",]$origin="DiveRS"
filled[filled$trait=="Outcrossingrate",]$origin="DiveRS"
filled[filled$trait=="Pollen.ovuleratio",]$origin="DiveRS"
filled[filled$trait=="Floralreward",]$origin="DiveRS"
#indicate whether they are DiveRS core traits
filled[filled$trait=="Dispersalsyndrome",]$retain="yes"
filled[filled$trait=="Self.incompatibilitysystem.genetic.",]$retain="yes"
filled[filled$trait=="Phenotypicmatingsystem",]$retain="yes"
filled[filled$trait=="Pollinationsyndrome",]$retain="yes"
filled[filled$trait=="Lifehistory",]$retain="yes"
filled[filled$trait=="Outcrossingrate",]$retain="yes"
filled[filled$trait=="Pollen.ovuleratio",]$retain="yes"
filled[filled$trait=="Floralreward",]$retain="yes"
filled[filled$trait=="Plantsexualsystem",]$retain="yes"
filled[filled$trait=="Floralstructuralsex",]$retain="yes"
filled[filled$trait=="Ovaryposition",]$retain="yes"
filled[filled$trait=="Symmetryofperianth",]$retain="yes"
filled[filled$trait=="Habit",]$retain="yes"
filled[filled$trait=="Numberoffertilestamens",]$retain="yes"
filled[filled$trait=="Numberofstructuralcarpels",]$retain="yes"
filled[filled$trait=="Numberofovulesperfunctionalcarpel",]$retain="yes"
filled[filled$trait=="Maincolorofperianthatanthesis",]$retain="yes"
filled[filled$trait=="Flowerdiameter",]$retain="yes"
filled[filled$trait=="Flowerlength",]$retain="yes"
#create and save summary graph
library(ggplot2)
ggplot(data=filled,aes(y=reorder(trait,n),x=100*n/nspecies,fill=origin,alpha=retain))+geom_bar(stat="identity")+labs(x="percentage",y="") + theme(legend.position="none")
ggsave("progress.pdf")

cleandata=data.frame(species=alldata$X)
datanames=data.frame(name=character(),st00=character(),st01=character(),st02=character(),st03=character(),st04=character(),st05=character(),st05=character())

#Plantsexualsystem
cleandata$sexmorphs <- "?"
#first monomorphic then dimorphic: suppose that when dimorphic has been reported, that refutes monomorphy
cleandata[grep("[01357]", alldata$Plantsexualsystem),]$sexmorphs="0"
cleandata[grep("[246]", alldata$Plantsexualsystem),]$sexmorphs="1"

newRow <- data.frame(name='sexmorphs',st00='monomorphic',st01='dimorphic',st02='',st03='',st04='',st05='') 
datanames=rbind(datanames,newRow)

#Floralstructuralsex
cleandata$flowersex <- "?"
#polymorphic
cleandata[!is.na(alldata$Floralstructuralsex),]$flowersex="0"
#bisexual
cleandata[!is.na(alldata$Floralstructuralsex) & alldata$Floralstructuralsex=="0",]$flowersex="2"
#unisexual
cleandata[!is.na(alldata$Floralstructuralsex) & (alldata$Floralstructuralsex=="1" | alldata$Floralstructuralsex=="2" | alldata$Floralstructuralsex=="12"),]$flowersex="1"
newRow <- data.frame(name='flowersex',st00='polymorphic',st01='unisexual',st02='bisexual',st03='',st04='',st05='') 
datanames=rbind(datanames,newRow)

#Symmetryofperianth
cleandata$symmetry <- "?"
cleandata[!is.na(alldata$Symmetryofperianth) & alldata$Symmetryofperianth=="2",]$symmetry="1"
cleandata[!is.na(alldata$Symmetryofperianth) & (alldata$Symmetryofperianth=="0"| alldata$Symmetryofperianth=="1" | alldata$Symmetryofperianth=="5" | alldata$Symmetryofperianth=="6") ,]$symmetry="0"
newRow <- data.frame(name='symmetry',st00='actinomorphic',st01='zygomorphic',st02='',st03='',st04='',st05='') 
datanames=rbind(datanames,newRow)

#Pollinationsyndrome
cleandata$pollination <- "?"
#first: suppose all data we have are biotically pollinated plants, we'll overwrite the other cases later
cleandata[!is.na(alldata$Pollinationsyndrome),]$pollination="1"
#only autonomous pollination: remove
cleandata[!is.na(alldata$Pollinationsyndrome) & alldata$Pollinationsyndrome=="H",]$pollination="?"
#abiotic
cleandata[grep("[01]",alldata$Pollinationsyndrome),]$pollination="0"
#mixed
cleandata[grep("([01])([2-9A-G])",alldata$Pollinationsyndrome),]$pollination="0&1"
newRow <- data.frame(name='pollination',st00='abiotic',st01='biotic',st02='',st03='',st04='',st05='') 
datanames=rbind(datanames,newRow)

#Habit
cleandata$habit <- "?"
cleandata[!is.na(alldata$Habit) & (alldata$Habit=="0" | alldata$Habit=="1" | alldata$Habit=="2" | alldata$Habit=="01" | alldata$Habit=="02" | alldata$Habit=="12" | alldata$Habit=="012"),]$habit="1" 
cleandata[!is.na(alldata$Habit) & (alldata$Habit=="3" | alldata$Habit=="4" | alldata$Habit=="5" | alldata$Habit=="34" | alldata$Habit=="35" | alldata$Habit=="45" | alldata$Habit=="345"),]$habit="0" 
newRow <- data.frame(name='habit',st00='herb',st01='woody',st02='',st03='',st04='',st05='') 
datanames=rbind(datanames,newRow)

#Lifehistory
cleandata$lifespan <- "?"
#species that can be annual or biennial are considered short-lived, no matter if they can be perennial as well
cleandata[!is.na(alldata$Lifehistory),]$lifespan="1"
cleandata[grep("[01]",alldata$Lifehistory),]$lifespan="0"
newRow <- data.frame(name='lifespan',st00='short',st01='long',st02='',st03='',st04='',st05='') 
datanames=rbind(datanames,newRow)

#Dispersalsyndrome
#not enough data yet

#Flower size 
alldata[is.na(alldata$Flowerdiameter),]$Flowerdiameter=0
alldata[is.na(alldata$Flowerlength),]$Flowerlength=0
flowersize=pmax(alldata$Flowerdiameter,alldata$Flowerlength)
cleandata[flowersize > 10,]$flowersize="2"
cleandata[flowersize <= 10,]$flowersize="1"
cleandata[flowersize <= 1,]$flowersize="0"
cleandata[flowersize < 0.001,]$flowersize="?"
newRow <- data.frame(name='flowersize',st00='less_than_1cm',st01='between_1_and_10cm',st02='more_than_10cm',st03='',st04='',st05='') 
datanames=rbind(datanames,newRow)

#Floralreward
#not enough data yet

#Maincolorofperianthatanthesis
#later

#Ovaryposition
cleandata$ovaryposition <- "?"
cleandata[!is.na(alldata$Ovaryposition),]$ovaryposition="1"
cleandata[!is.na(alldata$Ovaryposition) & alldata$Ovaryposition=="0",]$ovaryposition="2"
cleandata[!is.na(alldata$Ovaryposition) & alldata$Ovaryposition=="1",]$ovaryposition="0"
newRow <- data.frame(name='ovaryposition',st00='inferior',st01='intermediate/variable',st02='superior',st03='',st04='',st05='') 
datanames=rbind(datanames,newRow)

#Numberofstructuralcarpels
#Numberofovulesperfunctionalcarpel

#cleandata$numberofovules <- alldata$Numberofovulesperfunctionalcarpel*alldata$Numberofstructuralcarpels
#cleandata[!is.na(cleandata$numberofovules) & cleandata$numberofovules>9998,]$numberofovules=9999
cleandata$numberofovules <- "?"
numberofovules <- alldata$Numberofovulesperfunctionalcarpel*alldata$Numberofstructuralcarpels
cleandata[!is.na(numberofovules) & numberofovules > 50.5,]$numberofovules="5"
cleandata[!is.na(numberofovules) & numberofovules < 50.5,]$numberofovules="4"
cleandata[!is.na(numberofovules) & numberofovules < 10.5,]$numberofovules="3"
cleandata[!is.na(numberofovules) & numberofovules < 3.5,]$numberofovules="2"
cleandata[!is.na(numberofovules) & numberofovules < 2.5,]$numberofovules="1"
cleandata[!is.na(numberofovules) & numberofovules < 1.5,]$numberofovules="0"
#cleandata[!is.na(cleandata$numberofovules) & cleandata$numberofovules>9998,]$numberofovules=9999
newRow <- data.frame(name='numberofovules',st00='one',st01='two',st02='three',st03='four_to_ten',st04='ten_to_fifty',st05='more_than_fifty') 
datanames=rbind(datanames,newRow)

#Numberoffertilestamens
cleandata$numberofstamens="?"
cleandata[!is.na(alldata$Numberoffertilestamens) & alldata$Numberoffertilestamens > 50.5,]$numberofstamens="5"
cleandata[!is.na(alldata$Numberoffertilestamens) & alldata$Numberoffertilestamens < 50.5,]$numberofstamens="4"
cleandata[!is.na(alldata$Numberoffertilestamens) & alldata$Numberoffertilestamens < 10.5,]$numberofstamens="3"
cleandata[!is.na(alldata$Numberoffertilestamens) & alldata$Numberoffertilestamens < 3.5,]$numberofstamens="2"
cleandata[!is.na(alldata$Numberoffertilestamens) & alldata$Numberoffertilestamens < 2.5,]$numberofstamens="1"
cleandata[!is.na(alldata$Numberoffertilestamens) & alldata$Numberoffertilestamens < 1.5,]$numberofstamens="0"
newRow <- data.frame(name='numberofstamens',st00='one',st01='two',st02='three',st03='four_to_ten',st04='ten_to_fifty',st05='more_than_fifty') 
datanames=rbind(datanames,newRow)

#Outcrossingrate
#cleandata$outcrossingrate=alldata$Outcrossingrate

#Pollen.ovuleratio
#cleandata$poratio=alldata$Pollen.ovuleratio

#Phenotypicmatingsystem
#add information from other traits
cleandata$matingsystem <- "?"
cleandata[!is.na(alldata$Phenotypicmatingsystem) & alldata$Phenotypicmatingsystem=="0",]$matingsystem="0"
cleandata[!is.na(alldata$Phenotypicmatingsystem) & alldata$Phenotypicmatingsystem=="1",]$matingsystem="1"
cleandata[!is.na(alldata$Phenotypicmatingsystem) & alldata$Phenotypicmatingsystem=="2",]$matingsystem="2"
cleandata[!is.na(alldata$Self.incompatibilitysystem.genetic.) & (alldata$Self.incompatibilitysystem.genetic.=="1" |  alldata$Self.incompatibilitysystem.genetic.=="2" | alldata$Self.incompatibilitysystem.genetic.=="3"),]$matingsystem="2"
cleandata[!is.na(alldata$Plantsexualsystem) & alldata$Plantsexualsystem=="2",]$matingsystem="2"
cleandata[!is.na(alldata$Outcrossingrate),]$matingsystem="1"
cleandata[!is.na(alldata$Outcrossingrate) & alldata$Outcrossingrate < 0.2,]$matingsystem="0"
cleandata[!is.na(alldata$Outcrossingrate) & alldata$Outcrossingrate > 0.8,]$matingsystem="2"
newRow <- data.frame(name='matingsystem',st00='selfing',st01='intermediate',st02='outcrossing',st03='',st04='',st05='') 
datanames=rbind(datanames,newRow)

#Self.incompatibilitysystem.genetic.

cleandata$compatibility <- "?"
cleandata[!is.na(alldata$Self.incompatibilitysystem.genetic.) & alldata$Self.incompatibilitysystem.genetic.=="0",]$compatibility = "0"
cleandata[!is.na(alldata$Self.incompatibilitysystem.genetic.) & alldata$Self.incompatibilitysystem.genetic.=="1",]$compatibility = "1"
cleandata[!is.na(alldata$Self.incompatibilitysystem.genetic.) & alldata$Self.incompatibilitysystem.genetic.=="2",]$compatibility = "1"
cleandata[!is.na(alldata$Self.incompatibilitysystem.genetic.) & alldata$Self.incompatibilitysystem.genetic.=="3",]$compatibility = "1"
newRow <- data.frame(name='compatibility',st00='SC',st01='SI',st02='',st03='',st04='',st05='') 
datanames=rbind(datanames,newRow)

#Autonomousselfing
#cleandata$autonomouselfing <- "?"
#cleandata[!is.na(alldata$Autonomousselfing),]$autonomouselfing="1"
#cleandata[!is.na(alldata$Autonomousselfing) & alldata$Autonomousselfing=="0",]$autonomouselfing="2"
#cleandata[!is.na(alldata$Autonomousselfing) & alldata$Autonomousselfing=="2",]$autonomouselfing="1"
#newRow <- data.frame(name='autonomousselfing',st00='no',st01='variable',st02='yes',st03='',st04='',st05='') 
#datanames=rbind(datanames,newRow)

write.table(row.names=FALSE, sep=";",datanames,file="datanames.csv",quote=FALSE)
write.table(row.names=FALSE, sep=",",cleandata,file="cleandata.csv",quote=FALSE)


