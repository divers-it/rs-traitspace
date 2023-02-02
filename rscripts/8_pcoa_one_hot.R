#Modified script using one-hot encoding allowing variables to be plotted

rm(list=ls())
library(dplyr)
library(ggplot2)
library(vegan)

#load formatted data
df2<-readRDS(file = here::here("outputs/df_filt_trans_one_hot.rds"))
df2[sapply(df2, is.factor)] <- lapply(df2[sapply(df2, is.factor)],
                                     as.integer)


#Not necessary anymore
#ADDED: one-hot encoding (should be done before)

#library(caret)
#dummy <- dummyVars(" ~ .", data=df)
#df <- data.frame(predict(dummy, newdata=df))

#ADDED: resolve polymorphisms

#polytraits=names(df)[grep("_",names(df))]

#for(trait in polytraits) {
#	print(trait)
#	basetrait=strsplit(trait[1],"\\.")[[1]][1]
#	states=strsplit(trait[1],"\\.")[[1]][2]
#	for(state in strsplit(states,"_")[[1]]){
#		targettrait=paste(basetrait,state,sep=".")
#		print(targettrait)
#		if(targettrait %in% colnames(df)){ #sometimes the trait column doesn't exist if the state only has been observed as a polymorphism
#			df[rownames(df[df[ , eval(trait)]==1,]),eval(targettrait)]=1
#			df[rownames( df[ !is.na(df[ , eval(trait)]) & df[ , eval(trait)]==1,] ),eval(targettrait)]=1
#			print(nrow(df))
#		}
#	}
#	 df=df[,!names(df) %in% trait] #remove polymorphic trait
#	print(nrow(df))
#}

#dissimilarity matrix calculation
library(cluster)
#df3=df2[,c(seq(1:57))]
#df3=df2[,c(seq(1:36),seq(40,50))]
#gower_df3 <- daisy(df3, metric = "gower" )
gower_df2 <- daisy(df2, metric = "gower" )

summary(gower_df2)

dataset_dist2 <- stats::as.dist(gower_df2)
#dataset_dist3 <- stats::as.dist(gower_df3)

#ADDED: change pcoa method (using the method from the "vegan" package)
dataset_pcoa2 <- wcmdscale(d = dataset_dist2, eig = TRUE, add = "lingoes")
#dataset_pcoa3 <- wcmdscale(d = dataset_dist3, eig = TRUE, add = "lingoes")
#dataset_pcoa <- ape::pcoa(dataset_dist)

#ADDED: plot PCOA with variables (function envfit from vegan)

traitvectors12=scores(envfit(dataset_pcoa2, df2, na.rm = T, add = T),display="vectors")
traitvectors34=scores(envfit(dataset_pcoa2, df2, na.rm = T, add = T,choices=c(3,4)),display="vectors")
traitd=as.data.frame(traitvectors12)
traitd$Dim3=as.data.frame(traitvectors34)$Dim3
traitd$Dim4=as.data.frame(traitvectors34)$Dim4
traitd$trait=rownames(traitd)

speciesv=scores(dataset_pcoa2,display="species")
speciesd=as.data.frame(speciesv)

# use original traits (multistate) for visualization
dft<-readRDS(file = here::here("outputs/df_filt_trans_one_hot.rds"))
df_ord=as.data.frame(c(dft,speciesd))

library(ggrepel)

p1=ggplot() + geom_point(data=df_ord,aes(x=Dim1,y=Dim2))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))

p2=ggplot() + geom_point(data=df_ord,aes(x=Dim2,y=Dim3))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim2/2,yend=Dim3/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim2/2,y=Dim3/2,label=trait))

p3=ggplot() + geom_point(data=df_ord,aes(x=Dim3,y=Dim4))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim3/2,yend=Dim4/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim3/2,y=Dim4/2,label=trait))

library(patchwork)
p1 / p2 / p3

ggsave("figures/pcoa_loadings_one_hot.png",
       width = 20,
       height = 40,
       units = 'cm')


#ggplot() + geom_point(data=df_ord,aes(x=Dim1,y=Dim2,shape=as.factor(Woodiness),color=as.factor(Mating)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))

#ggsave("figures/pcoa_woody_mating_one_hot.pdf")

#ggplot() + geom_point(data=df_ord,aes(x=Dim1,y=Dim2,shape=as.factor(Woodiness),color=as.factor(FlowerSex)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))

#ggsave("figures/pcoa_woody_flowersex_one_hot.pdf")

#ggplot() + geom_point(data=df_ord,aes(x=Dim1,y=Dim2,shape=as.factor(Woodiness),color=as.factor(SexualSystem)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))

#ggsave("figures/pcoa_woody_sexualsystem_one_hot.pdf")

#vegan plot methods (without ggplot)
#plot(dataset_pcoa,type="p")
#plot(envfit(dataset_pcoa, df, na.rm = T, add = T),p.max=0.1,cex=0.7)
#modify p.max to get more or less variables


