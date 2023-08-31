library(ape)
library(phytools)
library(geiger)
library(castor)
library(cluster)
library(ggplot2)
library(ggpubr)
library(dplyr)



# Loading the tree and the trait matrix ####
angio_tree <- read.tree("outputs/pruned_tree.tre")
df_final <- readRDS("outputs/6_df_filt_trans.rds")
df_onehot <- readRDS("outputs/one_hot_6_df_filt_trans.rds")
df_recoded <- readRDS("outputs/imputed_recoded.RDS")

# Split the initial matrix into continuous and discrete trait matrices
species_names <- angio_tree$tip.label
continous_traits <- df_recoded[,c(1:7)]
discrete_traits <- df_recoded[,c(8:21)]

# Fitting evolution models for each trait ####

# Function to fit an OU model to a continuous trait and return the model parameters
fit_cont <- function(ttf) {
  # rates
  trait_to_fit <- ttf
  names(trait_to_fit) <- species_names
  fit <- fitContinuous(
    phy = angio_tree,
    dat = trait_to_fit,
    model = "OU",
    bounds = list(alpha=c(exp(-500),exp(5))) # max extended a bit / default
    )
  # ancestral state under the BM model
  # OU should be used but it's very unstable
  anc <- reconstruct(ttf,angio_tree,method = "ML")
  mean_root <- anc$ace[1] # The first node is the root
  sd_root <- sqrt((anc$CI95[1,2]-anc$CI95[1,1])/(2*1.96))
  param <- c(unlist(fit$opt[c(1:3)]),mean_root,sd_root)
  names(param) <- c("alpha","sigma","z0","mean_root","sd_root")
  return(param)
}
# creating and saving the list of OU parameters
ou_list <- lapply(continous_traits,fit_cont)
saveRDS(ou_list,"outputs/ou_list.rds")

# Function to transform the list of rates into the corresponding transition matrix
# For the fitDiscrete function, not useful for the asr_mk_model function
# rate2mat <- function(r) {
#   k <- length(r)
#   n <- as.integer((sqrt(4*k+1) + 1)/2)
#   for(i in c(0:(n-1))) r <- append(r,0,i*n+i)
#   mat <- matrix(r,n,n)
#   for(i in c(1:n)) mat[i,i] <- -colSums(mat)[i]
#   return(mat)
# }

# Function to fit a Mk model to a discrete trait and return the transition matrix
fit_disc <- function(ttf) {
  trait_to_fit <- as.numeric(ttf)
  names(trait_to_fit) <- species_names
  # Much more rapid than the fitDiscrete function
  fit <- asr_mk_model(angio_tree,as.numeric(ttf),rate_model = "ARD")
  #fit <- fitDiscrete(angio_tree,trait_to_fit,model = "ARD")
  #rates <- unlist(fit$opt[1:fit$opt$k])
  #mat <- rate2mat(rates)
  print("done")
  return(list("mat"=fit$transition_matrix,"anc_st"=fit$ancestral_likelihoods[1,]))
}
# creating and saving the list of transition matrices (take some time)
mk_list <- lapply(discrete_traits,fit_disc)
saveRDS(mk_list,"outputs/mk_list.rds")

# Simulating independent traits according to the fitted evolution model ####

# Function to simulate a dataset according to the evolutionary processes
# Take as in put:
# - A list of vectors with OU parameters + anc state for continuous characters
# - A list of matrices with transition rates + anc state proba for discrete characters
# - A tree
# Return a data.frame with simulated values

simul_dataset <- function(ou,mk,tree) {
  # continuous traits
  sim_continuous <- lapply(
    ou,
    function(x) rTraitCont(
      phy = tree,
      model = "OU",
      alpha = x[1],
      sigma = x[2],
      theta =  x[3],
      root.value = rnorm(1,x[4],x[5]))
  )
  # discrete traits
  sim_discrete <- lapply(
    mk,
    function(x) {
      st <- as.character(1:length(x$anc_st))
      root <- sample(st,1,prob = x$anc_st)
      as.factor(
        sim.history(
          tree = tree,
          Q = x$mat,
          anc = root,
          nsim = 1,
          message=F
        )$states
      )
    }
  )
  df <- data.frame(sim_continuous,sim_discrete)
  return(df)
}


#  Simulate n datasets ####
# Can be uncommented to resimulate
#ou_list <- readRDS("outputs/ou_list.rds")
#mk_list <- readRDS("outputs/mk_list.rds")
#n_sim <- 1000
#set.seed(123)
#simdata <- replicate(n_sim,simul_dataset(ou_list,mk_list,angio_tree),simplify = F)
#saveRDS(simdata,"outputs/phylo_simulated_datasets.rds")


# Analyzing simulations ####

sim_list <- readRDS("outputs/phylo_simulated_datasets.rds")

# Function to compute gower distance, run pcoa and return eigenvalues from one simulated data_frame
pc_ev <- function(df) {
  # We only keep the number of axes equal to the number of variables
  n_var <- dim(df)[2]
  # To ensure that discrete variables are factors
  df[,c(8:21)] <- lapply(df[,c(8:21)],as.factor)
  dist <- daisy(df,"gower")
  pc <- pcoa(dist)
  r_eig <- pc$values$Relative_eig
  r_eig <- r_eig[r_eig>0]/sum(r_eig[r_eig>0])
  return(r_eig[c(1:n_var)])
}

# Can be uncommented to recompute eigenvalues
# sim_ev <- sapply(sim_list,pc_ev)
# saveRDS(data.frame(t(sim_ev)),"outputs/phylo_simulated_ev.rds")


# Analysing the results

# Creating a dataset with observed an simulated ev
sim_ev <- readRDS("outputs/phylo_simulated_ev.rds")
sim_ev_mean <- apply(sim_ev,2,mean)
sim_ev_min <- apply(sim_ev,2,function(x) quantile(x,probs = c(0.05)))
sim_ev_max <- apply(sim_ev,2,function(x) quantile(x,probs = c(0.95)))
ef <- pcoa(daisy(df_final,"gower"))$values$Relative_eig
ef <- ef[ef>0]/sum(ef[ef>0])
final_ev <- ef[c(1:21)]
eo <- pcoa(daisy(df_onehot,"gower"))$values$Relative_eig
eo <- eo[eo>0]/sum(eo[eo>0])
onehot_ev <- eo[c(1:21)]
er <- pcoa(daisy(df_recoded,"gower"))$values$Relative_eig
er <- er[er>0]/sum(er[er>0])
recoded_ev <- er[c(1:21)]
rm(ef,eo,er)
df_ev <- data.frame(
  list(
    "sim_mean"=sim_ev_mean,
    "sim_min"=sim_ev_min,
    "sim_max"=sim_ev_max,
    "final"=final_ev,
    "onehot"=onehot_ev,
    "recoded"=recoded_ev)
)
# Plot
ggplot(data = df_ev,aes(x=c(1:21))) +
  geom_line(aes(y=final,col="final")) +
  geom_point(aes(y=final,col="final")) +
  geom_line(aes(y=onehot,col="onehot")) +
  geom_point(aes(y=onehot,col="onehot")) +
  geom_line(aes(y=recoded,col="recoded")) +
  geom_point(aes(y=recoded,col="recoded")) +
  geom_errorbar(aes(ymin=sim_min,ymax=sim_max)) +
  xlab("Rank of eigenvalues") + ylab("Relative eigenvalues") +
  labs(col = "Trait matrix")
ggsave("figures/observed_vs_simulated_eigenvalues.pdf")


# Looking at the traits in the simulated datasets
# As traits are independent, each axis should +/- to one trait
# In most simulations woodiness appears on Axis.1
# On other axis it's often the same other traits: 3 examples  are given
# TO INVESTIGATE MORE: it should be linked to the rate of evolution of the traits

df_sim <- sim_list[[1]] # pick one dataset
df_sim[,c(8:21)] <- lapply(df_sim[,c(8:21)],as.factor)
pp <- data.frame(pcoa(daisy(df_sim,"gower"))$vectors)
# Coloring the pcoa with different label
g1 <- ggplot(data=pp,aes(x=Axis.1,y=Axis.2,col=df_sim$woodiness)) + geom_point()
g2 <- ggplot(data=pp,aes(x=Axis.2,y=Axis.3,col=df_sim$dispersaldist)) + geom_point()
ggarrange(plotlist = list(g1,g2),ncol = 1,nrow = 2)
ggsave("figures/pcoa_sim_ex1.pdf")

df_sim <- sim_list[[2]] # pick one dataset
df_sim[,c(8:21)] <- lapply(df_sim[,c(8:21)],as.factor)
pp <- data.frame(pcoa(daisy(df_sim,"gower"))$vectors)
# Coloring the pcoa with different label
g1 <- ggplot(data=pp,aes(x=Axis.1,y=Axis.2,col=df_sim$woodiness)) + geom_point()
g2 <- ggplot(data=pp,aes(x=Axis.2,y=Axis.3,col=df_sim$showiness)) + geom_point()
ggarrange(plotlist = list(g1,g2),ncol = 1,nrow = 2)
ggsave("figures/pcoa_sim_ex2.pdf")

df_sim <- sim_list[[3]] # pick one dataset
df_sim[,c(8:21)] <- lapply(df_sim[,c(8:21)],as.factor)
pp <- data.frame(pcoa(daisy(df_sim,"gower"))$vectors)
# Coloring the pcoa with different label
g1 <- ggplot(data=pp,aes(x=Axis.1,y=Axis.2,col=df_sim$woodiness)) + geom_point()
g2 <- ggplot(data=pp,aes(x=Axis.2,y=Axis.3,col=df_sim$flowersex)) + geom_point()
ggarrange(plotlist = list(g1,g2),ncol = 1,nrow = 2)
ggsave("figures/pcoa_sim_ex3.pdf")



# df_rates <- data.frame(
#   trait = as.factor(names(mk_list)),
#   #rate = as.numeric(unlist(lapply(mk_list,function(x) mean(abs(x$mat)))))
#   rate = as.numeric(unlist(lapply(mk_list,function(x) max(abs(x$mat)))))
#   )
# 
# df_rates %>%
#   arrange(desc(rate)) %>%
#   mutate(trait=factor(trait, levels=trait)) %>%
#   ggplot(aes(x=trait, y=rate,fill=trait)) + 
#     geom_bar(stat = "identity",show.legend = F) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     xlab("Traits") + ylab("Evolutionary rates")
