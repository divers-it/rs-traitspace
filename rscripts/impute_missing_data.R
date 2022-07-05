#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

require(missMDA)
data(ozone)
res.impute <- imputeFAMD(df, method="EM", ncp=10,maxiter = 1000) 
res.afdm <- FAMD(ozone,tab.disj=res.impute$tab.disj) 
FactoMineR::FAMD(df,ncp=50)
