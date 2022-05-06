  # Read in trait data ----

  library(dplyr)

  #load formatted data
  df<-readRDS(file = here::here("outputs/df_filt.rds"))

  #numeric columns only
  nums <- unlist(lapply(df, is.numeric))
  facts <- unlist(lapply(df, is.factor))

  df2<-cbind(df[ , nums],df[ , facts])

  str(df2)



  # Correlation between Traits ----

  dataset<-df2

  dataset_cor <- matrix(0, ncol(dataset), ncol(dataset))

  for (i in 1:ncol(dataset)) {

    for (j in i:ncol(dataset)) {

      dataset_cor[i, j] <- stats::cor(
        x      = rank(dataset[ , i]),
        y      = rank(dataset[ , j]),
        method = "kendall"
      )
    }
  }

  dataset_cor[lower.tri(dataset_cor)] <- t(
    dataset_cor)[lower.tri(dataset_cor)]

  diag(dataset_cor) <- NA

  dataset_cor <- data.frame(
    mean_cor = mean(abs(dataset_cor), na.rm = TRUE),
    sd_cor   = sd(abs(dataset_cor),   na.rm = TRUE),
    max_cor  = max(abs(dataset_cor),  na.rm = TRUE),
    min_cor  = min(abs(dataset_cor),  na.rm = TRUE)
  )

  view(dataset_cor)


