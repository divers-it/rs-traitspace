myclprofiles <- function (object, x, vars = NULL, col = NULL) 
{
  if (length(object$cluster) != nrow(x)) 
    stop("Size of x does not match cluster result!")
  if (is.null(vars)) 
    vars <- 1:ncol(x)
  if (!is.numeric(vars)) 
    vars <- sapply(vars, function(z) return(which(colnames(x) == 
                                                    z)))
  if (length(vars) < 1) 
    stop("Specified variable names do not match x!")
  if (is.null(col)) {
    k <- max(unique(object$cluster))
    if (k > 2) 
      col <- brewer.pal(k, "Set3")
    if (k == 2) 
      col <- c("lightblue", "orange")
    if (k == 1) 
      col <- "lightblue"
  }
  clusids <- sort(unique(object$cluster))
  if (length(col) != max(clusids)) 
    warning("Length of col should match number of clusters!")
  #par(ask = TRUE)
  for (i in vars) {
    if (is.numeric(x[, i])) {
      boxplot(x[, i] ~ object$cluster, col = col, main = colnames(x)[i])
      legend("topright", legend = clusids, fill = col)
    }
    if (is.factor(x[, i])) {
      tab <- table(x[, i], object$cluster)
      for (j in 1:length(object$size)) tab[, j] <- tab[, 
                                                       j]/object$size[j]
      barplot(t(tab), beside = TRUE, main = colnames(x)[i], 
              col = col)
    }
  }
  par(ask = FALSE)
  invisible()
}
