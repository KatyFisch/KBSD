library(Rfast)


# kernel
rbf.normal <- function(x, point, disthalf) {
  gamma <- -log(0.5)/(disthalf)^2
  if (gamma == Inf){
    result <- ifelse(x==point, 1, 0)
  } else {
    result <- exp(-gamma * (x-point)^2 ) 
  }
  result
}

# calc diagnostic
singlediagnostic <- function(data_observed, obs_intervened, disthalf_vec, kernel=rbf.normal){
  foreach(i=1:nrow(data_observed), .combine='+') %:%
    foreach(coln=colnames(data_observed), .combine='*') %do% {
      if (coln %in% names(disthalf_vec)){
        res <- kernel(data_observed[i, coln], obs_intervened[[coln]], disthalf=unname(disthalf_vec[coln]))
      }
      else{1}
    }
}


#########



# DIAGNOSTIC
kbsd <- function(data, int_data_list, disthalf_vec, type="Rfast", kernel="rbf.normal", parallel=TRUE){
  
  # to handle division by zero issues
  disthalf_vec[disthalf_vec == 0] <- .Machine$double.xmin
  
  # set the kernel
  if (kernel=="rbf.normal"){
    kernel_used <- rbf.normal
  } else {
    kernel_used <- kernel
  }
  
  # only keep relevant variables for the observed data
  original_data <- data
  vartoconsider <- names(disthalf_vec)
  data <- data[vartoconsider]

  # apply the diagnostic to all interventions
  result_df <- data.frame(matrix(NA, nrow=0, ncol=3))
  names(result_df) <- c("intervention", "diagnostic", "observation")
  for (i in seq_along(int_data_list)){
    # intervened on data
    data_intervened <- int_data_list[[i]][vartoconsider]
    
    # calculation of the diagnostic
    if (type=="Rfast") {
      dat_int_scaled <- sweep(data_intervened, 2, 1/disthalf_vec, FUN='*')
      dat_scaled <- sweep(data, 2, 1/disthalf_vec, FUN='*')
      dia <- rowsums(exp(log(0.5)*dista(dat_int_scaled, dat_scaled, square = TRUE)), parallel = parallel)
    } else {
      dia <- apply(data_intervened, 1, function(x){singlediagnostic(data, x, disthalf_vec,  
                                                                    kernel=kernel_used)})
    }
    new_data <- data.frame(shift= rep(i, length(dia)), diagnostic=unname(dia), observation=1:length(dia) )
    result_df <- rbind(result_df, new_data)
  }
  
  # return result
  return(result_df)
}

# PLOTTING
kbsd_plot <- function(result_df){
  groups <- split(result_df$diagnostic, result_df$shift)
  
  # Quantile berechnen
  q <- lapply(groups, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
  
  # In Matrix umwandeln
  stats <- sapply(q, identity)
  
  # Boxplot zeichnen
  bx <- bxp(list(
    stats = stats,
    n = sapply(groups, length),
    names = names(groups)
  ),
  ylim = c(0, max(result_df$diagnostic, na.rm = TRUE)),
  xlab = "intervention",
  ylab = "EDP",
  main = "KBSD: Adjusted Boxplots"
  )
  abline(h = 0, col = "lightgray", lwd = 1, lty = 2)
  
  # Ausreißer (unter 5% oder über 95%) hinzufügen
  for(i in seq_along(groups)) {
    x <- groups[[i]]
    out <- x[x < stats[1,i] | x > stats[5,i]]
    points(rep(i, length(out)), out, pch = 16, cex=0.5)
  }
}

get_disthalf_vec <- function(data, interventions = NULL, remove = NULL, no_extrapolation = NULL){
  
  # standard deviations (numeric only)
  res <- sapply(data, function(x) {
    if (is.factor(x) || is.logical(x)) {
      0
    } else {
      sd(x, na.rm = TRUE)
    }
  })
  
  # no extrapolation, e.g. binary variables
  if (!is.null(no_extrapolation)) {
    res[names(res) %in% no_extrapolation] <- 0
  }
  
  # remove outcome etc.
  if (!is.null(remove)) {
    res <- res[setdiff(names(res), remove)]
  }
  
  # halve intervention variables
  if (!is.null(interventions)) {
    res[interventions] <- res[interventions] / 2
  }
  
  return(res)
}


#### Test Function

# simulated data

library(simcausal)
n=100

## CREATE DAG
D.base <- DAG.empty() + 
  node("L", distr = "rbern", prob = 0.3) +
  node("A",  distr = "rnorm", mean = (0.21)*L + 0.45, sd=0.14)
Dset <- simcausal::set.DAG(D.base)

# get data for plotting
Odat <- simcausal::sim(DAG = Dset, n = n, rndseed = 1, verbose=F) 
Odat$ID <- NULL

# 
Odat_int1 <- Odat
Odat_int1$A = 1

Odat_int2 <- Odat
Odat_int2$A = Odat_int2$A + 0.1

Odat_int3 <- Odat
Odat_int3$A = Odat_int3$A + 0.5



data_result <- kbsd(data=Odat, 
                             int_data_list=list(Odat_int1, Odat_int2, Odat_int3), 
                             disthalf_vec=c(L=0.5, A=0.1)
                             )

kbsd_plot(data_result)



# EFV single time point

library(CICI)

data(EFV)
efv_subset <- EFV[, c("log_age", "efv.0", "adherence.1")]

# intervention 1: static intervention = 1
efv_sub1 <- transform(efv_subset, efv.0 = 1)

# intervention 2: Modified Treatment Policy increment + 1
efv_sub2 <- transform(efv_subset, efv.0 = efv.0 + 1)

# apply diagnostic
data_result <- kbsd(data=efv_subset, 
                    int_data_list=list(efv_sub1, efv_sub2), 
                    disthalf_vec=c(log_age=1, efv.0=0.5, adherence.1=0)
)
# plot diagnostic
kbsd_plot(data_result)


# EFV multiple time points

# strategy 1
EFV_strat1 <- EFV; EFV_strat1[paste0("efv.", 0:4)] <- 3
# strategy 2
EFV_strat2 <- EFV; EFV_strat2[paste0("efv.", 0:4)] <- EFV_strat2[paste0("efv.", 0:4)] + 0.1


# apply diagnostic
data_result <- kbsd(data = EFV, 
                    int_data_list = list(EFV_strat1, EFV_strat2), 
                    disthalf_vec = get_disthalf_vec(EFV, interventions=paste0("efv.", 0:4), remove=c("VL.4", "metabolic", "NRTI"), no_extrapolation="sex")
)
# plot diagnostic
kbsd_plot(data_result)


