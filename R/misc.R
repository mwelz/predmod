#' Partition a vector into quantile groups
#'
#' Partitions a vector into quantile groups and returns a logical matrix indicating group membership.
#'
#' @param x A numeric vector to be partitioned.
#' @param cutoffs A numeric vector denoting the quantile cutoffs for the partition. Default are the quartiles: \code{c(0.25, 0.5, 0.75)}.
#'
#' @return
#' An object of type \code{quantile_group}, which is a logical matrix indicating group membership.
#'
#' @examples
#' set.seed(1)
#' x <- runif(100)
#' cutoffs <- c(0.25, 0.5, 0.75)
#' quantile_group(x, cutoffs)
#'
#' @export
quantile_group <- function(x,
                           cutoffs = c(0.25, 0.5, 0.75)){
  
  
  # input checks
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(cutoffs))
  stopifnot(0 < min(cutoffs) & max(cutoffs) < 1)
  
  # return
  quantile_group_NoChecks(x = x, cutoffs = cutoffs)
  
  
} # FUN


# same as above, just w/o input checks
quantile_group_NoChecks <- function(x = x,
                                    cutoffs = cutoffs){
  
  # get quantiles
  q         <- stats::quantile(x, cutoffs)
  q         <- c(-Inf, q, Inf)
  
  # check if breaks are unique: if x exhibits low variation, there might be empty quantile bins, which can cause an error in the cut() function. In this case, we add random noise to x to induce variation. NB: this bug has been spotted and fixed by Lucas Kitzmueller. All credits for this fix go to him!
  if(length(unique(q)) != length(q)){
    # specify standard deviation of the noise (x may have zero variation)
    sd <- ifelse(stats::var(x) == 0, 0.001, sqrt(stats::var(x) / 20))
    # add noise and updare quantiles
    x <- x + stats::rnorm(length(x), mean = 0, sd = sd)
    q <- stats::quantile(x, cutoffs)
    q <- c(-Inf, q, Inf)
  } # IF
  
  # get the grouping matrix
  group.mat <- group_matrix(x = x, breaks = q)
  
  # return
  return(structure(group.mat, type = "quantile_group"))
  
} # FUN


# returns a grouping matrix
group_matrix <- function(x, breaks)
{
  stopifnot(is.infinite(breaks[1L]) && is.infinite(breaks[length(breaks)]))
  
  groups    <- as.character(cut(x, breaks = breaks,
                                include.lowest = TRUE, right = FALSE, dig.lab = 3))
  group.nam <- unique(groups)
  group.nam <- group.nam[order(
    as.numeric(substr(sub("\\,.*", "", group.nam), 2, stop = 1e8L)),
    decreasing = FALSE)] # ensure the order is correct
  
  # get the grouping matrix
  group.mat <- sapply(1:length(group.nam), function(j) groups == group.nam[j])
  colnames(group.mat) <- gsub(",", ", ", gsub(" ", "", group.nam))
  return(group.mat)
}


intervals.quantile <- function(cutoffs){
  K   <- length(cutoffs) + 1
  nam <- rep(NA_character_, K)
  
  for(j in 1:K){
    if(j == 1){
      nam[j] <- paste0("<", 100*cutoffs[j], "%")
    } else if (j == K){
      nam[j] <- paste0(">=", 100*cutoffs[j-1], "%")
    } else{
      nam[j] <- paste0("[", 100*cutoffs[j-1], ",", 100*cutoffs[j], ")%")
    }
  } # FOR
  
  nam
  
} # FOR


# helper function to calculate the regression output. Accounts for the case of multicollinearity
# 'model' is typically a glm object
get_coefs <- function(model, cmprsk = FALSE)
{
  if(is.null(model)) return(NULL)
  
  smry <- summary(model)
  
  if(!cmprsk){
    cfs <- smry$coefficients
  } else{
    cfs <- smry$coef
  }
  
  aliased <- smry$aliased
  
  ## if there are aliased variables (i.e. collinear variables), add them to output
  if(any(aliased)){
    NA_mat <- matrix(NA_real_, nrow = sum(aliased), ncol = ncol(cfs))
    rownames(NA_mat) <- names(aliased)[aliased]
    cfs <- rbind(cfs, NA_mat)
  }
  cfs
} # FUN


# helper function for effect models
interacted_matrix <- function(X, w, interacted){
  
  X.nam <- colnames(X)
  
  if(is.character(interacted)){
    idx   <- which(X.nam %in% interacted)
    stopifnot(all(interacted %in% X.nam))
  } else if(is.numeric(interacted)){
    idx <- interacted
    stopifnot(max(interacted) <= ncol(X))
  } else stop("Indices must be either numeric or character")
  
  # make the matrix
  intmat <- sapply(idx, function(j) w * X[, j, drop = FALSE] )
  out <- cbind(X, w, intmat)
  colnames(out) <- c(X.nam, "w", paste0("w.", X.nam[idx]))
  out
  
} # FOR
