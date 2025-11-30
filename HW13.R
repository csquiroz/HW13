### Cody Quiroz ###
### Zoo 800 ###
### 30-Nov-2025 ###


# Objective 1 -------------------------------------------------------------
#LR param est: analytical solution

#load data
dragons <- read.csv("dragon_data.csv")
x <- dragons$size
y <- dragons$acres_on_fire

X <- cbind(1,x) #X

XtX <- t(X) %*% X #compute (X^T X)
XtY <- t(X) %*% y #compute (X^T y)

#bera hat = (X^T X)^(-1) * (X^T y)
beta_hat <- solve(XtX) %*% XtY
beta_hat

#get int and slope
int <- beta_hat[1]
slope <- beta_hat[2]

cat("Analytical Intercept =", int)
cat("Analytical Slope =", slope)


# Objective 2 -------------------------------------------------------------
#LR param est: ordinary least squares (OLS)

## part a.
#grid search (brute force)
m <- length(y)
#make SSE function (sum of squared errs)
SSE <- function(a, b) {
  sum((y - (a + b*x))^2)
}
#grid
int_grid <- seq(-20,100, by=0.1)
slope_grid <- seq(-5,5,by=0.1)

#matrix to store SSE
SSE_matrix <- matrix(NA, nrow = length(int_grid), ncol = length(slope_grid))
for (i in seq_along(int_grid)) {
  for (j in seq_along(slope_grid)) {
    SSE_matrix[i, j] <- SSE(int_grid[i], slope_grid[j])
  }
}

#find best (a,b)
min_idx <- which.min(SSE_matrix)         
idx_row <- arrayInd(min_idx, dim(SSE_matrix))[1]
idx_col <- arrayInd(min_idx, dim(SSE_matrix))[2]
best_int_grid <- int_grid[idx_row]
best_slope_grid     <- slope_grid[idx_col]

cat("OLS Intercept =", best_int_grid)
cat("OLS Slope =", best_slope_grid)

##part b.
#uing optim()

#SSE as function for optim
min_SSE <- function(par, x, y) {
  a <- par[1]
  b <- par[2]
  sum((y - (a + b*x))^2)
}

#initial guess
init <- c(0.2, 1.3)

#run optim
ols_result <- optim(par = init, fn = min_SSE, x = x, y = y)
ols_result$par #est slope and int

##part c.
#verify optimization converged and not sensitive to start vals

#try a bunch of diff initial (a,b)
starts <- list(c(0, 1),
  c(10, -2),
  c(-37, 5),
  c(100, -10),
  c(-50, 67))

#empty matrix for storage
results <- matrix(NA, ncol = 2, nrow = length(starts))
colnames(results) <- c("intercept", "slope") #make column names

#loop thru each pair using optim
for (i in seq_along(starts)) {
  fit <- optim(par = starts[[i]], fn = min_SSE, x = x, y = y)
  results[i, ] <- fit$par
}
results #all values within 0.01 of each other, therefore not sensitive to starting values


# Objective 3 -------------------------------------------------------------
#LR param est: maximum likelihood est (MLE)

##part a.
#grid search

#negative log likelihood function
NLL <- function(a,b){
  sigma <- sqrt(mean((y-(a+b*x))^2))
  n <- length(y)
  nll <- (n/2)*log(2*pi*sigma^2) + sum((y-(a+b*x))^2)/(2*sigma^2)
  return(nll)
}

#storage, use same grid as obj 2
NLL_matrix <- matrix(NA, nrow=length(int_grid),ncol=length(slope_grid))

for (i in seq_along(int_grid)) {
  for (j in seq_along(slope_grid)) {
    NLL_matrix[i,j] <- NLL(int_grid[i], slope_grid[j])
  }
}

#find minimum
min_idx2 <- which.min(NLL_matrix)
idx_row2 <- arrayInd(min_idx2, dim(NLL_matrix))[1]
idx_col2 <- arrayInd(min_idx2, dim(NLL_matrix))[2]
best_int_grid2 <- int_grid[idx_row2]
best_slope_grid2 <- slope_grid[idx_col2]

cat("MLE Intercept =", best_int_grid2)
cat("MLE Slope =", best_slope_grid2)

##part b.
#using optim()

#NLL function for optiim
NLL_optim <- function(par, x, y) {
  a <- par[1]
  b <- par[2]
  sigma <- sqrt(mean((y - (a + b*x))^2))
  n <- length(y)
  nll <- (n/2)*log(2*pi*sigma^2) + sum((y - (a + b*x))^2)/(2*sigma^2)
  return(nll)
}

#using same guess as obj 2 (init)
mle_result <- optim(par = init, fn = NLL_optim, x = x, y = y)
mle_result$par

##part c.
#check for sensitivity to starting values

#storage
results2 <- matrix(NA, ncol=2, nrow=length(starts))
colnames(results2) <- c("intercept","slope") #name columns of storage matrix

#for loop using same starts as obj 2
for (i in seq_along(starts)) {
  fit <- optim(par = starts[[i]], fn = NLL_optim, x = x, y = y)
  results2[i,] <- fit$par
}
#all values withing 0.001 of each other regardless of start vals


# Objective 4 -------------------------------------------------------------
#All intercept estimates were within 0.01 of -1.37 from all 3 approaches
#All slope estimates were 1.34 using all 3 approaches


