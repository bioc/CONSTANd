# CONSTANd
# Normalizes the data matrix <data> by raking the Nrows by Ncols matrix such that
# the row means and column means equal Ncols and Nrows, respectively. Missing 
# information needs to be presented as nan values and not as zero values because 
# CONSTANd employs the Matlab/NumPy functionality 'nanmean' that is able to ignore
# nan-values when calculating the mean. The variable <maxIterations> is an
# integer value that denotes the number of raking cycles. The variable <precision>
# defines the stopping criteria based on the L1-norm as defined by
# Friedrich Pukelsheim, Bruno Simeone in "On the Iterative Proportional
# Fitting Procedure: Structure of Accumulation Points and L1-Error
# Analysis"
# 
# © Dirk Valkenborg & Jef Hooyberghs, 2014
# Ported to R by Joris Van Houtven, 2020

CONSTANd <- function(data, precision=1e-5, maxIterations=50){
  # Return the normalized version of the input data (matrix) as an ndarray, as 
  # well as the convergence trail (residual error after each iteration) and the 
  # row and column multipliers R and S.
  
  Nrows = nrow(data)
  Ncols = ncol(data)
  TARGET = 1
  
  convergence_trail = rep(NA, 2*maxIterations)
  convergence = Inf
  R = rep(1, Nrows)
  S = rep(1, Ncols)
  RM = rowMeans(data, na.rm = TRUE)
  i = 0
  while (precision<convergence && i<maxIterations) {
    # fit the rows
    tempR = TARGET * 1/RM
    data = data * tempR
    R = R * tempR
    
    CM = colMeans(data, na.rm = TRUE)
    convergence_trail[2*i+1] = Nrows * sum(abs(CM - TARGET), na.rm = TRUE) / 2
    
    # fit the columns
    tempS = TARGET * 1/CM
    data = t(t(data) * tempS)
    S = S * tempS
    
    RM = rowMeans(data, na.rm = TRUE)
    convergence_trail[2*i+2] = Ncols * sum(abs(RM - TARGET), na.rm = TRUE) / 2
  
    convergence = convergence_trail[2*i+2]
    i = i+1
  }
  if (i == maxIterations) {
    warning(cat(sprintf("Max number of CONSTANd iterations (%i) reached. Attained precision: %f.5", maxIterations, convergence)))
  }
  result = list("normalized_data" = data, "convergence_trail" = convergence_trail[1:(2*i)], "R" = R, "S" = S)
  return(result)
}
