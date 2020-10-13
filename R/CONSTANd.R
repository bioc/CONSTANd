CONSTANd <- function(data, h=1e-6, maxIterations=50){
  
  m = nrow(data)
  n = ncol(data)
  
  f = vector(mode ="numeric", length=2*maxIterations)
  R = 1+vector(mode ="numeric", length=m)
  S = 1+vector(mode ="numeric", length=n)
  
  convergence = Inf
  eta = 0
  
  RM = rowMeans(data, na.rm = TRUE)
  while (h<convergence && eta<maxIterations) {
    
    tempR = 1/n * 1/RM
    data = data * tempR
    R = R * tempR
    
    CM = colMeans(data, na.rm = TRUE)
    f[2*eta+1] = m*0.5*sum(abs(CM - 1/n), na.rm = TRUE)
    
    
    tempS = 1/n * 1/CM
    data = t(t(data) * tempS)
    S = S * tempS
    
    RM = rowMeans(data, na.rm = TRUE)
    f[2*eta+2] = n*0.5*sum(abs(RM - 1/n), na.rm = TRUE)
    
  
    convergence = f[2*eta+2]
    eta = eta+1
    
  }
  
  result = list("data" = data, "f" = f[1:(2*eta)], "R" = R, "S" = S)
  return(result)
  
}