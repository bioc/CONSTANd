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
# (c) Dirk Valkenborg & Jef Hooyberghs, 2014
# Updated and ported to R by Joris Van Houtven, 2020

CONSTANd <- function(data, precision=1e-5, maxIterations=50, target=1){
    # Return the normalized version of the input data (matrix) as an ndarray, as
    # well as the convergence trail (residual error after each iteration) and the
    # row and column multipliers R and S.
    if (any(data[!is.na(data)]<0)) {
        stop("Negative values detected in quantification matrix.
             Are you using log-transformed ratio's?
             If so, use intensities instead.")
    }
    if (any(data[!is.na(data)]==0)) {
        warning("Zeros in quantification matrix detected; replacing with NA.")
        data[data==0] <- NA
    }

    Nrows = nrow(data)
    Ncols = ncol(data)

    convergence_trail = rep(NA, 2*maxIterations)
    convergence = Inf
    R = rep(1, Nrows)
    S = rep(1, Ncols)
    RM = rowMeans(data, na.rm = TRUE)
    i = 0
    while (precision<convergence && i<maxIterations) {
        # fit the rows
        tempR = target * 1/RM
        data = data * tempR
        R = R * tempR

        CM = colMeans(data, na.rm = TRUE)
        convergence_trail[2*i+1] = Nrows * sum(abs(CM - target), na.rm = TRUE) / 2

        # fit the columns
        tempS = target * 1/CM
        data = t(t(data) * tempS)
        S = S * tempS

        RM = rowMeans(data, na.rm = TRUE)
        convergence_trail[2*i+2] = Ncols * sum(abs(RM - target), na.rm = TRUE) / 2

        convergence = convergence_trail[2*i+2]
        i = i+1
    }
    if (i == maxIterations) {
        warning(sprintf("\nMax number of CONSTANd iterations (%i) reached. Attained precision: %.5f", maxIterations, convergence))
    }
    result = list("normalized_data" = data, "convergence_trail" = convergence_trail[seq(2*i)], "R" = R, "S" = S)
    return(result)
}
