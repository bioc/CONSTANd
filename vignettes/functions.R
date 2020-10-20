# Extra functions to use in the vignette.
# They automate all the boring stuff like data cleaning and rearranging.

clean_organs <- function(organs) {
    # remove high isolation interference cases
    organs <- lapply(organs, function(x) x[x$Isolation.Interference.... < 30,])
    # remove duplicate PSMs (peptide-spectrum matches): select only Mascot results
    organs <- lapply(organs, function(x) x[x$Identifying.Node.Type=='Mascot',])
    # identify quantification columns
    quanCols <- colnames(BR1)[(ncol(BR1)-7):ncol(BR1)]
    # summarize to peptide level (median of each channel's quantification value for each peptide)
    organs <- lapply(organs, function(x) as.data.frame(x %>% group_by(Sequence) %>% summarise_at(.vars = quanCols, .funs = median, na.rm=TRUE) %>% drop_na(., Sequence)))
    ## keep only quantification and peptide info; set rownames to peptide sequences
    for (i in seq_along(organs)) {
      rownames(organs[[i]]) <- organs[[i]]$Sequence
      organs[[i]]$Sequence <- NULL}
    # remove X from colnames
    organs <- lapply(organs, setNames, nm = unname(unlist(lapply(strsplit(quanCols, 'X'), function(x) x[[2]]))))
    quanCols
    # remove PSMs with less than 4 quantification values
    ## set NA in quantification columns to zero to calc median
    for (i in seq_along(organs)) { organs[[i]][is.na(organs[[i]])] <- 0 }
    ## remove rows with too many missing values by using median
    organs <- lapply(organs, function(x) x[apply(x, 1, median, na.rm = FALSE) > 0,])
    ## set NA values again for use with CONSTANd
    for (i in seq_along(organs)) { organs[[i]][organs[[i]]==0] <- NA }
    return(organs)
}

rearrange_organs_design <- function(study.design) {
    rownames(study.design) <- study.design$V1
    study.design$V1 <- NULL
    channels <- unname(unlist(lapply(sapply(as.character(study.design[1,]), strsplit, ':'), function(x) x[2])))
    tbl <- tibble(run=rep(rownames(study.design), each=8), channel=rep(channels,4), condition=rep("", 4*8))
    for (i in seq(nrow(study.design))) {
        row <- as.character(study.design[i,])
        ss.l <- sapply(row, strsplit, ':')
        ss.conditions <- unname(unlist(lapply(ss.l, function(x) x[1])))
        ss.channels <- unname(unlist(lapply(ss.l, function(x) x[2])))
        correct.order <- match(ss.channels, tbl$channel)
        tbl[tbl$run==paste0('BR',i),]$condition[correct.order] <- ss.conditions
    }
    tbl <- tbl %>% unite(uniqueChannel, c(run,channel), remove = FALSE)
    return(column_to_rownames(tbl, 'uniqueChannel'))
}

extract_legend <- function(my_ggp) {
    step1 <- ggplot_gtable(ggplot_build(my_ggp))
    step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
    step3 <- step1$grobs[[step2]]
    return(step3)
}
