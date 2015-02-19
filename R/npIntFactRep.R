#' @title Nonparametric Interaction Tests for Factorial Designs with Repeated Measures
#' @description Nonparametric aligned rank tests for interactions in two-way factorial designs with repeated measures.
#' @usage npIntFactRep(dat,rank)
#' @references
#' Higgins, J.J., & Tashtoush, S. (1994). "An aligned rank transform test for interaction". Nonlinear World, 1, 201-211.
#' Beasley, T.M., & Zumbo, B.D. (2009). "Aligned rank tests for interactions in split-plot designs: Distributional assumptions and stochastic homogeneity". Journal of Modern Applied Statistical Methods, 8, 16-50.
#' @details Returns an ANOVA table using the 'ezANOVA' function (from the 'ez' package by Michael A. Lawrence), which is well-suited for factorial design experiments with repeated measures. It yields ANOVA F- and p-values, generalized effect sizes and assumption checks. Type III sums of squares are calculated. ONLY the resulting values for the 'group:time' (between x within) INTERACTION are relevant. There is a choice between three types of aligned ranks for the tests: (1) REGULAR, (2) FRIEDMAN, and (3) KOCH ranks.
#' @author Jos Feys
#' @param dat  The R data set in wide (one row per subject) format. This should NOT have any missing values (NA). Missing data should be replaced before running. We suggest using the 'mi', 'mice' or 'Amelia' package for multiple imputations. The data set MUST have subj as name (header) of the (numeric or character) variable for units of observations and MUST have group as name of the (numeric or character) variable for between factor.
#' @param rank The numbers 1, 2 or 3; resp. for Regular, Friedman, or Koch ranks.
#' @examples
#' \dontrun{
#' dat1<-read.csv(file="c:/R/wide.csv",head=T)
#' #REGULAR
#' npIntFactRep(dat1,1)
#' }
#' dats<-read.table(header=TRUE,text="
#' subj group x1 x2 x3 x4
#' p1 a 1 2 3 4
#' p2 a 2 2 3 3
#' p3 a 1 3 3 4
#' p4 b 8 6 4 2
#' p5 b 6 6 4 4
#' p6 b 8 6 6 2
#' p7 c 4 4 4 3
#' p8 c 3 3 3 4
#' p9 c 3 4 4 2
#' ")
#' #FRIEDMAN
#' npIntFactRep(dats,2)
#' @import ez
#' @import plyr
#' @return Anova F-tests and p-values for the 'group:time' interaction
#' @export

npIntFactRep <- function (dat,rank){
  dep <- NULL # Setting the dependent variable to NULL first
  # make shure names are all lower case
  colnames(dat) <- tolower(colnames(dat))
  # sort by subj, group
  dat <- dat[order(dat$subj,dat$group),]
  # subset 2 column vectors
  group <- dat["group"]
  subj  <- dat["subj"]
  # remove "subj" & "group"
  # only repeated measures are selected
  # one does not need to know/specify the number of measures
  colsdel <- c("group","subj")
  y <- dat[,!(names (dat) %in% colsdel)]
  # number of columns (rep.meas.) & observations (subj's)
  k <- ncol(y)
  n <- nrow(y)
  # initializations (probably not necessary, but ...)
  ad <- matrix(0,n,k)
  ar <- matrix(0,n,k)
  fr <- matrix(0,n,k)
  q  <- matrix(0,n,k)
  # Higgins & Tashtoush formula (1994, p. 108)
  rmean <- colMeans(y)
  rmean <- matrix(rep(rmean, each=n), ncol=k)
  pmean <- rowMeans(y)
  gmean <- mean (pmean)
  #aligned data
  ad <- (y - rmean - pmean) + gmean
  
  if ((rank != 2  ) & (rank != 3)) {
    # default: aligned REGULAR ranks
    ar <- matrix(rank(ad, ties='average'), ncol=k)
    ar_full <- cbind(subj, group, ar)
    l <- ncol(ar_full)
    ar_long <- reshape(ar_full, direction="long", 
                       varying=list(names(ar_full)[3:l]), 
                       v.names="dep", 
                       idvar=c("subj","group"))
    ## Convert variables to factor
    ar_long <- within(ar_long, {
      group <- factor(group)
      time <- factor(time)
      subj <- factor(subj)
    })
  
    r1 <- ezANOVA(data=ar_long, dv = .(dep), wid=.(subj), between = .(group),
                      within=.(time), type=3, detailed=FALSE)
    cat(format("Regular Ranks", width=80, justify="centre"))
    return(r1)
  }
  
  else if (rank == 2) {
    #FRIEDMAN
    for(i in 1:n){
      fr[i,] <- matrix(rank(ad[i,]), ncol=k)
    }
    # voor adjustment
    # fr
    part1 = (k+1)/2
    part2 = ((k**2)-1)/12
    fr <- (fr-part1)/part2
    # na adjustment
    # fr
    fr_full <- cbind(subj, group, fr)
    l <- ncol(fr_full)
    fr_long <- reshape(fr_full, direction="long", 
                       varying=list(names(fr_full)[3:l]), 
                       v.names="dep", 
                       idvar=c("subj","group"))

    ## Convert variables to factor
    fr_long <- within(fr_long, {
      group <- factor(group)
      time <- factor(time)
      subj <- factor(subj)
    })
    
    r2 <- ezANOVA(data=fr_long, dv = .(dep), wid=.(subj), between = .(group),
                      within=.(time), type=3, detailed=FALSE)
    cat(format("Friedman Ranks", width=80, justify="centre"))
     return(r2)
  }
  
  else if (rank == 3) {
    #KOCH
    dif <- array(0, dim=c(n,k,k))
    rdif <- array(0, dim=c(n,k,k))
    for(j in 1:k){
      for(i in 1:k){
        dif[,i,j] <- y[,j] - y[,i]
        rdif[,i,j] <- matrix(rank(dif[,i,j]), nrow=n)
      }
    }
    som <- matrix(0,n,k)
    for(i in 1:k){
      som[,i] <- matrix(rowSums (rdif[,,i]))
    }
    # voor adjustment
    # som
    deel1 = (n+1)/2
    deel2 = (k-1)*(n+1)
    # na adjustment
    q <- (som-deel1)/deel2
    # q
    q_full <- cbind(subj, group, q)
    l <- ncol(q_full)
    q_long <- reshape(q_full, direction="long", 
                      varying=list(names(q_full)[3:l]), 
                      v.names="dep", 
                      idvar=c("subj","group"))
    
    ## Convert variables to factor
    q_long <- within(q_long, {
      group <- factor(group)
      time <- factor(time)
      subj <- factor(subj)
    })
    
    r3 <- ezANOVA(data=q_long, dv = .(dep), wid=.(subj), between = .(group),
                      within=.(time), type=3, detailed=FALSE)
    cat(format("Koch Ranks", width=80, justify="centre"))
    return(r3)
  }
}


