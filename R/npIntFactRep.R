#' @title Nonparametric Interaction Tests for Factorial Designs with Repeated Measures
#' @description Nonparametric aligned rank tests for interactions in two-way factorial designs with repeated measures (e.g., Beasley & Zumbo, 2006). There is a choice between 4 tests: (0) CHECK on the alignment for the interaction (1) REGULAR aligned ranks, (2) FRIEDMAN aligned ranks, and (3) KOCH ranks.
#' @usage npIntFactRep(dat,rank)
#' @references
#' Higgins, J.J., & Tashtoush, S. (1994). "An aligned rank transform test for interaction". Nonlinear World, 1, 201-211.
#' Beasley, T.M., & Zumbo, B.D. (2009). "Aligned rank tests for interactions in split-plot designs: Distributional assumptions and stochastic homogeneity". Journal of Modern Applied Statistical Methods, 8, 16-50.
#' @details Returns an ANOVA table using the 'ezANOVA' function (from the 'ez' package by Michael A. Lawrence), which is well-suited for factorial designs with repeated measures. It yields F- and p-values, along with generalized eta squared (ges) effect size statistics and the sphericity test. Type III sum of squares are calculated. ONLY the resulting values for the INTERACTION group:rep are relevant. There is a choice between 4 rank numbers: (0) Alignment of data (Higgins & Tashtoush, 1994) CHECK (p-values for group and for rep should be = 1.0), (1) REGULAR ranks of aligned data, (2) FRIEDMAN ranks of aligned data, (3) KOCH ranks.
#' @author Jos Feys
#' @param dat The R data set in wide (one row per subject or ID) format. This set should NOT have any missing values (NA). Missing data should be replaced (e.g., with 'mi', 'mice' or 'Amelia') before running. The data set MUST have subj as name (header) for the (numeric or character) variable for units of observations and MUST have group as name for the (numeric or character) variable for the between factor. The repeated measures (within) factor is named rep.
#' @param rank Numbers 0, 1, 2 or 3; resp. for Aligned data Check, Regular, Friedman, or Koch ranks.
#' @examples
#' \dontrun{
#' dat1 <- read.csv (file="c:/R/wide.csv", head=T)
#' #REGULAR
#' npIntFactRep(dat1,1)
#' }
#' dat2 <- read.table(header = TRUE, text = "
#' subj group x1 x2 x3 x4
#'  p1 a 1 2 3 4
#'  p2 a 2 2 3 3
#'  p3 a 1 3 3 4
#'  p4 b 8 6 4 2
#'  p5 b 6 6 4 4
#'  p6 b 8 6 6 2
#'  p7 c 4 4 4 3
#'  p8 c 3 3 3 4
#'  p9 c 3 4 4 2
#'  ")
#' #Aligned data check
#' npIntFactRep(dat2,0)
#' #FRIEDMAN
#' npIntFactRep(dat2,2)
#' dat3 <- read.table(header = TRUE, text = "
#' subj group 1 2 3 4 5
#' 1 1 1 2 3 4 5
#' 2 1 2 2 3 3 3
#' 3 1 1 3 3 4 2
#' 4 2 8 6 4 2 1
#' 5 2 6 6 4 4 4
#' 6 2 8 6 6 2 3
#' ")
#' #KOCH Ranks
#' npIntFactRep(dat3,3)
#' @import ez
#' @import plyr
#' @import stats
#' @return summary of Anova
#' @export

npIntFactRep <- function (dat,rank){
  dep <- NULL # Setting the variables to NULL first
  # make shure names are all lower case
  colnames(dat) <- tolower(colnames(dat))
  #sort by subj, group
  dat <- dat[order(dat$subj,dat$group),]

  group <- dat["group"]
  subj  <- dat["subj"]
  # remove columns "subj" & "group"
  # only  rep. meas. columns are left
  # one does not need the know the number
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

  if (rank == 0 ) {
    # aligned data
    ad_full <- cbind(subj, group, ad)
    # ad_full
    l <- ncol(ad_full)
    ad_long <- reshape(ad_full, direction="long",
                       varying=list(names(ad_full)[3:l]),
                       v.names="dep", timevar="rep",
                       idvar=c("subj","group"))
      ## Convert variables to factor
      ad_long <- within(ad_long, {
      group <- factor(group)
      rep <- factor(rep)
      subj <- factor(subj)
    })

    r0 <- ezANOVA(data=ad_long, dv = .(dep), wid=.(subj), between = .(group),
                  within=.(rep), type=3, detailed=FALSE)

    #  demo0.aov <- aov(dep ~ group * rep + Error(subj/rep), data = ad_long)
    #  summary(demo1.aov)

    cat(format("Aligned Data: p-values for group & for rep should be = 1.0)", width=80, justify="centre"))
    return(r0)
  }

  if ((rank != 0) && (rank != 2) && (rank != 3)) {
    # aligned ranks DEFAULT
    ar <- matrix(rank(ad, ties='average'), ncol=k)
    ar_full <- cbind(subj, group, ar)

    l <- ncol(ar_full)
    ar_long <- reshape(ar_full, direction="long",
                       varying=list(names(ar_full)[3:l]),
                       v.names="dep", timevar="rep",
                       idvar=c("subj","group"))

    ## Convert variables to factor
    ar_long <- within(ar_long, {
      group <- factor(group)
      rep <- factor(rep)
      subj <- factor(subj)
    })


    r1 <- ezANOVA(data=ar_long, dv = .(dep), wid=.(subj), between = .(group),
                      within=.(rep), type=3, detailed=FALSE)

  #  demo1.aov <- aov(dep ~ group * rep + Error(subj/rep), data = ar_long)
  #  summary(demo1.aov)

    cat(format("Regular Ranks", width=80, justify="centre"))
    return(r1)
  }

  if (rank == 2) {
    #FRIEDMAN
    for(i in 1:n){
      fr[i,] <- matrix(rank(ad[i,]), ncol=k)
    }
    #before adjustment
    #fr
    part1 = (k+1)/2
    part2 = ((k**2)-1)/12
    fr <- (fr-part1)/part2
    #after adjustment
    #fr
    fr_full <- cbind(subj, group, fr)

    l <- ncol(fr_full)
    fr_long <- reshape(fr_full, direction="long",
                       varying=list(names(fr_full)[3:l]),
                       v.names="dep", timevar="rep",
                       idvar=c("subj","group"))

    ## Convert variables to factor
    fr_long <- within(fr_long, {
      group <- factor(group)
      rep <- factor(rep)
      subj <- factor(subj)
    })

    r2 <- ezANOVA(data=fr_long, dv = .(dep), wid=.(subj), between = .(group),
                      within=.(rep), type=3, detailed=FALSE)

    # demo2.aov <- aov(dep ~ group * rep + Error(subj/rep), data = fr_long)
    # summary(demo2.aov)

    cat(format("Friedman Ranks", width=80, justify="centre"))

    return(r2)
  }

  if (rank == 3) {
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
    #voor adjustment
    #som
    deel1 = (n+1)/2
    deel2 = (k-1)*(n+1)
    #na adjustment
    q <- (som-deel1)/deel2
    #q
    q_full <- cbind(subj, group, q)

    l <- ncol(q_full)
    q_long <- reshape(q_full, direction="long",
                      varying=list(names(q_full)[3:l]),
                      v.names="dep", timevar="rep",
                      idvar=c("subj","group"))
    ## Convert variables to factor
    q_long <- within(q_long, {
      group <- factor(group)
      rep <- factor(rep)
      subj <- factor(subj)
    })

    r3 <- ezANOVA(data=q_long, dv = .(dep), wid=.(subj), between = .(group),
                      within=.(rep), type=3, detailed=FALSE)

    # demo3.aov <- aov(dep ~ group * rep + Error(subj/rep), data = q_long)
    # summary(demo3.aov)

    cat(format("Koch Ranks", width=80, justify="centre"))
    return(r3)
  }
}

