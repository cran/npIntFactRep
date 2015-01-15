#' @title Nonparametric Aligned Ranks Function
#' @description Nonparametric test for interactions in two-way factorial designs
#' with repeated measures - "summary" of full model "aov" on (for interaction) aligned ranks
#' @usage npIntFactRep(dat,rank)
#' @references
#' Higgins, J.J., & Tashtoush, S. (1994). "An aligned rank transform test for interaction". 
#' Nonlinear World, 1, 201-211.
#' Beasley, T.M., & Zumbo, B.D. (2009). "Aligned rank tests for interactions 
#' in split-plot designs: Distributional assumptions and stochastic homogeneity".  
#' Journal of Modern Applied Statistical Methods, 8, 16-50.
#' @details desired ranks numbers:
#' 1 Regular ranks of aligned data
#' 2 Friedman ranks of aligned data
#' 3 Koch ranks
#' Returns ANOVA F test and probabilities. 
#' Only the resulting values for 'group x time' ('between x within') are relevant
#' @author Jos Feys
#' @param dat  data set in wide (one row per subject) format, see examples 
#' - must have subj as name (header) of numeric var for units of obs id
#' - must have group as name of numeric var for between factor
#' @param rank numbers 1, 2 or 3; resp. for Regular, Friedman, or Koch ranks
#' @examples
#' \dontrun{
#' dat1 <- read.csv (file="c:/R/wide.csv", head=T)
#' #REGULAR
#' npIntFactRep(dat1,1)
#' }
#' dat2 <- read.table(header = TRUE, text = "
#' subj group m1 m2 m3 m4
#' 1 1 1 2 3 4
#' 2 1 2 2 3 3
#' 3 1 1 3 3 4
#' 4 2 8 6 4 2
#' 5 2 6 6 4 4
#' 6 2 8 6 6 2
#' ")
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
#' @return summary of (aov) Anova F-tests and p-value
#' @export  

npIntFactRep <- function (dat,rank){

# afzonderlijke matrices/vectoren: "subj", "group", "rep"
  
#sort by subj
dat <- dat[order(dat$subj),]

group <- dat["group"]
subj  <- dat["subj"]
#verwijder kolommen "subj" & "group"
#alleen  repeated measures blijven over
#men hoeft aantal niet te kennen
colsdel <- c("group","subj")
y <- dat[,!(names (dat) %in% colsdel)]

k <- ncol(y)
n <- nrow(y)
q <- matrix(0,n,k)
fr <- matrix(0,n,k)

rmean <- colMeans(y)
rmean <- matrix(rep(rmean, each=n), ncol=k)
pmean <- rowMeans(y)
gmean <- mean (pmean)
#aligned data
ad <- (y - rmean - pmean) + gmean

if (rank ==1 ) {
# aligned ranks
ar <- matrix(rank(ad, ties='average'), ncol=k)
ar_full <- cbind(subj, group, ar)

l <- ncol(ar_full)
ar_long <- reshape(ar_full, direction="long", 
                   varying=list(names(ar_full)[3:l]), 
                   v.names="value", 
                   new.row.names = 1:10000,
                   idvar=c("subj","group"))
ar_long
aov.out = aov(value ~ group*time + Error(subj/group), 
              data=ar_long)
cat(format("Regular Ranks", width=80, justify="centre"))
r1 <- summary(aov.out)
return(r1)
}

if (rank == 2) {
#FRIEDMAN
for(i in 1:n){
fr[i,] <- matrix(rank(ad[i,]), ncol=k)
}
#voor adjustment
#fr
part1 = (k+1)/2
part2 = ((k**2)-1)/12
fr <- (fr-part1)/part2
#na adjustment
#fr
fr_full <- cbind(subj, group, fr)

l <- ncol(fr_full)
fr_long <- reshape(fr_full, direction="long", 
                   varying=list(names(fr_full)[3:l]), 
                   v.names="value",
                   new.row.names = 1:10000,
                   idvar=c("subj","group"))

aov.out = aov(value ~ group*time + Error(subj/group), 
              data=fr_long)

cat(format("Friedman Ranks", width=80, justify="centre"))
r2 <- summary(aov.out)
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
q_long  <- reshape(q_full, direction="long", 
                  varying=list(names(q_full)[3:l]), 
                  v.names="value", 
                  new.row.names = 1:10000,
                  idvar=c("subj","group"))

aov.out = aov(value ~ group*time + Error(subj/group), 
              data=q_long)
cat(format("Koch Ranks", width=80, justify="centre"))
r3 <- summary(aov.out)
return(r3)
}
}


