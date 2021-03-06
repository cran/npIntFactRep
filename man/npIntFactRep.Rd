% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/npIntFactRep.R
\name{npIntFactRep}
\alias{npIntFactRep}
\title{Nonparametric Interaction Tests for Factorial Designs with Repeated Measures}
\usage{
npIntFactRep(dat)
}
\arguments{
\item{dat}{The name of the R data set in wide (one row per subject) format. This set should NOT have any missing values (NA). Missing data should be replaced (e.g., with 'mi', 'mice' or 'Amelia') before running. The data set MUST have subj as name (header) for the (numeric or character values of the) units of observations (subjects/blocks) factor and MUST have group as name for the (numeric or character values of the) between factor. Other variables must be excluded, because the number of levels of the repeated measures (within) factor is calculated based upon the remaining columns after eliminating the columns group and subj. The repeated measures factor is named rep.}
}
\value{
5 ANOVA summary tables
}
\description{
Nonparametric aligned rank tests for interactions in two-way factorial designs with repeated measures (e.g., Beasley & Zumbo, 2006; Higgings & Tashtoush, 1994). Five ANOVAs are reported: (1) PARAMETRIC on the original data, (2) on the interaction alignments as a CHECK, (3) on the aligned REGULAR, (4) FRIEDMAN, and (5) KOCH ranks. In the rank tests, only the resulting values for the interaction are relevant.
}
\details{
Returns five ANOVA tables using the 'ezANOVA' function (from the 'ez' package by Michael A. Lawrence), which is well-suited for factorial designs with repeated measures. The Higgins & Tashtoush formula for repeated measures (1994, p. 208) is used for the interaction alignments. The ANOVAs include F- and p-values, along with generalized eta squared (ges) effect size statistics and the sphericity test. Type II sum of squares are calculated, which are appropriate for the alignments check with unbalanced designs. The first ANOVA is the PARAMETRIC one on the original data. The second one is on the interaction alignments as a CHECK (the p-values for 'group' and for 'rep' should both be = 1). The three last ANOVAs are aligned rank tests: one on the aligned REGULAR, one on the FRIEDMAN, and one on the KOCH ranks. In these ranks tests, ONLY the resulting values for the INTERACTION group:rep are relevant
}
\examples{
\dontrun{
dat1 <- read.csv (file="c:/R/wide.csv", head=T)
npIntFactRep(dat1)
}
# example with character values for both group and subj variables (columns)
dat2 <- read.table(header = TRUE, text = "
subj group x1 x2 x3 x4
 p1 a 1 2 3 4
 p2 a 2 2 3 3
 p3 a 1 3 3 4
 p4 b 8 6 4 2
 p5 b 6 6 4 4
 p6 b 8 6 6 2
 p7 c 4 4 4 3
 p8 c 3 3 3 4
 p9 c 3 4 4 2
 ")
npIntFactRep(dat2)
# example with numeric values for the subj variable and charater values for the group variable
dat3 <- read.table(header = TRUE, text = "
subj group 1 2 3 4 5
1 a1 1 2 3 4 5
2 a1 2 2 3 3 3
3 a1 1 3 3 4 2
4 a2 8 6 4 2 1
5 a2 6 6 4 4 4
6 a2 8 6 6 2 3
")
npIntFactRep(dat3)
}
\author{
Jos Feys
}
\references{
Higgins, J.J., & Tashtoush, S. (1994). "An aligned rank transform test for interaction". Nonlinear World, 1, 201-211.
Beasley, T.M., & Zumbo, B.D. (2009). "Aligned rank tests for interactions in split-plot designs: Distributional assumptions and stochastic homogeneity". Journal of Modern Applied Statistical Methods, 8, 16-50.
}

