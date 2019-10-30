\name{nsinc.z}
\alias{nsinc.z}
\title{Colocalization index of z-type}
\description{
\code{nsinc.z} is used to calculate the Pearson's correlation coefficient of the signal proportions of two channels with a z-score normalization based on complete spatial randomness (CSR) in a specified proximity of all signals or all signals of interested type as the colocalization index for a whole image. If a range of proximity sizes are concerned, the \code{nsinc.z} will take the average of the index values over the range. In the case of multiple-species data, the average of index values of all pairs at each proximity size is taken as the index for the image at that size of neighborhood.
}
\usage{
nsinc.z(data, membership, dim = 2, r.min = NULL,
        r.max = NULL, r.count = NULL, r.adjust = NULL,
        box = NULL, edge.effect = TRUE, strata = FALSE,
        base.member = NULL, r.model = "full", ...)
}
\arguments{
  \item{data}{a data frame (or object coercible by as.data.frame to a data frame) containing at least the columns \code{membership} and \code{x (xc, X or Xc)}, \code{y (yc, Y or Yc)} if \code{dim = 2} and \code{x (xc, X or Xc)}, \code{y (yc, Y or Yc)}, \code{z (zc, Z or Zc)} if \code{dim = 3}.}
  \item{membership}{a string describing the column name in the \code{data} representing the membership of data points. There should be no less than 2 levels in the membership.}
  \item{dim}{an integer either \code{= 2} or \code{= 3}. If \code{dim = 2}, the data are treated as two-dimensional; if \code{dim = 3}, the data are treated as three-dimensional.}
  \item{r.min}{the minimum proximity size that the user identifies as colocalization of signals. It should be numeric. If \code{r.model = "full"}, the function will automatically choose the smallest inter-point distance as the \code{r.min}; if \code{r.model = "r.med"}, the function will use the median inter-point distance for both \code{r.min} and \code{r.max}; if \code{r.model = "other"}, the user must specify \code{r.min}, which should be no larger than \code{r.max}.}
  \item{r.max}{the maximum proximity size that the user identifies as colocalization of signals. It should be numeric. If \code{r.model = "full"}, the function will automatically choose \bold{half of} the largest inter-point distance as the \code{r.max}; if \code{r.model = "r.med"}, the function will use the median inter-point distance for both \code{r.min} and \code{r.max}; if \code{r.model = "other"}, the user must specify \code{r.max}, which should be between the smallest and the largest inter-point distances and no smaller than \code{r.min}.}
  \item{r.count}{the total count of the series of proximity sizes between \code{r.min} and \code{r.max}. If \code{r.max = r.min} or \code{r.adjust = (r.max - r.min)/2}, then \code{r.count = 1}, otherwise \code{r.count = 30} by default or is specified by the user.}
  \item{r.adjust}{a very small adjustment for \code{r.min} and \code{r.max} to get the series of proximity sizes between \code{r.min + r.adjust} and \code{r.max - r.adjust} to avoid zero standard deviation of normalized proportions of signals at extremely small and large r's. The values of \code{r.adjust} depends on the choice of \code{r.model} and values of \code{r.min} and \code{r.max}. For most scenarios, it is suggested to use \code{r.adjust = NULL} and let the function choose the default value for \code{r.adjust}. In general, by default either \code{r.adjust = 0} or \code{r.adjust = (r.max - r.min)/(r.count + 1)}; otherwise it is a positive number specified by the user satisfying \code{r.adjust} \eqn{\le} \code{(r.max - r.min)/2}.}
  \item{box}{a one-row data frame describing the study region which must contain columns \code{xmin, xmax, ymin, ymax} if \code{dim = 2} and additionally \code{zmin, zmax} if \code{dim = 3}. If \code{box = NULL}, the function will detect the smallest box containing all data points and add a buffer edge in each dimension which is equal to the median of nearest neighbor distances in that dimension. If \code{box} is specified by the user, only the data enclosed in the specified box will be considered in the analysis and signals outside the \code{box} will be ignored.}
  \item{edge.effect}{a logical value showing whether the edge effect should be corrected. By default it should be corrected otherwise the results are not accurate.}
  \item{strata}{a logical value showing whether the user wants to consider single-direction or bi-direction colocalization. By default \code{strata = FALSE} is for bi-direction colocalization. In this case, all proximity regions around all signals are considered. If \code{strata = TRUE}, then \code{base.member} must be specified or the first membership that R detects in the membership column will be used by default and only the circular regions around signals in the base membership are considered. Then, colocalization will be single-direction in this case.}
  \item{base.member}{one level of the memberships that is designated as the base. It works only when \code{strata = TRUE}. If \code{strata = TRUE} and no \code{base.member} is specified by the user, the first membership that R detects in the membership column will be used by default for \code{base.member}.}
  \item{r.model}{equals either \code{"full"}, \code{"r.med"} or \code{"other"}. The \code{r.model} will be used to choose the proximity size ranges that the user defines for colocalization. \code{"full"} or \code{"r.med"} can be used if the user has no specific sense of proximity sizes for colocalization. In \code{"full"} model, the colocalization proximity sizes will range from the smallest inter-point distance to \bold{half of} the largest inter-point distances; in \code{"r.med"} model, the fixed proximity size is the median of inter-point distances; in \code{"other"} model, the user can define their research driven proximity sizes by specifying \code{r.min} and \code{r.max}.}
  \item{...}{Parameters passed to \code{cor}. The user could choose methods other than \code{Pearson} for calculating correlation.}
}
\details{
The function calculates the proportion of two types of signals normalized by a z-score under CSR in a specified r neighborhood with edge effect corrected of all signals or all base signals if \code{strata = TRUE}is specified, then obtains the Pearson correlation coefficients of each pair of channels and average them among all pairs at each r in the r series between \code{r.min} to \code{r.max}. In the case of multiple-species data, the average of index values of all pairs at each proximity size is taken as the index for the image at that size of neighborhood. The index for the whole image is named as \ifelse{html}{\out{NSInC<sup>z</sup>}}{\eqn{NSInC^z}} or NSInC of type z. The index will be close to 1 if signals are colocalized, 0 if random and -1 if dispersed. The function can deal with 2D or 3D data.

If the users have their specific proximity size, then they are encouraged to specify \code{r.model = "other"}, and values of \code{r.min} and \code{r.max}.

The difference from \code{\link{nsinc.d}} is the normalization of the signal proportions. The z-type normalization has no heterogeneity under CSR caused by the edge effects related to the locations of signals. In many cases, \code{nsinc.d} and \code{nsinc.z} can give similar results. However, if the user's proximity of interest is larger than half of the largest inter-point distances, then \code{nsinc.d} is suggested.
}
\value{
\code{nsinc.z} returns colocalization index values at each separate proximity size r, and the average colocalization index across all r's, the data that the colocalization index is calculated from, the study region, i.e., the carrying box, the original and normalized proportions of each type of signals in an r neighborhodd of all (base) signals, the r series, and some summary information:

  \item{method}{"nsinc.z"}
  \item{input.data.summary}{a list containing the number of membership levels and the signal counts in each channel or membership of the input data.}
  \item{post.data.summary}{a list showing the number of membership levels and signal counts in each channel of the data after removal of signals located outside the specified box  by the user. If there is no signals excluded, then \code{post.data.summary} presents the same results as \code{input.data.summary}.}
  \item{r.summary}{a data frame listing the \code{r.min}, \code{r.max}, \code{r.count}, \code{r.adjust} used in the calculation and the \code{r.model} specified by the user or the default. \code{r.summary} also gives the r range for the default full model, i.e., the minimum and half of the maximum of the inter-point distance of all signals, and the median value in addition.}
  \item{strata}{a list showing the default setting of strata or the specified strata by the user. It also presents the base membership used in the function if \code{strata} is TRUE.}
  \item{edge.effect}{a data frame containing a logical value indicating whether edge effect is corrected or not.}
  \item{index.all}{a data frame showing the colocalization index of z-type at each r.}
  \item{index}{the averaged colocalization index of z-type across all r's.}
  \item{post.data}{a data frame representing the data after removal of signals located outside the specified box by the user. If there is no signal excluded, then \code{post.data} presents the same observations as data.}
  \item{study.region}{the carrying box with the size of buffer width in each dimension.}
  \item{P.all}{the data frame showing all original and normalized proportions of each type of signals in an r-neighborhood around every (base) signal. Rows are (base) signals and columns are all memberships and r's.}
  \item{r}{the r series for which the colocalization indices are calculated.}

}
\author{
Xueyan Liu, Jiahui Xu, Cheng Cheng, Hui Zhang.
}
\references{
Liu, X., Xu, J., Guy C., Romero E., Green D., Cheng, C., Zhang, H. (2019). \emph{Unbiased and Robust Analysis of Co-localization in Super-resolution Images}. Manuscript submitted for publication.
}
\examples{
## a simulated 2D example data.
set.seed(1234)
x <- runif(300, min = -1, max = 1)
y <- runif(300, min = -1, max = 1)
red <- data.frame(x,y, color = "red")
x <- runif(50, min = -1, max = 1)
y <- runif(50, min = -1, max = 1)
green <- data.frame(x,y, color = "green")

mydata <- rbind(red,green)
plot(mydata$x,mydata$y,col = mydata$color)

mydata.results <- nsinc.z(data = mydata, membership = "color", dim = 2,
                  r.model = "other", r.min = 0.01, r.max = 0.5, r.count = 5, r.adjust = 0)

mydata.results$index.all
mydata.results$index


## a simulated 3D example data.
data("twolines")

\donttest{
library("rgl")
plot3d(twolines[,c("x","y","z")], type='s', size=0.7, col = twolines$membership)
aspect3d("iso")

twolines.results <- nsinc.z(data = twolines, membership = "membership",
                            dim = 3, r.model = "full")

twolines.results$index
}
}
