\name{summary.colocal}
\alias{summary.colocal}
\title{
Summarizing the colocalization results
}
\description{
This function is used to summarize the results returned from the \code{\link{nsinc.d}} or \code{\link{nsinc.z}} function.
}
\usage{
\method{summary}{colocal}(object, \dots)
}
\arguments{
  \item{object}{an object of class \code{"colocal"}, containing the results returned from \code{\link{nsinc.d}} or \code{\link{nsinc.z}}.}
  \item{\dots}{further arguments to be passed from or to other methods.}
}
\details{
The results successfully returned from \code{\link{nsinc.d}} or \code{\link{nsinc.z}} functions give a list of length 12 encompassing the summarized information of the calculation of colocalization index and the detailed quantities used to calculate the index.

First the summary prints the strategic parameters for calculation of colocalization index, such as, the method, the strata, the edge effect, dimension, study region, membership levels and number of the observed signals in each channel.

The summary also prints the summarized information of proximity sizes, i.e., r, including the r model, the r range, the length of r series, etc.

Then the separate index results are listed at each r. The average colocalization index for the whole image is given at last.
}
\value{
\code{summary.colocal} does not return values. It only prints summarized results returned from \code{\link{nsinc.d}} or \code{\link{nsinc.z}} functions.
}
\author{
Xueyan Liu, Jiahui Xu, Cheng Cheng, Hui Zhang.
}
\references{
Liu, X., Xu, J., Guy C., Romero E., Green D., Cheng, C., Zhang, H. (2019). \emph{Unbiased and Robust Analysis of Co-localization in Super-resolution Images}. Manuscript submitted for publication.
}
\seealso{
\code{\link{nsinc.d}}, \code{\link{nsinc.z}}, \code{\link{summary}}
}
\examples{
## a simulated 3D example data.
data("twolines")

twolines.results <- nsinc.d(data = twolines, membership = "membership",
                            dim = 3, r.model = "r.med")

summary(twolines.results)
}

