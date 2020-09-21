\name{plot.colocal}
\alias{plot.colocal}
\title{Making scatter plots for signal proportions before and after d-type or z-type normalization}
\description{
This function is used to make scatter plots for signal proportions based on the results returned from the \code{nsinc.d} or \code{nsinc.z} function.
}
\usage{
\method{plot}{colocal}(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{"colocal"}, containing the results returned from \code{\link{nsinc.d}} or \code{\link{nsinc.z}}.}
  \item{\dots}{further arguments to be passed from or to other methods.}
}
\details{
The function currently works for results from \code{nsinc.d} or \code{nsinc.z} with bi-direction colocalization for dual-color images. At each proximity size r, the function makes two panels of scatter plots for signal proportions at all signals before and after d-type or z-type normalizations. Each signal in the original image contributes a point in the scatter plots whose x coordinate is the proportion of signals in one channel and y coordinate is the other channel. The scatter plots have the same color codes as signals in the original image. If the returned results contain colocalization results at multiple r's, then the scatter plots at each r are generated.
}
\value{
\code{plot.colocal} returns a list of plots which summarize the results returned from the \code{nsinc.d} or \code{nsinc.z} function.
}
\author{
Xueyan Liu, Jiahui Xu, Cheng Cheng, Hui Zhang.
}
\references{
Liu, X., Xu, J., Guy C., Romero E., Green D., Cheng, C., Zhang, H. (2019). \emph{Unbiased and Robust Analysis of Co-localization in Super-resolution Images}. Manuscript submitted for publication.
}
\examples{
## a simulated 3D example data.
data("twolines")

twolines.results <- nsinc.d(data = twolines, membership = "membership",
                            dim = 3, r.model = "r.med")

##plot(twolines.results)
}

