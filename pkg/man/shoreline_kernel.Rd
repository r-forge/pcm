\name{shoreline_kernel}
\alias{shoreline_kernel}
\title{
Computes a kernel
}
\description{
This function generates a squared kernel which is used to convolve with a land cover and AVHRR ch1,ch4 data.
}
\usage{
shoreline_kernel(x)
}
\arguments{
  \item{x}{
A numeric value with the kernel size.
}
}
\value{
x times x squared matrix
}
\author{
Jan Musial
}
\seealso{
\code{\link{PCM}},\code{\link{lc_prepare}}
}
\examples{
shoreline_kernel(3)
}
\keyword{shoreline_kernel}
