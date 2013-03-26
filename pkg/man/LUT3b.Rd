\name{LUT3b}
\alias{LUT3b}
\docType{data}
\title{
Look UP Table for AVHRR data with channel 3b configuration 
}
\description{
This data is supposed to be internaly used by the \code{\link{PCM}} function.
}
\usage{data(LUT3b)}
\format{
List containing the following entries:
  \describe{
    \item{\code{LUT}}{a multidimensional array with the PCM classification probability estimates for the 3b channel configuration}
    \item{\code{F1.vals}}{a numeric vector with threshold values apllied to 1st spectral feature}
    \item{\code{F2.vals}}{a numeric vector with threshold values apllied to 2nd spectral feature}
    \item{\code{F3.vals}}{a numeric vector with threshold values apllied to 3rd spectral feature}
    \item{\code{lc.vals}}{a numeric vector with land cover categories}
    \item{\code{dtn.vals}}{a numeric vector with daytime categories}
    \item{\code{azi.vals}}{a numeric vector with relative azimuth sectors}
    \item{\code{view.vals}}{a numeric vector with satellite zenith sectors}
    \item{\code{text.vals}}{a numeric vector with threshold values apllied to texture feature}
  }
}

\source{
Musial, J., Probabilistic approach to cloud and snow detection on AVHHR imagery, Atmos. Chem. Phys (in prep.)
}
\references{
\code{\link{PCM}}
}
\examples{
data(LUT3b)
}
\keyword{datasets}
