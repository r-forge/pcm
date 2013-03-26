\name{PCM-package}
\alias{PCM-package}
\docType{package}
\title{
Probabilistic cloud and snow detection on AVHRR imagery
}
\description{
This package consists of set of procedures used to: 
pre-process input ancillary data (land cover, DEM, surface skin temperature);
generate probabilistic cloud & snow mask for AVHRR data; 
post-classify the results into a binary form; 
compute cloud shadow mask
}
\details{
\tabular{ll}{
Package: \tab PCM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-03-19\cr
License: \tab GNU >= 2.0\cr
}
The data pre-processing is done in 3 steps. First the land cover data are remapped and interpolated (nearest neighbor) to desired projection and geographic extent using the \code{\link{lc_prepare}} function. Further (optionally) the Digital Elevation Model is remapped and interpolated (bilinear) to desired projection and geographic extent together with a temperature correction file which consists of a difference between high resolution DEM and DEM provided by the NWP model (geopotential at surface) multiply by the constant lapse rate of 0.6K/100m (\code{\link{dem_prepare}} function). Next (optionally) the NWP Skin Surface Temperature data is remapped and interpolated (bilinear) to desired projection and geographic extent and the temperature correction data are added (\code{link{skt_prepare}} function). Finally the cloud & snow detection is performed with optional DEM & SKT input and the classification probability is stored in the GeoTif file (\code{\link{PCM}} function). Optionally the post-classification can be performed which transforms the probability estimates to binary form and crude cloud shadow determination is added.
}
\author{
Jan Musial

Maintainer: Jan Musial <jmusial84@gmail.com>
}
\references{
Musial, J., Probabilistic approach to cloud and snow detection on AVHHR imagery, Atmos. Chem. Phys (in prep.)
}
\keyword{PCM}
