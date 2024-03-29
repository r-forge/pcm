\name{skt_prepare}
\alias{skt_prepare}
\title{
Prepare NWP Skin surface Temperature (SKT) data
}
\description{
This function remaps and interpolates (bilinear) to desired projection and geographic extent Skin surface Temperature data originating from the NWP model. Moreover the temperature correction provided by the \code{\link{dem_prepare}} function is added to the final results in order to improve temperature estimates over mountainous regions. 
}
\usage{
skt_prepare(NWP.skt.file = NWP.skt.file, out.dir = out.dir, out.proj = out.proj, out.extent = out.extent, out.dim = out.dim, correction.file = correction.file, NWP.skt.proj = NWP.skt.proj, NWP.skt.extent = NWP.skt.extent, NWP.skt.half.cell = NWP.skt.half.cell, NWP.skt.scale = NWP.skt.scale, NWP.skt.offset = NWP.skt.offset, NWP.skt.band = NWP.skt.band)
}
\arguments{
  \item{NWP.skt.file}{
Required! String with an absolute path to Numerical Weather Prediction Skin surface Temperature file. In can be stored in the .grib version 1 or 2 format but it has to be readable by \code{readGDAL}. If the input file consists of multiple bands please specify NWP.skt.band. 
}
  \item{out.dir}{
Required! Output directory where the NWP_skt.tif file will be stored.
}
  \item{out.proj}{
Required! String with the Proj4 library coordinate system definition for the output file. To retrieve it type gdalinfo -proj4 your_file in the command line.
}
  \item{out.extent}{
Required! Numeric vector with the coordinates of the CENTRE of upper-left pixel and the CENTRE of lower-right pixel in a form c(ULX,ULY,LRX,LRY).
}
  \item{out.dim}{
Required! Numeric vector with the output dimensions in the form c(ncols,nrows).
}
  \item{correction.file}{
Required! String with an absolute path to the NWP_correction.tif file generated by the \code{\link{dem_prepare}} function.
}
  \item{NWP.skt.proj}{
Optional. String with the Proj4 library coordinate system definition for the NWP.skt.file. This overwrites the information which is retrieved by the \code{readGDAL}.
}
  \item{NWP.skt.extent}{
Optional. Numeric vector with the coordinates of the CENTRE of upper-left pixel and the CENTRE of lower-right pixel in a form c(ULX,YLY,LRX,LRY) for the NWP.skt.file. This overwrites the information which is retrieved by the \code{readGDAL}.
}
  \item{NWP.skt.half.cell}{
Optional. Numeric vector with the tie point of upper-left pixel in the NWP.skt.file. By default c(0,1). For further details see \code{readGDAL}.
}
  \item{NWP.skt.scale}{
Optional. Numeric value with a scale factor which is applied to skt data.
}
  \item{NWP.skt.offset}{
Optional. Numeric value with an offset which is applied to skt data.
}
  \item{NWP.skt.band}{
Optional. The number of the skt layer in the multiband NWP.skt.file. By default 1.
}
}
\value{
0 on success or 1 on failure. Success indicates that the output NWP_skt.tif file in the GeoTiff format was written to out.dir.
}
\references{
Musial, J., Probabilistic approach to cloud and snow detection on AVHHR imagery, Atmos. Chem. Phys (in prep.)
}
\author{
Jan Musial
}
\note{
If you find one of the provided functions of this package useful please cite!
}

\seealso{
\code{\link{PCM}},\code{\link{dem_prepare}}
}
\examples{
\dontrun{
skt_prepare(NWP.skt.file='/path/to/skt_file',
out.dir='/path/to/outdir/',
out.proj="+proj=laea +a=6370997.0 +b=6370997.0 +lat_0=50.0 +lon_0=20.0",
out.extent=c(-2521000.0,2401000.0,980000.0,-2000000.0),
out.dim=c(3500,4400),
correction.file='/path/to/NWP_correction.tif')
}
}
\keyword{SKT}
\keyword{skt}
