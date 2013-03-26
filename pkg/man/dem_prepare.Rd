\name{dem_prepare}
\alias{dem_prepare}
\title{
Prepare Digital Elevation Model data and compute Skin Surface Temperature correction data
}
\description{
This function remaps and interpolates (bilinear) to desired projection and geographic extent the high resolution Digital Elevation Model (DEM) data. Moreover it creates a Skin surface Temperature (SKT) correction data which is derived as a difference between high resolution DEM and DEM provided by the NWP model (geopotential height at surface) multiply by the constant lapse rate of 0.6K/100m. This data is further used by the \code{\link{skt_prepare}} to improve the temperature estimates over the mountainous regions. The remapped DEM data is used by the PCM function together with the SKT data. If the SKT data are not available the DEM data are not needed.
}
\usage{
dem_prepare(NWP.dem.file = NWP.dem.file, dem.file = dem.file, out.dir = out.dir, out.proj = out.proj, out.extent = out.extent, out.dim = out.dim, NWP.dem.proj = NWP.dem.proj, NWP.dem.extent = NWP.dem.extent, dem.proj = dem.proj, dem.extent = dem.extent, NWP.dem.half.cell = NWP.dem.half.cell, dem.half.cell = dem.half.cell, NWP.dem.scale = NWP.dem.scale, NWP.dem.offset = NWP.dem.offset, dem.scale = dem.scale, dem.offset = dem.offset, NWP.dem.band = NWP.dem.band)
}

\arguments{
  \item{NWP.dem.file}{
Required! String with an absolute path to Numerical Weather Prediction DEM which is usually defined as geopotential height at surface. The file can be in the .grib version 1 or 2 format but it has to be readable by \code{readGDAL}. If the input file consists of multiple bands please specify NWP.dem.band. 
}
  \item{dem.file}{
Required! String with an absolute path to a single band "high" resolution DEM (preferably 1x1 km) file  which is readable by \code{readGDAL}.
}
  \item{out.dir}{
Required! Output directory where the dem.tif and NWP_correction.tif files will be stored.
}
  \item{out.proj}{
Required! String with the Proj4 library coordinate system definition for the output file. To retrieve it type gdalinfo -proj4 your_file in the command line.
}
  \item{out.extent}{
Required! Numeric vector with the coordinates of the CENTRE of upper-left pixel and the CENTRE of lower-right pixel in a form c(ULX,YLY,LRX,LRY).
}
  \item{out.dim}{
Required! Numeric vector with the output dimensions in the form c(ncols,nrows).
}
  \item{NWP.dem.proj}{
Optional. String with the Proj4 library coordinate system definition for the NWP.dem.file. This overwrites the information which is retrieved by the \code{readGDAL}.
}
  \item{NWP.dem.extent}{
Optional. Numeric vector with the coordinates of the CENTRE of upper-left pixel and the CENTRE of lower-right pixel in a form c(ULX,YLY,LRX,LRY) for the NWP.dem file. This overwrites the information which is retrieved by the \code{readGDAL}.
}
  \item{dem.proj}{
Optional. String with the Proj4 library coordinate system definition for the dem.file. This overwrites the information which is retrieved by the \code{readGDAL}.
}
  \item{dem.extent}{
Optional. Numeric vector with the coordinates of the CENTRE of upper-left pixel and the CENTRE of lower-right pixel in a form c(ULX,YLY,LRX,LRY) for the NWP.dem file. This overwrites the information which is retrieved by the \code{readGDAL}.
}
  \item{NWP.dem.half.cell}{
Optional. Numeric vector with the tie point of upper-left pixel in the NWP.dem.file. By default c(0,1). For further details see \code{readGDAL}.
}
  \item{dem.half.cell}{
Optional. Numeric vector with the tie point of upper-left pixel in the dem.file. By default c(0,1). For further details see \code{readGDAL}.
}
  \item{NWP.dem.scale}{
Optional. Numeric value with a scale factor which is applied to NWP.dem data.
}
  \item{NWP.dem.offset}{
Optional. Numeric value with an offset which is applied to NWP.dem data.
}
  \item{dem.scale}{
Optional. Numeric value with a scale factor which is applied to dem data.
}
  \item{dem.offset}{
Optional. Numeric value with an offset which is applied to dem data.
}
  \item{NWP.dem.band}{
Optional. The number of the DEM layer in the multiband NWP.dem.file. By default 1.
}
}

\value{
0 on success or 1 on failure. Success indicates that the output dem.tif and NWP_correction.tif files in the GeoTiff format were written to out.dir.
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
\code{\link{PCM}},\code{\link{skt_prepare}}
}
\examples{
\dontrun{
dem_prepare(NWP.dem.file='/path/to/NWP_dem.grib',
            dem.file='/path/to/DEM.tif',
            out.dir='/path/to/outdir/',
            out.proj="+proj=laea +a=6370997.0 +b=6370997.0 +lat_0=50.0 +lon_0=20.0",
            out.extent=c(-2521000.0,2401000.0,980000.0,-2000000.0),
            out.dim=c(3500,4400),
            NWP.dem.scale=0.1)
}
}
\keyword{DEM}
\keyword{dem}
