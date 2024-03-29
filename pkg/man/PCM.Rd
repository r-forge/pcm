\name{PCM}
\alias{PCM}
\title{
Probabilistic cloud & snow mask for AVHRR satellite data
}
\description{
This procedure computes probabilistic cloud & snow mask on the basis of 5/6 channel AVHRR data (with channel 3a and/or 3b), AVHRR acquisition angles, and land cover data derived from JRC GLC product (\url{http://bioval.jrc.ec.europa.eu/products/glc2000/products.php}). Optionally Skin Surface Temperature (SKT) data provided by a Numerical Weather Prediction (NWP) model (e.g. ECMWF data: \url{http://www.ecmwf.int/research/ifsdocs/CY25r1/Physics/Physics-04-06.html}) together with Digital Elevation Model can be utilize to improve classification results (highly recommended). The input data sets should be in one of the format directly readable by the rgdal package. Data in the hdf4/5 format should be converted to multiband tif file using command line utilities: gdal_translate -sds -of GTiff /your/file.hdf /your/output/dir; gdalbuildvrt -separate -input_file_list my_tif_files.txt my.vrt; gdal_translate my.vrt my_multiband.tif
}
\usage{
PCM(avhrr.file = avhrr.file, avhrr.band.order = avhrr.band.order, angles.file = angles.file, angles.band.order = angles.band.order, lc.file = lc.file, out.dir = out.dir, satellite = satellite, year = year, month = month, day = day, hour = hour, minute = minute, dem.file = dem.file, skt.file = skt.file, avhrr.scale = avhrr.scale, avhrr.offset = avhrr.offset, angles.scale = angles.scale, angles.offset = angles.offset, skt.scale = skt.scale, skt.offset = skt.offset, dem.scale = dem.scale, dem.offset = dem.offset, avhrr.noData = avhrr.noData, angles.noData = angles.noData, lc.noData = lc.noData, skt.noData = skt.noData, dem.noData = dem.noData, backingfile = backingfile, ahvrr.normalise = ahvrr.normalise, ahvrr.distance = ahvrr.distance,postclassification=postclassification,sunazi.satazi=sunazi.satazi,satazi.sunazi=satazi.sunazi)
}
\arguments{
  \item{avhrr.file}{
  Required! String with an absolute path to a multiband raster file with 5 (3a/3b) or 6 (3a&3b) AVHRR channels which could be read by \code{readGDAL}.
}
  \item{avhrr.band.order}{
  Required! Character vector of a form c('ch1','ch2','ch3b','ch4','ch5','ch3a') indicating the location of bands in the avhrr.file. It has to have exactly the same length as the number of bands in the avhrr.file. If a particular band does not contain any data put "empty" eg. c('ch1','ch2','empty','ch4','ch5','ch3a').
}
  \item{angles.file}{
  Required! String with an absolute path to a multiband raster file with at least 3 data sets containing sun zenith angle, satellite zenith angle, relative azimuth angle (values within 0-180 degree range,0 degree=forward scattering) which could be read by \code{readGDAL}.
}
  \item{angles.band.order}{
Required! Character vector of a form c('sunz','satz','razi') or c('sunz','satz','sunazi',satazi) indicating the location of bands in the angles.file. It has to have exactly the same length as the number of bands in the angles.file. If a particular band does not contain any data put "empty" eg. c('sunz','satz','razi','empty').
}
  \item{lc.file}{
Required! String with an absolute path to a single band file land cover file preferably generated by the \code{\link{lc_prepare}} which could be read by \code{readGDAL}. It HAS to contain a georeference information, because the final output file takes the projection information/geographic extent directly and exclusively from this file. 
}
  \item{out.dir}{
Required! Output directory where the final classification file will be stored as well as some temporary files.
}
  \item{satellite}{
Required! A lower-case string with a satellite name eg:'noaa07','noaa18','metopA','metopB'.
}
  \item{year}{
Required! Numeric value with a year of an image acquisition eg: 2009
}
  \item{month}{
Required! Numeric value with a month of an image acquisition eg: 5
}
  \item{day}{
Required! Numeric value with a day of an image acquisition eg: 22
}
  \item{hour}{
Required! Numeric value with an hour of an image acquisition (0-23) in the UTC time-zone.
}
  \item{minute}{
Required! Numeric value with a minute of an image acquisition in the UTC time-zone.
}
  \item{dem.file}{
Optional. String with an absolute path to a single band digital elevation model file preferably generated by the \code{\link{dem_prepare}} which could be read by \code{readGDAL}. Required if the skt.file is provided.
}
  \item{skt.file}{
Optional. String with an absolute path to a single band skin surface temperature file preferably generated by the \code{\link{skt_prepare}} which could be read by \code{readGDAL}.
}
  \item{avhrr.scale}{
Optional. Numeric vector with scale factors which are applied to bands specified by the avhrr.band.order. It has to have exactly the same length as the number of bands in the avhrr.file eg. c(0.0001,0.0001,0.01,0.01,0.01,0.0001).
}
  \item{avhrr.offset}{
Optional. Numeric vector with offset values which are applied to bands specified by the avhrr.band.order. It has to have exactly the same length as the number of bands in the avhrr.file eg. c(0,0,0,0,0,0).
}
  \item{angles.scale}{
Optional. Numeric vector with scale factors which are applied to bands specified by the angles.band.order. It has to have exactly the same length as the number of bands in the angles.file eg. c(0.01,0.01,0.01).
}
  \item{angles.offset}{
Optional. Numeric vector with offset values which are applied to bands specified by the angles.band.order. It has to have exactly the same length as the number of bands in the angles.file eg. c(0,0,0).
}
  \item{skt.scale}{
Optional. Numeric value with a scale factor which is applied to skt data.
}
  \item{skt.offset}{
Optional. Numeric value with an offset which is applied to skt data.
}
  \item{dem.scale}{
Optional. Numeric value with a scale factor which is applied to dem data.
}
  \item{dem.offset}{
Optional. Numeric value with an offset which is applied to dem data.
}
  \item{avhrr.noData}{
Optional. Numeric vector with noData values which are applied to bands specified by the avhrr.band.order. It has to have exactly the same length as the number of bands in the avhrr.file. NA values are allowed.
}
  \item{angles.noData}{
Optional. Numeric vector with noData values which are applied to bands specified by the angles.band.order. It has to have exactly the same length as the number of bands in the angles.file. NA values are allowed.
}
  \item{lc.noData}{
noData value for the land cover data. NA value is allowed.
}
  \item{skt.noData}{
Optional. noData value for the skin surface temperature data. NA value is allowed.
}
  \item{dem.noData}{
Optional. noData value for the dem data. NA value is allowed.
}
  \item{backingfile}{
Optional. String with an absolute path to a temporary file of a class \code{big.matrix} which is written to a hard drive to safe some RAM during the computations over vast areas. It can be used if you run out of RAM but this will slow down the processing. 
}
  \item{ahvrr.normalise}{
Optional. Logical vector indicating which one of bands specified by the avhrr.band.order should be normalised by the cosine of the sun zenith angle. It has to have exactly the same length as the number of bands in the avhrr.file eg. c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE).
}
  \item{ahvrr.distance}{
Optional. Logical vector indicating which one of bands specified by the avhrr.band.order should be corrected for sun/earth distance variations. It has to have exactly the same length as the number of bands in the avhrr.file eg. c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE).
}
  \item{postclassification}{
Optional. Logical value indicating if the classification probability estimates should be recoded to binary form and the cloud shadow mask to be computed. By default FALSE.
}
\item{sunazi.satazi}{
Optional. Logical value indicating if relative azimuth angle should be computed as an absolute difference between sun azimuth and satellite azimuth. By default FALSE. Acquired values should be within 0-180 degree range and 0 degree should indicate forward scattering.
}
\item{satazi.sunazi}{
Optional. Logical value indicating if relative azimuth angle should be computed as an absolute difference between satellite azimuth and sun azimuth. By default FALSE. Acquired values should be within 0-180 degree range and 0 degree should indicate forward scattering.
}

}
\details{
The output classification probability file is coded according to the following convention: 
  0 value means 0 percent probability of snow class and 100 percent probability of clear-sky class
 50 value means 50 percent probability of snow class and 50 percent probability of clear-sky class
100 value means 100 percent probability of snow class and 0 percent probability of clear-sky/cloudy classes
150 value means 50 percent probability of snow class and 50 percent probability of cloudy class
200 value means 0 percent probability of snow/clear-sky classes and 100 percent probability of cloudy class
250 value means 50 percent probability of clear-sky class and 50 percent probability of cloudy class
300 value means 1-0 percent probability of clear-sky class and 0 percent probability of cloudy class

The output binary classification is coded according to the following convention:
0 <- no data
1 <- water class taken from land cover
2 <- land class taken from land cover
3 <- snow class
4 <- cloud shadow over water
5 <- cloud shadow over land
6 <- cloud shadow over snow
7 <- pixel adjacent to cloud over water
8 <- pixel adjacent to cloud over land
9 <- pixel adjacent to cloud over snow
10<- cloudy class

So many categories allow user to ignore the adjacent and cloud shadow discrimination and to reclassify the results into water,land,snow,cloudy categories.

It the SKT data are not provided the cloud shadow are computed assuming the contrast cloud top height of around 3.5 km.

The original implementation of the PCM algorithm utilizes way more LUTs which are suited separately for different hours of the satellite overpass as well as different seasons. Here the classification process is based on two generalized LUTs separately for 3a/3b channel combinations. Therefore quality of the classification is slightly lower. 
}
\section{Warning}{
All the input data should be stacked i.e. have the same number of columns and lines and the same geographic extent. The land cover file has to have the georeference information included and visible to GDAL library because this information is directly copied to the output files. This is always the case if the routine lc_prepare was used. Moreover the algorithm was applied only to LAC AVHRR data covering Europe os application to other regions involves preparation of land cover data with exactly the same codes as you can find in GLC2000 codes (\url{http://bioval.jrc.ec.europa.eu/products/glc2000/legend.php}). Land cover codes for deserted areas and agriculture class are changed, thus please generate use function \code{\link{lc_prepare}} to see final coding.

The cloud shadow mask is generated on a basis of simple geometrical relationships and has nothing to do with spectral analysis of an image.

The algorithm was tested for NOAA16,17,18,19,MetopA satellites over Europe and northern parts of Africa and the applicability to other regions is still to be determined.
}
\value{
0 on success or 1 on failure. Success indicates that the output classification probability file and/or binary classification file in the GeoTiff format was/were written to out.dir
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
\code{\link{lc_prepare}},\code{\link{dem_prepare}},\code{\link{skt_prepare}}
}
\examples{
\dontrun{
##NOAA17 with 3a channel
res<-PCM(avhrr.file=system.file("examples/noaa17_20080101_1000_28693.eu1km.avhrr.tif", package = "PCM"),
         angles.file=system.file("examples/noaa17_20080101_1000_28693.eu1km.sunsatangles.tif", package = "PCM"),
         lc.file=system.file("examples/lc.tif", package = "PCM"),
         dem.file=system.file("examples/dem.tif", package = "PCM"),
         skt.file=system.file("examples/NWP_skt.tif", package = "PCM")[1],
         out.dir='/tmp/',
         avhrr.band.order=c('ch1','ch2','empty','ch4','ch5','ch3a'),
         ahvrr.normalise=c(TRUE,TRUE,FALSE,FALSE,FALSE,TRUE),
         ahvrr.distance=c(TRUE,TRUE,FALSE,FALSE,FALSE,TRUE),
         angles.band.order=c('sunz','satz','razi'),
         avhrr.scale=c(0.0001,0.0001,0.01,0.01,0.01,0.0001),
         avhrr.offset=c(0,0,273.15,273.15,273.15,0),
         angles.scale=c(0.01,0.01,0.01),
         angles.offset=c(0,0,0),
         skt.scale=0.1,
         skt.offset=0,
         dem.scale=1,
         dem.offset=0,
         avhrr.noData=c(-32000,-32000,-32000,-32000,-32000,-32000),
         angles.noData=c(0,0,0),
         dem.noData=65535,
         lc.noData=65535,
         skt.noData=65535,
         satellite='noaa17',
         year=2008,
         month=1,
         day=1,
         hour=10,
         minute=0
)

#NOAA18 with 3b channel
res<-PCM(avhrr.file=system.file("examples/noaa18_20080101_1216_13483.eu1km.avhrr.tif", package = "PCM"),
         angles.file=system.file("examples/noaa18_20080101_1216_13483.eu1km.sunsatangles.tif", package = "PCM"),
         lc.file=system.file("examples/lc.tif", package = "PCM"),
         dem.file=system.file("examples/dem.tif", package = "PCM"),
         skt.file=system.file("examples/NWP_skt.tif", package = "PCM")[1],
         out.dir='/tmp/',
         avhrr.band.order=c('ch1','ch2','ch3b','ch4','ch5','empty'),
         ahvrr.normalise=c(TRUE,TRUE,FALSE,FALSE,FALSE,TRUE),
         ahvrr.distance=c(TRUE,TRUE,FALSE,FALSE,FALSE,TRUE),
         angles.band.order=c('sunz','satz','razi'),
         avhrr.scale=c(0.0001,0.0001,0.01,0.01,0.01,0.0001),
         avhrr.offset=c(0,0,273.15,273.15,273.15,0),
         angles.scale=c(0.01,0.01,0.01),
         angles.offset=c(0,0,0),
         skt.scale=0.1,
         skt.offset=0,
         dem.scale=1,
         dem.offset=0,
         avhrr.noData=c(-32000,-32000,-32000,-32000,-32000,-32000),
         angles.noData=c(0,0,0),
         dem.noData=65535,
         lc.noData=65535,
         skt.noData=65535,
         satellite='noaa18',
         year=2008,
         month=1,
         day=1,
         hour=12,
         minute=16
)
}
}
\keyword{PCM}
\keyword{AVHRR}
\keyword{pcm}
\keyword{satellite}
\keyword{cloud}
\keyword{snow}
\keyword{detection}
