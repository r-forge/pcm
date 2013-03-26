\name{lc_prepare}
\alias{lc_prepare}
\title{
Prepare land cover data
}
\description{
This function remaps and interpolates (nearest neighbor) to desired projection and geographic extent the GLC2000 land cover data for Europe and Africa which could be obtained from \url{http://bioval.jrc.ec.europa.eu/products/glc2000/products.php} free of charge for scientific usage. The data should be in binary form and correspond to Africa v5.0 and Europe v2.0 products. 
}
\usage{
lc_prepare(lc.africa.file = lc.africa.file, lc.europe.file = lc.europe.file, out.dir = out.dir, out.proj = out.proj, out.extent = out.extent, out.dim = out.dim, Africa = Africa)
}

\arguments{
  \item{lc.africa.file}{
Required! String with an absolute path to the africa_v5.bil file.
}
  \item{lc.europe.file}{
Required! String with an absolute path to the GLC_EU_V2.bil file.
}
  \item{out.dir}{
Required! Output directory where the lc.tif will be stored.
}
  \item{out.proj}{
Required! String with the Proj4 library coordinate system definition. To retrieve it type gdalinfo -proj4 your_file in the command line.
}
  \item{out.extent}{
Required! Numeric vector with the coordinates of the CENTRE of upper-left pixel and the CENTRE of lower-right pixel in a form c(ULX,YLY,LRX,LRY).
}
  \item{out.dim}{
Required! Numeric vector with the output dimensions in the form c(ncols,nrows).
}
  \item{Africa}{
Optional. Logical value indicating if Africa data should be utilized. In case your region of interest is only in the Europe you can set this keyword to FALSE. By default TRUE.
}
}
\value{
0 on success or 1 on failure. Success indicates that the output lc.tif file in the GeoTiff format was written to out.dir.
}
\references{
The Land Cover Map for Africa in the Year 2000. P. Mayaux, E. Bartholome, A. Cabral, M. Cherlet, P. Defourny, A. Di Gregorio, O. Diallo, M. Massart, A. Nonguierma, J.-F. Pekel, C. Pretorius, C. Vancutsem, M. Vasconcelos. GLC2000 database, European Commision Joint Research Centre, 2003.
The Land Cover Map for Europe in the Year 2000. A. Hartley, J-F.Pekel, M.Ledwith, J-L.Champeaux, E.De Badts, and S.A.Bartalev. GLC2000 database, European Commision Joint Research Centre, 2006. 
}
\author{
Jan Musial
}
\note{
Please read JRC disclaimer information: \url{http://bioval.jrc.ec.europa.eu/products/glc2000/disclaimer.php}
If you use the GLC data in your publication please cite:"Global Land Cover 2000 database. European Commission, Joint Research Center, 2003, http://www-gem.jrc.it/glc2000"
}

\seealso{
\code{\link{PCM}}
}
\examples{
\dontrun{
lc_prepare(lc.africa.file='/path/to/africa_v5.bil',
               lc.europe.file='/path/to/GLC_EU_V2.bil',
               out.dir='/path/to/outdir/',
               out.proj="+proj=laea +a=6370997.0 +b=6370997.0 +lat_0=50.0 +lon_0=20.0",
               out.extent=c(-2521000.0,2401000.0,980000.0,-2000000.0),
               out.dim=c(3500,4400),
               Africa=T)
}
}
\keyword{lc_prepare}
\keyword{GLC2000}
\keyword{landcover}
