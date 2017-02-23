\name{EU_RB}
\alias{EU_RB}
\docType{data}
\title{
Daily mean PM10 concentrations over Europe in June and July 2005
}
\description{
Daily mean PM10 concentrations over Europe in June and July 2005
}
\usage{data("EU_RB")}
\format{
  The format is:
Formal class 'STSDF' [package "spacetime"] with 5 slots
  ..@ data   :'data.frame':	11834 obs. of  2 variables:
  .. ..$ PM10: num [1:11834] 14 9.7 7.8 21.9 11.2 9 11 6.1 7.4 7.4 ...
  .. ..$ EMEP: num [1:11834] 11.62 5.02 3.94 3.82 7.01 ...
  ..@ index  : int [1:11834, 1:2] 1 2 3 4 5 6 7 8 9 10 ...
  ..@ sp     :Formal class 'SpatialPointsDataFrame' [package "sp"] with 5 slots
  .. .. ..@ data       :'data.frame':	194 obs. of  1 variable:
  .. .. .. ..$ station_altitude: int [1:194] 525 581 918 560 172 117 665 1137 330 330 ...
  .. .. ..@ coords.nrs : num(0) 
  .. .. ..@ coords     : num [1:194, 1:2] 4592866 4761515 4658756 4690954 4799839 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr [1:2] "coords.x1" "coords.x2"
  .. .. ..@ bbox       : num [1:2, 1:2] 2749697 1647732 6412269 4604814
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:2] "coords.x1" "coords.x2"
  .. .. .. .. ..$ : chr [1:2] "min" "max"
  .. .. ..@ proj4string:Formal class 'CRS' [package "sp"] with 1 slots
  .. .. .. .. ..@ projargs: chr " +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
  ..@ time   :An 'xts' object on 2005-06-01/2005-07-31 containing:
  Data: int [1:61, 1] 5266 5267 5268 5269 5270 5271 5272 5273 5274 5275 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr "..1"
  Indexed by objects of class: [POSIXct,POSIXt] TZ: GMT
  xts Attributes:  
 NULL
  ..@ endTime: POSIXct[1:61], format: "2005-06-02 02:00:00" "2005-06-03 02:00:00" "2005-06-04 02:00:00" ...
}
\source{
Obtained from the european Air Qualtiy airbase http://acm.eionet.europa.eu/databases/airbase/.
}
\references{
http://acm.eionet.europa.eu/databases/airbase/

Graeler, B., L. E. Gerharz, & E. Pebesma (2012): Spatio-temporal analysis and interpolation of PM10 measurements in Europe. ETC/ACM Technical Paper 2011/10, January 2012.
http://acm.eionet.europa.eu/reports/ETCACM_TP_2011_10_spatio-temp_AQinterpolation
}
\examples{
data("EU_RB")
str(EU_RB)
}
\keyword{datasets}
