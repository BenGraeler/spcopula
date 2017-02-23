\name{EU_RB_2005}
\alias{EU_RB_2005}
\docType{data}
\title{
Daily mean PM10 concentrations over Europe in 2005 as used in the JSS manuscript
}
\description{
Daily mean PM10 concentrations over Europe in 2005 as used in the JSS manuscript
}
\usage{data("EU_RB_2005")}
\format{
  The format is:
Formal class 'STFDF' [package "spacetime"] with 4 slots
  ..@ data   :'data.frame':  70810 obs. of  3 variables:
  .. ..$ PM10         : num [1:70810] 28 7 11.9 12.9 14.6 30 31.1 8.4 37.8 37.8 ...
  .. ..$ EMEP         : num [1:70810] 6.36 4.13 5.84 4.93 5.86 ...
  .. ..$ logResidKrige: num [1:70810] 12.8 12.4 10.6 11.6 17.1 ...
  ..@ sp     :Formal class 'SpatialPointsDataFrame' [package "sp"] with 5 slots
  .. .. ..@ data       :'data.frame':	194 obs. of  8 variables:
  .. .. .. ..$ station_altitude     : int [1:194] 525 581 918 560 172 117 665 1137 330 330 ...
  .. .. .. ..$ station_european_code: Factor w/ 7734 levels "AD0942A","AD0944A",..: 12 61 112 69 73 14 194 184 23 25 ...
  .. .. .. ..$ country_iso_code     : Factor w/ 39 levels "AD","AL","AT",..: 3 3 3 3 3 3 3 3 3 3 ...
  .. .. .. ..$ station_start_date   : Factor w/ 2344 levels "1900-01-01","1951-04-01",..: 1117 377 296 411 649 134 658 429 672 684 ...
  .. .. .. ..$ station_end_date     : Factor w/ 811 levels "","1900-01-01",..: 1 1 1 1 1 1 1 1 1 736 ...
  .. .. .. ..$ type_of_station      : Factor w/ 5 levels "","Background",..: 2 2 2 2 2 2 2 2 2 2 ...
  .. .. .. ..$ station_type_of_area : Factor w/ 5 levels "","rural","suburban",..: 2 2 2 2 2 2 2 2 2 2 ...
  .. .. .. ..$ street_type          : Factor w/ 5 levels "","Canyon street: L/H < 1.5",..: 1 1 5 4 4 1 4 1 2 1 ...
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
  ..@ time   :An 'xts' object on 2005-01-01/2005-12-31 containing:
  Data: int [1:365, 1] 5115 5116 5117 5118 5119 5120 5121 5122 5123 5124 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr "..1"
  Indexed by objects of class: [POSIXct,POSIXt] TZ: GMT
  xts Attributes:  
 NULL
  ..@ endTime: POSIXct[1:365], format: "2005-01-02 01:00:00" "2005-01-03 01:00:00" "2005-01-04 01:00:00" "2005-01-05 01:00:00" ...
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
data("EU_RB_2005")
str(EU_RB_2005)
}
\keyword{datasets}
