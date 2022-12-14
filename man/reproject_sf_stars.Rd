% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reproject_sf_stars.R
\name{reproject_sf_stars}
\alias{reproject_sf_stars}
\title{Reproject shapefile, stars or terra object with right EPSG}
\usage{
reproject_sf_stars(EPSG = 3163, sf_stars_object = NULL, terra_object = NULL)
}
\arguments{
\item{EPSG}{target EPSG, default 3163}

\item{sf_stars_object}{sf or stars object to reproject}

\item{terra_object}{terra object to reproject}
}
\value{
return same type then input with `EPSG` projection
}
\description{
Only one variable between, `sf_stars_object` and `terra_object`.
}
