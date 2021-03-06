% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SG.R
\name{snapPointsToLines}
\alias{snapPointsToLines}
\title{Updated version of maptools snapPointsToLines}
\usage{
snapPointsToLines(points, lines, maxDist = NA, withAttrs = TRUE)
}
\arguments{
\item{points}{A SpatialPointsDataFrame in a planar (i.e. UTM) projection}

\item{lines}{A SpatialLinesDataFrame in a planar (i.e. UTM) projection}

\item{maxDist}{Numeric value for establishing a maximum distance to avoid snapping points that are farther apart; its default value is NA.}

\item{withAttrs}{Boolean value for preserving (TRUE) or getting rid (FALSE) of the original point attributes. Default: TRUE. This parameter is optional.}
}
\description{
Updated version of maptools snapPointsToLines
}
\details{
see snapPointsToLines function from maptools.
Updated version from 	http://rstudio-pubs-static.s3.amazonaws.com/12524_7de6eb887f2945389c5d12869388be14.html
}
