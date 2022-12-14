% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_var_JSDM.R
\name{extract_var_JSDM}
\alias{extract_var_JSDM}
\title{Extract layers and create power wanted of each variable}
\usage{
extract_var_JSDM(
  stars_object,
  variables_names,
  power_variable = rep(1, length(variables_names)),
  scale = rep(TRUE, length(variables_names))
)
}
\arguments{
\item{stars_object}{stars object. output of "read_stars" of "stars" library with layer's names.}

\item{variables_names}{a character vector. of layer's name needed}

\item{power_variable}{int vector same length then `variables_names` . maximal power needed for each variable}

\item{scale}{boolean vector same length then `variables_names` . scale variables and its power}
}
\value{
stars object with only variables and power asked.
}
\description{
Extract layers and create power wanted of each variable
}
