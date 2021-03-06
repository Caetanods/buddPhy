% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_phylip_data.R
\name{make_phylip_data}
\alias{make_phylip_data}
\title{Creates BioGeoBEARS phylip file}
\usage{
make_phylip_data(data, states, file = "biogeobears_data.txt")
}
\arguments{
\item{data}{a named vector with the tip data.}

\item{states}{a vector with the states of the data.}

\item{file}{the name for the file. The default is "biogeobears_data.txt".}
}
\value{
Writes a phylip format data file compatible with BioGeoBEARS.
}
\description{
Creates BioGeoBEARS phylip file
}
