% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{concat}
\alias{concat}
\title{Concatenate metrosamp objects from continuation runs into a single object}
\usage{
concat(mslist1, mslist2)
}
\arguments{
\item{mslist1}{First metrosamp object or list of metrosamp objects to concatenate.}

\item{mslist2}{Second metrosamp object or list of metrosamp objects to concatenate.}
}
\value{
A metrosamp object, or a list of metrosamp objects, representing the concatenated runs.
}
\description{
This function takes \code{\link{metrosamp}} objects from an initial run and
metrosamp objects from a continuation run and combines them into a single
metrosamp object for the entire collection of iterations.  The arguments can
be either two individual metrosamp objects, or two lists of metrosamp
objects.  In the latter case, the objects are concatenated pairwise to
produce a list of objects for the concatenated runs.
}
\details{
Debuggging information from the runs is kept if it is present in both
objects; otherwise it is dropped.

The intent of this function is to concatenate runs with subsequent runs that
started where the first run left off, but no check is made to ensure that the
second run really is a continuation of the first.  Concatenating two runs
that aren't continuations may result in an error, or in bogus results.
}
