% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\docType{class}
\name{CNV.anno-class}
\alias{CNV.anno-class}
\alias{show,CNV.anno-method}
\title{CNV.anno class}
\usage{
\S4method{show}{CNV.anno}(object)
}
\arguments{
\item{object}{\code{CNV.anno} object}
}
\value{
\code{CNV.anno} class.
}
\description{
Annotations required for CNV analysis are stored in this class.
}
\details{
This class does not contain any sample data. Use \code{CNV.create_anno} to create.
}
\examples{
# create object
anno <- CNV.create_anno()

# general information
anno
show(anno)
}
\author{
Volker Hovestadt, Bjarne Daenekas \email{conumee@hovestadt.bio}
}
