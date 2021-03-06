% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DIGREsyn.R
\docType{data}
\name{KEGGnet.mat}
\alias{KEGGnet.mat}
\title{Universal KEGG pathway Gene Network}
\format{An object of class \code{matrix} with 2709 rows and 2709 columns.}
\usage{
data(KEGGnet.mat)
}
\description{
A matrix contains the global pathway (GP) gene-gene interaction network merged from KEGG pathway database.
}
\details{
We selected KEGG pathways belonging to genetic information processing, environmental 
information processing, cellular processing, and cancer disease. From this set of selected pathways we removed 
any pathway with fewer than 10 edges. Finally we merged the remaining 32 KEGG pathways into a global pathway (GP) 
which included 11642 interactions among 2322 genes.
}
\keyword{datasets}
