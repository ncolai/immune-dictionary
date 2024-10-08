% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/irea-function.R
\name{GetEnrichmentScoreProjection}
\alias{GetEnrichmentScoreProjection}
\title{IREA analysis for transcriptome matrix input}
\usage{
GetEnrichmentScoreProjection(
  input_profile,
  input_celltype,
  genediff_cutoff = 0.25,
  species = "mouse"
)
}
\arguments{
\item{input_profile}{Gene expression matrix}

\item{input_celltype}{Choose from one of the listed cell types that most resemble the input}

\item{genediff_cutoff}{Only include the genes that are differentially expressed above this threshold
between cytokine-treated samples and PBS samples to speed up computation}
}
\description{
\code{GetEnrichmentScoreProjection} Compute enrichment score using the Wilcoxon rank sum test between
cosine similarity scores with cytokine-treated samples and cosine similarity scores with control samples
}
