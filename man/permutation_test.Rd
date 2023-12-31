% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistical_tests.R
\name{permutation_test}
\alias{permutation_test}
\title{Permutation test}
\usage{
permutation_test(label, amino_acid, aaoi, proteome_count, n)
}
\arguments{
\item{label}{a vector with immunogenicity labels ("immunogenic-CD8" or "non-immunogenic-CD8")}

\item{amino_acid}{vector of amino acids (has to have same size as label)}

\item{aaoi}{amino acid of interest. For this amino acid, the p value is calculated}

\item{proteome_count}{a tibble with the column name for amino acid and value with the count of this amino acid (used for background frequencies)}

\item{n}{number of permutations that should be conducted}
}
\value{
p-value

differences of random frequencies between groups for all permutations

difference (between immunogenic and non immunogenic) found in original data
}
\description{
Calculate the statistical significance of the difference of occurences of an amino acid
}
\details{
Calculation conducted on the description of Calis et al. (2013).
}
