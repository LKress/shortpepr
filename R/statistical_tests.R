#' @title Permutation test
#' @description Calculate the statistical significance of the difference of occurences of an amino acid
#' @details Calculation conducted on the description of Calis et al. (2013).
#'
#' @param label a vector with immunogenicity labels ("immunogenic-CD8" or "non-immunogenic-CD8")
#' @param amino_acid vector of amino acids (has to have same size as label)
#' @param aaoi amino acid of interest. For this amino acid, the p value is calculated
#' @param proteome_count a tibble with the column name for amino acid and value with the count of this amino acid (used for background frequencies)
#' @param n number of permutations that should be conducted
#'
#' @return p-value
#' @return differences of random frequencies between groups for all permutations
#' @return difference (between immunogenic and non immunogenic) found in original data
#'
#'
#'
permutation_test <- function(label, amino_acid, aaoi, proteome_count, n){

  distribution=c()
  result=0

  # get the amino acids for immunogenic and non immunogenic condition
  aa_immu <- amino_acid[label == "immunogenic-CD8"]
  aa_non_immu <- amino_acid[label == "non-immunogenic-CD8"]
  # calculate the frequencies
  freq_immu <- (table(aa_immu)/length(aa_immu))[aaoi]
  freq_non_immu <- (table(aa_non_immu)/length(aa_non_immu))[aaoi]
  # if no amino acid is given in one of the sets set the frequency to 0
  if (is.na(freq_immu)){freq_immu <- 0}
  if (is.na(freq_non_immu)){freq_non_immu <- 0}
  # the absolute difference is used to compare with the absolute differences of
  # the random sample
  diff_true <- abs(freq_non_immu - freq_immu)

  for(i in 1:n){
    # randomly draw amino acids with the given background distribution
    random_background <- sample(proteome_count$name,
                                length(label),
                                prob=proteome_count$value,
                                replace=TRUE)
    # divide into immunogenic and non immunogenic set
    random_immu <- random_background[label == "immunogenic-CD8"]
    random_non_immu <- random_background[label == "non-immunogenic-CD8"]
    # calculate the frequencies
    random_freq_immu <- (table(random_immu)/length(random_immu))[aaoi]
    random_freq_non_immu <- (table(random_non_immu)/length(random_non_immu))[aaoi]
    # if no amino acid is given in one of the sets set the frequency to 0
    if (is.na(random_freq_immu)){random_freq_immu <- 0}
    if (is.na(random_freq_non_immu)){random_freq_non_immu <- 0}
    # store the frequency differences between the two conditions in a vector
    distribution[i]=abs(random_freq_non_immu - random_freq_immu)
  }
  # the p value is the frequency of random samples, where the difference between
  # the two conditions is greater than in the original data
  pvalue=sum(distribution >= diff_true)/(n)
  return(list(pvalue, distribution, diff_true))
}
