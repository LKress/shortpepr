amino_acids = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                "M", "N", "P", "Q", "R", "S", "T", "W", "V", "Y")

#' @title Sequences to PWM
#' @description Create position weight matrix (PWM) form list of equal length strings
#' @details ...
#'
#' @param sequences vector of input sequences
#' @param alphabet list of characters that can occure in the sequences (default: 20 amino acids)
#' @param pc pseudo count (default 1)
#' @param background A named list with the character as name and backgroud frequency for that character as value
#'
#' @return A tibble containing the position weight matrix (columns: letters; rows: positions)
#'
#' @import tidyverse
#' @export
seqs2pwm <- function(sequences, background, alphabet = amino_acids, pc = 1){
  library(tidyverse)
  # do all steps from pcm creation to pwm generation
  pcm <- seqs2pcm(sequences, alphabet)
  pcm_pseudo <- add_pseudocount(pcm, pc)
  pfm <- pcm2pfm(pcm_pseudo, background=NULL)
  pwm <- pfm2pwm(pfm, background)
  return(pwm)
}

#' @title Sequences to PFM
#' @description Create position frequency matrix (PFM) form list of equal length strings
#' @details The frequency can also be normalised by a background frequency.
#'
#' @param sequences vector of input sequences
#' @param alphabet list of characters that can occure in the sequences (default: 20 amino acids)
#' @param background A named list with the character as name and backgroud frequency for that character as value
#'
#' @return A tibble containing the position frequency matrix (columns: letters; rows: positions)
#'
#' @import tidyverse
#' @export
seqs2pfm <- function(sequences, background = NULL, alphabet = amino_acids, pc=1){
  library(tidyverse)
  # do all steps from pcm creation to pwm generation
  pcm <- seqs2pcm(sequences, alphabet)
  pcm_pseudo <- add_pseudocount(pcm, pc)
  pfm <- pcm2pfm(pcm_pseudo, background)
  return(pfm)
}

#' @title Sequences to PCM
#' @description Create position count matrix (PCM) form list of equal length strings
#' @details For every position the different characters are counted and returned
#'     as a matrix.
#'
#' @param sequences vector of input sequences
#' @param alphabet list of characters that can occure in the sequences (default: 20 amino acids)
#'
#' @return A tibble containing the position count matrix (columns: letters; rows: positions)
#'
#' @import tidyverse
#' @export
seqs2pcm <- function(sequences, alphabet = amino_acids){
  library(tidyverse)
  stopifnot("Elements of your input vector do not have the same length!" =
              all(str_length(sequences[1])== str_length(sequences)))
  seq_length <- str_length(sequences[1])
  seq_tibble <- tibble(sequences = sequences)
  pcm <- seq_tibble %>%
    # split the sequences into different columns
    separate(sequences, into = as.character(0:seq_length), sep="") %>%
    # separate with "" will produce an empty string before the first character
    # which is removed with the following line
    dplyr::select(-"0") %>%
    pivot_longer(as.character(1:seq_length), names_to = "position", values_to = "alphabet") %>%
    # the position has to be sorted numerically, otherwise the positions are
    # in incorrect alphabetical order
    mutate(position = as.integer(position)) %>%
    arrange(position) %>%
    count(position, alphabet) %>%
    pivot_wider(names_from = alphabet, values_from = n) %>%
    dplyr::select(-position) %>%
    replace(is.na(.), 0)
  # if the provided sequences do not contain all characters from the alphabet,
  # the missing characters are added as empty columns
  pcm[,setdiff(alphabet, colnames(pcm))] = 0
  return(pcm)
}

#' @title Add pseudocount
#' @description Add a pseudo count to a PCM or a PFM
#' @details ...
#'
#' @param matrix PCM or PFM matrix
#' @param pc The pseoudocount that should be added (default 1)
#'
#' @return A tibble containing the input matrix with the added pseudo count.
#'
#' @import tidyverse
#' @export
add_pseudocount <- function(matrix, pc = 1){
  library(tidyverse)
  pc_matrix <- tibble(matrix + pc)
  return(pc_matrix)
}

#' @title PCM to PFM
#' @description Create a position frequency matrix (PFM) from a position count matrix (PCM)
#' @details The frequencies of the characters in the PCM for each position are computed
#'     and returned as PFM tibble. The frequencies can be normalised by a background
#'     frequency if needed.
#'
#' @param pcm The position count matrix
#' @param background Normalize the frequencies by a background frequency
#'
#' @return A tibble containing position frequency matrix (PFM).
#'
#' @import tidyverse
#' @export
pcm2pfm <- function(pcm, background = NULL){
  library(tidyverse)
  # if no background is given, all values are set to 1 so the normalisation using
  # the resulting named vector has no effect
  if (is.null(background)){
    len <- dim(pcm)[2]
    background <- setNames(rep(1, len), amino_acids)
  }

  pfm <- pcm %>%
    mutate(position = 1:nrow(.)) %>%
    pivot_longer(-position, names_to = "chr", values_to = "counts") %>%
    group_by(position) %>%
    mutate(
      # normalise by the background frequencies
      freq = (counts/sum(counts))/background[chr],
      # scale back to proportions (range 0-1)
      freq = freq/sum(freq)
    ) %>%
    ungroup() %>%
    dplyr::select(-counts) %>%
    pivot_wider(names_from = "chr", values_from = "freq") %>%
    dplyr::select(-position)
}

#' @title PFM to PWM
#' @description Create a position weight matrix (PWM) from a position frequency matrix (PFM)
#' @details For each position and character, the log likelihood is calculated, to
#'     get the position weight matrix. This matrix can be used to score other sequences.
#'
#' @param pfm The position frequency matrix
#' @param background A named list with the character as name and backgroud frequency for that character as value
#'
#' @return A tibble containing position frequency matrix (PFM).
#'
#' @import tidyverse
#' @export
pfm2pwm <- function(pfm, background){
  library(tidyverse)
  stopifnot("Background frequencies do not match PWM" =
              nrow(background) == ncol(pfm))
  cols <- c()
  # calculate log likelihood for each column that contains a character
  for (chr in colnames(pfm)) {
    cols <- c(cols, log2(pfm[chr]/background[chr]))
  }
  pwm <- bind_cols(cols)
}

#' @title Score sequence
#' @description Score a sequence with a position weight matrix (PWM)
#' @details For each character in the sequence, the score from the PWM is summed
#'     up to get the final score
#'
#' @param seq List of character sequences matching the PWM
#' @param pwm Position weight matrix tibble
#'
#' @return Vector of scores
#'
#' @import tidyverse
#' @export
score_seq <- function(seq, pwm){
  library(tidyverse)
  score <- c()
  # iterate over all given sequences
  for (sequence in seq){
    # sequence and pwm have to be compatible (= same length)
    stopifnot("PWM and sequence do not match" =
      dim(pwm)[1] == str_length(sequence))

    chr_list <- str_split(sequence, "")
    score_sequence <- 0
    for (pos in 1:str_length(sequence)){
      chr <- chr_list[[1]][pos]
      stopifnot("Character in sequence that is not given in PWM" =
                  chr %in% colnames(pwm))
      pwm_value <- pwm[chr][[1]][pos]
      score_sequence <- score_sequence + pwm_value
    }
    score <- c(score, score_sequence)
  }
  return(score)
}

#' @title Score sequence allele combination
#' @description Score a sequence allele length combination with a position weight matrix (PWM)
#' @details For each character in the sequence, the score from the PWM is summed
#'     up to get the final score
#'
#' @param seq List of character sequences matching the PWM
#' @param allele List of alleles that have the same length as seq
#' @param pwm tibble containing columns allele, peptide length, pwm
#' @param normalise If TRUE, the score is normalised by the maximum score of the given PWM
#'
#' @return Vector of scores
#'
#' @import tidyverse
#' @export
score_seq_allele <- function(seq, allele, pwm, normalise = T){
  library(tidyverse)
  stopifnot("Sequence and allele vector lengths have to match" =
              length(seq) == length(allele))
  score <- c()
  # iterate over all given sequences
  for (i in length(seq)){
    normalisation_factor <- 1
    seq_len <- str_length(seq[i])
    current_allele <- allele[i]
    current_pwm <- pwm %>%
      filter(len == seq_len, allele == current_allele) %>% pull(pwm)

    # if there is no PWM given in the PWM-tibble for this allele length combination
    # set score to na
    if (length(current_pwm) < 1){
      return(NA)
    }
    current_pwm <- current_pwm[[1]]
    if (normalise){
      normalisation_factor <- pwm_max(current_pwm)
    }
    # sequence and pwm have to be compatible (= same length)
    stopifnot("PWM and sequence do not match" =
                dim(current_pwm)[1] == seq_len)

    chr_list <- str_split(seq[i], "")
    score_sequence <- 0
    for (pos in 1:str_length(seq[i])){
      chr <- chr_list[[1]][pos]
      stopifnot("Character in sequence that is not given in PWM" =
                  chr %in% colnames(current_pwm))
      pwm_value <- current_pwm[chr][[1]][pos]
      score_sequence <- score_sequence + pwm_value
    }
    # normalise the score by the maximum possible score of the given PWM if normalise
    # is TRUE, otherwise normalisation_factor is 1 and this operation has no effect
    score_sequence <- score_sequence/normalisation_factor
    score <- c(score, score_sequence)
  }
  return(score)
}

#' @title Score sequence allele combination for specific positions
#' @description Score a sequence allele length combination with a position weight matrix (PWM) for one position in the given sequence
#' @details ...
#'
#' @param seq List of character sequences matching the PWM
#' @param allele List of alleles that have the same length as seq
#' @param pos List of positions
#' @param pwm tibble containing columns allele, peptide length, pwm
#' @param normalise If TRUE, the score is normalised by the maximum score of the given PWM
#'
#' @return List of PWM scores at given position
#'
#' @import tidyverse
#' @export
score_seq_allele_pos <- function(seq, mhc, pos, pwm_tibble, normalise = T){
  library(tidyverse)

  chr <- str_sub(seq, pos, pos)
  seq_len <- str_length(seq)

  current_pwm <- pwm_tibble %>%
    filter(len == seq_len, allele == mhc) %>% pull(pwm)

  # if there is no PWM given in the PWM-tibble for this allele length combination
  # set score to na
  if (length(current_pwm) < 1){
    return(NA)
  }

  current_pwm <- current_pwm[[1]]


  # if no normalisation should be applied, the score is divided by 1 which has no effect
  normalisation_factor <- 1
  if (normalise){
    normalisation_factor <- pwm_max(current_pwm)
  }

  # sequence and pwm have to be compatible (= same length)
  stopifnot("PWM and sequence do not match" =
              dim(current_pwm)[1] == seq_len)

  stopifnot("Character in sequence that is not given in PWM" =
              chr %in% colnames(current_pwm))

  score <- current_pwm[chr][[1]][pos]

  # normalise the score by the maximum possible score of the given PWM if normalise
  # is TRUE, otherwise normalisation_factor is 1 and this operation has no effect
  score <- score/normalisation_factor

  return(score)
}

score_seq_allele_pos <- Vectorize(
  score_seq_allele_pos,
  vectorize.args = c("seq", "mhc", "pos")
)


#' @title Score sequence allele combination position wise
#' @description Score a sequence allele length combination with a position weight matrix (PWM)
#' @details For each character in the sequence, the PWM score is added as an extra column
#'
#' @param seq List of character sequences matching the PWM
#' @param allele List of alleles that has the same length as seq
#' @param pwm tibble containing columns allele, peptide length, pwm
#' @param normalise If TRUE, the score is normalised by the maximum score of the given PWM
#'
#' @return Vector of scores
#'
#' @import tidyverse
#' @export
score_seq_allele_position_wise <- function(aa, pos, allele, len, pwm, normalise =T){
  library(tidyverse)

  current_pwm <- pwm %>%
    filter(len == len, allele == allele) %>% pull(pwm)

  # if there is no PWM given in the PWM-tibble for this allele length combination
  # set score to na
  if (length(current_pwm) < 1){
    return(NA)
  }
  current_pwm <- current_pwm[[1]]
  if (normalise){
    normalisation_factor <- pwm_max(current_pwm)
  }
  # sequence and pwm have to be compatible (= same length)
  stopifnot("PWM and sequence do not match" =
              dim(current_pwm)[1] == len)

  stopifnot(
    "Character in sequence that is not given in PWM" = aa %in% colnames(current_pwm))
  pwm_value <- current_pwm[aa][[1]][pos]

  # normalise the score by the maximum possible score of the given PWM if normalise
  # is TRUE, otherwise normalisation_factor is 1 and this operation has no effect
  pwm_value <- pwm_value/normalisation_factor

  return(pwm_value)
}

score_seq_allele_position_wise <- Vectorize(
  score_seq_allele_position_wise,
  vectorize.args = c('aa','pos', 'allele', 'len')
)



#' @title Calculate information content
#' @description Calculate the information content at each position of a PFM
#' @details The information content for a given PFM is calculated as described in ggplot.
#'
#' @param pfm Position frequency matrix
#'
#' @return Vector of information content values (same length as input PFM)
#'
#' @import tidyverse
#' @export
calculate_information_content <- function(pfm){
  library(tidyverse)
  # the alphabet size defines the maximum amount of information content
  alphabet_size <- ncol(pfm)
  entropy <- pfm %>%
    rowwise() %>%
    # calculate information content as described in https://doi.org/10.1016/j.xcrm.2021.100194
    #mutate(ic = 1 + sum(across(everything())*(log10(across(everything()))/log10(alphabet_size)))) %>%
    # calculate entropy as described in https://bioinformaticshome.com/bioinformatics_tutorials/sequence_alignment/introduction_to_information_theory_page3.html
    # and used in ggseqlogo
    mutate(e = - sum(across(everything())*log2(across(everything())))) %>%
    pull(e)
  ## calculate information_content as described in https://bioinformaticshome.com/bioinformatics_tutorials/sequence_alignment/introduction_to_information_theory_page3.html
  return(log2(alphabet_size) - entropy)
}

#' @title Get max PWM score
#' @description Calculate the maximum possible score of an PWM
#' @details ...
#'
#' @param pwm Position frequency matrix
#'
#' @return Maximum score that can be reached for given PWM
#'
#' @import tidyverse
#' @export
pwm_max <- function(pwm){
  library(tidyverse)
  pwm_max <- pwm %>%
    pivot_longer(cols = everything(), names_to = "chr", values_to = "pwm_score") %>%
    mutate(pos = rep(1:nrow(pwm), each=ncol(pwm))) %>%
    group_by(pos) %>%
    summarise(max_chr = max(pwm_score)) %>%
    pull(max_chr) %>%
    sum()
  return(pwm_max)
}


#' @title Create PWM tibble
#' @description Create a tibble, that contains for every allele, peptide length a PWM
#' @details ...
#'
#' @param df Input dataframe containing the "peptide" that should be converted into a PWM,
#'    the "allele" and the "len" which denotes the length of the peptide
#' @param min_datapoints Minimum number of datapoints for this length allele combination
#'    to create a PWM
#' @param background The background that should be used to create the PWM
#'
#' @return Tibble that contains the PWM for each length allele combination
#'
#' @import tidyverse
#' @export
create_pwm_df <- function(df, min_datapoints, background){
  library(tidyverse)
  stopifnot("Input dataframe has to contain the columns: peptide, allele, len" =
              all(c("peptide", "allele", "len") %in% colnames(df)))
  pwm_tibble <- df %>%
    count(len, allele) %>%
    filter(n >= min_datapoints) %>%
    mutate(
      pwm = df %>%
        dplyr::select(allele, len, peptide) %>%
        group_by(len,allele) %>%
        filter(n() >= min_datapoints) %>%
        # returns a list of PWMs for each allele length group
        group_map(~ shortpepr::seqs2pwm(.x$peptide, background=background))
    )
  return(pwm_tibble)
}
