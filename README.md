
# shortpepr

The shortpepr tool can be used to generate PWMs from peptides and to
calculate the information content. It is intended to use with peptides
that are bound to MHCI molecules, but can be applied to all gap free
peptide alignments.

## Installation

You can install the development version of shortpepr like so:

``` r
if(!require(devtools)){
    install.packages("devtools")
    library(devtools)
}

devtools::build()

install.packages("<path_to_repo>/shortpepr_0.1.1.tar.gz", repos = NULL, type="source")
```

## Examples

Shortpepr is designed to work with tibbles from tidyverse. As example
epitopes from the IEDB are loaded.

``` r
suppressMessages({
library(shortpepr)
library(tidyverse)
})

# dataframe containing the peptide, allele and the length of the peptide
example_epitopes <- system.file(
  "extdata",
  "test_epitopes.tsv", 
  package = "shortpepr")
example_epitopes_df <- read_tsv(example_epitopes)
#> Rows: 888 Columns: 3
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (2): peptide, allele
#> dbl (1): len
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
head(example_epitopes_df)
#> # A tibble: 6 × 3
#>   peptide    allele   len
#>   <chr>      <chr>  <dbl>
#> 1 AIVNYANL   H2Kb       8
#> 2 TSYKFESV   H2Kb       8
#> 3 VSLDYINTM  H2Kb       9
#> 4 AAFEFINSL  H2Kb       9
#> 5 YAPVSPIVI  H2Db       9
#> 6 YSLPNAGDVI H2Db      10
```

To generate PWMs background frequencies are required.

``` r
# the background count and frequencies of amino acids (e.g. from the proteome)
example_background <- system.file(
  "extdata",
  "test_background.tsv", 
  package = "shortpepr")
example_background_df <- read_tsv(example_background) %>% select(name, freq) %>% deframe()
#> Rows: 20 Columns: 3
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr (1): name
#> dbl (2): value, freq
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
example_background_df
#>          A          C          D          E          F          G          H 
#> 0.06851810 0.02221453 0.04853844 0.07083450 0.03628429 0.06333933 0.02598077 
#>          I          K          L          M          N          P          Q 
#> 0.04353847 0.05723104 0.09947218 0.02223334 0.03581281 0.06203675 0.04845813 
#>          R          S          T          V          W          Y 
#> 0.05617254 0.08586854 0.05443012 0.06068323 0.01188343 0.02646946
```

A PWM can be generated from a vector of peptides.

``` r
# a PWM can be directly generated from a list of peptides and a given background frequency
test_pwm <- example_epitopes_df %>% 
  filter(allele == "H2Db", len==9) %>% 
  pull(peptide) %>% 
  seqs2pwm(background = example_background_df)

test_pwm
#> # A tibble: 9 × 20
#>        A       C       F      G       H      I      K      L       M      N
#>    <dbl>   <dbl>   <dbl>  <dbl>   <dbl>  <dbl>  <dbl>  <dbl>   <dbl>  <dbl>
#> 1  0.796  0.251   1.54   -0.261 -0.238   0.280 -0.115 -1.18  -0.0135  0.884
#> 2  1.74   1.25   -2.04    0.962 -1.56   -1.30  -1.11  -1.33   1.25   -1.02 
#> 3 -1.96  -1.33   -0.235  -1.26  -0.238   1.70  -1.70   0.503  1.25    1.22 
#> 4 -0.637 -0.0123 -1.46   -0.261  0.0248  1.15  -1.11  -0.590  0.665   1.06 
#> 5 -2.37  -0.334  -0.720  -0.846 -0.560  -1.30  -2.70  -1.04  -1.34    4.23 
#> 6  0.626 -0.0123  0.417   0.154 -0.975   0.280 -0.530 -0.105  0.987  -0.701
#> 7  0.363 -0.334  -0.0420 -0.846  0.247   0.695 -1.38  -1.33  -1.34   -0.216
#> 8 -0.152  0.836   1.48   -2.85   0.440  -0.135 -0.699 -0.327 -0.335  -0.701
#> 9 -1.37  -0.0123  0.128  -3.85  -2.56    2.55  -2.70   1.77   2.57   -2.02 
#> # ℹ 10 more variables: Q <dbl>, R <dbl>, S <dbl>, T <dbl>, V <dbl>, W <dbl>,
#> #   Y <dbl>, E <dbl>, P <dbl>, D <dbl>
```

The funciton `create_pwm_df` generates a PWM for all allele length
combinations for which at least `min_datapoints` are available.

``` r
pwm_df <- example_epitopes_df %>%
  create_pwm_df(min_datapoints=50, background=example_background_df)

pwm_df
#> # A tibble: 4 × 4
#>     len allele     n pwm              
#>   <dbl> <chr>  <int> <list>           
#> 1     8 H2Kb     217 <tibble [8 × 20]>
#> 2     9 H2Db     207 <tibble [9 × 20]>
#> 3     9 H2Kb     110 <tibble [9 × 20]>
#> 4     9 H2Kd     116 <tibble [9 × 20]>
```

Peptides can be scored based on a given PWM dataframe:

``` r
example_peptide = "ASNENAETM"
score_seq_allele(example_peptide, "H2Db", pwm_df)
#> [1] 0.7596006
```

For a given frequency matrix, the information content at the different
positions can be determined by using the function
`calculate_information_content`.

``` r
example_epitopes_df %>%
  filter(allele == "H2Db", len==9) %>%
  pull(peptide) %>% 
  seqs2pfm() %>% 
  calculate_information_content()
#> [1] 0.4061437 0.8142406 0.6516111 0.2584148 2.1069134 0.1866548 0.3894873
#> [8] 0.3018932 1.4782398
```

# Developer information

Before commiting changes update the documentation and this README:

``` r
devtools::document()
devtools::build_readme()
```
