# mashr-paper

This repository contains R source code implementing Empirical Bayes
methods for simultaneously testing many effects in many conditions or
outcomes. For more information on these methods and their application
to the the
[Genotype-Tissue Expression (GTEx) study](https://www.gtexportal.org),
please see [Urbut *et al*, 2017](https://doi.org/10.1101/096552).

*Note:* The primary purpose of this repository is to implement the
analyses presented in the Urbut *et al* manuscript. We are also
developing an R package for broader application, and this R package
can be found [here](https://github.com/stephenslab/mashr).

## Using the code

Although this repository has the standard structure of an
[R package](http://r-pkgs.had.co.nz/package.html), this package is in
development, so it is recommended that you load the function
definitions directly into your R environment rather than attempting to
install this repository as a package. For example, this can be done by
cloning or downloading this repository, setting your working directory
to this repository on your computer, then running the following
commands in R:

```r
library(devtools)
load_all(export_all = TRUE)
```

This will load all the mashr functions into your environment without
actually installing the package.

To start, we recommend walking through
[a small demonstration of mashr on simulated
data](inst/demos/Reprosim/newsim.Rmd), as well as a
[more advanced demonstration of mashr for estimating tissue-specific
effects on gene expression](inst/demos/Advanced/TissueSpecificVignette.Rmd).

## Overview

A typical mashr analysis will include the following steps:

+ Generate a list of covariance matrices to capture patterns of
sharing in the data. This may be done with predetermined matrices, or
using a method that adapts the matrices to the data (e.g, using
[Extreme Deconvoluion](https://github.com/jobovy/extreme-deconvolution).

+ Estimate weights reflecting the relative frequency of each pattern
of sharing in the data.

+ (Optionally) Assess the model fit on a held-out test set, which is
useful for comparing across models.

+ Compute the J x P conditional likelihood matrix, where J is the
number of estimated effects and P is the number of mixture
components.
   
+ Compute intermediate posterior quantities; e.g., means, tail
probabilities, marginal subgroup-specific posterior variances.

+ Compute posterior quantities of interest; e.g., posterior mean,
local false sign rate, marginal variance for the effect size in
each condition.

## License

Copyright (C) 2017, Sarah Urbut.

This program is free software: you can redistribute it and/or modify
it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html) as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
**without any warranty**; without even the implied warranty of
**merchantability** or **fitness for a particular purchase**. See the
[GNU General Public License](LICENSE) for more details.
