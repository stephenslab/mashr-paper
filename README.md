Welcome to Multivariate Adaptive Shrinkage (R package)!

The package implements methods to estimate and test many effects in many conditions (or many effects on many outcomes).

The methods use Empirical Bayes methods to estimate patterns of similarity among conditions, and then exploit
those patterns of similarity among conditions to improve accuracy of effect estimates.
See (preprint)[preprint link to appear] for further details.

To install:

```
devtools::install_github("stephenslab/mashr")
library("mashr")
```

See the Manual and vignettes: [Repo Sim](Vignettes/Reprosim/newsim.html) and [Tissue Specific](Vignettes/Advanced/TissueSpecificVignette.html) in the Vignettes for package execution and more details in the Advanced branch on comparing between methods.

A typical `mashr` analysis will involve the following steps:

* Generate a list of covariance matrices which aim to capture all the patterns of sharing in the data. This may be done in a fixed or iteratively 'learned' fashion, with or without the `ExtremeDeconvolution` package.

* Estimate weights reflecting the relative frequency of each pattern of sharing in the data. 

* (Optionally) Assess the model fit on a held-out test set (useful for comparing across models).
	* Compute its likelihood as a *J effects x P components* matrix of likelihoods at each component
   
   * Posterior quantities of interest (means, tail probabilities and marginal subgroup specific posterior variances) stored in a *J effects x P components x R subgroup* array.

* Compute the posterior quantities of interest: posterior mean, local false sign rate and marginal variance for the effect size in each condition.

