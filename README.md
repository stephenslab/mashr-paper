Welcome to Matrix Ash!


To install:

```
devtools::install_github("surbut/matrix_ash")
library("mash")

```

The intention of this package is to produce an estimate of the posterior effect of a SNP on gene expression across multiple subgroups by modeling these effects as distributed according to a mixture of multivariate normals where each component of the mixture is specified by its prior covariance matrix specifying the relationship in effect sizes between tissues. 
 

```
 g(\cdot; \bm{\pi},\bm{U}) = \sum_{k=1}^K \sum_{l=1}^L \pi_{k,l} \; N_R(\cdot;\bm{0}, \omega_l U_{k}).

```


The novelty of our approach is to use a list of data-sensitive covariance matrices which aim to recapitulate these patterns of sharing. Please see the Manual and vignettes: [Repo Sim](Vignettes/Reprosim/newsim.html) and [Tissue Specific](Vignettes/Advanced/TissueSpecificVignette.html) in the Vignettes for package execution and more details in the Advanced branch on comparing between methods.

Essentially, matrix ash (`mash`) works in a modular fashion 

* Generate a list of covariance matrices which aim to capture all the patterns of sharing in the data. This may be done in a fixed or iteratively 'learned' fashion, with or without the `ExtremeDeconvolution` package.

* Compute the likelihood of a training data set to produce a set of hierarchical weights reflecting the relative frequency of each pattern of sharing in the data. The choice of data will influence the appropriate `shrinkage` of the data set.

* Compute the mixture distribution on the test set:
	* Compute its likelihood as a *J effects x P components* matrix of likelihoods at each component
   
   * Posterior quantities of interest (means, tail probabilities and marginal subgroup specific posterior variances) stored in a *J effects x P components x R subgroup* array.

* Using the prior weights computed in step (1) and the  likelihood computed in (2a) compute the posterior weighted quantities to understand the posterior mean, local false sign rate and marginal variance for the effect size in each tissue.

Please see extensive documentation in the manual and Urbut _et al 2016_ for details on specific functions and derivation.