Welcome to Matrix Ash!


To install:

```
devtools::install_github("surbut/matrix_ash")
library("mash")

```

The intention of this package is to produce an estimate of the posterior effect of a SNP on gene expression across multiple subgroups by modeling these effects as distributed according to a mixture of multivariate normals where each component of the mixture is specified by its prior covariance matrix specifying the relationship in effect sizes between tissues. 
 

```
 $\bm{b} | \bm{\pi},\bf{U} \sim \sum_{k,l} \pi_{k,l} \;{\it N}_R(\bm{0}, \omega_l U_{k})$
```


The novelty of our approach is to use a list of data-sensitive covariance matrices which aim to recapitulate these patterns of sharing. Please see the Manual and vignettes (compute.per.snp and separatingtrainingandtesting) for package execution and information on comparing between methods.

Essentially, matrix ash works in a modular fashion 

1) Generate a list of covariance matrices which aim to catpure all the patterns of sharing in the data. This may be done in a fixed or iteratively 'learned' fashion.

2) Compute the likelihood of a training data set to produce a set of hierarchical weights reflecting the relative frequency of each pattern of sharing in the data. The choice of largely null data will influence the shrinkage of the data set.

2) Compute the mixture distribution on the test set:
   a) Compute its likelihood as a J x K matrix of likelihoods at each component
   b) Posterior quantities of interest (means, tail probabilities and covariance matrices) stored in a J snps x K components x R subgroup array.

3) Using the prior weights computed in step (1) and the  likelihood computed in (2a) compute the posterior weighted quantities to understand the posterior mean, local false sign rate and marginal variance for the effect size in each tissue.

Please see extensive documentation in the manual and Documentation section for details on separate functions and derivation.
