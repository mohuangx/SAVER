# SAVER

SAVER (Single-cell Analysis Via Expression Recovery) implements a regularized regression prediction and empirical Bayes method to recover the true gene expression profile in noisy and sparse single-cell RNA-seq data.

## News and Updates
April 25, 2018
* Version 1.0.0 released
  * First official release

February 28, 2018
* Version 0.4.0 released
  * Implemented faster version of SAVER which is 20-30x faster
  * No longer necessary to specify ```parallel=TRUE```
  * Implemented progress bar

January 29, 2018
* Version 0.3.1 released.

September 19, 2017
* Version 0.3.0 released.

August 13, 2017
* Version 0.2.2 released.

## Installation

You can install the most recent updates of SAVER from github with:

```R
# install.packages("devtools")
devtools::install_github("mohuangx/SAVER")
```
To install the stable release of SAVER:

```R
# install.packages("devtools")
devtools::install_github("mohuangx/SAVER@*release")
```

## saver Example

We will demonstrate the ```saver``` function using the ```linnarsson``` dataset.

```R
data("linnarsson")

# predictions for top 5 highly expressed genes
saver1 <- saver(linnarsson, npred = 5)

# predictions for certain genes
genes <- c("Thy1", "Mbp", "Stim2", "Psmc6", "Rps19")
genes.ind <- which(rownames(linnarsson) %in% genes)
saver2 <- saver(linnarsson, pred.genes = genes.ind)

# only return output for predicted genes
saver3 <- saver(linnarsson, pred.genes = genes.ind, pred.genes.only = TRUE)

# run predictions for all genes
saver4 <- saver(linnarsson)

# run in parallel
library(doParallel)
cl <- makeCluster(4, outfile = "")
registerDoParallel(cl)
saver5 <- saver(linnarsson)
stopCluster(cl)
```

Recovered normalized expression estimates are stored in the ```estimate``` element of the returned ```saver``` object.

Library-size normalization is performed by default. If the input expression matrix is already normalized or normalization is not desired, use the option ```size.factor = 1```. A vector of other size factors can also be assigned to the ```size.factor``` argument.

## Combining SAVER objects

Running SAVER on large datasets may run into memory issues, especially when parallelized across cores where each core has limited memory. One solution would be to register a parallel backend using SOCK or MPI to parallelize across nodes. Another solution would be to split the genes into multiple runs using `pred.genes` and `pred.genes.only`. For example, say we have a dataset `x` with 10,000 genes. We can split this up into 4 runs on different machines and combine using the `combine.saver` function:

```R
saver1 <- saver(x, pred.genes = 1:2500, pred.genes.only = TRUE)
saver2 <- saver(x, pred.genes = 2501:5000, pred.genes.only = TRUE)
saver3 <- saver(x, pred.genes = 5001:7500, pred.genes.only = TRUE)
saver4 <- saver(x, pred.genes = 7501:10000, pred.genes.only = TRUE)

saver.all <- combine.saver(list(saver1, saver2, saver3, saver4))
```

## Sampling example

Oftentimes for downstream analysis, you may want to sample from the posterior distributions to account for estimation uncertainty. This is similar to performing multiple imputation to obtain multiple imputed datasets. To sample from the posterior distributions, we use the function ```sample.saver```, which takes in the ```saver``` result and outputs sampled datasets. For example,

```R
samp1 <- sample.saver(saver1, rep = 1, seed = 50)
samp5 <- sample.saver(saver1, rep = 5, seed = 50)
```

## Correlation example

Because the SAVER estimates contain uncertainty, correlations between genes and cells cannot be directly calculated using the SAVER estimates. To adjust for the uncertainty, we provide the functions ```cor.genes``` and ```cor.cells``` to construct gene-to-gene and cell-to-cell correlation matrices respectively for the SAVER output. These functions take in the ```saver``` result and outputs the gene-to-gene or cell-to-cell correlation matrix. For example,

```R
saver1.cor.gene <- cor.genes(saver1)
saver1.cor.cell <- cor.cells(saver1)
```




