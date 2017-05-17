# SAVER

SAVER (Single-cell Analysis Via Expression Recovery) implements a regularized regression prediction and empirical Bayes method to recover the true gene expression profile in noisy and sparse single-cell RNA-seq data.

## Installation

You can install SAVER from github with:

```R
# install.packages("devtools")
devtools::install_github("mohuangx/SAVER")
```

## Example

We will demonstrate the ```saver``` function using the ```linnarsson``` dataset.

```R
data("linnarsson")

# predictions for top 5 highly expressed genes
saver1 <- saver(linnarsson, npred = 5)

# predictions for certain genes
genes <- c("Thy1", "Mbp", "Stim2", "Psmc6", "Rps19")
genes.ind <- which(rownames(linnarsson) %in% genes)
saver2 <- saver(linnarsson, pred.genes = genes.ind)

# run predictions for all genes
saver3 <- saver(linnarsson)

# run in parallel
library(doParallel)
registerDoParallel(cores = 4)
saver4 <- saver(linnarsson, parallel = TRUE)
```

Recovered normalized expression estimates are stored in the ```estimate``` element of the list returned by ```saver```.

Library-size normalization is performed by default. If the input expression matrix is already normalized or normalization is not desired, use the option ```size.factor = 1```. A vector of other size factors can also be assigned to the ```size.factor``` argument.
