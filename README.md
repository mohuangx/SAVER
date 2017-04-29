# SAVER

SAVER (Single-cell Analysis Via Expression Recovery) implements a regularized regression prediction and empirical Bayes method to recover the true gene expression profile in noisy and sparse single-cell RNA-seq data.

## Installation

You can install SAVER from github with:

```R
# install.packages("devtools")
devtools::install_github("SAVER/mohuangx")
```

## Example

We will demonstrate the ```saver``` function using the ```linnarsson``` dataset.

```R
data("linnarsson")
linnarsson.saver <- saver(linnarsson)
```

SAVER estimates are stored in ```linnarsson.saver$estimate```.
