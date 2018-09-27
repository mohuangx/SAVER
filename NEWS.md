## SAVER 1.1.1 (September 26, 2018)

* Fix bug when inputted mu values are far from observed
* Add `get.mu` function to get prior mean predictions.

## SAVER 1.1.0 (September 2, 2018)

* Improve finish time estimation
* Allow input of prior predictions
* Can set `ncores` to set up parallel backend
* Update documentation
* Make SAVER tutorial vignette
* Allow input of sparse matrices
* Allow return of estimates matrix only

## SAVER 1.0.0 (April 25, 2018)

* First official release

## SAVER 0.4.0 (February 28, 2018)
  
* Implemented faster version of SAVER which is 20-30x faster
* No longer necessary to specify ```parallel=TRUE```
* Implemented progress bar

## SAVER 0.3.1 (January 29, 2018)

* Incorporate `iterators` package
* Implement mean expression cutoff for prediction matrix
* Fix minor bugs

## SAVER 0.3.0 (September 19, 2017)

* Combine parallel output using `unlist` instead of `Reduce`. 
  * `Reduce` is computationally intensive. Might cause issue #1.
* Add nlambda option to `saver`
* Incorporate `bigmemory` package
* Export `cor.cells`


## SAVER 0.2.2 (August 13, 2017)

* Rewrote combine function of parallel implementation of `saver`
* Add null model option to `saver`

## SAVER 0.2.1 (August 1, 2017)

* Rewrote parallel implementation of SAVER
* Added verbose, pred.cells, and predict.time options to `saver`

## SAVER 0.2.0 (July 23, 2017)

* Addressed issue #1
* Add option to filter genes out for prediction in `saver`  
* Reduced `saver` output
* Created correlation adjustment functions `cor.genes` and `cor.cells`


## SAVER 0.1.3 (July 18, 2017)

* Addressed issue #1 
* Remove requirement of `saver` object
* Created `sample.saver` function
* More accurate calculation of finish time
* Update version number

## SAVER 0.1.2 (May 29, 2017)

* Added `sample.saver` to sample datasets from the posterior distribution.
* Documented dataset and package.
* Included size factor, number of genes and cells, and names to `saver` output.
* Created `saver` class.

## SAVER 0.1.1 (May 9, 2017)

Prerelease version with documentation.



