## Test environments
* local Windows 7 install, R 3.5.1
* ubuntu 14.04.5 (on travis-ci), R 3.5.1

## R CMD check results

There were no ERRORs, WARNINGs

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Mo Huang <mohuangx@gmail.com>'

New submission

* checking examples ... NOTE
Examples with CPU or elapsed time > 5s
               user system elapsed
saver         41.70   1.11   42.87
cor_adjust    20.47   0.55   21.03
combine.saver 18.23   0.81   19.08
sample.saver  17.03   0.43   17.49
get.mu         8.82   0.25    9.08

Practical examples were difficult to make shorter.
