## Test environments

* local OS X install, R 3.5.0
* ubuntu 14.04 (on travis-ci), R 3.5.0 and R-devel
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs or NOTES.

## Reverse Dependencies

None

## Features

* added S3 method `variable.names` to extract variable names in the highest probability model, median probability
model, and best probability model for objects created by `predict`.

## Bugs

Fixed incorrect documentation in `predict.basglm` which had that  `type = "link"` was the default for prediction [issue #18](https://github.com/merliseclyde/BAS/issues/18)


## Comments

added example in vignette may take time to build  (total check time < 15 minutes)

