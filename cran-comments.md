This is a submission (first release)

---

## Test environments

* Mac (travis-ci): R-release 3.4.3 [12/18/17]
* Mac OS X El Capitan (local machine): R 3.4.0 [12/18/17]
* Linux (travis-ci): R-release 3.4.2 [12/18/17]
* ubuntu 17.10 (local machine): R 3.3.2 [12/18/17]
* win-builder: R-release 3.4.3, R-dev [12/18/17]

## R CMD check results

Mac, Linux:

0 errors | 0 warnings | 0 notes

Win:

0 errors | 0 warnings | 0 notes

Dev versions:

* Seemingly have issues with igraph and phyloseq installations. Both of these
  dependencies are well maintained and shouldn't be an issue when new R versions
  release.
  
## Responses to Previous Submission [12/18/17]

* Removed Akima depends due to licensing issues; it's not longer used in the pkg
* Vignette takes ~15s to knit; tests take ~47s to complete

## Responses to Previous Submission [12/20/17]

* Link in description now works
* No code is now run in vignette builds -- to improve check time
* Removes set of tests that generate figures -- to improve check time
* Linux check time: 247s
* Travis: linux job time 396s, mac job time 1018s
* Win-builder: install time 25s, check time 387s

## License 

* License components with restrictions and base license permitting such:
  MIT + file LICENSE
  
## Downstream dependencies

* None (first release)