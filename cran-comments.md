This is a submission (first release)

---

## Test environments
* Mac (travis-ci): R-release 3.4.3 [12/18/17]
* Mac OS X El Capitan (local machine): R 3.4.0 [12.18.17]
* Linux (travis-ci): R-release 3.4.2 [12/18/17]
* ubuntu 17.10 (local machine): R 3.3.2 [12.18.17]
* win-builder: R 3.4.3 [12/18/17]

## R CMD check results

Mac, Linux:

0 errors | 0 warnings | 0 notes

Win:

0 errors | 0 warnings | 1 notes

* Package has a FOSS license but eventually depends on the following package 
  which restricts use:
    + akima
* Response: akima has an ACM license which grants noncommercial use:
    + https://www.acm.org/publications/policies/software-copyright-notice
    
Dev version:

* Seemingly have issues with igraph and phyloseq installations. Both of these
  dependencies are well maintained and shouldn't be an issue when new R versions
  release.

## License 

* License components with restrictions and base license permitting such:
  MIT + file LICENSE
  
## Downstream dependencies

* None (first release)