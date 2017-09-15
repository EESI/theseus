# test qualitycontour (old qual)
#library(theseus)
# fns<-c(
#   file.path("/Users/jprice/Desktop","sample1_R1.fastq.gz"),
#   file.path("/Users/jprice/Desktop","sample1_R2.fastq.gz"),
#   file.path("/Users/jprice/Desktop","sample2_R1.fastq.gz"),
#   file.path("/Users/jprice/Desktop","sample2_R2.fastq.gz")
# )
library(theseus)
library(ggplot2)
getwd()
fns<-sort(
  list.files(
    # not sure how you want to access the file list in this directory. 
    file.path(getwd(), "inst/extdata"), 
    #system.file("extdata", package="theseus"), 
  full.names = TRUE)
)

fns
# added ".fastq.gz" in the event we have more than just sequences in this directory
f_path <- fns[grepl('R1.fastq.gz', fns)] 
r_path <- fns[grepl('R2.fastq.gz', fns)]
# p1 <- theseus::qual(f_path,r_path,n_seqs=2,verbose=TRUE,percentile=.25,nc=1)
p1 <- theseus::qualcontour(f_path,r_path,n_samples=2,verbose=TRUE,percentile=.25,nc=1)
p1 + geom_hline(yintercept=175) + geom_vline(xintercept=275)
