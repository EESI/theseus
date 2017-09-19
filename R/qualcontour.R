#' @import ggplot2
NULL

#' Plots quality assessment
#'
#' This function generates a 2-D contour map representing the location of quality scores for a designated percentile. It is intended to assist with deciding where trimming should be performed.
#'
#' @param f_path (required) A character vector locating the forward read (Read 1) .fastq files
#' @param r_path (required) A character vector locating the reverse read (Read 2) .fastq files
#' @param idx Indexes (within f_path and r_path) identifying specific .fastq files for analysis
#' @param percentile The percentile to be targeted . Defaults to .25 (i.e. the first quartile).
#' @param n_trim asdf. Defaults to 100.
#' @param n_samples Integer indicating the number of samples to include in the visualization. Defaults to 12.
#' @param q A numeric vector designating Phred quality scores to be represented on the plot. Defaults to 25, 30, and 35.
#' @param bins Integer designating the number of bins each read should be seperated into. For example, visualizing a 250 bp read with 50 bins would imply that each bin represents 5 cycles/bp. Increasing the number of bins improves granularity at the cost of memory and processing speed. Defaults to 50.
#' @param nc The number of cores to use when multithreading. Defaults to 1.
#' @param seed An integer value to be used when randomly selecting the subset of samples to be visualized.
#' @param verbose If set to TRUE, provides verbose output. Defaults to FALSE.
#'
#' @return A ggplot object with the following attributes:
#' \describe{
#' \item{idx}{Samples used to generate the plot.}
#' \item{seed}{Seed used to select the samples used to generate the plot.}
#' }
#'
#' @references
#' X
#' Y
#' Z
#'
#' @seealso \code{\link[ShortRead]{qa}} \code{\link[dada2]{plotQualityProfile}}
#'
#' @examples
#' \dontrun{
#' fns <- sort(list.files(data_path,full.names=TRUE))
#' f_path <- fns[grepl('R1', fns)]
#' r_path <- fns[grepl('R2', fns)]
#' p1 <- qualcontour(f_path,r_path,n_samples=12,verbose=TRUE,percentile=.25,nc=1)
#' p1 + geom_hline(yintercept=175) + geom_vline(xintercept=275)
#' }
#'
#' @export

qualcontour <- function(f_path,r_path,idx,percentile=.25,n_trim=100,n_samples=12,q=c(25,30,35),bins=50,nc=1,
                 seed=sample.int(.Machine$integer.max,1),verbose=FALSE){

  set.seed(seed)

  if (missing(idx)){
    if (verbose) message(sprintf('No sequence file indexes provided. Randomly selecting %s sequence pairs.',n_samples))
    idx <- sample(seq_along(f_path),n_samples)
  }

  if (verbose) message(sprintf('Calculating scores for %s sequence pairs: ',n_samples))
  scores <- parallel::mclapply(seq_along(idx),function(i){
    q_f <- qual_scores(f_path[i],percentile=percentile)
    q_r <- qual_scores(r_path[i],percentile=percentile,forward=FALSE)
    if (verbose) cat('.')
    (q_f + q_r)/2
  },mc.cores=nc)

  len <- nrow(scores[[1]])

  contours <- lapply(scores,function(s){
    cont <- grDevices::contourLines(1:nrow(s),1:nrow(s),s,levels=q)
    cont <- do.call(rbind,lapply(cont,function(x) data.frame(x=x$x,y=x$y,level=x$level)))
  })
  contours <- do.call(rbind,contours)


  df_trim <- data.frame(x=0:n_trim,y=(0:n_trim)*-1 + n_trim)

  suppressWarnings(
    p1 <- ggplot(contours,aes_(~x,~y,alpha=~..count..,fill=~as.factor(level))) +
      geom_bin2d(bins=bins) +
      scale_fill_brewer(palette='Set1') +
      geom_abline(slope=-1,intercept=n_trim,linetype=3) +
      geom_ribbon(data=df_trim,aes(x,ymin=1,ymax=y,count=0),fill='black',alpha=.1) +
      theme_classic() +
      theme(aspect.ratio=1,axis.text.x=element_text(angle=90,vjust=.5,hjust=1)) +
      scale_x_continuous(expand=c(0,0),limits=c(-1,len),breaks=seq(0,len,50)) +
      scale_y_continuous(expand=c(0,0),limits=c(-1,len),breaks=seq(0,len,50)) +
      guides(alpha=FALSE) +
      labs(x='Forward',y='Reverse',fill='Q')
  )

  attr(p1,'seed') <- seed
  attr(p1,'idx') <- idx

  p1

}