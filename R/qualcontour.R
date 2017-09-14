#' @import ggplot2
NULL

#' Plots quality assessment
#'
#' Description paragraph
#'
#' @param f_path (required) asdf
#' @param r_path (required) asdf
#' @param idx asdf
#' @param percentile asdf. Defaults to .25.
#' @param n_trim asdf. Defaults to 100.
#' @param n_seqs asdf. Defaults to 12.
#' @param q asdf. Defaults to 25, 30, and 35.
#' @param bins adsf. Defaults to 50.
#' @param nc asdf. Defaults to 1.
#' @param seed asdf.
#' @param verbose asdf. Defaults to FALSE.
#'
#' @return A ggplot object with the following attributes:
#' \describe{
#' \item{idx}{asdf}
#' \item{seed}{asdf}
#' }
#'
#' @references
#' X
#' Y
#' Z
#'
#' @seealso \code{\link[ShortRead]{qa}}
#'
#' @examples
#' \dontrun{
#' fns <- sort(list.files(data_path,full.names=TRUE))
#' f_path <- fns[grepl('R1', fns)]
#' r_path <- fns[grepl('R2', fns)]
#' p1 <- qual(f_path,r_path,n_seqs=12,verbose=TRUE,percentile=.25,nc=12)
#' p1 + geom_hline(yintercept=175) + geom_vline(xintercept=275)
#' }
#'
#' @export

qual <- function(f_path,r_path,idx,percentile=.25,n_trim=100,n_seqs=12,q=c(25,30,35),bins=50,nc=1,
                 seed=sample.int(.Machine$integer.max,1),verbose=FALSE){

  set.seed(seed)

  if (missing(idx)){
    if (verbose) message(sprintf('No sequence file indexes provided. Randomly selecting %s sequence pairs.',n_seqs))
    idx <- sample(seq_along(f_path),n_seqs)
  }

  if (verbose) message(sprintf('Calculating scores for %s sequence pairs: ',n_seqs))
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
    p1 <- ggplot(contours,aes(x,y,alpha=..count..,fill=as.factor(level))) +
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