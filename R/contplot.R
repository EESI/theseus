#' @import ggplot2
NULL

#' Plots a color/contour plot of a three dimensional surface
#'
#' Function contplot fits a surface through a 3-dimensional irregularly-spaced dataset and presents the fitted surface as a 2-dimensional color/contour plot. 'linear', 'spline', and 'loess' fitting methods are currently supported for this function. Please note that missing values in 'data' are not tolerated (see Details). 
#'
#' @param x (required) Vector of x-coordinates (within 'data').
#' @param y (required) Vector of y-coordinates (within 'data').
#' @param z (required) Vector of z-coordinates (within 'data').
#' @param data (required) A dataframe-like object containing x, y, and z.
#' @param method Surface fitting method. Options are 'linear', 'spline', and 'loess' (see Details). Defaults to 'linear'.
#' @param removeMissing Remove entries with missing values (see Details). Defaults to FALSE.
#' @param ... Additional arguments for loess.
#'
#' @details
#' This function can be used for both regularly and irregularly spaced data.  
#' \subsection{'method' - surface fitting methods}{
#' \itemize{
#'   \item 'linear' uses the \code{\link[akima]{interp}} function in the akima R package (Akima and Gebhardt, 2016) to carry out bivariate linear intrpolation. 
#'   \item 'spline' also uses the \code{\link[akima]{interp}} function in the akima R package (Akima and Gebhardt, 2016), but the 'linear' argument passed to \code{\link[akima]{interp}} is set to 'FALSE', resulting in bicubic spline interpolation. 
#'   \item 'loess' uses the \code{\link[stats]{loess}} (R Core Team, 2017) and \code{\link[splancs]{inpip}} (Rowlingson and Diggle, 2017) functions to carry out locally weighted scatterplot smoothing (LOESS).
#' }
#' }
#' 
#' \subsection{Missing Values}{
#' Missing values in 'x', 'y', or 'z' are not tolerated by these functions and will result in an error. Setting 'removeMissing' to TRUE will remove any and all entries that contain missing values before interpolation is carried out.
#' }
#' 
#' @return A ggplot object.
#'
#' @references
#' Akima, H., and Gebhardt, A. (2016). akima: Interpolation of Irregularly and Regularly Spaced Data. R package version 0.6-2. \url{https://CRAN.R-project.org/package=akima}.  
#' R Core Team (2017). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. \url{https://www.R-project.org/}.
#' Rowlingson, B., and Diggle, P. (2017). splancs: Spatial and Space-Time Point Pattern Analysis. R package version 2.01-40. \url{https://CRAN.R-project.org/package=splancs}. 
#'
#' @seealso \code{\link[stats]{loess}} \code{\link[splancs]{inpip}} \code{\link[akima]{interp}}
#'
#' @examples
#' \dontrun{
#' library(theseus)
#' data(akima, package='akima')
#' str(akima)
#' dat <- data.frame(cbind(a=akima$x, b=akima$y, c=akima$z))
#' p.cp.l <- contplot(x='a', y='b', z='c', data=dat, method='linear')
#' p.cp.l + ggtitle("contour plot with linear surface fitting")
#' p.cp.s <- contplot(x='a', y='b', z='c', data=dat, method='spline')
#' p.cp.s + ggtitle("contour plot with spline surface fitting")
#' p.cp.loess <- contplot(x='a', y='b', z='c', data=dat, method='loess')
#' p.cp.loess + ggtitle("contour plot with loess surface fitting")
#' }
#'
#' @export

contplot <- function(x,y,z,data,method=c('linear','spline','loess'),removeMissing=FALSE,...){

  method <- match.arg(method)

  if (!is.data.frame(data))
    data <- data.frame(data,stringsAsFactors=FALSE)

  if (removeMissing==TRUE){
    data <- subset(data,!(is.na(x) | is.na(y) | is.na(z)))
  }

  if (method == 'loess'){
    form <- stats::as.formula(sprintf('%s ~ %s * %s',z,x,y))
    lo <- stats::loess(form,data,...)
    lo_x <- seq(min(data[,x]),max(data[,x]),length=length(data[,x]))
    lo_y <- seq(min(data[,y]),max(data[,y]),length=length(data[,y]))
    lo_grid <- expand.grid(x=lo_x,y=lo_y)
    colnames(lo_grid) <- c(x,y)
    lo_z <- stats::predict(lo,newdata=lo_grid)
    df <- data.frame(lo_grid,z=matrix(lo_z))
    colnames(df) <- c('x','y','z')

    df <- df[suppressWarnings(splancs::inpip(df[,1:2],data[grDevices::chull(data[,x],data[,y]),1:2])),]
  }

  if (method == 'linear' | method == 'spline'){

    if (method == 'linear') linear <- TRUE else linear <- FALSE

    S <- akima::interp(x=data[,x],y=data[,y],z=data[,z],
                       xo=seq(min(data[,x]),max(data[,x]),length=length(data[,x])),
                       yo=seq(min(data[,y]),max(data[,y]),length=length(data[,y])),
                       linear=linear)

    df <- with(S,data.frame(x=x,y=rep(y,each=ncol(z)),z=matrix(z))) %>%
      tidyr::drop_na()

  }

  ggplot() +
    geom_tile(data=df,aes(x=x,y=y,fill=z),alpha=.7) +
    geom_point(data=data,aes_string(x=x,y=y)) +
    viridis::scale_fill_viridis() +
    theme_classic() +
    theme(aspect.ratio=1) +
    labs(x=x,y=y,fill=z)

}
