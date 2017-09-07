#' @import ggplot2
NULL

#' Plots quality assessment
#'
#' Description paragraph
#'
#' @param x (required) asdf
#' @param y (required) asdf
#' @param z (required) asdf
#' @param data (required) asdf
#' @param method asdf. Defaults to linear.
#' @param ... Additional arguments for loess.
#'
#' @return A ggplot object.
#'
#' @references
#' X
#' Y
#' Z
#'
#' @seealso \code{\link[stats]{loess}} \code{\link[splancs]{inpip}} \code{\link[akima]{interp}}
#'
#' @examples
#' \dontrun{
#' cont('xvar','yvar','zvar',dat,method='loess',degree = 2, span = 0.1)
#' }
#'
#' @export

cont <- function(x,y,z,data,method=c('linear','spline','loess'),...){

  method <- match.arg(method)

  if (!is.data.frame(data))
    data <- data.frame(data,stringsAsFactors=FALSE)

  if (method == 'loess'){
    form <- as.formula(sprintf('%s ~ %s * %s',z,x,y))
    lo <- stats::loess(form,data,...)
    lo_x <- seq(min(data$x),max(data$x),length=length(data$x))
    lo_y <- seq(min(data$y),max(data$y),length=length(data$y))
    lo_grid <- expand.grid(x=lo_x,y=lo_y)
    colnames(lo_grid) <- c(x,y)
    lo_z <- predict(lo,newdata=lo_grid)
    df <- data.frame(lo_grid,z=matrix(lo_z))
    colnames(df) <- c('x','y','z')

    df <- df[suppressWarnings(splancs::inpip(df[,1:2],data[chull(data$xvar,data$yvar),1:2])),]
  }

  if (method == 'linear' | method == 'spline'){

    if (method == 'linear') linear <- TRUE else linear <- FALSE

    S <- akima::interp(x=data$x,y=data$y,z=data$z,
                       xo=seq(min(data$x),max(data$x),length=length(data$x)),
                       yo=seq(min(data$y),max(data$y),length=length(data$y)),
                       linear=linear)

    df <- with(S,data.frame(x=x,y=rep(y,each=ncol(z)),z=matrix(z))) %>%
      drop_na()

  }

  ggplot() +
    geom_tile(data=df,aes(x=x,y=y,fill=z),alpha=.7) +
    geom_point(data=data,aes_string(x=x,y=y)) +
    viridis::scale_fill_viridis() +
    theme_classic() +
    theme(aspect.ratio=1) +
    labs(x=x,y=y,fill=z)

}