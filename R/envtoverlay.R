#' @import ggplot2
NULL

#' Environmental variable fitting to unconstrained ordination diagrams
#'
#' Fits environmental variables as vectors (via \code{\link[vegan]{envfit}}) and
#' smooth surfaces (via \code{\link[vegan]{ordisurf}}) to an ordination diagram.
#' The figure is faceted if multiple variables are specified.
#'
#' @param PS (required) A phyloseq object.
#' @param covariates (required) A character vector of covariates present in the
#'   phyloseq objects sample_data().
#' @param ordmet Ordination method. Options are Principal Component Analysis
#'   ("PCA") or Correspondence Analysis ("CA"). Defaults to "PCA".
#'
#' @return A ggplot object.
#'
#' @seealso \code{\link[vegan]{rda}} \code{\link[vegan]{cca}}
#'   \code{\link[vegan]{envfit}}
#'   \code{\link[vegan]{ordisurf}}
#'
#' @examples
#' \dontrun{
#' library(theseus)
#' library(phyloseq)
#' library(ggplot2)
#' data('WWTP_Impact')
#' cv <- c('log_NO3N', 'log_PO4')
#' p.eo <- envtoverlay(WWTP_Impact, covariates=cv)
#' p.eo
#' }
#'
#' @export

envtoverlay <- function(PS, covariates, ordmet='PCA'){

  graphics::plot.new()

  OTU <- PS@otu_table@.Data
  SAMP <- suppressWarnings(dplyr::as_tibble(PS@sam_data))

  if (!all(covariates %in% colnames(SAMP)))
    stop('Some covaraites not found in sample data!')

  if (PS@otu_table@taxa_are_rows) OTU <- t(OTU)

  if (ordmet=='PCA') {
    R <- vegan::rda(OTU)
  } else if (ordmet=='CA') {
    R <- vegan::cca(OTU)
  }

  form <- stats::as.formula(sprintf('%s ~ %s','R',paste(covariates,collapse=' + ')))
  env <- vegan::envfit(form,SAMP,na.rm=TRUE)

  v <- (sqrt(env$vectors$r) * env$vectors$arrows)[,1:2,drop=FALSE]

  df <- lapply(covariates,function(k){

    form <- stats::as.formula(sprintf('%s ~ %s','R',k))
    O <- vegan::ordisurf(form,SAMP,plot=FALSE)

    df1 <- with(O,data.frame(x=grid$x,y=rep(grid$y,each=ncol(grid$z)),z=matrix(grid$z))) %>%
      dplyr::filter_(~!is.na(z)) %>%
      dplyr::mutate_(covariate=~k)

    df2 <- with(O,data.frame(x=model$x1,y=model$x2,z=NA)) %>%
      dplyr::mutate_(covariate=~k)

    list(df1=df1,df2=df2)

  })

  df1 <- do.call(rbind,lapply(df,function(x) x[[1]]))
  df2 <- do.call(rbind,lapply(df,function(x) x[[2]]))

  v <- v * max(df1[,1:2])
  df_arrow <- data.frame(rownames(v),0,v,stringsAsFactors=FALSE)
  colnames(df_arrow) <- c('covariate','origin','PC1','PC2')

  ggplot() +
    stat_contour(data=df1,aes_(~x,~y,z=~z,color=~..level..)) +
    geom_point(data=df2,aes_(~x,~y)) +
    geom_segment(data=df_arrow,aes_(x=~origin,y=~origin,xend=~PC1,yend=~PC2),
                 arrow=arrow(length=unit(0.03,'npc'))) +
    facet_wrap(~covariate) +
    viridis::scale_color_viridis() +
    theme_classic() +
    theme(aspect.ratio=1) +
    labs(x='Axis 1',y='Axis 2')

}
