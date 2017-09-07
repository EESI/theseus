#' @import ggplot2
NULL

#' Plots quality assessment
#'
#' Description paragraph
#'
#' @param PS (required) asdf
#' @param covariates (required) asdf
#'
#' @return A ggplot object.
#'
#' @references
#' X
#' Y
#' Z
#'
#' @seealso \code{\link[vegan]{rda}} \code{\link[vegan]{envfit}} \code{\link[vegan]{ordisurf}}
#'
#' @examples
#' \dontrun{
#' covariates <- c('SL_NPOC','SL_NO3','SL_NH4')
#' surf(PS,covariates)
#' }
#'
#' @export

surf <- function(PS,covariates){

  plot.new()

  OTU <- PS@otu_table@.Data
  SAMP <- suppressWarnings(as_tibble(PS@sam_data))

  if (!all(covariates %in% colnames(SAMP)))
    stop('Some covaraites not found in sample data!')

  if (PS@otu_table@taxa_are_rows) OTU <- t(OTU)

  R <- vegan::rda(OTU)

  form <- as.formula(sprintf('%s ~ %s','R',paste(covariates,collapse=' + ')))
  env <- vegan::envfit(form,SAMP,na.rm=TRUE)

  v <- (sqrt(env$vectors$r) * env$vectors$arrows)[,1:2,drop=FALSE]

  df <- lapply(covariates,function(k){

    form <- as.formula(sprintf('%s ~ %s','R',k))
    O <- vegan::ordisurf(form,SAMP,plot=FALSE)

    df1 <- with(O,data.frame(x=grid$x,y=rep(grid$y,each=ncol(grid$z)),z=matrix(grid$z))) %>%
      filter(!is.na(z)) %>%
      mutate(covariate=k)

    df2 <- with(O,data.frame(x=model$x1,y=model$x2,z=NA)) %>%
      mutate(covariate=k)

    list(df1=df1,df2=df2)

  })

  df1 <- do.call(rbind,lapply(df,function(x) x[[1]]))
  df2 <- do.call(rbind,lapply(df,function(x) x[[2]]))

  v <- v * max(df1[,1:2])
  df_arrow <- data.frame(rownames(v),0,v,stringsAsFactors=FALSE)
  colnames(df_arrow) <- c('covariate','origin','PC1','PC2')

  ggplot() +
    stat_contour(data=df1,aes(x,y,z=z,color=..level..)) +
    geom_point(data=df2,aes(x,y)) +
    geom_segment(data=df_arrow,aes(x=origin,y=origin,xend=PC1,yend=PC2),
                 arrow=arrow(length=unit(0.03,'npc'))) +
    facet_wrap(~covariate) +
    viridis::scale_color_viridis() +
    theme_classic() +
    theme(aspect.ratio=1) +
    labs(x='Axis 1',y='Axis 2')

}