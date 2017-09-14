# test envtoverlay (old "surf")


# load wwtp data
#data(WWTP_Impact, package="theseus")
load(file.path(getwd(),"data/WWTP_Impact.RData"))
WWTP_Impact
covariates <-c("log_Cl", "log_Si", "log_Sr", "log_F")
psurf<-surf(WWTP_Impact, covariates)
psurf
