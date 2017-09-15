# test constord (old conord)
library(theseus)
#library(phyloseq)
#library(ggplot2)


# let's try to prepare the dune dataset for use with theseus
?data
data(dune, package="vegan")
dune
rownames(dune)<-paste("Sample", sprintf("%02d",1:20), sep="")

data(dune.env, package="vegan")
dune.env
rownames(dune.env)<-paste("Sample", sprintf("%02d",1:20), sep="")

names(dune)
dim(dune)
tx<-data.frame(
  Kingdom=rep(as.character(NA),30), 
  Phylum=rep(as.character(NA),30), 
  Class=rep(as.character(NA),30), 
  Order=rep(as.character(NA),30), 
  Family=rep(as.character(NA),30), 
  Genus=rep(as.character(NA),30), 
  Species=rep(as.character(NA),30), 
  stringsAsFactors = FALSE
  )
str(tx)
tx$Species<-names(dune)
names(tx)
rownames(tx)<-names(dune)

PS<-phyloseq(
  otu_table(dune, taxa_are_rows = FALSE),
  sample_data(dune.env),
  tax_table(as.matrix(tx))
)
PS

# trouble with theseus function
# Error in overscope_eval_next(overscope, expr) : 
#  object 'Kingdom' not found
conord(PS, ~Management,tax_level = "Species")

# it looks like the constructed phyloseq object is fine. 
plot_ordination(PS,ordinate(PS), type="taxa", color="Species")

# save the dune phyloseq object as an .RData file. 
save(PS, 
     file=file.path(getwd(),"data/psdune.RData")
     )

# So the above dataset doesn't work. I've added another dataset WWTP_Impact from a paper we have currently in revision. 
# let's test constord with that. 
load(file.path(getwd(),"data/WWTP_Impact.RData"))
WWTP_Impact
pconord<- conord(WWTP_Impact, ~Position+Location)
pconord
# it works! 
# There's something wrong with the way I'm preparing the dune dataset. Any ideas why? 

