library(gage)
library(pathview)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
}

run <- function() {}

output <- function(outputfile) {
ref.idx=5:8
samp.idx=1:4
#data(kegg.gs)
cnts.norm <- read.csv(paste(pfix, parameters["norm", 2], sep="/"))

#knockdown and control samples are unpaired
cnts.kegg.p <- readRDS(paste(pfix, parameters["gage", 2], sep="/"))

pdf(outputfile)
#gage(cnts.norm, gsets = kegg.gs, ref = ref.idx,
               #     samp = samp.idx, compare ="unpaired")


###################################################
### code chunk number 11: pathview
###################################################
#differential expression: log2 ratio or fold change, uppaired samples
cnts.d= cnts.norm[, samp.idx]-rowMeans(cnts.norm[, ref.idx])

#up-regulated pathways (top 3) visualized by pathview
sel <- cnts.kegg.p$greater[, "q.val"] < 0.1 &
         !is.na(cnts.kegg.p$greater[,"q.val"])
path.ids <- rownames(cnts.kegg.p$greater)[sel]
print(cnts.kegg.p)
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(
                      gene.data = cnts.d, pathway.id = pid,
                      species = "hsa"))

#down-regulated pathways  (top 3) visualized by pathview
sel.l <- cnts.kegg.p$less[, "q.val"] < 0.1 &
           !is.na(cnts.kegg.p$less[,"q.val"])
path.ids.l <- rownames(cnts.kegg.p$less)[sel.l]
path.ids.l2 <- substr(path.ids.l, 1, 8)
pv.out.list.l <- sapply(path.ids.l2[1:3], function(pid) pathview(
                      gene.data = cnts.d, pathway.id = pid,
                      species = "hsa"))
}
