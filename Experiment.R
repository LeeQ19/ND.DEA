#########################################################################################################################
### Project  : Network DEA with malmquist for research institutes
#########################################################################################################################

# Load data
df.f.2d <- read.csv(url("https://docs.google.com/spreadsheets/d/e/2PACX-1vSaNq2LrKyvSWG2pisX4QnJw8ui7lj2lfQ4SVzwfFY5tl2BWf1AS5ORIfy1544dCNvfpAr8McUMiJk_/pub?output=csv"), header = T)
df.f.3d <- simplify2array(by(df.f.2d[, -c(1)], df.f.2d$YEAR, as.matrix))

# source("nd.dea.nlp.neos.R")
# Engine = c("ANTIGONE", "BARON", "CONOPT", "Couenne", "DICOPT", "Ipopt", "Knitro", "LINDOGlobal", "MINOS", "PATHNLP", "SBB", "SNOPT")
# res.nd.nlp <- nd.dea.nlp(df.f.3d[, id.xr, ], df.f.3d[, id.xt, ], df.f.3d[, id.yr, ], df.f.3d[, id.yt, ], df.f.3d[, id.z, ], engine = "DICOPT")
load("res.nlp.Rdata")

res.nd.nlp.DICOPT$eff.t
