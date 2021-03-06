#########################################################################################################################
### Project  : Network DEA with malmquist for research institutes
#########################################################################################################################
# Load library
source("nd.dea.nlp.neos_zsm.R")


# Load data
df.f.2d <- read.csv(url("https://docs.google.com/spreadsheets/d/e/2PACX-1vSaNq2LrKyvSWG2pisX4QnJw8ui7lj2lfQ4SVzwfFY5tl2BWf1AS5ORIfy1544dCNvfpAr8McUMiJk_/pub?output=csv"), header = T)
df.f.3d <- simplify2array(by(df.f.2d[, -c(1)], df.f.2d$YEAR, as.matrix))

# Set parameter
id.out <- c(2)
id.t <- c(1)
id.x1 <- c(2, 13)
id.x2 <- c(4, 11)
id.y1 <- c(5, 6)
id.y2 <- c(9)
id.z <- c(7)
id.zl <- c(14)
rts  <- "crs"
orientation  <- "i"
wv <- NULL
engine <- "Ipopt"

df.f.3d <- df.f.3d[-id.out, , ]

# Measure efficiency
# Engine = c("ANTIGONE", "BARON", "CONOPT", "Ipopt", "Knitro", "LINDOGlobal", "MINOS", "PATHNLP", "SBB", "SNOPT")
res.nd.nlp <- nd.dea.nlp(df.f.3d[, id.x1, ], df.f.3d[, id.x2, ], df.f.3d[, id.y1, ], df.f.3d[, id.y2, ], df.f.3d[, id.z, ], df.f.3d[, id.zl, ], engine = engine)
# load("res.nlp.Rdata")

write.csv(res.nd.nlp, file = paste0(engine, ".csv"))
