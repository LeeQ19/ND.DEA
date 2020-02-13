#################################################################################
#################################################################################
library(ROI)
mat <- matrix(c(3, 4, 2,
                2, 1, 2,
                1, 3, 2), nrow=3, byrow=TRUE)
x <- OP(objective = c(2, 4, 3),
        constraints = L_constraint(L = mat,
                                   dir = c("<=", "<=", "<="),
                                   rhs = c(60, 40, 80)),
        maximum = TRUE)
opt <- ROI_solve(x, solver = "neos", method = "Bonmin")
opt
solution(opt)

#################################################################################
#################################################################################
pkgs <- c("dplyr", "ROI", "ROI.plugin.glpk", "ompr", "ompr.roi", "CVXR", "slam")
sapply(pkgs, require, character.only = T)
full <- Variable(1)
compact <- Variable(1)
p_full <- Variable(1)
p_compact <- Variable(1)
objective <- Maximize(full * (p_full - 150) + compact * (p_compact - 100))
constr <- list(2 * full + compact <= 500,
               2 * full + 3 * compact <= 800, full <= 490 - p_full,
               compact <= 640 - 2 * p_compact, full >= 0, compact >= 0,
               p_full >= 0, p_compact >= 0)
magnetron <- Problem(objective, constr)
cvxr_sol <- solve(magnetron)

Q <- simple_triplet_matrix(i = 1:4, j = c(2, 1, 4, 3), rep(1, 4))
as.matrix(Q)
var_names <- c("full", "price_full", "compact", "price_compact")
o <- OP(Q_objective(Q = Q, L = c(-150, 0, -100, 0), names = var_names), 
        L_constraint(rbind(c(2, 0, 1, 0), c(2, 0, 3, 0), c(1, 1, 0, 0), c(0, 0, 1, 2)), 
                     dir = leq(4), rhs = c(500, 800, 490, 640)), maximum = TRUE)
sol <- ROI_solve(o, solver = "neos", method = "BARON")
solution(sol)
#################################################################################
#################################################################################
pkgs <- c("ROI", "CVXR", "slam")
sapply(pkgs, require, character.only = T)
obj  <- c(1, 1)
cnst <- Q_constraint(Q = simple_triplet_matrix(i = c(1, 2), j = c(1, 2), c(-2, 0)),
                     L = c(4, 1), 
                     dir = leq(1), 
                     rhs = c(5))
o <- OP(obj, cnst, maximum = F)
sol <- ROI_solve(o, solver = "neos", method = "Bonmin")
sol
solution(sol)
#################################################################################
#################################################################################
pkgs <- c("ROI", "CVXR", "slam")
sapply(pkgs, require, character.only = T)
obj  <- c(1, 1)
cnst <- Q_constraint(Q = list(simple_triplet_matrix(i = rep(1:2, each = 2), j = rep(1:2, times = 2), c(2, 0, 0, 0)),
                              simple_triplet_matrix(i = rep(1:2, each = 2), j = rep(1:2, times = 2), c(0, 0, 0, 2))),
                     L = rbind(c(0, 0), c(0, 0)), 
                     dir = leq(2), 
                     rhs = c(4, 9))
o <- OP(obj, cnst, maximum = F)
sol <- ROI_solve(o, solver = "neos", method = "Bonmin")
sol
solution(sol)
#################################################################################
#################################################################################
pkgs <- c("ROI", "CVXR", "slam")
sapply(pkgs, require, character.only = T)
obj  <- c(1, 1)
cnst <- Q_constraint(Q = list(simple_triplet_matrix(i = rep(1:2, each = 2), j = rep(1:2, times = 2), c(0, 0, 0, 0)),
                              simple_triplet_matrix(i = rep(1:2, each = 2), j = rep(1:2, times = 2), c(0, 0, 0, 0))),
                     L = rbind(c(1, 0), c(0, 1)), 
                     dir = leq(2), 
                     rhs = c(2, 3))
o <- OP(obj, cnst, maximum = T)
sol <- ROI_solve(o, solver = "neos", method = "Bonmin")
sol
solution(sol)
#################################################################################
#################################################################################
pkgs <- c("ROI", "CVXR", "slam")
sapply(pkgs, require, character.only = T)
obj  <- c(1, 1)
cnst <- Q_constraint(Q = list(simple_triplet_matrix(i = rep(1:2, each = 2), j = rep(1:2, times = 2), c(0, 1, 1, 0)), 
                              simple_triplet_matrix(i = rep(1:2, each = 2), j = rep(1:2, times = 2), c(2, 0, 0, 0)), 
                              simple_triplet_matrix(i = rep(1:2, each = 2), j = rep(1:2, times = 2), c(0, 0, 0, 2))),
                     L = rbind(c(0, 0), 
                               c(0, 0), 
                               c(0, 0)), 
                     dir = leq(3), 
                     rhs = c(4, 16, 16))
o <- OP(obj, cnst, maximum = T)
sol <- ROI_solve(o, solver = "neos", method = "Bonmin")
sol
solution(sol)
#################################################################################
## Do not run
#################################################################################
pkgs <- c("ROI", "CVXR", "slam")
sapply(pkgs, require, character.only = T)
obj  <- c(1, 1)
cnst <- list(Q_constraint(Q = list(simple_triplet_matrix(i = rep(1:2, each = 2), j = rep(1:2, times = 2), c(2, 0, 0, 0)),
                                   simple_triplet_matrix(i = rep(1:2, each = 2), j = rep(1:2, times = 2), c(0, 0, 0, 2))),
                          L = rbind(c(0, 0), c(0, 0)), 
                          dir = leq(2), 
                          rhs = c(4, 9)), 
             L_constraint(L = rbind(c(1, 0)), 
                          dir = leq(1), 
                          rhs = c(1)))
o <- OP(obj, cnst, maximum = T)
sol <- ROI_solve(o, solver = "neos", method = "Bonmin")
sol
solution(sol)

#################################################################################
#################################################################################
pkgs <- c("ROI", "ROI.plugin.glpk", "CVXR", "slam")
sapply(pkgs, require, character.only = T)

# Load data
df.f.2d <- read.csv("Dataset_ToyData.csv", header = T)
df.f.3d <- simplify2array(by(df.f.2d[, -c(1)], df.f.2d$t, as.matrix))

# Parameter
id.t <- c(1)
id.x1 <- c(2)
id.x2 <- c(5)
id.y1 <- c(3)
id.y2 <- c(6)
id.z <- c(4)
id.zl <- c(7)
rts  <- "crs"
orientation  <- "i"
wv <- NULL
engine <- "IPOPT"
j <- 3

x1data <- df.f.3d[, id.x1, ]
x2data <- df.f.3d[, id.x2, ]
y1data <- df.f.3d[, id.y1, ]
y2data <- df.f.3d[, id.y2, ]
zdata  <- df.f.3d[, id.z, ]
zlower <- df.f.3d[, id.zl, ]
