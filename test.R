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
df.f.2d <- read.csv(url("https://docs.google.com/spreadsheets/d/e/2PACX-1vSaNq2LrKyvSWG2pisX4QnJw8ui7lj2lfQ4SVzwfFY5tl2BWf1AS5ORIfy1544dCNvfpAr8McUMiJk_/pub?output=csv"), header = T)
df.f.3d <- simplify2array(by(df.f.2d[, -c(1)], df.f.2d$YEAR, as.matrix))

# Parameter
id.t <- c(1)
id.xr <- c(2, 13)
id.xt <- c(4, 11)
id.yr <- c(5, 6)
id.yt <- c(9)
id.z <- c(7)
rts  <- "crs"
orientation  <- "i"
wv <- NULL

x1data <- df.f.3d[, id.xr, ]
x2data <- df.f.3d[, id.xt, ]
y1data <- df.f.3d[, id.yr, ]
y2data <- df.f.3d[, id.yt, ]
zdata  <- df.f.3d[, id.z, ]

make.l <- function (xt, indices, len = p.end - 1) {
  l <- rep(0, len)
  l[indices] <- xt
  return(l)
}
{
  names <- if (is.null(rownames(x1data))) {
    1:n
  } else {
    rownames(x1data)
  }
  x1data <- if (length(dim(x1data)) != 3) {
    array(x1data, c(dim(x1data)[1], 1, dim(x1data)[2]))
  } else {
    as.array(x1data)
  }
  x2data <- if (length(dim(x2data)) != 3) {
    array(x2data, c(dim(x2data)[1], 1, dim(x2data)[2]))
  } else {
    as.array(x2data)
  }
  y1data <- if (length(dim(y1data)) != 3) {
    array(y1data, c(dim(y1data)[1], 1, dim(y1data)[2]))
  } else {
    as.array(y1data)
  }
  y2data <- if (length(dim(y2data)) != 3) {
    array(y2data, c(dim(y2data)[1], 1, dim(y2data)[2]))
  } else {
    as.array(y2data)
  }
  zdata <- if (length(dim(zdata)) != 3) {
    array(zdata, c(dim(zdata)[1], 1, dim(zdata)[2]))
  } else {
    as.array(zdata)
  }
  
  n  <- dim(x1data)[1]
  m1 <- dim(x1data)[2]
  m2 <- dim(x2data)[2]
  s1 <- dim(y1data)[2]
  s2 <- dim(y2data)[2]
  b  <- dim(zdata)[2]
  t  <- dim(x1data)[3]
  
  wv <- if (is.null(wv)) {
    rep(1, t)
  } else {
    as.vector(wv)
  }
  
  results.lambda1 <- array(NA, dim = c(n, n, t), dimnames = list(names, names))
  results.lambda2 <- array(NA, dim = c(n, n, t), dimnames = list(names, names))
  results.efficiency.s <- array(NA, dim = c(n, 1), dimnames = list(names, "Eff.Sys"))
  results.efficiency.t <- array(NA, dim = c(n, t), dimnames = list(names, paste0("Eff.T.", 1:t)))
  results.z21data <- array(NA, dim = c(n, n, b, t), dimnames = list(names, names, paste0("z21data.", 1:b)))
  results.z22data <- array(NA, dim = c(n, n, b, t), dimnames = list(names, names, paste0("z22data.", 1:b)))
  results.x1slack <- array(NA, dim = c(n, m1, t), dimnames = list(names, paste0("x1slack.", 1:m1)))
  results.x2slack <- array(NA, dim = c(n, m2, t), dimnames = list(names, paste0("x2slack.", 1:m2)))
  results.y1slack <- array(NA, dim = c(n, s1, t), dimnames = list(names, paste0("y1slack.", 1:s1)))
  results.y2slack <- array(NA, dim = c(n, s2, t), dimnames = list(names, paste0("y2slack.", 1:s2)))
  results.z1slack <- array(NA, dim = c(n, b, t), dimnames = list(names, paste0("z1slack.", 1:b)))
  results.z2slack <- array(NA, dim = c(n, b, t), dimnames = list(names, paste0("z2slack.", 1:b)))
  
  p.lm2 <- n * t + 1
  p.eff <- n * t + n * t + 1
  p.z21 <- n * t + n * t + t + 1
  p.z22 <- n * t + n * t + t + n * b * t + 1
  p.x1s <- n * t + n * t + t + n * b * t + n * b * t + 1
  p.x2s <- n * t + n * t + t + n * b * t + n * b * t + m1 * t + 1
  p.y1s <- n * t + n * t + t + n * b * t + n * b * t + m1 * t + m2 * t + 1
  p.y2s <- n * t + n * t + t + n * b * t + n * b * t + m1 * t + m2 * t + s1 * t + 1
  p.z1s <- n * t + n * t + t + n * b * t + n * b * t + m1 * t + m2 * t + s1 * t + s2 * t + 1
  p.z2s <- n * t + n * t + t + n * b * t + n * b * t + m1 * t + m2 * t + s1 * t + s2 * t + b * t + 1
  p.end <- n * t + n * t + t + n * b * t + n * b * t + m1 * t + m2 * t + s1 * t + s2 * t + b * t + b * t + 1
  
  if (orientation == "i") 
    obj <- c(rep(0, p.eff - 1), wv, rep(0, p.end - p.z21))
  if (orientation == "o") 
    obj <- c(rep(0, p.eff - 1), -wv, rep(0, p.end - p.z21))
  
  Q <- NULL
  L <- NULL
  rhs <- NULL
  
  for (k in 1:t) {
    if (orientation == "i") {
      for (i in 1:m1) {
        Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
        L <- rbind(L, make.l(c(x1data[, i, k], -x1data[j, i, k], 1), indices = c(((k - 1) * n + 1):(k * n), p.eff - 1 + k, p.x1s - 1 + m1 * (k - 1) + i)))
        rhs <- c(rhs, 0)
      }
      for (i in 1:m2) {
        Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
        L <- rbind(L, make.l(c(x2data[, i, k], -x2data[j, i, k], 1), indices = c((p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), p.eff - 1 + k, p.x2s - 1 + m2 * (k - 1) + i)))
        rhs <- c(rhs, 0)
      }
      for (r in 1:s1) {
        Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
        L <- rbind(L, make.l(c(y1data[, r, k], -1), indices = c(((k - 1) * n + 1):(k * n), p.y1s - 1 + s1 * (k - 1) + r)))
        rhs <- c(rhs, y1data[j, r, k])
      }
      for (r in 1:s2) {
        Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
        L <- rbind(L, make.l(c(y2data[, r, k], -1), indices = c((p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), p.y2s - 1 + s2 * (k - 1) + r)))
        rhs <- c(rhs, y2data[j, r, k])
      }
    }
    else {
      for (i in 1:m1) {
        Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
        L <- rbind(L, make.l(c(x1data[, i, k], 1), indices = c(((k - 1) * n + 1):(k * n), p.x1s - 1 + m1 * (k - 1) + i)))
        rhs <- c(rhs, x1data[j, i, k])
      }
      for (i in 1:m2) {
        Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
        L <- rbind(L, make.l(c(x2data[, i, k], 1), indices = c((p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), p.x2s - 1 + m2 * (k - 1) + i)))
        rhs <- c(rhs, x2data[j, i, k])
      }
      for (r in 1:s1) {
        Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
        L <- rbind(L, make.l(c(y1data[, r, k], -y1data[j, r, k], -1), indices = c(((k - 1) * n + 1):(k * n), p.eff - 1 + k, p.y1s - 1 + s1 * (k - 1) + r)))
        rhs <- c(rhs, 0)
      }
      for (r in 1:s2) {
        Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
        L <- rbind(L, make.l(c(y2data[, r, k], -y2data[j, r, k], -1), indices = c((p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), p.eff - 1 + k, p.y2s - 1 + s2 * (k - 1) + r)))
        rhs <- c(rhs, 0)
      }
    }
    for (l in 1:b) {
      Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
      L <- rbind(L, make.l(c(zdata[, l, k], -1), indices = c(((k - 1) * n + 1):(k * n), p.z1s - 1 + b * (k - 1) + l)))
      rhs <- c(rhs, zdata[j, l, k])
      
      if (k == 1) {
        Q <- append(Q, list(simple_triplet_matrix(i = c(1:(p.end - 1), (p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n)), 
                                                  j = c(1:(p.end - 1), (p.z21 + n * b * (k - 1) + n * (l - 1)):(p.z21 - 1 + n * b * (k - 1) + n * l)), 
                                                  c(rep(0, p.end - 1), rep(2, n)))))
        L <- rbind(L, make.l(c(-1, 1), indices = c(p.z21 - 1 + n * b * (k - 1) + n * (l - 1) + j, p.z2s - 1 + b * (k - 1) + l)))
        rhs <- c(rhs, 0)
      }
      else {
        Q <- append(Q, list(simple_triplet_matrix(i = c(1:(p.end - 1), (p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), (p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n)), 
                                                  j = c(1:(p.end - 1), (p.z21 + n * b * (k - 1) + n * (l - 1)):(p.z21 - 1 + n * b * (k - 1) + n * l), (p.z22 + n * b * (k - 2) + n * (l - 1)):(p.z22 - 1 + n * b * (k - 2) + n * l)), 
                                                  c(rep(0, p.end - 1), rep(2, n), rep(2, n)))))
        L <- rbind(L, make.l(c(-1, -1, 1), indices = c(p.z21 - 1 + n * b * (k - 1) + n * (l - 1) + j, p.z22 + n * b * (k - 2) + n * (l - 1) + j, p.z2s - 1 + b * (k - 1) + l)))
        rhs <- c(rhs, 0)
      }
      
      for (h in 1:n) {
        Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
        L <- rbind(L, make.l(c(1, 1), indices = c(c(p.z21, p.z22) - 1 + n * b * (k - 1) + n * (l - 1) + h)))
        rhs <- c(rhs, zdata[h, l, k])
      }
    }
  }
  for (h in 1:n) {
    for (l in 1:b) {
      Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
      L <- rbind(L, make.l(c(1), indices = c(p.z22 - 1 + n * b * (t - 1) + n * (l - 1) + h)))
      rhs <- c(rhs, 0)
    }
  }
  cnst <- Q_constraint(Q = Q, L = L, dir = eq(length(rhs)), rhs = rhs)
  op <- OP(obj, cnst)
  res <- ROI_solve(op, solver = "neos", method = "ANTIGONE")
  results.efficiency.s[j] <- abs(res$objval)/sum(wv)
  temp.p <- res$solution
  results.lambda1[j, , ]    <- array(temp.p[1:(p.lm2 - 1)], c(n, t))
  results.lambda2[j, , ]    <- array(temp.p[p.lm2:(p.eff - 1)], c(n, t))
  results.efficiency.t[j, ] <- temp.p[p.eff:(p.z21 - 1)]
  results.z21data[j, , , ]  <- array(temp.p[p.z21:(p.z22 - 1)], c(n, b, t))
  results.z22data[j, , , ]  <- array(temp.p[p.z22:(p.x1s - 1)], c(n, b, t))
  results.x1slack[j, , ]    <- array(temp.p[p.x1s:(p.x2s - 1)], c(m1, t))
  results.x2slack[j, , ]    <- array(temp.p[p.x2s:(p.y1s - 1)], c(m2, t))
  results.y1slack[j, , ]    <- array(temp.p[p.y1s:(p.y2s - 1)], c(s1, t))
  results.y2slack[j, , ]    <- array(temp.p[p.y2s:(p.z1s - 1)], c(s2, t))
  results.z1slack[j, , ]    <- array(temp.p[p.z1s:(p.z2s - 1)], c(b, t))
  results.z2slack[j, , ]    <- array(temp.p[p.z2s:(p.end - 1)], c(b, t))
}
