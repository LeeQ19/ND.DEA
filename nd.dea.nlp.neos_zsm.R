pkgs <- c("ROI", "ROI.plugin.glpk", "CVXR", "slam")
sapply(pkgs, require, character.only = T)

make.l <- function (xt, indices, len) {
  l <- rep(0, len)
  l[indices] <- xt
  return(l)
}

nd.dea.nlp <- function (x1data, x2data, y1data, y2data, zdata, zlower = array(0, dim(zdata)), rts = "crs", orientation = "i", wv = NULL, se = FALSE, engine = "BARON") {
  if (length(unique(dim(x1data)[1], dim(y1data)[1], dim(x2data)[1], dim(y2data)[1], dim(zdata)[1])) != 1)
    stop("Data must be balanced.")
  if (length(unique(rev(dim(x1data))[1], rev(dim(y1data))[1], rev(dim(x2data))[1], rev(dim(y2data))[1], rev(dim(zdata))[1])) != 1)
    stop("Data must be balanced.")
  if (!(rts %in% c("crs", "vrs", "irs", "drs")))
    stop("rts must be \"crs\", \"vrs\", \"irs\", or \"drs\".")
  if (!(orientation %in% c("i", "o")))
    stop("orientation must be \"i\", \"o\", or \"n\".")
  if (!all(wv >= 0))
    stop("wv must be >= 0.")
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
  zlower <- if (length(dim(zlower)) != 3) {
    array(zlower, c(dim(zlower)[1], 1, dim(zlower)[2]))
  } else {
    as.array(zlower)
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
  p.wrs <- n * t + n * t + t + n * b * t + n * b * t + m1 * t + m2 * t + s1 * t + s2 * t + b * t + b * t + 1
  p.end <- n * t + n * t + t + n * b * t + n * b * t + m1 * t + m2 * t + s1 * t + s2 * t + b * t + b * t + n * b * t + 1
  
  for (j in 1:n) {
    if (orientation == "i") 
      obj <- make.l(c(wv), indices = c(p.eff:(p.z21 - 1)), p.end - 1)
    if (orientation == "o") 
      obj <- make.l(c(-wv), indices = c(p.eff:(p.z21 - 1)), p.end - 1)
    
    Q <- NULL
    L <- NULL
    dir <- NULL
    rhs <- NULL
    
    for (k in 1:t) {
      if (orientation == "i") {
        for (i in 1:m1) {
          Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
          L <- rbind(L, make.l(c(x1data[, i, k], -x1data[j, i, k], 1), indices = c((1 + (k - 1) * n):(k * n), p.eff - 1 + k, p.x1s - 1 + m1 * (k - 1) + i), p.end - 1))
          dir <- append(dir, "==")
          rhs <- c(rhs, 0)
        }
        for (i in 1:m2) {
          Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
          L <- rbind(L, make.l(c(x2data[, i, k], -x2data[j, i, k], 1), indices = c((p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), p.eff - 1 + k, p.x2s - 1 + m2 * (k - 1) + i), p.end - 1))
          dir <- append(dir, "==")
          rhs <- c(rhs, 0)
        }
        for (r in 1:s1) {
          Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
          L <- rbind(L, make.l(c(y1data[, r, k], -1), indices = c((1 + (k - 1) * n):(k * n), p.y1s - 1 + s1 * (k - 1) + r), p.end - 1))
          dir <- append(dir, "==")
          rhs <- c(rhs, y1data[j, r, k])
        }
        for (r in 1:s2) {
          Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
          L <- rbind(L, make.l(c(y2data[, r, k], -1), indices = c((p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), p.y2s - 1 + s2 * (k - 1) + r), p.end - 1))
          dir <- append(dir, "==")
          rhs <- c(rhs, y2data[j, r, k])
        }
      }
      else {
        for (i in 1:m1) {
          Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
          L <- rbind(L, make.l(c(x1data[, i, k], 1), indices = c((1 + (k - 1) * n):(k * n), p.x1s - 1 + m1 * (k - 1) + i), p.end - 1))
          dir <- append(dir, "==")
          rhs <- c(rhs, x1data[j, i, k])
        }
        for (i in 1:m2) {
          Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
          L <- rbind(L, make.l(c(x2data[, i, k], 1), indices = c((p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), p.x2s - 1 + m2 * (k - 1) + i), p.end - 1))
          dir <- append(dir, "==")
          rhs <- c(rhs, x2data[j, i, k])
        }
        for (r in 1:s1) {
          Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
          L <- rbind(L, make.l(c(y1data[, r, k], -y1data[j, r, k], -1), indices = c((1 + (k - 1) * n):(k * n), p.eff - 1 + k, p.y1s - 1 + s1 * (k - 1) + r), p.end - 1))
          dir <- append(dir, "==")
          rhs <- c(rhs, 0)
        }
        for (r in 1:s2) {
          Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
          L <- rbind(L, make.l(c(y2data[, r, k], -y2data[j, r, k], -1), indices = c((p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), p.eff - 1 + k, p.y2s - 1 + s2 * (k - 1) + r), p.end - 1))
          dir <- append(dir, "==")
          rhs <- c(rhs, 0)
        }
      }
      for (l in 1:b) {
        Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
        L <- rbind(L, make.l(c(zdata[, l, k], -1), indices = c((1 + (k - 1) * n):(k * n), p.z1s - 1 + b * (k - 1) + l), p.end - 1))
        dir <- append(dir, "==")
        rhs <- c(rhs, zdata[j, l, k])
        
        if (k == 1) {
          Q <- append(Q, list(simple_triplet_matrix(i = c(1:(p.end - 1), (p.lm2):(p.lm2 - 1 + n)), 
                                                    j = c(1:(p.end - 1), (p.z21 + n * (l - 1)):(p.z21 - 1 + n * l)), 
                                                    c(rep(0, p.end - 1), rep(2, n)))))
          L <- rbind(L, make.l(c(-1, 1), indices = c(p.z21 - 1 + n * (l - 1) + j, p.z2s - 1 + l), p.end - 1))
          dir <- append(dir, "==")
          rhs <- c(rhs, 0)
        }
        else {
          Q <- append(Q, list(simple_triplet_matrix(i = c(1:(p.end - 1), (p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), (p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n)), 
                                                    j = c(1:(p.end - 1), (p.z21 + n * b * (k - 1) + n * (l - 1)):(p.z21 - 1 + n * b * (k - 1) + n * l), (p.z22 + n * b * (k - 2) + n * (l - 1)):(p.z22 - 1 + n * b * (k - 2) + n * l)), 
                                                    c(rep(0, p.end - 1), rep(2, n), rep(2, n)))))
          L <- rbind(L, make.l(c(-1, -1, 1), indices = c(p.z21 - 1 + n * b * (k - 1) + n * (l - 1) + j, p.z22 -1 + n * b * (k - 2) + n * (l - 1) + j, p.z2s - 1 + b * (k - 1) + l), p.end - 1))
          dir <- append(dir, "==")
          rhs <- c(rhs, 0)
        }
        for (h in 1:n) {
          Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
          L <- rbind(L, make.l(c(1, 1), indices = c(c(p.z21, p.z22) - 1 + n * b * (k - 1) + n * (l - 1) + h), p.end - 1))
          dir <- append(dir, "==")
          rhs <- c(rhs, zdata[h, l, k])
          
          if (k == 1) {
            Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
            L <- rbind(L, make.l(c(1, -1), indices = c(p.z21 - 1 + n * (l - 1) + h, p.wrs - 1 +  n * (l - 1) + h), p.end - 1))
            dir <- append(dir, "==")
            rhs <- c(rhs, zlower[h, l, 1])
          }
          else {
            Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
            L <- rbind(L, make.l(c(1, 1, -1), indices = c(p.z21 - 1 + n * b * (k - 1) + n * (l - 1) + h, p.z22 - 1 + n * b * (k - 2) + n * (l - 1) + h, p.wrs - 1 + n * b * (k - 1) + n * (l - 1) + h), p.end - 1))
            dir <- append(dir, "==")
            rhs <- c(rhs, zlower[h, l, k])
          }
        }
      }
    }
    for (h in 1:n) {
      for (l in 1:b) {
        Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
        L <- rbind(L, make.l(c(1), indices = c(p.z22 - 1 + n * b * (t - 1) + n * (l - 1) + h), p.end - 1))
        dir <- append(dir, "==")
        rhs <- c(rhs, 0)
      }
    }
    cnst <- Q_constraint(Q = Q, L = L, dir = dir, rhs = rhs)
    op <- OP(obj, cnst)
    res <- ROI_solve(op, solver = "neos", method = engine)
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
    
    # Stage II
    obj <- make.l(c(rep(-1, p.x1s - p.z22)), indices = c(p.z22:(p.x1s - 1)), p.end - 1)

    for (k in 1:t) {
      Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
      L <- rbind(L, make.l(c(1), indices = c(p.eff - 1 + k), p.end - 1))
      dir <- append(dir, "==")
      rhs <- c(rhs, results.efficiency.t[j, k])
    }
    cnst <- Q_constraint(Q = Q, L = L, dir = dir, rhs = rhs)
    op <- OP(obj, cnst)
    res <- ROI_solve(op, solver = "neos", method = engine)
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
  results <- list(eff.s = results.efficiency.s, 
                  eff.t = results.efficiency.t, 
                  lambda1 = results.lambda1, 
                  lambda2 = results.lambda2,
                  z21data = results.z21data, 
                  z22data = results.z22data, 
                  x1slack = results.x1slack, 
                  x2slack = results.x2slack, 
                  y1slack = results.y1slack, 
                  y2slack = results.y2slack, 
                  z1slack = results.z1slack, 
                  z2slack = results.z2slack)
  return(results)
}
