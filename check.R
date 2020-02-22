pkgs <- c("ROI", "ROI.plugin.glpk", "CVXR", "slam")
sapply(pkgs, require, character.only = T)

make.l <- function (xt, indices, len) {
  l <- rep(0, len)
  l[indices] <- xt
  return(l)
}

check <- function (op, v) {
  obj.val <- as.vector(op$objective$L) %*% v
  Q.val <- sapply(op$constraints$Q, function (temp) {0.5 * v %*% as.matrix(temp) %*% v})
  L.val <- as.matrix(op$constraints$L) %*% v
  lhs <- c(Q.val + L.val)
  dif <- lhs - op$constraints$rhs
  
  return(list(obj.val = obj.val, difference = dif))
}


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


if (orientation == "i") 
  obj <- make.l(c(wv), indices = c(p.eff:(p.z21 - 1)), p.end - 1)
if (orientation == "o") 
  obj <- make.l(c(-wv), indices = c(p.eff:(p.z21 - 1)), p.end - 1)

Q <- NULL
L <- NULL
dir <- NULL
rhs <- NULL
cat <- NULL

for (k in 1:t) {
  if (orientation == "i") {
    for (i in 1:m1) {
      Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
      L <- rbind(L, make.l(c(x1data[, i, k], -x1data[j, i, k], 1), indices = c(((k - 1) * n + 1):(k * n), p.eff - 1 + k, p.x1s - 1 + m1 * (k - 1) + i), p.end - 1))
      dir <- append(dir, "==")
      rhs <- c(rhs, 0)
      cat <- append(cat, paste0("x1^", k, "_", i))
    }
    for (i in 1:m2) {
      Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
      L <- rbind(L, make.l(c(x2data[, i, k], -x2data[j, i, k], 1), indices = c((p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), p.eff - 1 + k, p.x2s - 1 + m2 * (k - 1) + i), p.end - 1))
      dir <- append(dir, "==")
      rhs <- c(rhs, 0)
      cat <- append(cat, paste0("x2^", k, "_", i))
    }
    for (r in 1:s1) {
      Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
      L <- rbind(L, make.l(c(y1data[, r, k], -1), indices = c(((k - 1) * n + 1):(k * n), p.y1s - 1 + s1 * (k - 1) + r), p.end - 1))
      dir <- append(dir, "==")
      rhs <- c(rhs, y1data[j, r, k])
      cat <- append(cat, paste0("y1^", k, "_", r))
    }
    for (r in 1:s2) {
      Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
      L <- rbind(L, make.l(c(y2data[, r, k], -1), indices = c((p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), p.y2s - 1 + s2 * (k - 1) + r), p.end - 1))
      dir <- append(dir, "==")
      rhs <- c(rhs, y2data[j, r, k])
      cat <- append(cat, paste0("y2^", k, "_", r))
    }
  }
  else {
    for (i in 1:m1) {
      Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
      L <- rbind(L, make.l(c(x1data[, i, k], 1), indices = c(((k - 1) * n + 1):(k * n), p.x1s - 1 + m1 * (k - 1) + i), p.end - 1))
      dir <- append(dir, "==")
      rhs <- c(rhs, x1data[j, i, k])
      cat <- append(cat, paste0("x1^", k, "_", i))
    }
    for (i in 1:m2) {
      Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
      L <- rbind(L, make.l(c(x2data[, i, k], 1), indices = c((p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), p.x2s - 1 + m2 * (k - 1) + i), p.end - 1))
      dir <- append(dir, "==")
      rhs <- c(rhs, x2data[j, i, k])
      cat <- append(cat, paste0("x2^", k, "_", i))
    }
    for (r in 1:s1) {
      Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
      L <- rbind(L, make.l(c(y1data[, r, k], -y1data[j, r, k], -1), indices = c(((k - 1) * n + 1):(k * n), p.eff - 1 + k, p.y1s - 1 + s1 * (k - 1) + r), p.end - 1))
      dir <- append(dir, "==")
      rhs <- c(rhs, 0)
      cat <- append(cat, paste0("y1^", k, "_", r))
    }
    for (r in 1:s2) {
      Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
      L <- rbind(L, make.l(c(y2data[, r, k], -y2data[j, r, k], -1), indices = c((p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), p.eff - 1 + k, p.y2s - 1 + s2 * (k - 1) + r), p.end - 1))
      dir <- append(dir, "==")
      rhs <- c(rhs, 0)
      cat <- append(cat, paste0("y2^", k, "_", r))
    }
  }
  for (l in 1:b) {
    Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
    L <- rbind(L, make.l(c(zdata[, l, k], -1), indices = c(((k - 1) * n + 1):(k * n), p.z1s - 1 + b * (k - 1) + l), p.end - 1))
    dir <- append(dir, "==")
    rhs <- c(rhs, zdata[j, l, k])
    cat <- append(cat, paste0("z1^", k, "_", l))

    if (k == 1) {
      Q <- append(Q, list(simple_triplet_matrix(i = c(1:(p.end - 1), (p.lm2):(p.lm2 - 1 + n)),
                                                j = c(1:(p.end - 1), (p.z21 + n * (l - 1)):(p.z21 - 1 + n * l)),
                                                c(rep(0, p.end - 1), rep(2, n)))))
      L <- rbind(L, make.l(c(-1, 1), indices = c(p.z21 - 1 + n * (l - 1) + j, p.z2s - 1 + l), p.end - 1))
      dir <- append(dir, "==")
      rhs <- c(rhs, 0)
      cat <- append(cat, paste0("z2^", k, "_", l))
    }
    else {
      Q <- append(Q, list(simple_triplet_matrix(i = c(1:(p.end - 1), (p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n), (p.lm2 + (k - 1) * n):(p.lm2 - 1 + k * n)),
                                                j = c(1:(p.end - 1), (p.z21 + n * b * (k - 1) + n * (l - 1)):(p.z21 - 1 + n * b * (k - 1) + n * l), (p.z22 + n * b * (k - 2) + n * (l - 1)):(p.z22 - 1 + n * b * (k - 2) + n * l)),
                                                c(rep(0, p.end - 1), rep(2, n), rep(2, n)))))
      L <- rbind(L, make.l(c(-1, -1, 1), indices = c(p.z21 - 1 + n * b * (k - 1) + n * (l - 1) + j, p.z22 - 1 + n * b * (k - 2) + n * (l - 1) + j, p.z2s - 1 + b * (k - 1) + l), p.end - 1))
      dir <- append(dir, "==")
      rhs <- c(rhs, 0)
      cat <- append(cat, paste0("z2^", k, "_", l))
    }
    for (h in 1:n) {
      Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
      L <- rbind(L, make.l(c(1, 1), indices = c(c(p.z21, p.z22) - 1 + n * b * (k - 1) + n * (l - 1) + h), p.end - 1))
      dir <- append(dir, "==")
      rhs <- c(rhs, zdata[h, l, k])
      cat <- append(cat, paste0("zsum^", k, "_", l, ",", h))

      if (k == 1) {
        Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
        L <- rbind(L, make.l(c(1, -1), indices = c(p.z21 - 1 + n * (l - 1) + h, p.wrs - 1 +  n * (l - 1) + h), p.end - 1))
        dir <- append(dir, "==")
        rhs <- c(rhs, zlower[h, l, 1])
        cat <- append(cat, paste0("zlower^", k, "_", l, ",", h))
      }
      else {
        Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
        L <- rbind(L, make.l(c(1, 1, -1), indices = c(p.z21 - 1 + n * b * (k - 1) + n * (l - 1) + h, p.z22 - 1 + n * b * (k - 2) + n * (l - 1) + h, p.wrs - 1 + n * b * (k - 1) + n * (l - 1) + h), p.end - 1))
        dir <- append(dir, "==")
        rhs <- c(rhs, zlower[h, l, k])
        cat <- append(cat, paste0("zlower^", k, "_", l, ",", h))
      }
    }
  }
}
for (l in 1:b) {
  for (h in 1:n) {
    Q <- append(Q, list(simple_triplet_matrix(i = 1:(p.end - 1), j = 1:(p.end - 1), rep(0, p.end - 1))))
    L <- rbind(L, make.l(c(1), indices = c(p.z22 - 1 + n * b * (t - 1) + n * (l - 1) + h), p.end - 1))
    dir <- append(dir, "==")
    rhs <- c(rhs, 0)
    cat <- append(cat, paste0("zfinal", "_", l, ",", h))
  }
}
cnst <- Q_constraint(Q = Q, L = L, dir = dir, rhs = rhs)
op <- OP(obj, cnst)
res <- ROI_solve(op, solver = "neos", method = "IPOPT")

ch <- check(op, res$solution)
round(ch$difference, 4)
cat[which(round(ch$difference, 4) != 0)]
