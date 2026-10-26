library(ggplot2)
library(patchwork)
theme_set(theme_bw())

###Matérn nu = 3/2----
#sigma^2 (1 + sqrt(3) d/phi) exp(-sqrt(3) d/phi), phi = range
mrate <- function(phi){
  if(getOption("matern_param", "range") == "range"){
      sqrt(3) / phi 
  }
  else phi
}
cor32 <- function(d, phi) { 
  a <- mrate(phi)
  (1 + a * d) * exp(-a * d)
}
# parámetro identificable (Zhang 2004): sigma^2 a^3, expresado como delta
delta_fun <- function(phi, tau){
  if(getOption("matern_param", "range") == "range"){
    tau * phi^3
  }
  else phi^3 / tau
}
Rmat  <- function(D, phi){
  cor32(D, phi) + diag(1e-6, nrow(D))
}

grid_centers <- function(n) {
  g <- (seq_len(n) - 0.5) / n
  as.matrix(expand.grid(x = g, y = g))
}
dist_mat <- function(A, B = A){
  sqrt(outer(A[, 1], B[, 1], "-")^2 + outer(A[, 2], B[, 2], "-")^2)
}

rgp <- function(D, phi, tau){ 
  drop(crossprod(chol(Rmat(D, phi)), rnorm(nrow(D)))) / sqrt(tau)
}

gp_ld <- function(th, phi, tau, D) {
  U <- chol(Rmat(D, phi))
  z <- backsolve(U, th, transpose = TRUE)
  -sum(log(diag(U))) + 0.5 * length(th) * log(tau) - 0.5 * tau * sum(z^2)
}

# ---- Samplers ----
ess <- function(f, L, loglik, ll) {
  nu <- drop(L %*% rnorm(length(f)))
  ly <- ll + log(runif(1))
  th <- runif(1, 0, 2 * pi)
  lo <- th - 2 * pi
  hi <- th
  repeat {
    fp <- f * cos(th) + nu * sin(th)
    llp <- loglik(fp)
    if (llp > ly) return(list(f = fp, ll = llp))
    if (th < 0) lo <- th else hi <- th
    th <- runif(1, lo, hi)
  }
}

slice1 <- function(x, logf, w, lo = -Inf, hi = Inf) {
  ly <- logf(x) + log(runif(1))
  a <- x - w * runif(1)
  b <- a + w
  while (a > lo && logf(a) > ly) a <- a - w
  while (b < hi && logf(b) > ly) b <- b + w
  a <- max(a, lo)
  b <- min(b, hi)
  repeat {
    x1 <- runif(1, a, b)
    if (logf(x1) > ly) return(x1)
    if (x1 < x) a <- x1 else b <- x1
  }
}

upd_hyper <- function(th, phi, tau, D, prior) {
  phi <- slice1(phi, function(p) gp_ld(th, p, tau, D), 0.05, prior$phi[1], prior$phi[2])
  U <- chol(Rmat(D, phi))
  q <- sum(backsolve(U, th, transpose = TRUE)^2)
  tau <- rgamma(1, prior$tau[1] + length(th) / 2, prior$tau[2] + q / 2)
  list(phi = phi, tau = tau, L = t(U) / sqrt(tau), U = U)
}

shift_step <- function(mu, th, U, tau, sd_mu) {
  a <- backsolve(U, rep(1, length(th)), transpose = TRUE)
  b <- backsolve(U, th, transpose = TRUE)
  P <- 1 / sd_mu^2 + tau * sum(a^2)
  c <- rnorm(1, (-mu / sd_mu^2 + tau * sum(a * b)) / P, 1 / sqrt(P))
  list(mu = mu + c, th = th - c)
}

# llfun(log r, log lambda0) -> log-verosimilitud; s = estado inicial
# true_th = list(r = theta_r real, l = theta_lambda real): añade la comparación de campos en los gráficos en vivo
run_mcmc <- function(llfun, D, s, prior, n_iter, burn, thin, every = 10000, true = NULL, n_last = 500,
                     shift = TRUE, true_th = NULL, nc = FALSE) {
  nc_update <- function(th, phi, tau, ll_of) {
    z <- sqrt(tau) * backsolve(chol(Rmat(D, phi)), th, transpose = TRUE)
    fld <- function(p, t) drop(crossprod(chol(Rmat(D, p)), z)) / sqrt(t)
    phi <- slice1(phi, function(p) ll_of(fld(p, tau)), 0.05, prior$phi[1], prior$phi[2])
    lt <- slice1(log(tau), function(x) ll_of(fld(phi, exp(x))) +
                   dgamma(exp(x), prior$tau[1], prior$tau[2], log = TRUE) + x, 0.3)
    list(th = fld(phi, exp(lt)), phi = phi, tau = exp(lt))
  }
  BR <- BL <- matrix(NA, n_last, nrow(D))
  vb <- vario_bins(D)
  last <- matrix(NA, n_last, 8, dimnames = list(NULL, c("mu_r", "mu_l", "tau_r", "tau_l", "phi_r", "phi_l", "delta_r", "delta_l")))
  m <- nrow(D)
  Lr <- t(chol(Rmat(D, s$phi_r))) / sqrt(s$tau_r)
  Ll <- t(chol(Rmat(D, s$phi_l))) / sqrt(s$tau_l)
  ll <- llfun(s$mu_r + s$th_r, s$mu_l + s$th_l)
  nk <- (n_iter - burn) %/% thin
  par <- matrix(NA, nk, 6, dimnames = list(NULL, c("mu_r", "mu_l", "tau_r", "tau_l", "phi_r", "phi_l")))
  TR <- TL <- matrix(NA, nk, m); j <- 0
  lpm <- function(x) dnorm(x, 0, prior$sd_mu, log = TRUE)
  for (it in seq_len(n_iter)) {
    e <- ess(s$th_r, Lr, function(f) llfun(s$mu_r + f, s$mu_l + s$th_l), ll); s$th_r <- e$f
    e <- ess(s$th_l, Ll, function(f) llfun(s$mu_r + s$th_r, s$mu_l + f), e$ll); s$th_l <- e$f
    s$mu_r <- slice1(s$mu_r, function(x) llfun(x + s$th_r, s$mu_l + s$th_l) + lpm(x), 0.1)
    s$mu_l <- slice1(s$mu_l, function(x) llfun(s$mu_r + s$th_r, x + s$th_l) + lpm(x), 0.1)
    ll <- llfun(s$mu_r + s$th_r, s$mu_l + s$th_l)
    h <- upd_hyper(s$th_r, s$phi_r, s$tau_r, D, prior)
    s$phi_r <- h$phi
    s$tau_r <- h$tau
    Lr <- h$L
    if(shift){ z <- shift_step(s$mu_r, s$th_r, h$U, s$tau_r, prior$sd_mu)
                 s$mu_r <- z$mu
                 s$th_r <- z$th
    }
    if(nc){
      z <- nc_update(s$th_r, s$phi_r, s$tau_r, function(f) llfun(s$mu_r + f, s$mu_l + s$th_l))
      s$th_r <- z$th
      s$phi_r <- z$phi
      s$tau_r <- z$tau
      Lr <- t(chol(Rmat(D, s$phi_r))) / sqrt(s$tau_r)
    }
    h <- upd_hyper(s$th_l, s$phi_l, s$tau_l, D, prior)
    s$phi_l <- h$phi
    s$tau_l <- h$tau
    Ll <- h$L
    if(shift){z <- shift_step(s$mu_l, s$th_l, h$U, s$tau_l, prior$sd_mu)
              s$mu_l <- z$mu
              s$th_l <- z$th 
    }
    if(nc){
      z <- nc_update(s$th_l, s$phi_l, s$tau_l, function(f) llfun(s$mu_r + s$th_r, s$mu_l + f))
      s$th_l <- z$th
      s$phi_l <- z$phi
      s$tau_l <- z$tau
      Ll <- t(chol(Rmat(D, s$phi_l))) / sqrt(s$tau_l)
      ll <- llfun(s$mu_r + s$th_r, s$mu_l + s$th_l)
    }
    if (it > burn && (it - burn) %% thin == 0) {
      j <- j + 1
      par[j, ] <- c(s$mu_r, s$mu_l, s$tau_r, s$tau_l, s$phi_r, s$phi_l)
      TR[j, ] <- s$th_r
      TL[j, ] <- s$th_l
    }
    last[(it - 1) %% n_last + 1, ] <- c(s$mu_r, s$mu_l, s$tau_r, s$tau_l, s$phi_r, s$phi_l,
                                        delta_fun(s$phi_r, s$tau_r), delta_fun(s$phi_l, s$tau_l))
    BR[(it - 1) %% n_last + 1, ] <- s$th_r
    BL[(it - 1) %% n_last + 1, ] <- s$th_l
    if (it %% every == 0) {
      message(sprintf("it %d/%d | loglik %.2f | mu_r %.3f mu_l %.3f | tau_r %.2f tau_l %.2f | phi_r %.3f phi_l %.3f | delta_r %.2e delta_l %.2e",
                      it, n_iter, ll, s$mu_r, s$mu_l, s$tau_r, s$tau_l, s$phi_r, s$phi_l,
                      delta_fun(s$phi_r, s$tau_r), delta_fun(s$phi_l, s$tau_l)))
      nb <- seq_len(min(it, n_last))
      live_plot(last[nb, , drop = FALSE], true, it,
                list(r = BR[nb, , drop = FALSE], l = BL[nb, , drop = FALSE]), true_th, vb)
    }
  }
  list(par = par, th_r = TR, th_l = TL)
}

ID_PARS <- c("mu_r", "mu_l", "delta_r", "delta_l")

TH_LAB <- c(r = "theta_r", l = "theta_lambda")

vario_bins <- function(D, nb = 12) {
  ut <- upper.tri(D)
  d <- D[ut]
  b <- cut(d, seq(0, max(D) / 2, length.out = nb + 1))
  keep <- !is.na(b)
  ut[ut] <- keep; d <- d[keep]
  b <- droplevels(b[keep])
  list(ut = ut, b = b, h = as.vector(tapply(d, b, mean)))
}
emp_vario <- function(th, vb) as.vector(tapply(0.5 * outer(th, th, "-")[vb$ut]^2, vb$b, mean))

live_plot <- function(X, true, it, B = NULL, true_th = NULL, vb = NULL) {
  full <- !is.null(true_th)
  op <- par(mfrow = if (full) c(4, 4) else c(2, 2), mar = c(4, 3.5, 2.5, 1), mgp = c(2, 0.7, 0))
  on.exit(par(op))
  live_densities(X, true, it)
  if (!full) return(invisible())
  cen <- function(x) x - mean(x)
  ns <- sqrt(ncol(B$r))
  pal <- hcl.colors(64, "viridis")
  for (k in c("r", "l")) {                       
    tr <- cen(true_th[[k]]); es <- cen(colMeans(B[[k]])); zl <- range(tr, es)
    image(matrix(tr, ns), zlim = zl, col = pal, axes = FALSE, main = paste(TH_LAB[k], "true"))
    image(matrix(es, ns), zlim = zl, col = pal, axes = FALSE, main = paste(TH_LAB[k], "estimated (mean)"))
  }
  for (k in c("r", "l")) {                      
    tr <- cen(true_th[[k]]); es <- cen(colMeans(B[[k]]))
    plot(tr, es, pch = 16, cex = 0.5, xlab = "True (centred)", ylab = "Estimated (centred)",
         main = sprintf("%s  cor = %.2f", TH_LAB[k], cor(tr, es)))
    abline(0, 1, col = "red", lwd = 2)
  }
  for (k in c("r", "l")) {                       
    vt <- emp_vario(true_th[[k]], vb)
    sel <- unique(round(seq(1, nrow(B[[k]]), length.out = min(50, nrow(B[[k]])))))
    ve <- sapply(sel, function(i) emp_vario(B[[k]][i, ], vb))
    q <- apply(ve, 1, quantile, c(.05, .5, .95), na.rm = TRUE)
    plot(vb$h, vt, type = "n", ylim = range(0, vt, q, na.rm = TRUE), xlab = "Distance", ylab = "Semivariance",
         main = paste(TH_LAB[k], "semivariogram"))
    polygon(c(vb$h, rev(vb$h)), c(q[1, ], rev(q[3, ])), col = "mistyrose", border = NA)
    lines(vb$h, q[2, ], col = "red", lwd = 2)
    lines(vb$h, vt, type = "b", pch = 16, lwd = 2)
    legend("bottomright", c("True", "Draws (median, 90%)"), col = c("black", "red"), lwd = 2, bty = "n", cex = 0.8)
  }
  for (k in c("r", "l")) {                       # error relativo de la varianza de cada iteración
    v <- apply(B[[k]], 1, var)
    vt <- var(true_th[[k]])
    hist(v, breaks = 30, col = "grey85", border = "grey60", xlim = range(v, vt), xlab = "Sample variance of field",
         main = paste(TH_LAB[k], "variance"))
    abline(v = vt, col = "red", lwd = 2)
  }
}

live_densities <- function(X, true, it) {
  for (p in ID_PARS) {
    x <- X[, p]
    x <- x[is.finite(x)]
    if (length(unique(x)) < 2) { plot.new(); title(p); next}
    d <- density(x)
    xr <- range(d$x, if (!is.null(true) && p %in% names(true)) true[[p]])
    plot(d, xlim = xr, main = p, xlab = sprintf("Last %d iterations (it %d)", nrow(X), it), ylab = "")
    if (!is.null(true) && p %in% names(true)) abline(v = true[[p]], col = "red", lwd = 2)
  }
}

post_summary <- function(fit) {
  p <- fit$par
  cbind(p, delta_r = delta_fun(p[, "phi_r"], p[, "tau_r"]), delta_l = delta_fun(p[, "phi_l"], p[, "tau_l"]))
}
post_table <- function(post, true) {
  v <- c("mu_r", "mu_l", "delta_r", "delta_l")
  q <- apply(post[, v], 2, quantile, c(.5, .025, .975))
  data.frame(true = true[v], median = q[1, ], low = q[2, ], high = q[3, ])
}

# ---- Proceso predictivo: theta~ y derivadas en G a partir de theta* en los knots K ----
pp_pred <- function(th, phi, D, DX, DY, esc) {
  w <- solve(Rmat(D, phi), th)
  a <- mrate(phi)
  d <- esc * sqrt(DX^2 + DY^2)
  e <- exp(-a * d)
  g <- -a^2 * esc^2 * e
  list(f = drop(((1 + a * d) * e) %*% w), fx = drop((g * DX) %*% w), fy = drop((g * DY) %*% w))
}

vel_summary <- function(fit, K, G, tvel, esc, dfun, ndraw = nrow(fit$par), sc = c(1, 1, 1)) {
  D <- esc * dist_mat(K)
  DX <- outer(G[, 1], K[, 1], "-"); DY <- outer(G[, 2], K[, 2], "-")
  vars <- c("lam", "lt", "lx", "ly", "v")
  S <- array(0, c(nrow(G), length(vars), length(tvel)), dimnames = list(NULL, vars, NULL))
  idx <- unique(round(seq(1, nrow(fit$par), length.out = ndraw)))
  for (i in idx) {
    p <- fit$par[i, ]
    a <- pp_pred(fit$th_r[i, ], p[["phi_r"]], D, DX, DY, esc)
    b <- pp_pred(fit$th_l[i, ], p[["phi_l"]], D, DX, DY, esc)
    r <- exp(p[["mu_r"]] + a$f); l0 <- exp(p[["mu_l"]] + b$f)
    for (k in seq_along(tvel)) {
      d <- dfun(r, r * a$fx, r * a$fy, l0, l0 * b$fx, l0 * b$fy, tvel[k])
      lt <- d$lt / sc[3]
      lx <- d$lx / sc[1]
      ly <- d$ly / sc[2]
      S[, , k] <- S[, , k] + cbind(d$lam, lt, lx, ly, lt / sqrt(lx^2 + ly^2))
    }
  }
  S <- S / length(idx)
  out <- do.call(rbind, lapply(seq_along(tvel), function(k)
    data.frame(x = G[, 1], y = G[, 2], t = tvel[k], S[, , k])))
  nrm <- sqrt(out$lx^2 + out$ly^2)
  out$ux <- sign(out$lt) * out$lx / nrm
  out$uy <- sign(out$lt) * out$ly / nrm
  out
}

#######GRAPHS
tl <- labeller(t = function(v) paste("t =", v))

sim_plot <- function(Y, tt, K, title) {
  p1 <- ggplot(data.frame(t = tt, n = colSums(Y)), aes(t, n)) + geom_col(width = diff(tt)[1]) +
    labs(x = "Time", y = "Number of events")
  p2 <- ggplot(data.frame(x = K[, 1], y = K[, 2], n = rowSums(Y)), aes(x, y, fill = n)) +
    geom_tile() + scale_fill_viridis_c(name = "Events") + coord_equal()
  p1 + p2 + plot_annotation(title = title)
}

field_plot <- function(K, fields, title) {
  d <- do.call(rbind, lapply(names(fields), function(n)
    data.frame(x = K[, 1], y = K[, 2], value = fields[[n]], panel = n)))
  d$panel <- factor(d$panel, names(fields))
  ggplot(d, aes(x, y, fill = value)) + geom_tile() + facet_wrap(~ panel) +
    scale_fill_viridis_c() + coord_equal() + labs(title = title)
}

post_plot <- function(post, true) {
  post <- post[, ID_PARS]
  d <- data.frame(par = rep(colnames(post), each = nrow(post)), value = c(post))
  p <- ggplot(d, aes(value)) + geom_density(fill = "grey80") +
    facet_wrap(~ par, scales = "free") + labs(x = NULL, y = "Posterior density")
  if (length(true)) {
    tr <- data.frame(par = names(true), value = true)
    p <- p + geom_vline(data = tr[tr$par %in% d$par, ], aes(xintercept = value), colour = "red")
  }
  p
}

tlab_default <- function(t) paste("t =", t)
per_t <- function(df, panel, title, ncol = 2)
  wrap_plots(lapply(split(df, df$t), panel), ncol = ncol) + plot_annotation(title = title)

map_plot <- function(df, var, title, lab = var, tlab = tlab_default)
  per_t(df, function(d) ggplot(d, aes(x, y, fill = .data[[var]])) + geom_raster() +
          scale_fill_viridis_c(name = lab) + coord_equal() +
          labs(subtitle = tlab(d$t[1]), x = "x", y = "y"), title)

cmp_plot <- function(df, title, tlab = tlab_default)
  per_t(df, function(d) ggplot(d, aes(x, y, fill = lam)) + geom_tile(aes(width = w, height = w)) +
          facet_wrap(~ type, ncol = 1) + scale_fill_viridis_c(name = expression(lambda)) +
          coord_equal() + labs(subtitle = tlab(d$t[1]), x = "x", y = "y"),
        title, ncol = length(unique(df$t)))

arrow_plot <- function(df, title, step = 4, len = 0.025, tlab = tlab_default) {
  xs <- sort(unique(df$x))
  s <- xs[seq(2, length(xs), by = step)]
  per_t(df, function(d) {
    a <- d[d$x %in% s & d$y %in% s, ]
    ggplot(d, aes(x, y)) + geom_raster(aes(fill = lam)) +
      geom_segment(data = a, aes(xend = x + len * ux, yend = y + len * uy),
                   arrow = arrow(length = unit(0.08, "cm")), linewidth = 0.3) +
      scale_fill_viridis_c(name = expression(lambda)) + coord_equal() +
      labs(subtitle = tlab(d$t[1]), x = "x", y = "y")
  }, title)
}

series_plot <- function(est, truth, tt, cells, ylab, title) {
  q <- apply(est, c(2, 3), quantile, c(.025, .5, .975))
  d <- data.frame(cell = factor(rep(cells, times = length(tt))), t = rep(tt, each = length(cells)),
                  low = c(q[1, , ]), med = c(q[2, , ]), high = c(q[3, , ]), true = c(truth))
  ggplot(d, aes(t)) + geom_line(aes(y = low), colour = "green4", linetype = 2) +
    geom_line(aes(y = high), colour = "green4", linetype = 2) +
    geom_line(aes(y = true)) + geom_line(aes(y = med), colour = "red") +
    facet_wrap(~ cell, labeller = label_both) + labs(x = "Time", y = ylab, title = title)
}

pick_typical <- function(S) {
  s <- apply(S, 2, mad); s[s == 0] <- 1
  Z <- sweep(sweep(S, 2, apply(S, 2, median)), 2, s, "/")
  which.min(rowSums(Z^2))
}

gp_stats <- function(th, D, lags)
  c(mean = mean(th), var = var(th),
    setNames(sapply(lags, function(h) {
      ij <- which(abs(D - h) < 1e-8, arr.ind = TRUE)
      0.5 * mean((th[ij[, 1]] - th[ij[, 2]])^2)
    }), paste0("semivar_lag", seq_along(lags))))

pick_gp <- function(n, D, phi, tau, lags) {
  C <- replicate(n, rgp(D, phi, tau))
  S <- t(apply(C, 2, gp_stats, D = D, lags = lags))
  k <- pick_typical(S)
  list(th = C[, k], S = S, k = k)
}

pick_counts <- function(n, rsim, stats) {
  cand <- lapply(seq_len(n), function(i) rsim())
  S <- t(sapply(cand, stats))
  k <- pick_typical(S)
  list(x = cand[[k]], S = S, k = k)
}

cand_plot <- function(S, k, title, ref = NULL) {
  d  <- data.frame(stat = rep(colnames(S), each = nrow(S)), value = c(S))
  ch <- data.frame(stat = colnames(S), value = S[k, ])
  p <- ggplot(d, aes(value)) + geom_histogram(bins = 30, fill = "grey80", colour = "grey50") +
    geom_vline(data = ch, aes(xintercept = value), colour = "red") +
    facet_wrap(~ stat, scales = "free") + labs(title = title, x = NULL, y = "Count")
  if (!is.null(ref))
    p <- p + geom_vline(data = data.frame(stat = names(ref), value = ref), aes(xintercept = value), linetype = 2)
  p
}

gp_ml <- function(th, D, lo = 0.005, hi = 1) {
  n <- length(th)
  q <- function(phi){ 
    U <- chol(Rmat(D, phi))
    list(U = U, q = sum(backsolve(U, th, transpose = TRUE)^2)) 
  }
  prof <- function(lp){
    z <- q(exp(lp))
    -sum(log(diag(z$U))) - 0.5 * n * log(z$q / n)
  }
  phi <- exp(optimize(prof, log(c(lo, hi)), maximum = TRUE)$maximum)
  tau <- n / q(phi)$q
  c(tau = tau, phi = phi, delta = delta_fun(phi, tau))
}

pick_gp_ml <- function(n, D, phi, tau) {
  C <- replicate(n, rgp(D, phi, tau))
  S <- t(apply(C, 2, gp_ml, D = D))
  err <- sweep(log(S), 2, log(c(tau, phi, delta_fun(phi, tau))))
  k <- which.min(rowSums(err^2))
  list(th = C[, k], S = S, k = k, ref = c(tau = tau, phi = phi, delta = delta_fun(phi, tau)))
}

vario_plot <- function(fit, true_th, D, title) {
  vb <- vario_bins(D)
  d <- do.call(rbind, lapply(c("r", "l"), function(k) {
    th <- if (k == "r") fit$th_r else fit$th_l
    ve <- apply(th, 1, emp_vario, vb = vb)
    q <- apply(ve, 1, quantile, c(.025, .5, .975), na.rm = TRUE)
    data.frame(field = TH_LAB[[k]], h = vb$h, low = q[1, ], med = q[2, ], high = q[3, ],
               true = emp_vario(true_th[[k]], vb))
  }))
  ggplot(d, aes(h)) + geom_ribbon(aes(ymin = low, ymax = high), fill = "mistyrose") +
    geom_line(aes(y = med), colour = "red") + geom_line(aes(y = true)) + geom_point(aes(y = true)) +
    facet_wrap(~ field, scales = "free_y") +
    labs(x = "Distance", y = "Semivariance", title = title,
         subtitle = "Black: true field; red: posterior median with 95% band")
}