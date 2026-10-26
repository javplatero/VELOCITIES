setwd("~/Alan")
source("funciones_comunes.R")
options(matern_param = "range")   # Matérn 3/2, phi = rango; delta = tau * phi^3
set.seed(2025)
dir.create("Graphs/Section4", recursive = TRUE, showWarnings = FALSE)

NS <- 20
NT <- 200
NG <- 80
GAM <- 1
ESC <- 1

BERN_VOL <- TRUE
N_ITER <- 200000
BURN <- 100000
THIN <- 100

true <- c(mu_r = 0.84, mu_l = 7.75, tau_r = 9.4, tau_l = 11.4, phi_r = 0.18, phi_l = 0.17)
true <- c(true, delta_r = delta_fun(true[["phi_r"]], true[["tau_r"]]), delta_l = delta_fun(true[["phi_l"]], true[["tau_l"]]))
prior <- list(sd_mu = 10, tau = c(10, 1), phi = c(0.12, 0.2))

K <- grid_centers(NS)
D <- ESC * dist_mat(K)
G <- grid_centers(NG)
tt <- (seq_len(NT) - 0.5) / NT
V <- (1 / NS^2) * (1 / NT)
tvel <- c(0.2525, 0.5025, 0.7525, 0.9975)
VB <- if (BERN_VOL) V else 1

lam_lin <- function(lr, ll, tt) {
  r <- exp(lr)
  E <- exp(outer(r, tt))
  E * (exp(ll) + GAM / r) - GAM / r
}
d_lin <- function(r, rx, ry, l0, l0x, l0y, t) {
  ert <- exp(r * t); lam <- ert * (l0 + GAM / r) - GAM / r
  c1 <- t * (lam + GAM / r) + GAM * (1 - ert) / r^2
  list(lam = lam, lt = r * lam + GAM, lx = rx * c1 + ert * l0x, ly = ry * c1 + ert * l0y)
}

#######SImulation
N_CAND <- 500

gr <- pick_gp_ml(N_CAND, D, true[["phi_r"]], true[["tau_r"]])
th_r <- gr$th
gl <- pick_gp_ml(N_CAND, D, true[["phi_l"]], true[["tau_l"]])
th_l <- gl$th
print(rbind(theta_r = gr$S[gr$k, ], true_r = gr$ref, theta_l = gl$S[gl$k, ], true_l = gl$ref), digits = 4)

lam_true <- lam_lin(true[["mu_r"]] + th_r, true[["mu_l"]] + th_l, tt)
P <- lam_true * VB / (1 + lam_true * VB)
blk <- cut(seq_len(NT), 4, labels = FALSE)
y_stats <- function(Y)
  c(total = sum(Y), setNames(as.vector(tapply(colSums(Y), blk, sum)), paste0("block", 1:4)),
    spatial_cor = cor(rowSums(Y), rowSums(P)))
gy <- pick_counts(N_CAND, function() matrix(rbinom(length(P), 1, P), nrow(P)), y_stats); Y <- gy$x
print(paste0("Events: ", sum(Y), " / ", length(Y), " cells"))

####Ajuste
llfun <- function(lr, ll){ 
  m <- lam_lin(lr, ll, tt) * VB
  sum(Y * log(m) - log1p(m))
}
init <- list(mu_r = 1, mu_l = log((sum(Y[, 1:20]) + 1) / (NS^2 * 20 * VB)),
             th_r = rep(0, NS^2), th_l = rep(0, NS^2),
             tau_r = 10, tau_l = 10, phi_r = 0.15, phi_l = 0.15)
fit <- run_mcmc(llfun, D, init, prior, N_ITER, BURN, THIN, every = 200, true = true,
                true_th = list(r = th_r, l = th_l))
post <- post_summary(fit)
print(post_table(post, true), digits = 4)

####Velocities
vel <- vel_summary(fit, K, G, tvel, ESC, d_lin)
vel$vcap <- pmin(vel$v, quantile(vel$v, .975, na.rm = TRUE))

######Probabilities column 9
col9 <- which(round(K[, 1] * NS + 0.5) == 9)
pest <- array(NA, c(nrow(post), length(col9), NT))
for (i in seq_len(nrow(post))) {
  m <- lam_lin(post[i, "mu_r"] + fit$th_r[i, col9], post[i, "mu_l"] + fit$th_l[i, col9], tt) * VB
  pest[i, , ] <- m / (1 + m)
}

####PLOTS
sv <- function(p, f, w = 8, h = 6) ggsave(file.path("Graphs/Section4", f), p, width = w, height = h, dpi = 300)
sv(cand_plot(gr$S, gr$k, "theta_r candidates: ML estimates from the field (red = selected, dashed = true)", gr$ref), "CandidatesThetaR.png", 9, 6)
sv(cand_plot(gl$S, gl$k, "theta_lambda candidates: ML estimates from the field (red = selected, dashed = true)", gl$ref), "CandidatesThetaL.png", 9, 6)
sv(cand_plot(gy$S, gy$k, "Simulated data candidates (red = selected)"), "CandidatesData.png", 9, 6)
sv(sim_plot(Y, tt, K, "Simulated presences"), "SimulacionGEOS.png", 10, 4.5)
sv(field_plot(K, list(theta_r = th_r, theta_lambda = th_l), "Simulated spatial fields"), "SimulacionCampos.png", 8, 4)
sv(post_plot(post, true), "MarginalPosteriors.png", 10, 7)
sv(field_plot(K, list(`theta_r (true)` = th_r, `theta_r (estimated)` = colMeans(fit$th_r),
                      `theta_lambda (true)` = th_l, `theta_lambda (estimated)` = colMeans(fit$th_l)),
              "True vs estimated fields"), "ComparacionCampos.png", 8, 8)
sv(vario_plot(fit, list(r = th_r, l = th_l), D, "Empirical semivariogram of the fields"), "SemivariogramFields.png", 10, 4.5)
sv(series_plot(pest, P[col9, ], tt, seq_along(col9), "Probability",
               "True (black) vs estimated (red) probabilities, column 9"), "EstimacionProbabilidades.png", 10, 8)
lt_true <- do.call(rbind, lapply(tvel, function(t)
  data.frame(x = K[, 1], y = K[, 2], t = t, w = 1 / NS, type = "True",
             lam = drop(lam_lin(true[["mu_r"]] + th_r, true[["mu_l"]] + th_l, t)))))
sv(cmp_plot(rbind(lt_true, data.frame(vel[, c("x", "y", "t", "lam")], w = 1 / NG, type = "Estimated")),
            "True (15x15) vs estimated (80x80) intensity"), "ComparisonRealEstimate.png", 12, 6)
sv(map_plot(vel, "lt", "Temporal derivative", "d lambda/dt"), "DERT.png")
sv(map_plot(vel, "lx", "Partial derivative in x", "d lambda/dx"), "DERX.png")
sv(map_plot(vel, "ly", "Partial derivative in y", "d lambda/dy"), "DERY.png")
sv(map_plot(vel, "vcap", "Minimum velocity (capped at 97.5th percentile)", "Velocity"), "Velocities.png")
sv(arrow_plot(vel, "Estimated intensity and directions of minimum velocity"), "IntensityAndDirections.png")

saveRDS(list(th_r = th_r, th_l = th_l, Y = Y, P = P, fit = fit, vel = vel,
             post = post, true = true,
             cand_stats = list(th_r = gr$S, th_l = gl$S, Y = gy$S), chosen = c(th_r = gr$k, th_l = gl$k, Y = gy$k)),
        "res_geoestadistico.rds")