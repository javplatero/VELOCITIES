setwd("~/Alan")
source("funciones_comunes.R")
options(matern_param = "range")   # Matérn 3/2 del artículo, phi = rango; delta = tau * phi^3
set.seed(123456)
dir.create("Graphs/Section5", recursive = TRUE, showWarnings = FALSE)

##SETTINGS SIM
NS <- 20
NT <- 200
NG <- 80
T0 <- 0.25
ESC <- 1

#MCMC SETTINGS
N_ITER <- 200000
BURN <- 100000
THIN <- 100

true <- c(mu_r = 0.41, mu_l = 8.14, tau_r = 7.5, tau_l = 12.7, phi_r = 0.17, phi_l = 0.19)
true <- c(true, delta_r = delta_fun(true[["phi_r"]], true[["tau_r"]]), delta_l = delta_fun(true[["phi_l"]], true[["tau_l"]]))
prior <- list(sd_mu = 10, tau = c(10, 1), phi = c(0.12, 0.2))

K <- grid_centers(NS)
D <- ESC * dist_mat(K)
G <- grid_centers(NG)
tt <- T0 + (seq_len(NT) - 0.5) * (1 - T0) / NT
Am <- 1 / NS^2
V <- Am * (1 - T0) / NT
tvel <- c(0.321875, 0.509375, 0.753125, 0.996875)   # las velocidades se evalúan en cualquier t

loglam_weib <- function(lr, ll, tt) {
  r <- exp(lr)
  ll + outer(r, log(tt / T0)) - outer(r + 1, tt, function(a, b) b^a) + T0^(r + 1)
}
d_weib <- function(r, rx, ry, l0, l0x, l0y, t) {
  lam <- l0 * (t / T0)^r * exp(-(t^(r + 1) - T0^(r + 1)))
  c1 <- (T0^(r + 1) - 1) * log(T0) - (t^(r + 1) - 1) * log(t)
  list(lam = lam, lt = (r / t - (r + 1) * t^r) * lam,
       lx = lam * (l0x / l0 + rx * c1), ly = lam * (l0y / l0 + ry * c1))
}

###Simulation
N_CAND <- 500
lags <- ESC * c(1, 2, 4, 8) / NS          # 1, 2, 4 y 8 celdas
gr <- pick_gp(N_CAND, D, true[["phi_r"]], true[["tau_r"]], lags)
th_r <- gr$th
gl <- pick_gp(N_CAND, D, true[["phi_l"]], true[["tau_l"]], lags)
th_l <- gl$th
lam_true <- exp(loglam_weib(true[["mu_r"]] + th_r, true[["mu_l"]] + th_l, tt))
l0_true <- exp(true[["mu_l"]] + th_l)
blk <- cut(seq_len(NT), 4, labels = FALSE)
n_stats <- function(z)
  c(initial = sum(z$N0), total = sum(z$N),
    setNames(as.vector(tapply(colSums(z$N), blk, sum)), paste0("block", 1:4)),
    spatial_cor = cor(rowSums(z$N), rowSums(lam_true)), initial_cor = cor(z$N0, l0_true))
gn <- pick_counts(N_CAND, function()
  list(N0 = rpois(NS^2, l0_true * Am), N = matrix(rpois(length(lam_true), lam_true * V), nrow(lam_true))),
  n_stats)
N0 <- gn$x$N0
N <- gn$x$N
print(paste0("Points: ", sum(N0) + sum(N), " (initial pattern: ", sum(N0), ")"))
# ---- Serie temporal simulada ----
ts_df <- data.frame(t = tt, observed = colSums(N), expected = colSums(lam_true) * V)
p_ts <- ggplot(ts_df, aes(t)) +
  geom_col(aes(y = observed), width = diff(tt)[1], fill = "grey75") +
  geom_line(aes(y = expected), colour = "red", linewidth = 0.8) +
  labs(x = "Time", y = "Number of events",
       title = "Simulated events per time bin (grey) vs expected (red)")
print(p_ts)
ggsave("Graphs/Section5/SerieTemporalSimulada.png", p_ts, width = 9, height = 5, dpi = 300)

#######Ajuste
llfun <- function(lr, ll) {
  lm <- loglam_weib(lr, ll, tt)
  sum(N * lm - exp(lm) * V) + sum(N0 * ll - exp(ll) * Am)
}
init <- list(mu_r = 0, mu_l = log(mean(N0) / Am),
             th_r = rep(0, NS^2), th_l = rep(0, NS^2),
             tau_r = 10, tau_l = 10, phi_r = 0.15, phi_l = 0.15)
fit <- run_mcmc(llfun, D, init, prior, N_ITER, BURN, THIN, every = 200, true = true,
                true_th = list(r = th_r, l = th_l))
post <- post_summary(fit)
print(post_table(post, true), digits = 4)

###VELOCITIES
vel <- vel_summary(fit, K, G, tvel, ESC, d_weib)
q <- quantile(vel$v, c(.025, .975), na.rm = TRUE)
vel$vcap <- pmax(pmin(vel$v, q[2]), q[1])

####INTENSITIES ROW 14
row14 <- which(round(K[, 2] * NS + 0.5) == 14)
lest <- array(NA, c(nrow(post), length(row14), NT))
for (i in seq_len(nrow(post)))
  lest[i, , ] <- exp(loglam_weib(post[i, "mu_r"] + fit$th_r[i, row14], post[i, "mu_l"] + fit$th_l[i, row14], tt))

##PLOTS
sv <- function(p, f, w = 8, h = 6) ggsave(file.path("Graphs/Section5", f), p, width = w, height = h, dpi = 300)
sv(cand_plot(gr$S, gr$k, "theta_r candidates (red = selected)"), "CandidatesThetaR.png", 9, 6)
sv(cand_plot(gl$S, gl$k, "theta_lambda candidates (red = selected)"), "CandidatesThetaL.png", 9, 6)
sv(cand_plot(gn$S, gn$k, "Simulated data candidates (red = selected)"), "CandidatesData.png", 9, 6)
sv(sim_plot(N, tt, K, "Simulated point pattern"), "SimulacionWeibullPP.png", 10, 4.5)
sv(field_plot(K, list(theta_r = th_r, theta_lambda = th_l), "Simulated spatial fields"), "SimulacionCamposW.png", 8, 4)
sv(post_plot(post, true), "PosteriorsWeibull.png", 10, 7)
sv(field_plot(K, list(`theta_r (true)` = th_r, `theta_r (estimated)` = colMeans(fit$th_r),
                      `theta_lambda (true)` = th_l, `theta_lambda (estimated)` = colMeans(fit$th_l)),
              "True vs estimated fields"), "ComparacionCamposW.png", 8, 8)
sv(vario_plot(fit, list(r = th_r, l = th_l), D, "Empirical semivariogram of the fields"), "SemivariogramFieldsW.png", 10, 4.5)
sv(series_plot(lest, lam_true[row14, ], tt, seq_along(row14), expression(lambda),
               "True (black) vs estimated (red) intensity, row 14"), "SeriesTemporalesFil14.png", 10, 8)
lt_true <- do.call(rbind, lapply(tvel, function(t)
  data.frame(x = K[, 1], y = K[, 2], t = t, w = 1 / NS, type = "True",
             lam = exp(drop(loglam_weib(true[["mu_r"]] + th_r, true[["mu_l"]] + th_l, t))))))
sv(cmp_plot(rbind(lt_true, data.frame(vel[, c("x", "y", "t", "lam")], w = 1 / NG, type = "Estimated")),
            "True (15x15) vs estimated (80x80) intensity"), "ComparacionIntensidades.png", 12, 6)
sv(map_plot(vel, "lt", "Temporal derivative", "d lambda/dt"), "DERT.png")
sv(map_plot(vel, "lx", "Partial derivative in x", "d lambda/dx"), "DERX.png")
sv(map_plot(vel, "ly", "Partial derivative in y", "d lambda/dy"), "DERY.png")
sv(map_plot(vel, "vcap", "Minimum velocity (2.5th-97.5th percentiles)", "Velocity"), "Velocities.png")
sv(arrow_plot(vel, "Estimated intensity and directions of minimum velocity"), "IntensitiesDirectionsW.png")

saveRDS(list(th_r = th_r, th_l = th_l, N0 = N0, N = N, lam = lam_true, fit = fit, vel = vel,
             post = post, true = true,
             cand_stats = list(th_r = gr$S, th_l = gl$S, N = gn$S), chosen = c(th_r = gr$k, th_l = gl$k, N = gn$k)),
        "res_patron_puntual.rds")


###SAMES times AS GEOSTATISTIC
tgeo <- c(0.2525, 0.5025, 0.7525, 0.9975)
brk <- T0 + (0:NT) * (1 - T0) / NT
tvel <- tt[findInterval(tgeo, brk, rightmost.closed = TRUE)]
print(rbind(tgeo, tvel))
vel <- vel_summary(fit, K, G, tvel, ESC, d_weib)
q <- quantile(vel$v, c(.025, .975), na.rm = TRUE)
vel$vcap <- pmax(pmin(vel$v, q[2]), q[1])
vel$t <- tgeo[match(vel$t, tvel)]                      # etiqueta t = valor del geoestadístico
lt_true <- do.call(rbind, lapply(seq_along(tvel), function(k)
  data.frame(x = K[, 1], y = K[, 2], t = tgeo[k], w = 1 / NS, type = "True",
             lam = exp(drop(loglam_weib(true[["mu_r"]] + th_r, true[["mu_l"]] + th_l, tvel[k]))))))
sv(cmp_plot(rbind(lt_true, data.frame(vel[, c("x", "y", "t", "lam")], w = 1 / NG, type = "Estimated")),
            "True vs estimated intensity"), "ComparacionIntensidades.png", 12, 6)
sv(map_plot(vel, "lt", "Temporal derivative", "d lambda/dt"), "DERT.png")
sv(map_plot(vel, "lx", "Partial derivative in x", "d lambda/dx"), "DERX.png")
sv(map_plot(vel, "ly", "Partial derivative in y", "d lambda/dy"), "DERY.png")
sv(map_plot(vel, "vcap", "Minimum velocity (2.5th-97.5th percentiles)", "Velocity"), "Velocities.png")
sv(arrow_plot(vel, "Estimated intensity and directions of minimum velocity"), "IntensitiesDirectionsW.png")

