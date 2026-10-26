setwd("~/Alan")
source("funciones_comunes.R")
options(matern_param = "range")   # Matérn 3/2 del artículo, phi = rango; delta = tau * phi^3

library(sf)
sf_use_s2(FALSE)
set.seed(35649)
dir.create("Graphs/Section6", recursive = TRUE, showWarnings = FALSE)

NS <- 20                           
T0 <- 1/3
FECHA_INI <- as.Date("2020-03-02"); FECHA_FIN <- as.Date("2020-07-25")   # primer día t = 0, último t = 1
dias <- as.Date(c("2020-04-19", "2020-05-29", "2020-07-08", "2020-07-23"))
N_ITER <- 250000
BURN <- 125000
THIN <- 100
prior <- list(sd_mu = 10, tau = c(10, 1), phi = c(0.02, 0.3))
NG <- 80
ESC <- 1
COL_FECHA <- NULL                  # nombre de la columna de fecha en DatosCali; si es NULL se usa TiempoC tal cual

# ---- Ventana ----
borde <- st_union(st_read("borde_cali.shp", quiet = TRUE))
bb <- st_bbox(borde)
cc <- st_coordinates(borde)
poly01 <- data.frame(x = (cc[, 1] - bb[["xmin"]]) / (bb[["xmax"]] - bb[["xmin"]]),
                     y = (cc[, 2] - bb[["ymin"]]) / (bb[["ymax"]] - bb[["ymin"]]),
                     g = interaction(as.data.frame(cc[, grepl("^L", colnames(cc)), drop = FALSE]), drop = TRUE))
K <- grid_centers(NS)
D <- ESC * dist_mat(K)
G <- grid_centers(NG)
Gll <- st_as_sf(data.frame(lon = bb[["xmin"]] + G[, 1] * (bb[["xmax"]] - bb[["xmin"]]),
                           lat = bb[["ymin"]] + G[, 2] * (bb[["ymax"]] - bb[["ymin"]])),
                coords = c("lon", "lat"), crs = st_crs(borde))
G <- G[lengths(st_intersects(Gll, borde)) > 0, ]   

#DATA
load("dataset-covid-Cali-modificado.Rdata")
SPAN <- as.numeric(FECHA_FIN - FECHA_INI)          
tcase <- if (is.null(COL_FECHA)) DatosCali$TiempoC else
  as.numeric(as.Date(DatosCali[[COL_FECHA]]) - FECHA_INI) / SPAN
dcase <- round(tcase * SPAN)                       
bin <- function(u, n) pmin(pmax(ceiling(u * n), 1), n)
ic <- bin(DatosCali$xmod, NS) + NS * (bin(DatosCali$ymod, NS) - 1)

td <- (0:SPAN) / SPAN                              
ini <- td <= T0
N0 <- tabulate(ic[ini[dcase + 1]], NS^2)           
dsel <- which(!ini) - 1
N  <- sapply(dsel, function(d) tabulate(ic[dcase == d], NS^2))
tu <- td[!ini]
Am <- 1 / NS^2; V <- Am / SPAN
print(paste0("Records: ", nrow(DatosCali), " | initial pattern (t <= t0): ", sum(N0), " | t > t0: ", sum(N)))

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

####AJUSTE
llfun <- function(lr, ll) {
  lm <- loglam_weib(lr, ll, tu)
  sum(N * lm - exp(lm) * V) + sum(N0 * ll - exp(ll) * Am)
}
init <- list(mu_r = 0, mu_l = log((mean(N0) + 0.1) / Am),
             th_r = rep(0, NS^2), th_l = rep(0, NS^2),
             tau_r = 10, tau_l = 10, phi_r = 0.15, phi_l = 0.15)
fit <- run_mcmc(llfun, D, init, prior, N_ITER, BURN, THIN, every = 5000)
post <- post_summary(fit)
v <- c("mu_r", "mu_l", "delta_r", "delta_l")
q <- apply(post[, v], 2, quantile, c(.5, .05, .95, .025, .975))
tab <- data.frame(median = q[1, ], low90 = q[2, ], high90 = q[3, ], low95 = q[4, ], high95 = q[5, ])
print(tab, digits = 4)
write.csv(tab, "Graphs/Section6/PosteriorTable.csv")

#VELOCITIES DAYS SELECTED
tvel <- as.numeric(dias - FECHA_INI) / SPAN
vel <- vel_summary(fit, K, G, tvel, ESC, d_weib)
vel$vabs <- pmin(abs(vel$v), quantile(abs(vel$v), .975, na.rm = TRUE))

#PLOS
day_lab <- function(t){
  d <- FECHA_INI + round(t * SPAN)
  paste(month.name[as.integer(format(d, "%m"))], as.integer(format(d, "%d")))
}
outline <- geom_path(data = poly01, aes(x, y, group = g), inherit.aes = FALSE, linewidth = 0.3)
sv <- function(p, f, w = 8, h = 7) ggsave(file.path("Graphs/Section6", f), p, width = w, height = h, dpi = 300)

daily <- data.frame(day = FECHA_INI + 0:SPAN, n = tabulate(dcase + 1, SPAN + 1))
sv(ggplot(daily, aes(day, n)) + geom_col(fill = "grey60") +
     geom_vline(xintercept = FECHA_INI + T0 * SPAN, linetype = 2) +
     geom_vline(xintercept = dias, colour = "red") +
     labs(x = "Date", y = "Daily cases", title = "Covid-19 cases in Cali (dashed: t0; red: selected days)"),
   "SerieTemporalCali.png", 9, 4)
sv(post_plot(post, NULL), "PosteriorParameters.png", 10, 6)
sv(map_plot(vel, "lam", "Estimated intensity", "lambda", tlab = day_lab) & outline, "Intensity.png")
sv(map_plot(vel, "lt", "Temporal derivative", "d lambda/dt", tlab = day_lab) & outline, "DERT.png")
sv(map_plot(vel, "lx", "Partial derivative in x", "d lambda/dx", tlab = day_lab) & outline, "DERX.png")
sv(map_plot(vel, "ly", "Partial derivative in y", "d lambda/dy", tlab = day_lab) & outline, "DERY.png")
sv(map_plot(vel, "vabs", "Posterior mean of minimum velocity (capped at 97.5th percentile)", "|Velocity|",
            tlab = day_lab) & outline, "Velocities.png")
sv(arrow_plot(vel, "Directions of minimum rate of change", tlab = day_lab) & outline, "Arrows.png")

saveRDS(list(N0 = N0, N = N, tu = tu, fit = fit, vel = vel, tab = tab), "res_cali.rds")




#FINAL DAYS PAPER
sm <- as.numeric(stats::filter(daily$n, rep(1 / 7, 7), sides = 2))

d0 <- min(dsel)                                                      
dd <- d0:SPAN
s <- sm[dd + 1]
ip <- which.max(s)
dpk <- dd[ip]

dmid <- dd[which(s[1:ip] >= (s[1] + s[ip]) / 2)[1]]
sp <- s[ip:length(s)]
lo <- min(sp, na.rm = TRUE)

dfin <- dd[ip - 1 + which(sp <= lo + 0.25 * (s[ip] - lo))[1]]
if (is.na(dfin) || dfin == dpk) { warning("No clear decline after the peak"); dfin <- SPAN - 3 }
dias <- FECHA_INI + c(d0, dmid, dpk, dfin)
print(data.frame(phase = c("start", "mid-rise", "peak", "decline"), day = dias))

tvel <- as.numeric(dias - FECHA_INI) / SPAN
vel <- vel_summary(fit, K, G, tvel, ESC, d_weib)
vel$vabs <- pmin(abs(vel$v), quantile(abs(vel$v), .975, na.rm = TRUE))
sv(ggplot(daily, aes(day, n)) + geom_col(fill = "grey60") +
     geom_line(aes(y = sm), colour = "black") +
     geom_vline(xintercept = FECHA_INI + T0 * SPAN, linetype = 2) +
     geom_vline(xintercept = dias, colour = "red") +
     labs(x = "Date", y = "Daily cases", title = "Covid-19 cases in Cali (dashed: t0; red: selected days)"),
   "SerieTemporalCali.png", 9, 4)
sv(map_plot(vel, "lam", "Estimated intensity", "lambda", tlab = day_lab) & outline, "Intensity.png")
sv(map_plot(vel, "lt", "Temporal derivative", "d lambda/dt", tlab = day_lab) & outline, "DERT.png")
sv(map_plot(vel, "lx", "Partial derivative in x", "d lambda/dx", tlab = day_lab) & outline, "DERX.png")
sv(map_plot(vel, "ly", "Partial derivative in y", "d lambda/dy", tlab = day_lab) & outline, "DERY.png")
sv(map_plot(vel, "vabs", "Posterior mean of minimum velocity (capped at 97.5th percentile)", "|Velocity|",
            tlab = day_lab) & outline, "Velocities.png")
sv(arrow_plot(vel, "Directions of minimum rate of change", tlab = day_lab) & outline, "Arrows.png")
saveRDS(list(N0 = N0, N = N, tu = tu, fit = fit, post= post, vel = vel, tab = tab, dias = dias), "res_cali.rds")
