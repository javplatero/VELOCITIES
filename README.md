# Velocities driven by spatial growth curves

R code accompanying the paper

> Platero, J., Mateu, J. and Gelfand, A. E. *Velocities driven by spatial growth curves for geostatistical and spatial point pattern models.* Submitted to *Biometrics*.

The code fits the two models of the paper with MCMC: a geostatistical model with binary responses and a spatio-temporal point pattern model. It then computes posterior velocities and directions of minimum velocity, and produces every table and figure in the paper and its Supplementary Material (Web Appendices A–D).

## Repository contents

| File | Description |
|---|---|
| `funciones_comunes.R` | Shared functions: Matérn ν = 3/2 covariance, Gaussian process simulation, MCMC sampler, predictive process gradients and velocities, plotting functions, and the selection of representative simulated realizations. Sourced by the three scripts below. |
| `GEOS_SIM_DEF.R` | Simulation study for the geostatistical model with the linear growth curve (Section 4.1, Web Appendix B). |
| `PointPattern_SIMULATION_DEF.R` | Simulation study for the point pattern model with the Weibull growth curve (Section 4.2, Web Appendix C). |
| `CALI_DEF.R` | Covid-19 case study in Cali, Colombia (Section 5, Web Appendix D). Needs data that are not distributed here (see [Data](#data)). |

## Requirements

- R (≥ 4.1)
- CRAN packages: `ggplot2` and `patchwork` for all scripts, `sf` for `CALI_DEF.R`, and optionally `coda` for MCMC diagnostics.

```r
install.packages(c("ggplot2", "patchwork", "sf", "coda"))
```

## How to run

Set the working directory to the repository folder and run each script from there. Every script sources `funciones_comunes.R`, fixes its own random seed, and writes its figures to `Graphs/Section4`, `Graphs/Section5` or `Graphs/Section6`. Those folders are created if they do not exist.

```r
source("GEOS_SIM_DEF.R")                 # Section 4.1 + Web Appendix B
source("PointPattern_SIMULATION_DEF.R")  # Section 4.2 + Web Appendix C
source("CALI_DEF.R")                     # Section 5   + Web Appendix D (needs the data)
```

Each script runs a long MCMC: 200,000 iterations for the simulations and 250,000 for Cali, on a 20 × 20 grid. Expect several hours per script. Each script also saves its results (`res_geoestadistico.rds`, `res_patron_puntual.rds`, `res_cali.rds`), so tables and figures can be rebuilt later without running the MCMC again.

### Settings used in the paper

| | Geostatistical | Point pattern | Cali |
|---|---|---|---|
| Seed | 2025 | 123456 | 35649 |
| Grid (space × time) | 20 × 20 × 200 | 20 × 20 × 200 | 20 × 20 × daily |
| Growth curve | Linear (γ = 1) | Weibull (t₀ = 1/4) | Weibull (t₀ = 1/3) |
| Iterations / burn-in / thinning | 200,000 / 100,000 / 100 | 200,000 / 100,000 / 100 | 250,000 / 125,000 / 100 |
| Prior on φ | U(0.12, 0.2) | U(0.12, 0.2) | U(0.02, 0.3) |

In all cases μ ~ N(0, 10²) and τ ~ Gamma(10, 1). The Matérn covariance uses φ as the range parameter, `options(matern_param = "range")`, and results are reported for the identifiable parameter δ = τφ³.

## Output: correspondence with the paper

### Main text

| Paper | Script | Output |
|---|---|---|
| Table 1 (geostatistical rows) | `GEOS_SIM_DEF.R` | printed by `post_table(post, true)` |
| Table 1 (point pattern rows) | `PointPattern_SIMULATION_DEF.R` | printed by `post_table(post, true)` |
| Table 2 | `CALI_DEF.R` | `Graphs/Section6/PosteriorTable.csv` |
| Figure 1 | `GEOS_SIM_DEF.R` | `Graphs/Section4/EstimacionProbabilidades.png` |
| Figure 2 | `GEOS_SIM_DEF.R` | `Graphs/Section4/IntensityAndDirections.png` |
| Figure 3 | `PointPattern_SIMULATION_DEF.R` | `Graphs/Section5/SeriesTemporalesFil14.png` |
| Figure 4 | `PointPattern_SIMULATION_DEF.R` | `Graphs/Section5/IntensitiesDirectionsW.png` |
| Figure 5 | `CALI_DEF.R` | `Graphs/Section6/SerieTemporalCali.png` |
| Figure 6 | `CALI_DEF.R` | `Graphs/Section6/Arrows.png` |

### Supplementary Material

| Web Appendix | Script | Files |
|---|---|---|
| B (geostatistical simulation) | `GEOS_SIM_DEF.R` | `Graphs/Section4/`: `SimulacionGEOS`, `CandidatesThetaR`, `CandidatesThetaL`, `CandidatesData`, `SimulacionCampos`, `MarginalPosteriors`, `ComparacionCampos`, `SemivariogramFields`, `ComparisonRealEstimate`, `DERT`, `DERX`, `DERY`, `Velocities` |
| C (point pattern simulation) | `PointPattern_SIMULATION_DEF.R` | `Graphs/Section5/`: `SimulacionWeibullPP`, `CandidatesThetaR`, `CandidatesThetaL`, `CandidatesData`, `SerieTemporalSimulada`, `SimulacionCamposW`, `PosteriorsWeibull`, `ComparacionCamposW`, `SemivariogramFieldsW`, `ComparacionIntensidades`, `DERT`, `DERX`, `DERY`, `Velocities` |
| D (Cali) | `CALI_DEF.R` | `Graphs/Section6/`: `Intensity`, `PosteriorParameters`, `DERT`, `DERX`, `DERY`, `Velocities` |

All files are `.png`. Web Appendix A (predictive process gradients) is analytical and is implemented in `pp_pred()` and `vel_summary()` in `funciones_comunes.R`.

## Data

The simulated data are fully reproducible with the scripts and seeds above.

The Covid-19 case data for Cali were provided by the Municipal Public Health Secretary of Cali (Secretaría de Salud Pública Municipal de Cali) under a data use agreement. They contain the residential locations of individual patients and therefore **are not publicly available and are not included in this repository**. Requests for access should be addressed to the Municipal Public Health Secretary of Cali.

`CALI_DEF.R` expects the following files in the working directory:

- `dataset-covid-Cali-modificado.Rdata`: a data frame `DatosCali` with one row per confirmed case. It needs the columns `xmod` and `ymod`, the residence coordinates rescaled to [0, 1], and `TiempoC`, the diagnosis date rescaled to [0, 1] so that 0 is 2 March 2020 and 1 is 25 July 2020.
- `borde_cali.shp` and its companion files: the urban boundary of Cali, used only to clip the maps.

## Citation

If you use this code, please cite the paper above.

## Contact

Javier Platero, Department of Mathematics, Universitat Jaume I, Castellón, Spain — platero@uji.es
