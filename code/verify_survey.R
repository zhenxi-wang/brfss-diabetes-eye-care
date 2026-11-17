# Optional independent R cross-check. Run at project root after code/analyse.py.
# Requires the survey and arrow packages, plus nnet from the recommended R packages.
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = FALSE)
suppressPackageStartupMessages(library(survey))

contrasts <- read.csv("results/adjusted_contrasts.csv")
profiles <- read.csv("results/profiles.csv")

jacobian <- function(f, estimate, h = 1e-5) {
  vapply(seq_along(estimate), function(j) {
    delta <- numeric(length(estimate))
    delta[j] <- h
    (f(estimate + delta) - f(estimate - delta)) / (2 * h)
  }, numeric(length(f(estimate))))
}

for (year in c(2025, 2022)) {
  d <- as.data.frame(arrow::read_parquet(sprintf("data/derived/%d_analysis.parquet", year)))
  d$stratum <- interaction(d[["_STATE"]], d[["_STSTR"]], drop = TRUE)
  d$unit <- d[["_PSU"]]
  d$weight <- d[["_LLCPWT"]]
  d$both <- as.numeric(d$dil & d$photo)
  d$dilation_only <- as.numeric(d$dil & !d$photo)
  d$photography_only <- as.numeric(!d$dil & d$photo)
  d$neither <- as.numeric(!d$dil & !d$photo)
  design <- svydesign(ids = ~unit, strata = ~stratum, weights = ~weight, nest = TRUE, data = d)
  paired <- subset(design, pair)
  syntax <- ~both + dilation_only + photography_only + neither
  raw <- svymean(syntax, paired)
  summary <- profiles[profiles$year == year & profiles$group == "All", ]
  known <- summary[match(c("Both", "Dilation only", "Photography only", "Neither"), summary$measure), ]
  stopifnot(max(abs(coef(raw)*100 - known$estimate)) < 1e-7,
            max(abs(SE(raw)*100 - known$se)) < 1e-5)
  raw_limits <- plogis(qlogis(coef(raw)) + outer(SE(raw)/(coef(raw)*(1-coef(raw))),
                         qt(c(.025, .975), degf(design))))
  stopifnot(max(abs(raw_limits*100 - as.matrix(known[c("lo", "hi")]))) < 1e-5)

  rhs <- ~ nonmetro + cost + age_group + sex + race + education + insulin + state_name
  keep <- d$pair & complete.cases(d[all.vars(rhs)])
  s <- d[keep, ]
  s$pattern <- factor(max.col(cbind(s$neither, s$both, s$dilation_only, s$photography_only)),
                      levels = 1:4)
  w <- s$weight / mean(s$weight)
  fit <- nnet::multinom(update(rhs, pattern ~ .), data = s, weights = w, Hess = TRUE,
                        trace = FALSE, maxit = 2000, reltol = 1e-14)
  stopifnot(fit$convergence == 0)
  x <- model.matrix(rhs, s)
  beta <- as.vector(t(coef(fit)))
  predict_prob <- function(design_matrix, b) {
    logits <- cbind(0, design_matrix %*% matrix(b, nrow = ncol(design_matrix), ncol = 3))
    mass <- exp(logits - apply(logits, 1, max))
    mass / rowSums(mass)
  }
  indicators <- model.matrix(~ pattern - 1, s)[, -1]
  residual <- indicators - predict_prob(x, beta)[, -1]
  scores <- do.call(cbind, lapply(seq_len(ncol(residual)), function(j) x * (w*residual[,j])))
  categories <- rbind(c(0, 0, 0), c(1, 1, 0), c(1, 1, 0), c(0, 1, 1))
  for (exposure in c("nonmetro", "cost")) {
    off <- on <- x
    off[, exposure] <- 0
    on[, exposure] <- 1
    individual <- function(b) (predict_prob(on, b) - predict_prob(off, b)) %*% categories
    marginal <- function(b) colSums(individual(b)*w)/sum(w)
    estimate <- marginal(beta)
    influence <- scores %*% solve(fit$Hessian, jacobian(marginal, beta))
    influence <- influence + sweep(individual(beta), 2, estimate) * w/sum(w)
    for (j in 1:3) {
      full <- numeric(nrow(d))
      full[keep] <- influence[, j]/s$weight
      design$variables[[paste0("contrast", j)]] <- full
    }
    errors <- SE(svytotal(~contrast1 + contrast2 + contrast3, design))
    expected <- contrasts[contrasts$year == year & contrasts$model == "M1" &
                          contrasts$exposure == exposure, ]
    expected <- expected[match(c("dil", "either", "Change in disparity"), expected$outcome), ]
    stopifnot(max(abs(estimate*100 - expected$difference_pp)) < .001,
              max(abs(errors*100 - expected$se)) < .001)
  }
  cat(year, "survey proportions, confidence limits, joint contrasts and gap-change SEs passed\n")
}

