library(data.table)
library(capybara)
library(aly18)

sigma <- 7
ref_country <- "0-DEU"

# tight tolerances, otherwise the fixed effects are only accurate to ~1e-4 and the error is amplified in the
# multilateral resistances
ctrl <- list(return_fe = TRUE, dev_tol = 1e-12, center_tol = 1e-12)

# step 1: baseline ----

aly18 <- copy(aly18)

# the border dummy is split into the border between the US and Canada and the border for the rest of the world
aly18[, BRDR_USACAN := as.integer(
  (exporter == "USA" & importer == "CAN") | (exporter == "CAN" & importer == "USA")
)]
aly18[BRDR_USACAN == 1, brdr := 0L]

fit_bln <- fepoisson(
  trade ~ ln_dist + cntg + brdr + BRDR_USACAN | exporter + importer,
  data = aly18,
  control = ctrl
)

summary(fit_bln)

# estibrated trade costs, equation (14): the gravity error (observed / predicted trade) is treated as a
# component of the trade costs, so the baseline trade costs match the observed trade flows perfectly
b <- coef(fit_bln)

aly18[, error_ij := trade / predict(fit_bln)]

# a zero flow has zero estibrated trade cost and is zero in all the scenarios, so it is dropped
aly18_full <- copy(aly18)
aly18 <- aly18[trade > 0]

aly18[, t_ij_bln := exp(
  b["ln_dist"] * ln_dist + b["cntg"] * cntg + b["brdr"] * brdr + b["BRDR_USACAN"] * BRDR_USACAN
) * error_ij]

# the baseline trade flows are the observed trade flows
aly18[, trade_bln := trade]

# baseline indexes, equations (7) and (8)
add_fes(aly18, fit_bln, "exp_fes_bln", "imp_fes_bln")

aly18[, omr_bln := output * expndr_deu / exp_fes_bln]
aly18[, imr_bln := expndr / (imp_fes_bln * expndr_deu)]

aly18[, rGDP_bln := sum((output / imr_bln^(1 / (1 - sigma)))[exporter == importer]), by = exporter]
aly18[, tot_exp_bln := sum(trade_bln[exporter != importer]), by = exporter]

# expenditure is a constant share of output (phi)
aly18[, phi := (expndr / output)[exporter == importer], by = importer]

# steps 2 and 3: counterfactual ----

# border_removed is TRUE for the observations in which the US-Canada border is abolished
full_ge <- function(aly18, border_removed) {
  aly18 <- copy(aly18)

  # step 2: conditional, the trade costs are fixed with an offset and only the fixed effects are estimated
  aly18[, t_ij_ctrf := t_ij_bln * exp(-b["BRDR_USACAN"] * BRDR_USACAN * border_removed)]
  aly18[, t_ij_ctrf_1 := log(t_ij_ctrf)]

  fit_cdl <- fepoisson(
    trade ~ 0 | exporter + importer,
    data = aly18,
    offset = ~t_ij_ctrf_1,
    control = ctrl
  )

  aly18[, trade_cdl := predict(fit_cdl)]
  add_fes(aly18, fit_cdl, "exp_fes_cdl", "imp_fes_cdl")

  aly18[, omr_cdl := output * expndr_deu / exp_fes_cdl]
  aly18[, imr_cdl := expndr / (imp_fes_cdl * expndr_deu)]

  # step 3: full endowment, output and expenditure change with the factory-gate prices, we iterate until
  # the prices converge
  # the columns ending in _old store the values from the previous iteration, initially those of the
  # conditional step
  aly18[, trade_pred := trade_cdl]
  aly18[, c("fes_old", "expndr_deu_old", "omr_old", "imr_old") := .(exp_fes_cdl, expndr_deu, omr_cdl, imr_cdl)]
  aly18[, fes_dom_old := exp_fes_cdl[exporter == importer], by = importer]

  # factory-gate prices of each exporter (p_exp) and of each importer as a producer (p_imp), and changes in
  # the multilateral resistances
  aly18[, p_exp := (exp_fes_cdl / exp_fes_bln)^(1 / (1 - sigma))]
  aly18[, p_imp := (fes_dom_old / exp_fes_bln[exporter == importer])^(1 / (1 - sigma)), by = importer]
  aly18[, c("omr_ch", "imr_ch") := 1]

  # in this scenario the iteration oscillates around the solution (the standard deviation of the change in
  # p_exp does not fall below 0.001), so it stops when the maximum of the change in p_exp is below 0.001
  p_exp_old <- 0
  sd_diff <- 1
  max_diff <- 1
  iter <- 1

  while (max_diff > 1e-3) {
    iter <- iter + 1

    # equation (14), rescale the trade flows and estimate equation (11)
    aly18[, trade_iter := trade_pred * p_exp * p_imp / (omr_ch * imr_ch)]

    fit_iter <- fepoisson(
      trade_iter ~ 0 | exporter + importer,
      data = aly18,
      offset = ~t_ij_ctrf_1,
      control = ctrl
    )

    aly18[, trade_pred := predict(fit_iter)]
    add_fes(aly18, fit_iter, "exp_fes", "imp_fes")
    aly18[, fes_dom := exp_fes[exporter == importer], by = importer]

    # update output and expenditure
    aly18[, output_iter := sum(trade_pred), by = exporter]
    aly18[, expndr_iter := (phi * output_iter)[exporter == importer], by = importer]
    aly18[, expndr_deu_iter := sum(trade_pred[importer == ref_country])]

    # convergence criteria with the change in p_exp from the previous iteration
    diffs <- aly18$p_exp - p_exp_old
    sd_diff <- sd(diffs)
    max_diff <- abs(max(diffs))
    p_exp_old <- aly18$p_exp

    # update prices
    aly18[, p_exp := ((exp_fes / fes_old) / (expndr_deu_iter / expndr_deu_old))^(1 / (1 - sigma))]
    aly18[, p_imp := ((fes_dom / fes_dom_old) / (expndr_deu_iter / expndr_deu_old))^(1 / (1 - sigma))]

    # update the multilateral resistances, equations (7) and (8)
    aly18[, omr_iter := output_iter * expndr_deu_iter / exp_fes]
    aly18[, imr_iter := expndr_iter / (imp_fes * expndr_deu_iter)]

    aly18[, c("omr_ch", "imr_ch") := .(omr_iter / omr_old, imr_iter / imr_old)]

    aly18[, c("fes_old", "fes_dom_old", "expndr_deu_old", "omr_old", "imr_old") := .(
      exp_fes, fes_dom, expndr_deu_iter, omr_iter, imr_iter
    )]

    cat("Iteration", iter, "- sd:", sd_diff, "- max:", max_diff, "\n")
  }

  # full endowment indexes, as in the original Stata code the expenditure of the reference country
  # is phi * output from the last iteration, and not the sum of the predicted flows used in the loop
  aly18[, expndr_deu_full := expndr_iter[importer == ref_country][1]]

  aly18[, p_full := ((exp_fes / exp_fes_bln) / (expndr_deu_full / expndr_deu))^(1 / (1 - sigma))]
  aly18[, output_full := p_full * output]
  aly18[, expndr_full := (phi * output_full)[exporter == importer], by = importer]

  aly18[, omr_full := output_full * expndr_deu_full / exp_fes]
  aly18[, imr_full := expndr_iter / (imp_fes * expndr_deu_full)]

  aly18[, rGDP_full := sum((output_full / imr_full^(1 / (1 - sigma)))[exporter == importer]), by = exporter]

  aly18[, trade_full := output_full * expndr_full * t_ij_ctrf / (imr_full * omr_full)]
  aly18[, tot_exp_full := sum(trade_full[exporter != importer]), by = exporter]

  # the exporter indexes are constant for each exporter and the importer indexes for each importer, so we
  # keep one observation of each
  indexes_exp <- aly18[, .(
    tot_exp_bln = tot_exp_bln[1], tot_exp_full = tot_exp_full[1],
    rGDP_bln = rGDP_bln[1], rGDP_full = rGDP_full[1],
    omr_bln = omr_bln[1], omr_full = omr_full[1],
    p_full = p_full[1]
  ), by = .(country = exporter)]

  indexes_imp <- aly18[, .(imr_bln = imr_bln[1], imr_full = imr_full[1]), by = .(country = importer)]

  indexes <- merge(indexes_exp, indexes_imp, by = "country")
  indexes[, country := fifelse(country == ref_country, "DEU", country)]
  setorder(indexes, country)

  # changes in percentage as reported in the paper
  indexes[, .(
    country,
    exports = (tot_exp_full / tot_exp_bln - 1) * 100,
    rGDP = (rGDP_full / rGDP_bln - 1) * 100,
    imr = price_change(imr_full, imr_bln),
    omr = price_change(omr_full, omr_bln),
    prices = (p_full - 1) * 100
  )]
}

# the border between the US and Canada is abolished in both directions
indexes_both <- full_ge(aly18, border_removed = TRUE)

# only the border for the exports of Canada to the US is abolished
indexes_asym <- full_ge(aly18, border_removed = aly18$exporter == "CAN" & aly18$importer == "USA")

# step 4: indexes ----

setnames(indexes_both, names(indexes_both)[-1], paste0(names(indexes_both)[-1], "_both"))
setnames(indexes_asym, names(indexes_asym)[-1], paste0(names(indexes_asym)[-1], "_asym"))

indexes_final <- merge(indexes_both, indexes_asym, by = "country")

cols <- setdiff(names(indexes_final), "country")
indexes_final[, (cols) := lapply(.SD, round, 2), .SDcols = cols]

indexes_final

fwrite(indexes_final, "results/table3-effect-of-removing-us-canada-border-based-on-estibrated-trade-costs.csv")

# step 5: validation with the non-linear solver ----

# solves the systems of the original MATLAB code (ge_ppml_indexes.m) directly, without iterating on regressions
nls_formula <- trade ~ ln_dist + cntg + brdr + BRDR_USACAN | exporter + importer

nls_scenario <- function(removed, suffix) {
  r <- nls_geppml(
    aly18_full, formula = nls_formula, border = "BRDR_USACAN", removed = removed,
    sigma = sigma, ref_country = ref_country, estibrated = TRUE
  )
  r <- r[, .(country, exports = exports_full, rGDP = rGDP_full, imr = imr_full, omr = omr_full,
             prices = prices_full)]
  setnames(r, names(r)[-1], paste0(names(r)[-1], suffix))
  r
}

indexes_nls <- merge(
  nls_scenario(TRUE, "_both"),
  nls_scenario(aly18_full$exporter == "CAN" & aly18_full$importer == "USA", "_asym"),
  by = "country"
)

indexes_nls <- indexes_nls[, lapply(.SD, round, 2), by = country]

cols <- setdiff(names(indexes_final), "country")

m <- merge(indexes_final, indexes_nls, by = "country", suffixes = c("_it", "_nls"))

diff_dt <- m[, c(list(country = country),
                 setNames(lapply(cols, function(v) abs(get(paste0(v, "_it")) - get(paste0(v, "_nls")))), cols))]

# max gap per index
diff_dt[, lapply(.SD, max), .SDcols = cols]
