source("00-helpers.r")

# step 1: baseline ----

fit_bln <- fepoisson(
  trade ~ LN_DIST + CNTG + BRDR | exporter + importer,
  data = dt,
  control = ctrl
)

summary(fit_bln)

# estibrated trade costs, equation (14): the gravity error (observed / predicted trade) is treated as a
# component of the trade costs, so the baseline trade costs match the observed trade flows perfectly
b <- coef(fit_bln)

dt[, error_ij := trade / predict(fit_bln)]

# a zero flow has zero estibrated trade cost and is zero in all the scenarios, so it is dropped
dt <- dt[trade > 0]

# in the counterfactual there are no borders and the internal trade costs remain the same
dt[, t_ij_bln := exp(b["LN_DIST"] * LN_DIST + b["CNTG"] * CNTG + b["BRDR"] * BRDR) * error_ij]
dt[, t_ij_ctrf := fifelse(
  exporter == importer, t_ij_bln, exp(b["LN_DIST"] * LN_DIST + b["CNTG"] * CNTG) * error_ij
)]
dt[, t_ij_ctrf_1 := log(t_ij_ctrf)]

# the baseline trade flows are the observed trade flows
dt[, trade_bln := trade]

# baseline indexes, equations (7) and (8)
add_fes(dt, fit_bln, "exp_fes_bln", "imp_fes_bln")

dt[, omr_bln := output * expndr_deu / exp_fes_bln]
dt[, imr_bln := expndr / (imp_fes_bln * expndr_deu)]

dt[, rGDP_bln := sum((output / imr_bln^(1 / (1 - sigma)))[exporter == importer]), by = exporter]
dt[, tot_exp_bln := sum(trade_bln[exporter != importer]), by = exporter]

# step 2: conditional ----

# equation (11), the trade costs are fixed with an offset and only the fixed effects are estimated
fit_cdl <- fepoisson(
  trade ~ 0 | exporter + importer,
  data = dt,
  offset = ~t_ij_ctrf_1,
  control = ctrl
)

dt[, trade_cdl := predict(fit_cdl)]
add_fes(dt, fit_cdl, "exp_fes_cdl", "imp_fes_cdl")

dt[, omr_cdl := output * expndr_deu / exp_fes_cdl]
dt[, imr_cdl := expndr / (imp_fes_cdl * expndr_deu)]

dt[, rGDP_cdl := sum((output / imr_cdl^(1 / (1 - sigma)))[exporter == importer]), by = exporter]
dt[, tot_exp_cdl := sum(trade_cdl[exporter != importer]), by = exporter]

# step 3: full endowment ----

# output and expenditure change with the factory-gate prices, we iterate until the prices converge
# expenditure is a constant share of output (phi)
# the columns ending in _old store the values from the previous iteration, initially those of the conditional step

dt[, phi := (expndr / output)[exporter == importer], by = importer]

dt[, trade_pred := trade_cdl]
dt[, c("fes_old", "expndr_deu_old", "omr_old", "imr_old") := .(exp_fes_cdl, expndr_deu, omr_cdl, imr_cdl)]
dt[, fes_dom_old := exp_fes_cdl[exporter == importer], by = importer]

# factory-gate prices of each exporter (p_exp) and of each importer as a producer (p_imp), and changes in the
# multilateral resistances
dt[, p_exp := (exp_fes_cdl / exp_fes_bln)^(1 / (1 - sigma))]
dt[, p_imp := (fes_dom_old / exp_fes_bln[exporter == importer])^(1 / (1 - sigma)), by = importer]
dt[, c("omr_ch", "imr_ch") := 1]

# converged when the standard deviation and the maximum of the change in p_exp are below 0.001
p_exp_old <- 0
sd_diff <- 1
max_diff <- 1
iter <- 1

while (sd_diff > 1e-3 | max_diff > 1e-3) {
  iter <- iter + 1

  # equation (14), rescale the trade flows and estimate equation (11)
  dt[, trade_iter := trade_pred * p_exp * p_imp / (omr_ch * imr_ch)]

  fit_iter <- fepoisson(
    trade_iter ~ 0 | exporter + importer,
    data = dt,
    offset = ~t_ij_ctrf_1,
    control = ctrl
  )

  dt[, trade_pred := predict(fit_iter)]
  add_fes(dt, fit_iter, "exp_fes", "imp_fes")
  dt[, fes_dom := exp_fes[exporter == importer], by = importer]

  # update output and expenditure
  dt[, output_iter := sum(trade_pred), by = exporter]
  dt[, expndr_iter := (phi * output_iter)[exporter == importer], by = importer]
  dt[, expndr_deu_iter := sum(trade_pred[importer == ref_country])]

  # convergence criteria with the change in p_exp from the previous iteration
  diffs <- dt$p_exp - p_exp_old
  sd_diff <- sd(diffs)
  max_diff <- abs(max(diffs))
  p_exp_old <- dt$p_exp

  # update prices
  dt[, p_exp := ((exp_fes / fes_old) / (expndr_deu_iter / expndr_deu_old))^(1 / (1 - sigma))]
  dt[, p_imp := ((fes_dom / fes_dom_old) / (expndr_deu_iter / expndr_deu_old))^(1 / (1 - sigma))]

  # update the multilateral resistances, equations (7) and (8)
  dt[, omr_iter := output_iter * expndr_deu_iter / exp_fes]
  dt[, imr_iter := expndr_iter / (imp_fes * expndr_deu_iter)]

  dt[, c("omr_ch", "imr_ch") := .(omr_iter / omr_old, imr_iter / imr_old)]

  dt[, c("fes_old", "fes_dom_old", "expndr_deu_old", "omr_old", "imr_old") := .(
    exp_fes, fes_dom, expndr_deu_iter, omr_iter, imr_iter
  )]

  cat("Iteration", iter, "- sd:", sd_diff, "- max:", max_diff, "\n")
}

# full endowment indexes, as in the original Stata code the expenditure of the reference country
# is phi * output from the last iteration, and not the sum of the predicted flows used in the loop
dt[, expndr_deu_full := expndr_iter[importer == ref_country][1]]

dt[, p_full := ((exp_fes / exp_fes_bln) / (expndr_deu_full / expndr_deu))^(1 / (1 - sigma))]
dt[, output_full := p_full * output]
dt[, expndr_full := (phi * output_full)[exporter == importer], by = importer]

dt[, omr_full := output_full * expndr_deu_full / exp_fes]
dt[, imr_full := expndr_iter / (imp_fes * expndr_deu_full)]

dt[, rGDP_full := sum((output_full / imr_full^(1 / (1 - sigma)))[exporter == importer]), by = exporter]

dt[, trade_full := output_full * expndr_full * t_ij_ctrf / (imr_full * omr_full)]
dt[, tot_exp_full := sum(trade_full[exporter != importer]), by = exporter]

# step 4: indexes ----

# the exporter indexes are constant for each exporter and the importer indexes for each importer, so we keep
# one observation of each
indexes_exp <- dt[, .(
  tot_exp_bln = tot_exp_bln[1], tot_exp_cdl = tot_exp_cdl[1], tot_exp_full = tot_exp_full[1],
  rGDP_bln = rGDP_bln[1], rGDP_cdl = rGDP_cdl[1], rGDP_full = rGDP_full[1],
  omr_bln = omr_bln[1], omr_cdl = omr_cdl[1], omr_full = omr_full[1],
  p_full = p_full[1]
), by = .(country = exporter)]

indexes_imp <- dt[, .(
  imr_bln = imr_bln[1], imr_cdl = imr_cdl[1], imr_full = imr_full[1]
), by = .(country = importer)]

indexes <- merge(indexes_exp, indexes_imp, by = "country")
indexes[, country := fifelse(country == ref_country, "DEU", country)]
setorder(indexes, country)

# changes in percentage as reported in the paper
indexes_final <- indexes[, .(
  country,
  exports_cdl = (tot_exp_cdl / tot_exp_bln - 1) * 100,
  exports_full = (tot_exp_full / tot_exp_bln - 1) * 100,
  rGDP_cdl = (rGDP_cdl / rGDP_bln - 1) * 100,
  rGDP_full = (rGDP_full / rGDP_bln - 1) * 100,
  imr_cdl = price_change(imr_cdl, imr_bln),
  imr_full = price_change(imr_full, imr_bln),
  omr_cdl = price_change(omr_cdl, omr_bln),
  omr_full = price_change(omr_full, omr_bln),
  prices_full = (p_full - 1) * 100
)]

cols <- setdiff(names(indexes_final), "country")
indexes_final[, (cols) := lapply(.SD, round, 2), .SDcols = cols]

indexes_final

fwrite(indexes_final, "results/table2-effect-of-removing-trade-borders-based-on-estibrated-trade-costs.csv")
