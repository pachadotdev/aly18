local({
  solve <- getFromNamespace("nls_geppml_", "aly2018")
  t_bln <- matrix(c(1, 0.2, 0.3, 0.4, 1, 0.25, 0.15, 0.35, 1), 3, 3,
                  byrow = TRUE)
  t_ctrf <- t_bln
  t_ctrf[1, 2] <- 0.6
  q <- c(2, 3, 4)
  phi <- c(0.9, 1.1, 1)
  inputs <- list(t_bln, t_ctrf, q, phi)
  before <- serialize(inputs, NULL)

  solution <- solve(t_bln, t_ctrf, q, phi, 7, 1e-12, 1400L, FALSE)
  expect_identical(serialize(inputs, NULL), before)
  expect_true(all(vapply(solution, is.double, logical(1))))
  expect_identical(
    solve(t_bln, t_ctrf, q, phi, 7, 1e-12, 1400L, FALSE), solution
  )

  residuals <- function(costs, imr, omr, output) {
    y1 <- output / sum(output)
    y2 <- phi * output / sum(phi * output)
    inward <- 1 - imr * drop(t(costs) %*% (y1 * omr))
    inward[length(q)] <- 1 - tail(imr, 1L)
    outward <- 1 - omr * drop(costs %*% (y2 * imr))
    c(inward, outward)
  }
  expect_true(max(abs(residuals(
    t_bln, solution$imr_bln, solution$omr_bln, q
  ))) < 1e-10)
  expect_true(max(abs(residuals(
    t_ctrf, solution$imr_cdl, solution$omr_cdl, q
  ))) < 1e-10)
  expect_true(max(abs(residuals(
    t_ctrf, solution$imr_full, solution$omr_full, solution$pstar * q
  ))) < 1e-10)
  y1 <- solution$pstar * q / sum(solution$pstar * q)
  expect_true(max(abs(
    1 - y1 * solution$omr_full / (solution$pstar^(1 - 7) * solution$beta)
  )) < 1e-10)
  expect_equal(solution$beta, q / sum(q) * solution$omr_bln)
  expect_true(all(solution$norm < 1e-12))

  unchanged <- solve(t_bln, t_bln, q, phi, 7, 1e-12, 1400L, FALSE)
  expect_equal(unchanged$imr_bln, unchanged$imr_cdl)
  expect_equal(unchanged$omr_bln, unchanged$omr_cdl)
  expect_equal(unchanged$pstar, rep(1, length(q)))

  expect_error(
    solve(t_bln, t_ctrf, q, phi[-1L], 7, 1e-12, 1400L, FALSE),
    "same length"
  )
  expect_error(
    solve(t_bln, t_ctrf[, -1L], q, phi, 7, 1e-12, 1400L, FALSE),
    "square matrices"
  )
})
