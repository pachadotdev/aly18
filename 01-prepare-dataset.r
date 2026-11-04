# Data preparation for "General equilibrium analysis with PPML"
# (Anderson, Larch, and Yotov, 2018).
# This creates output ($Y_i$) and expenditure ($E_j$) for each country

if (!require("data.table")) install.packages("data.table", repos = "https://cran.r-project.org")
if (!require("haven")) install.packages("haven", repos = "https://cran.r-project.org")

library(data.table)
library(haven)

dt <- as.data.table(read_dta("original-data/ge_ppml_data_orig.dta"))
dt <- dt[, .(exporter, importer, trade, LN_DIST, CNTG, BRDR)]

# Aggregate variables: output by exporter and expenditure by importer
dt[, output := sum(trade), by = exporter]
dt[, expndr := sum(trade), by = importer]

# Germany is the reference country. It is renamed to "ZZZ" to be the last
# category in Stata, so here we use "0-DEU" using the R convention for the
# reference level
dt[exporter == "DEU", exporter := "0-DEU"]
dt[importer == "DEU", importer := "0-DEU"]
dt[, expndr_deu := expndr[importer == "0-DEU"][1]]
dt[, pair := paste(exporter, importer, sep = "-")]

try(dir.create("results"))
fwrite(dt, "results/ge_ppml_data.csv")
