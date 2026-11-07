# Data preparation for "General equilibrium analysis with PPML"
# (Anderson, Larch, and Yotov, 2018).
# This creates output ($Y_i$) and expenditure ($E_j$) for each country

if (!require("data.table")) install.packages("data.table", repos = "https://cran.r-project.org")
if (!require("haven")) install.packages("haven", repos = "https://cran.r-project.org")

library(data.table)
library(haven)

aly18 <- as.data.table(read_dta("original-code-and-data/ge_ppml_data_orig.dta"))
aly18 <- aly18[, .(exporter, importer, trade, LN_DIST, CNTG, BRDR)]

# Aggregate variables: output by exporter and expenditure by importer
aly18[, output := sum(trade), by = exporter]
aly18[, expndr := sum(trade), by = importer]

# Germany is the reference country. It is renamed to "ZZZ" to be the last
# category in Stata, so here we use "0-DEU" using the R convention for the
# reference level
aly18[exporter == "DEU", exporter := "0-DEU"]
aly18[importer == "DEU", importer := "0-DEU"]
aly18[, expndr_deu := expndr[importer == "0-DEU"][1]]
aly18[, pair := paste(exporter, importer, sep = "-")]

colnames(aly18) <- tolower(colnames(aly18))

try(dir.create("results"))
fwrite(aly18, "results/ge_ppml_data.csv")
usethis::use_data(aly18, overwrite = TRUE)
