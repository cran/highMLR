## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
# library(highMLR)
# data(hnscc)
# 
# fit <- highmlr(
#   hnscc,
#   time   = "OS",
#   status = "Death",
#   method = "coxnet",
#   resampling = "cv",
#   folds = 5
# )
# 
# print(fit)
# plot(fit, top_n = 20)

## ----eval = FALSE-------------------------------------------------------------
# cmp <- highmlr_compare(
#   hnscc, "OS", "Death",
#   methods = c("coxnet", "rsf", "univariate")
# )
# cmp$summary

## ----eval = FALSE-------------------------------------------------------------
# data(srdata)
# keep <- highmlr_screen(srdata, "OS", "event",
#                        filter = "variance", keep = 500)
# fit  <- highmlr(srdata, "OS", "event",
#                 features = keep, method = "coxnet")

## ----eval = FALSE-------------------------------------------------------------
# ex <- highmlr_explain(fit, new_data = hnscc, method = "survshap")
# print(ex)
# plot(ex)

## -----------------------------------------------------------------------------
sessionInfo()

