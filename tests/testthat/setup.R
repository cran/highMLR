# Limit multithreading during R CMD check so total CPU time stays within
# the CRAN two-core policy. These settings affect the test session only and
# do not change behaviour for users running the package normally.

Sys.setenv(OMP_NUM_THREADS = "2")        # OpenMP (xgboost, aorsf, glmnet)
Sys.setenv(OMP_THREAD_LIMIT = "2")
Sys.setenv(RANGER_NUM_THREADS = "2")     # ranger
Sys.setenv(XGBOOST_NTHREAD = "2")        # xgboost

options(Ncpus = 2L)

# data.table-style guard in case any dependency consults it
Sys.setenv(R_DATATABLE_NUM_THREADS = "2")
