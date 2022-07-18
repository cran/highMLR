#' mlhighCox
#'
#' Performs feature Selection using Cox PH on high-dimensional data
#'
#' @description This function extracts desired number of features based on minimum log-Loss function
#' using Cox proportional hazard model as learner method on a high dimensional survival data.
#'
#' @details Using the Cox proportional hazard model on the given survival data, this function selects
#' the most significant feature based on a performance measure. The performance measure is considered as logarithmic loss function. It is defined as,
#' \deqn{L(f,t)=-log(f(t))}. The features with minimum log-loss function are extracted.
#'
#' @param cols A numeric vector of column numbers indicating the features for which the log Loss functions are to be computed
#' @param idSurv The name of the survival time variable
#' @param idEvent The name of the survival event variable
#' @param per Percentage of total features to be selected, default value 20
#' @param fold An integer denoting number of folds in cross validation, default value 3
#' @param data A data frame that contains the survival and covariate information for the subjects
#'
#' @import mlr3
#' @import mlr3learners
#' @import survival
#' @import utils
#' @import gtools
#' @import dplyr
#' @importFrom stats coef as.formula quantile BIC complete.cases
#' @return A dataframe containing desired number of features and the corresponding log Loss function.
#' @references Sonabend, R., Király, F. J., Bender, A., Bernd Bischl B. and Lang M. mlr3proba: An R Package for Machine Learning in Survival Analysis, 2021, Bioinformatics, <https://doi.org/10.1093/bioinformatics/btab039>
#'
#' @examples
#' \dontrun{
#' data(hnscc)
#' mlhighCox(cols=c(6:15), idSurv="OS", idEvent="Death", per=20, fold = 3, data=hnscc)
#' }
#' @export
#' @author Atanu Bhattacharjee, Gajendra K. Vishwakarma & Souvik Banerjee
#' @seealso mlhighKap, mlhighFrail

mlhighCox=function(cols, idSurv, idEvent, per=20, fold=3, data)
{
  format_range = function(range) {
    l = min(range)
    u = max(range)

    str = sprintf(
      "%s%s, %s%s",
      if (is.finite(l)) "[" else "(",
      if (is.finite(l)) c(l, l) else c("-\\infty", "-Inf"),
      if (is.finite(u)) c(u, u) else c("\\infty", "Inf"),
      if (is.finite(u)) "]" else ")")
    paste0("\\eqn{", str[1L], "}{", str[2L], "}")
  }

  format_types = function(types) {
    if (length(types) == 0) {
      return("-")
    } else {
      return(paste0(types, collapse = ", "))
    }
  }

  toproper = function(str, split = " ", fixed = TRUE) {
    str = strsplit(str, split, fixed)
    str = lapply(str, function(x) {
      paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, 1000)), collapse = split)
    })
    return(unlist(str))
  }

  check_subsetpattern = function(x, choices, empty.ok = TRUE) { # nolint
    if (all(grepl(paste0(choices, collapse = "|"), x))) {
      return(TRUE)
    } else {
      return(sprintf(
        "Must be a subset of %s, but is %s",
        paste0("{", paste0(choices, collapse = ", "), "}"),
        paste0("{", paste0(x, collapse = ", "), "}")))
    }
  }

  get_akritas_learner = function() {
    require_namespaces("mlr3extralearners")
    utils::getFromNamespace("LearnerSurvAkritas", "mlr3extralearners")
  }

  r6_private = function(x) {
    x$.__enclos_env__$private
  }


  Survnew = R6::R6Class("Survnew",
                        inherit = TaskSupervised,
                        public = list(
                          initialize = function(id, backend, time = "time", event = "event", time2,
                                                type = c("right", "left", "interval", "counting", "interval2", "mstate")) {

                            type = match.arg(type)

                            backend = as_data_backend(backend)

                            if (type != "interval2") {
                              c_ev = r6_private(backend)$.data[, event, with = FALSE][[1]]
                              if (type == "mstate") {
                                assert_factor(c_ev)
                              } else if (type == "interval") {
                                assert_integerish(c_ev, lower = 0, upper = 3)
                              } else if (!is.logical(c_ev)) {
                                assert_integerish(c_ev, lower = 0, upper = 2)
                              }
                            }

                            private$.censtype = type

                            if (type %in% c("right", "left", "mstate")) {
                              super$initialize(
                                id = id, task_type = "surv", backend = backend,
                                target = c(time, event))
                            } else if (type %in% c("interval", "counting")) {
                              super$initialize(
                                id = id, task_type = "surv", backend = backend,
                                target = c(time, time2, event))
                            } else {
                              super$initialize(
                                id = id, task_type = "surv", backend = backend,
                                target = c(time, time2))
                            }
                          },

                          truth = function(rows = NULL) {
                            # truth is defined as the survival outcome as a Survival object
                            tn = self$target_names
                            ct = self$censtype
                            d = self$data(rows, cols = self$target_names)
                            args = list(time = d[[tn[1L]]], type = self$censtype)

                            if (ct %in% c("right", "left", "mstate")) {
                              args$event = as.integer(d[[tn[2L]]])
                            } else if (ct %in% c("interval", "counting")) {
                              args$event = as.integer(d[[tn[3L]]])
                              args$time2 = d[[tn[2L]]]
                            } else {
                              args$time2 = d[[tn[2L]]]
                            }

                            if (allMissing(args$event) & allMissing(args$time)) {
                              return(suppressWarnings(invoke(Surv, .args = args)))
                            } else {
                              return(invoke(Surv, .args = args))
                            }
                          },

                          formula = function(rhs = NULL) {
                            # formula appends the rhs argument to Surv(time, event)~
                            tn = self$target_names
                            if (length(tn) == 2) {
                              lhs = sprintf("Surv(%s, %s, type = '%s')", tn[1L], tn[2L], self$censtype)
                            } else {
                              lhs = sprintf("Surv(%s, %s, %s, type = '%s')", tn[1L], tn[2L], tn[3L], self$censtype)
                            }
                            formulate(lhs, rhs %??% ".", env = getNamespace("survival"))
                          },

                          times = function(rows = NULL) {
                            truth = self$truth(rows)
                            if (self$censtype %in% c("interval", "counting", "interval2")) {
                              return(truth[, 1:2])
                            } else {
                              return(truth[, 1L])
                            }
                          },

                          f_status = function(rows = NULL) {
                            truth = self$truth(rows)
                            if (self$censtype %in% c("interval", "counting", "interval2")) {
                              f_status = truth[, 3L]
                            } else {
                              f_status = truth[, 2L]
                            }

                            as.integer(f_status)
                          },

                          unq_times = function(rows = NULL) {
                            if (self$censtype %in% c("interval", "counting", "interval2")) {
                              stop("Not implemented for 'interval', 'interval2', or 'counting', 'censtype'.")
                            }

                            sort(unique(self$times(rows)))
                          },

                          unq_event_times = function(rows = NULL) {
                            if (self$censtype %in% c("interval", "counting", "interval2")) {
                              stop("Not implemented for 'interval', 'interval2', or 'counting', 'censtype'.")
                            }

                            sort(unique(self$times(rows)[self$f_status(rows) != 0]))
                          },

                          risk_set = function(time = NULL) {
                            if (self$censtype %in% c("interval", "counting", "interval2")) {
                              stop("Not implemented for 'interval', 'interval2', or 'counting', 'censtype'.")
                            }

                            if (is.null(time)) {
                              self$row_ids
                            } else {
                              self$row_ids[self$times() >= time]
                            }
                          }
                        ),

                        active = list(
                          censtype = function() {
                            return(private$.censtype)
                          }
                        ),

                        private = list(
                          .censtype = character()
                        )
  )

  learn_method=mlr3::lrn("surv.coxph") #making the learner function,
  #options=kaplan, coxph
  learners=list(learn_method)
  if(per <= 0){
    stop("Wrong percentage value")
  }
  s=NULL
  for(i in cols) #cols=column numbers containing genes e.g. c(1,3,5:9,10)
  {print(i)
    f=0
    S_test=coxph(Surv(get(idSurv),get(idEvent))~data[,i],na.action=NULL,data=data)
    #if(is.na(coef(S_test))==TRUE) break
    task = Survnew$new(id = "data",
                        backend = data[,c(which(colnames(data)==idSurv),
                                          which(colnames(data)==idEvent),i)],
                        time = idSurv, event = idEvent)
    resample = rsmp("cv", folds = fold)
    design = benchmark_grid(task, learners, resample)
    sou=design$resampling[[1]]
    for(a in 1:fold)
    {
      cdata=data[sou$test_set(a),]
      if(sum(cdata[,idEvent])==0)
      { f=1
      }
    }
    if(f==1)
    {
      print("No events or all survival times are identical. Consider decreasing number of fold.")
      next
    }
    invisible(capture.output(bm <- benchmark(design)))
    c=cbind(colnames(data)[i],bm$aggregate(msr("surv.logloss")))
    s=rbind(s,c)
  }
  if(is.null(s) == TRUE){
    print("No possible gene")
  }
  else{
    colnames(s)[1]="gene"
  }
  s=head(s[order(s$surv.logloss),],(dim(s)[1]*per)/100)
  return(s)
}

utils::globalVariables(c("assert_factor","require_namespaces","assert_integerish","%??%","super","self","allMissing","invoke","formulate","private"))


