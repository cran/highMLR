#' highMLR: Machine Learning Feature Selection for High Dimensional Survival Data
#'
#' A unified, flexible framework for high dimensional feature selection in the
#' presence of a survival outcome. Provides multiple machine learning approaches
#' under a single interface: Cox elastic net, random survival forest,
#' accelerated oblique RSF, gradient-boosted Cox, stability selection,
#' classical univariate Cox screening, pseudo-observation bridging to any
#' regression learner, and Fine-Gray competing risks selection. Adds causal
#' survival forest estimation of heterogeneous treatment effects, conformal
#' survival prediction intervals, and time-dependent SHAP explanations via
#' SurvSHAP(t).
#'
#' @section Main functions:
#' \describe{
#'   \item{[highmlr()]}{Main entry point. Fit one of eight ML methods.}
#'   \item{[highmlr_compare()]}{Compare multiple methods side by side.}
#'   \item{[highmlr_stability()]}{Stability selection wrapper.}
#'   \item{[highmlr_explain()]}{Time-dependent SHAP via SurvSHAP(t).}
#'   \item{[highmlr_screen()]}{Pre-screening for very high p.}
#'   \item{[highmlr_report()]}{Generate a Quarto/Rmd report.}
#'   \item{[highmlr_causal()]}{Causal survival forest (experimental).}
#'   \item{[highmlr_conformal()]}{Conformal prediction intervals.}
#' }
#'
#' @section Bundled datasets:
#' \describe{
#'   \item{[hnscc]}{High dimensional head and neck cancer survival data.}
#'   \item{[srdata]}{High dimensional protein gene expression data.}
#' }
#'
#' @author Atanu Bhattacharjee \email{atanustat@gmail.com}
#' @keywords internal
#' @importFrom stats predict coef as.formula complete.cases sd quantile setNames
#' @importFrom utils head tail
#' @importFrom rlang .data abort warn inform
"_PACKAGE"

#' High dimensional head and neck cancer survival and gene expression data
#'
#' Survival and gene expression measurements for head and neck squamous
#' cell carcinoma patients, used to demonstrate high-dimensional feature
#' selection.
#'
#' @format A data frame with 565 rows (one per patient) and 104 columns.
#'   The first five columns are the identifier and outcome variables:
#'   \code{ID} (patient identifier), \code{Death} (overall survival event
#'   indicator, 1 = death, 0 = censored), \code{OS} (overall survival
#'   time), \code{PFS} (progression-free survival time), and \code{Prog}
#'   (progression event indicator, 1 = progression, 0 = none). The
#'   remaining 99 columns are numeric gene expression features named by
#'   gene symbol (for example \code{GJB1}, \code{HPN}, \code{PROM1}).
#' @source Bundled with the package since highMLR v0.1.1.
"hnscc"

#' High dimensional protein gene expression survival data
#'
#' Protein expression measurements with a survival outcome, used to
#' demonstrate high-dimensional feature selection.
#'
#' @format A data frame with 288 rows and 250 columns. The first four
#'   columns are the identifier and outcome variables: \code{ID} (subject
#'   identifier), \code{Visit} (visit number), \code{OS} (overall survival
#'   time), and \code{event} (survival event indicator, 1 = event,
#'   0 = censored). The remaining 246 columns are numeric protein
#'   expression features named by protein or marker (for example
#'   \code{C6kine}, \code{ActivinA}, \code{Adiponectin}).
#' @source Bundled with the package since highMLR v0.1.1.
"srdata"
