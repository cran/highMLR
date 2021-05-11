#' mlhighKap
#'
#' Performs feature selection using Kaplan Meier method
#'
#' @description This function extracts desired number of features based on minimum log-Loss function using Kaplan Meier model as learner
#' method on a high dimensional survival data.
#'
#' @details Using the Kaplan Meier method on the given survival data, this function selects the most significant feature based on a performance measure. The performance measure is considered as logarithmic loss function. It is defined as,
#' \deqn{L(f,t)=-log(f(t))}. The features with minimum log-loss function are extracted.
#'
#' @param cols A numeric vector of column numbers indicating the features for which the log Loss functions are to be computed
#' @param idSurv The name of the survival time variable
#' @param idEvent The name of the survival event variable
#' @param per Percentage of features to be selected, default value 20
#' @param fold An integer denoting number of folds in cross validation, default value 3
#' @param data A data frame that contains the survival and covariate information for the subjects
#'
#' @import mlr3
#' @import mlr3proba
#' @import mlr3learners
#' @import survival
#' @import utils
#' @import gtools
#' @import dplyr
#' @importFrom stats coef as.formula quantile BIC complete.cases
#' @return A dataframe containing desired number of features based on minimum log Loss function
#'
#' @examples
#' \donttest{
#' data(hnscc)
#' mlhighKap(cols=c(6:15), idSurv="OS", idEvent="Death", per=20, fold = 3, data=hnscc)
#' }
#' @export
#' @author Atanu Bhattacharjee, Gajendra K. Vishwakarma & Souvik Banerjee
#' @seealso mlhighCox

mlhighKap=function(cols, idSurv, idEvent, per=20, fold=3, data)
{
  learn_method=mlr3::lrn("surv.kaplan") #making the learner function,
  #options=kaplan
  learners=list(learn_method)
  if(per <= 0)
  {stop("Wrong percentage value")
  }
  s=NULL
  for(i in cols) #cols=column numbers containing genes e.g. c(1,3,5:9,10)
  {print(i)
    f=0
    S_test=coxph(Surv(get(idSurv),get(idEvent))~data[,i],na.action=NULL,data=data)
    #if(is.na(coef(S_test))==TRUE) break
    task = TaskSurv$new(id = "data",
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
