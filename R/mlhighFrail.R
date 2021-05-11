#' mlhighFrail
#'
#' Performs CoxPH frailty on high doimensional survival data
#'
#' @description This function extracts  features based on minimum log-Loss function using Cox proportional hazard model as learner method on a
#' high dimensional survival data. For those genes, we obtain frailty variances using CoxPH.
#'
#' @details Using the Cox proportional hazard model on the given survival data, this function selects the most significant feature based on minimum logarithmic loss function. The logarithmic loss function is defined as,
#' \deqn{L(f,t)=-log(f(t))}
#' After selcting the most significant features, a Cox proportional hazard frailty model is fitted on the selected features. The CoxPH frailty model is defined as,
#' \deqn{\lambda(t)=\lambda 0(t)\nu exp{X'\beta}} where \eqn{\nu} is called the frailty component. The variance of the
#' frailty term is considered as the heterogeneity among the subjects or patients. The distribution of frailty component is considered as either Gaussian, Gamma or t distribution.
#'
#' @param cols A numeric vector of column numbers indicating the features for which the log Loss functions are to be computed
#' @param idSurv The name of the survival time variable
#' @param idEvent The name of the survival event variable
#' @param idFrail The name of the frailty variable
#' @param dist The name of the frailty distribution. Options are "gamma", "gaussian" or "t", default is "gaussian"
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
#' @return A dataframe containing desired number of features with corresponding frailty variances.
#'
#' @examples
#' \donttest{
#' data(hnscc)
#' mlhighFrail(cols=c(10:20), idSurv="OS", idEvent="Death", idFrail="ID", dist="gaussian",
#' per=20, fold = 3, data=hnscc)
#' }
#' @export
#' @author Atanu Bhattacharjee, Gajendra K. Vishwakarma & Souvik Banerjee
#' @seealso mlhighHet, mlhighCox

mlhighFrail=function(cols, idSurv, idEvent, idFrail, dist="gaussian", per=20, fold=3, data)
{
  learn_method=mlr3::lrn("surv.coxph") #making the learner function,
  #options=kaplan, coxph
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
  m=NULL
  for(i in 1:dim(s)[1])
  {
    surv_cox=coxph(Surv(get(idSurv),get(idEvent))~get(s$gene[i])+frailty(get(idFrail),distribution = dist),na.action=NULL,data=data)
    f = cbind(s$gene[i],surv_cox$history$f$theta)
    m=rbind(m,f)
  }
  colnames(m)=c("gene","frailty variance")
  m=as.data.frame(m)
  m[,2]=as.numeric(m[,2])
  return(m)
}

