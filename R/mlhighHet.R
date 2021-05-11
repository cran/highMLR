#' @title mlhighHet
#'
#' Performs heterogeneity analysis in gene expression
#'
#' @description This function extracts features based on ML method, finds optimal cut-off values of features using sequencial
#' Cox PH model and obtain the most consistent level according to the cut-offs.
#'
#' @details This function extracts features based on minimum log-Loss function using Cox proportional hazard model as learner method on a high dimensional survival data.
#' For those selected genes, we obtain optimal cutoff values using minimum p-value in a Cox PH model. The Cox PH model is used sequencially for each combination of genes
#' and all possible gene combinations are tested to obtain best possible combination with minimum BIC value. The subjects are classified according to different levels of those genes. Using a Cox PH frailty model, we obtain the most consistent level for which the frailty variance is minimum.
#' The data is splited using cross validation technique. The performance measure is considered as logarithmic loss function. It is defined as,
#' \deqn{L(f,t)=-log(f(t))}
#' The CoxPH frailty model is defined as,
#' \deqn{\lambda(t)=\lambda 0(t)\nu exp{X'\beta}} where \eqn{\nu} is called the frailty. The variance of the
#' frailty term is considered as the heterogeneity among the subjects or patients. Gaussian distribution with mean 0 is considered for the distribution of frailty component.
#'
#' @param cols A numeric vector of column numbers indicating the features for which the log Loss functions are to be computed
#' @param idSurv The name of the survival time variable
#' @param idEvent The name of the survival event variable
#' @param idFrail The name of the frailty variable
#' @param num Number of features to be selected
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
#' @return dataframes containing optimal gene cutoff values and most consistent level according to those cut-offs with frailty variance.
#' @examples
#' \donttest{
#' data(hnscc)
#' mlhighHet(cols=c(27:32), idSurv="OS", idEvent="Death", idFrail="ID", num=2, fold = 3, data=hnscc)
#' }
#' @export
#' @author Atanu Bhattacharjee, Gajendra K. Vishwakarma & Souvik Banerjee
#' @seealso mlhighCox, mlhighFrail

mlhighHet=function(cols, idSurv, idEvent, idFrail, num, fold=3, data)
{
  learn_method=mlr3::lrn("surv.coxph") #making the learner function,
  learners=list(learn_method)
  if(num <= 0)
  {stop("Wrong percentage value")
  }
  s=NULL
  for(i in cols)
  {print(i)
    f=0
    S_test=coxph(Surv(get(idSurv),get(idEvent))~data[,i],na.action=NULL,data=data)
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
    stop("No possible gene")
  }else{
    colnames(s)[1]="gene"
  }
  s=head(s[order(s$surv.logloss),],num)

  ##################################################################################
  print("sequential cut-off selection started")
  permute=permutations(num,num,s$gene)
  opt_cut=data.frame()
  pb <- txtProgressBar(min = 0, max = dim(permute)[1], style = 3)
  for(q in 1:dim(permute)[1])
  {
    data[,colnames(data) %in% permute[q,]] = scale(data[,colnames(data) %in% permute[q,]], center = TRUE, scale = TRUE)
    for(r in 1:num)
    {
      assign("cuts",as.numeric(quantile(data[,colnames(data)==permute[q,r]],probs=seq(0,1,0.01))))
      xcuts=cuts[-c(which.min(cuts),which.max(cuts))]
      coxmod=data.frame()
      for(i in 1:length(xcuts))
      {
        assign(paste("gene",r,sep=""),ifelse(data[,colnames(data)==permute[q,r]]<xcuts[i],0,1))
        for(z in (r+1):num)
        {
          if(r==num) break
          assign(paste("gene",z,sep=""),ifelse(data[,colnames(data)==permute[q,z]]<quantile(data[,colnames(data)==permute[q,z]],0.5),0,1))
        }
        genename=paste("gene",1:num,sep="")
        model<-as.formula(paste("Surv(get(idSurv),get(idEvent)) ~ ", paste(genename, collapse= "+")))
        coxmodel = coxph(model,data=data)
        d=data.frame(i,xcuts[i], t(coef(coxmodel)), t(as.numeric(summary(coxmodel)[["coefficients"]][,"Pr(>|z|)"])), BIC(coxmodel))
        coxmod=rbind(coxmod,d)
      }
      colnames(coxmod)=c("i","cutoff",paste("coef_gene",1:num,sep=""),paste("p_gene",1:num,sep=""),"BIC")
      coxmod=coxmod[complete.cases(coxmod),]
      coxfix=coxmod[coxmod[,colnames(coxmod)==paste("p_gene",r,sep="")]==min(coxmod[,colnames(coxmod)==paste("p_gene",r,sep="")]),][1,]
      assign(paste("gene",r,sep=""),ifelse(data[,colnames(data)==permute[q,r]]<coxfix$cutoff,0,1))
      assign(paste("opt_cutoff_gene",r,sep=""),coxfix$cutoff)
    }
    fun<-function(x,y)
    {
      y<-get(paste(y,x,sep=""))
      return(y)
    }
    opcut = data.frame(t(permute[q,]),t(mapply(fun,1:num,"opt_cutoff_gene")),coxfix$BIC)
    colnames(opcut)=c(paste("gene",1:num,sep=""),paste("opt_cutoff_gene",1:num,sep=""),"BIC")
    opt_cut=rbind(opt_cut,opcut)
    setTxtProgressBar(pb, q)
  }
  close(pb)
  print("sequential cut-off selection finished")
  final=opt_cut[opt_cut$BIC==min(opt_cut$BIC),][1,]
  ###################################################################################
  name=colnames(data)
  data1=data
  for(i in 1:num)
  {
    assign(paste("gene_class",i,sep=""),ifelse(data1[,colnames(data1)==final[,colnames(final)==paste("gene",i,sep="")]]<final[,colnames(final)==paste("opt_cutoff_gene",i,sep="")],0,1))
    data1=cbind(data1,get(paste("gene_class",i,sep="")))
  }
  colnames(data1)<-c(name,paste("gene_class",1:num,sep=""))
  frailID=data1 %>% group_by(mapply(fun,1:num,"gene_class")) %>% group_indices()
  data1=cbind(data1,frailID)
  len=unique(data1$frailID)
  fid=paste("frailty.gaussian(",idFrail,")",sep="")
  frail_data=data.frame()
  for(l in 1:length(len))
  {
    if (dim(data1[data1$frailID==len[l],])[1]>5)
    {
      coxfrail=coxph(as.formula(paste(paste("Surv(get(idSurv),get(idEvent)) ~ ",
                                            paste(final[,colnames(final) %in% paste("gene",1:num,sep="")],
                                                  collapse="+")),fid,sep="+")),data=data1[data1$frailID==len[l],])
      u=cbind(len[l],coxfrail$history$f$theta)
      colnames(u)=c("frail_id","frail_variance")
      frail_data=rbind(frail_data,u)
    }
  }
  opt_frail=frail_data[frail_data$frail_variance==min(frail_data$frail_variance),]
  opt1=data1[data1$frailID %in% opt_frail$frail_id,]
  opt_f=cbind(final[,colnames(final) %in% paste("gene",1:num,sep="")],select(final,paste("opt_cutoff_gene",1:num,sep="")),
              head(select(opt1,paste("gene_class",1:num,sep="")),1),opt_frail$frail_variance)[1,]
  return(list(optimal_cutoff=final, consistency_level=opt_f))
}
