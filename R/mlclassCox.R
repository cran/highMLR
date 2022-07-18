#' Applications of machine learning in survival analysis by prognostic classification of genes by CoxPH model.
#'
#' @param m Starting column number from where high dimensional variates to be selected.
#' @param n Ending column number till where high dimensional variates to be selected.
#' @param idSurv "Column/Variable name" consisting duration of survival.
#' @param idEvent "Column/Variable name" consisting survival event.
#' @param Time "Column/Variable name" consisting Times of repeated observations.
#' @param s_ID "Column/Variable name" consisting unique identification for each subject.
#' @param per Percentage value for ordering, default=20.
#' @param fold Number of folds for re-sampling, default=3.
#' @param data High dimensional data containing survival observations with multiple covariates.
#' @return A list of genes as per their classifications
#' \describe{
#'   \item{GeneClassification}{List of genes classified using Cox proportional hazard model}
#'   \item{GeneClassification$Positive_Gene}{Sublist of genes classified as positive genes}
#'   \item{GeneClassification$Negative_Gene}{Sublist of genes classified as negative genes}
#'   \item{GeneClassification$Volatile_Gene}{Sublist of genes classified as volatile genes}
#'   \item{Result}{A dataframe consisting threshold values with corresponding coefficients and p-values.}
#' }
#' @import tibble
#' @import survival
#' @import coxme
#' @import mlr3
#' @import missForest
#' @export
#'
#' @examples
#' \dontrun{
#' data(srdata)
#' mlclassCox(m=50,n=59,idSurv="OS",idEvent="event",Time="Visit",s_ID="ID",per=20,fold=3,data=srdata)
#' }
mlclassCox<-function(m,n,idSurv,idEvent,Time,s_ID,per=20,fold=3,data)
{
  data1 <- data
  data2 <- data
  w=m
  t1<-c(unique(data1[,Time]))

  tlast <- t1[length(t1)]
  data3 <- subset(data1, get(Time) == tlast)
  smr<-summary(data3[,idEvent])

  if(smr[6]!=1){
    print("Event should be categorized as 0 and 1 only")
    }
  data303 <- subset(data3,select = c(get(s_ID),get(idSurv),get(idEvent),get(Time),m:n))
  m=5
  n=ncol(data303)

  if(sum(is.na(data303))>0){
    data.imp <- missForest(data303) #imputed tmc longit
    data3330 <- data.frame(data.imp$ximp)
    data3 <- data3330
  }
  if(sum(is.na(data303))==0){
    data3 <- data303
  }

  if(m<n){
  data <- data3
    m=5
    n=ncol(data)
  cols<-c(m:n)
  learn_method=mlr3::lrn("surv.coxph") #making the learner function,
  learners=list(learn_method)
  if(per <= 0){
    stop("Wrong percentage value")
  }
  s=NULL
  for(i in cols) #cols=column numbers containing genes e.g. c(1,3,5:9,10)
  {
    print(w-m+i)
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
    print("No possible gene")
  }else{
    colnames(s)[1]="gene"
  }
  s=head(s[order(s$surv.logloss),],(dim(s)[1]*per)/100)
  s1<-c(s$gene)
  }
  newd1 <- subset(data1,select=c(get(s_ID),get(idSurv),get(idEvent),get(Time)))
  newd2 <- subset(data1,select=c(s1))
  newdata <- data.frame(newd1,newd2)

  if(sum(is.na(newdata))>0){
    data.imp <- missForest(newdata)
    ndata <- data.frame(data.imp$ximp)
    newdata1 <- ndata
  }
  if(sum(is.na(newdata))==0){
    newdata1 <- newdata
  }

  result1=c()
  min_result1=data.frame(dummy=1)

  positive_gene_frail=c()
  negative_gene_frail=c()
  volatile_gene_frail=c()

  m1=5;n1=ncol(newdata1)
  names(newdata1)[names(newdata1)==Time]<-"time"
  for(i in m1:n1)
  {
    min_output1=data.frame()

    output_b=data.frame()
    coxdata=newdata1[newdata1[,i]>0,]
    t=unique(newdata1$time)

    for(l in 1:length(t))
    {
      coxdata1=coxdata[coxdata$time==t[l],]
      p=range(coxdata1[,i])
      q=quantile(coxdata1[,i],probs = seq(0,1,0.1))
      output=data.frame()
      output1=data.frame()
      for(j in 2:10)
      {
        indicator=c()
        for(k in 1:length(coxdata1[,i]))
        {
          if(coxdata1[k,i]<=q[j])
          {
            indicator[k]=0
          }else
          {
            indicator[k]= 1
          }
        }
        coxdata2=cbind(coxdata1,indicator)
        non_zero=table(coxdata2$indicator,coxdata2[,idEvent])
        if(sum(non_zero==0)==0)
        {
          fit_frail=coxme(Surv(get(idSurv),get(idEvent))~factor(indicator)+(1|get(s_ID)),data=coxdata2,x=TRUE)
          extract_coxme_table <- function (mod){
            beta <- mod$coefficients
            nvar <- length(beta)
            nfrail <- nrow(mod$var) - nvar
            se <- sqrt(diag(as.matrix(mod$var))[nfrail + 1:nvar])
            z<- round(beta/se, 10)
            p<- signif(1 - pchisq((beta/se)^2, 1), 10)
            return(p)
          }
          n=cbind(j,q[j],as.numeric(fit_frail$coefficients),as.numeric(extract_coxme_table(fit_frail)),sqrt(as.numeric(fit_frail$vcoef)))
          output1=rbind(output1,n)
        }
      }
      if(dim(output1)[1]==0) next
      f=cbind(t[l],output1)
      output_b=rbind(output_b,f)

      min_f=cbind(t[l],na.omit(output1[output1[,4]==min(output1[,4],na.rm=T),])[1,])
      min_output1=rbind(min_output1,min_f)
    }
    colnames(min_output1)<-c("time","partition","threshold",paste("min",colnames(newdata1)[i],"coef",sep='_'),paste("min",colnames(newdata1)[i],"p_value",sep='_'), paste("min",colnames(newdata1)[i],"random_part_sd",sep='_'))
    min_result1=cbind(min_result1,min_output1)

    if(all(as.numeric(as.character(min_output1[,4]))>0))
    {
      positive_gene_frail=rbind(positive_gene_frail,colnames(newdata1)[i])
    }else if (all(as.numeric(as.character(min_output1[,4]))<0))
    {
      negative_gene_frail=rbind(negative_gene_frail,colnames(newdata1)[i])
    }else
    {
      volatile_gene_frail=rbind(volatile_gene_frail,colnames(newdata1)[i])
    }

  }
  min_result1=min_result1[,-1]
  min_res1=list(result=min_result1)

  classify_by_Gene<-list(Positive_Gene=positive_gene_frail,Negative_Gene=negative_gene_frail,Volatile_Gene=volatile_gene_frail)

  Classified<-list(Result=min_res1, GeneClassification=classify_by_Gene)
  return(Classified)
}
utils::globalVariables(c("capture.output","head","quantile","pchisq","na.omit","TaskSurv"))


