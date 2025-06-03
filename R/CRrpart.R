#' @import cmprsk
NULL
#' @import casebase
NULL
#' @import gtools
NULL

#' Recursive partitioning with hypothesis testing at each split for a competing risks endpoint
#'
#' Recursively partition a dataset with competing risks endpoint using the maximum pseudo log-likelihood among potential splits to determine optimal split. The partitioning will stop when there are not more statistically significant splits. The package uses the competing risk analysis formulation and functions from the 'cmprsk' package.
#' @param x Predictor variables
#' @param ftime Failure/censoring time variable
#' @param fstatus Failure/censoring event variable
#' @param sig.level Significance level (default = 0.05)
#' @param n.splits Maximum number of splits (default =500)
#' @param minbucket Minimum number of observations needed in each daughter node. Either minbucket or support needs to be specified. (default = NULL)
#' @param support Minimum percent of observations needed in each daughter node. Either minbucket or support needs to be specified. (default = NULL)
#' @param failcode Code of fstatus corresponding to the event of interest. (default = 1)
#' @param cencode Code of fstatus corresponding to censored observations. (default = 0)
#' @param iter Number of iterations for p-value calculations. (default = 2000)
#' @param p.adj Adjust p-value for multiple testing? (default = FALSE). The adjustment is based upon the the distribution of the maximum log-likelihood statistic. It adjusts based on the number of input variables, not the number of valid data splits.
#' @return def: Terms used to create the recursive partitions
#' @return pvalue: pvalues for the significant terms
#' @return call: Input function call parameters and data
#' @return class.id: Subgroup identifier
#' @seealso \link[cmprsk]{crr}, \link[cmprsk]{cuminc}
#' @examples
#' \dontrun{
#' library(CRrpart)
#' g=CRrpart(minbucket=NULL,support=.2,ftime=bmtcrr$ftime,fstatus=bmtcrr$Status,
#' cencode=0,failcode=1,x=bmtcrr[,c('Sex','D','Phase','Age','Source')],
#' n.splits=500, sig.level=0.05,iter=2000,p.adj=FALSE)
#' }
#' @export

CRrpart<-function(x,ftime,fstatus,sig.level=0.05,n.splits=500,minbucket=NULL,support=NULL,
                  failcode=1,cencode=0,iter=2000,p.adj=FALSE)
{
  if(is.null(minbucket) & is.null(support)) stop('Either minbucket or support must be specified, not neither')
  if(!is.null(minbucket) & !is.null(support)) stop('Either minbucket or support must be specified, not both')

  if (p.adj==FALSE) p.use=1
  if (p.adj==TRUE) p.use=2

  y=cbind(ftime,fstatus)

  if(sum(is.na(x)>0)) stop('Missing values in predictor variables')
  if(sum(is.na(y)>0)) stop('Missing values in response variable')

  call<-list("x"=deparse(substitute(x)),"y"=deparse(substitute(y)),"support"=support,
             "sig.level"=sig.level,'minbucket'=minbucket,'failcode'=failcode,
             'cencode'=cencode,'p.adj'=p.adj)
  # number and levels of predictor variables
  x.vars<-ncol(x)
  x.var.def<-rep(FALSE,x.vars)
  for (i in 1:x.vars) x.var.def[i]<-is.factor(x[,i])
  y.org<-y
  x.org<-x

  iidd<-rep("0.",nrow(x.org))
  status=rep('open',nrow(x.org))
  rpart.def<-matrix("",n.splits,4)
  pvalue.out=matrix(NA,n.splits,2)
  colnames(pvalue.out)=c('Global_unadj','Global_adj')

  # loop through the maximum number of partitions to be created
  for (jj in 1:n.splits)
  {
    lev.used=min(as.numeric(iidd[status=='open']))
    ids.to.use=which(as.numeric(iidd)==lev.used & status=='open')
    if (length(ids.to.use)==2) break

    x<-x.org[ids.to.use,,drop=FALSE]
    y<-y.org[ids.to.use,,drop=FALSE]

    result=NULL
    #need at least 10 observations, with 2 failures and needing at least 5 observations per daughter node
    minnumber=ifelse(is.null(support),minbucket,ceiling(support*nrow(y)))
    if (nrow(y)<minnumber*2 | minnumber<5|nrow(y)<10 | sum(y[,2]==failcode)<=1) result$pvalue=1
    if (nrow(y)>=minnumber*2 & minnumber>=5 & nrow(y)>=10 & sum(y[,2]==failcode)>1) result<-splits.CRrpart(x.=x,y.=y,x.var.def=x.var.def,minbuck=minnumber,failcode=failcode,cencode=cencode,iter=iter)

    pvalue.out[jj,]=result$pvalue

    if (pvalue.out[jj,p.use]<sig.level)
    {
      rpart.def[jj,1]=unique(iidd[ids.to.use])
      iidd[ids.to.use]=paste(iidd[ids.to.use],result$iidd,sep='')
      rpart.def[jj,2]=paste(sort(unique(iidd[ids.to.use])),collapse='_')
      rpart.def[jj,3]=result$best.def
      rpart.def[jj,4]=result$best
    }
    if (pvalue.out[jj,p.use]>=sig.level) status[ids.to.use]='closed'

    if (all(status=='closed')) break
  }
  sig.node=which(rpart.def[,1]!='')
  if (length(sig.node)>0) ans=list('def'=rpart.def[sig.node,-4,drop=F],'pvalue'=pvalue.out[sig.node,],'call'=call,'class.id'=iidd)
  if (length(sig.node)==0) ans=list('def'="None",'pvalue'=pvalue.out[1],'call'=call,'class.id'=iidd)
  class(ans)<-"CRrpart"
  return(ans)
}

#' Find the best split  (function not to be called)
#'
#' Find the best split among potential splits
#' @param x. Predictor variables
#' @param y. Response variables
#' @param x.var.def identify categorical predictor variable
#' @param minbuck Minimum number of observations needed in each daughter node
#' @param failcode Code of fstatus corresponding to the event of interest
#' @param cencode Code of fstatus corresponding to censored observations
#' @param iter Number of iterations for p-value calculations. (default = 2000)
#' @export

splits.CRrpart=function(x.,y.,x.var.def,minbuck,failcode,cencode,iter)
{
    best=0
    best.def=''
    iidd=rep(0,nrow(y.))
    if (sum(!x.var.def)>=1)
    {
      for (i in which(!x.var.def)) #continuous variables
      {
        cont.tab=t(t(table(x.[,i])))
        cont.tab.cs=cumsum(cont.tab)
        min.stt=min(1+which(cont.tab.cs>=minbuck))
        if (sum(nrow(x.)-cont.tab.cs>=minbuck)==0) max.end=1
        if (sum(nrow(x.)-cont.tab.cs>=minbuck)>0) max.end=max(which(nrow(x.)-cont.tab.cs>=minbuck))+1
        if(min.stt<=nrow(cont.tab))
        {
          starts=as.numeric(names(cont.tab[min.stt,]))
          ends=as.numeric(names(cont.tab[max.end,]))
          cut.pts=sort(unique(x.[which(x.[,i]>=starts & x.[,i]<=ends),i]))
          if (length(cut.pts)>0)
          {
            for (j in 1:length(cut.pts))
            {
              g=crr(ftime=y.[,1],fstatus=y.[,2],cov1=as.numeric(x.[,i]>=cut.pts[j]),failcode=failcode,cencode=cencode)
              if (2*(g$loglik-g$loglik.null)>best)
              {
                best=2*(g$loglik-g$loglik.null)
                if (g$coef[1]>0)
                {
                  best.def=paste(names(x.)[i],'>=',cut.pts[j])
                  iidd=as.numeric(x.[,i]>=cut.pts[j])
                }
                if (g$coef[1]<0)
                {
                  best.def=paste(names(x.)[i],'<',cut.pts[j])
                  iidd=1-as.numeric(x.[,i]>=cut.pts[j])
                }
              }
            }
          }
        }
      }
    }
    if (sum(x.var.def)>=1)
    {
      for (i in which(x.var.def)) #categorical variables
      {
        n=nlevels(x.[,i])
        lst=list()
        if (n==2) lst[[1]]=1
        if (n>2 & n%%2==0)
        {
          for (j in 1:(n/2-1))
          {
            lst[[j]]=combinations(n,j)
          }
          tmp=combinations(n,n/2)
          tmp=tmp[which(tmp[,1]==1),]
          lst[[n/2]]=tmp
        }
        if (n>2 & n%%2==1)
        {
          for (j in 1:((n-1)/2))
          {
            lst[[j]]=combinations(n,j)
          }
        }
        for (j in 1:length(lst))
        {
          for (k in 1:length(lst[[j]]))
          {
            test.var=rep(0,nrow(x.))
            test.var[which(x.[,i] %in% levels(x.[,i])[lst[[j]][k]])]=1
            if(sum(test.var)>=minbuck & sum(1-test.var)>=minbuck)
            {
              g=crr(ftime=y.[,1],fstatus=y.[,2],cov1=test.var,failcode=failcode,cencode=cencode)
              if (2*(g$loglik-g$loglik.null)>best)
              {
                best=2*(g$loglik-g$loglik.null)
                if (g$coef[1]>0)
                {
                  best.def=paste(names(x.)[i],'%in%',paste(levels(x.[,i])[lst[[j]][k]],collapse=';'))
                  iidd=test.var
                }
                if (g$coef[1]<0)
                {
                  best.def=paste(names(x.)[i],'%in%',paste(levels(x.[,i])[-lst[[j]][k]],collapse=';'))
                  iidd=1-test.var
                }
              }
            }
          }
        }
      }
    }

    out=rep(NA,iter+1)
    out[1]=best
    for (j in 2:(iter+1))
    {
      m=nrow(x.)
      n=sample(minbuck:(m-minbuck),1)
      t.var=rep(0,m)
      t.var[sample(m,n)]=1
      g=crr(ftime=y.[,1],fstatus=y.[,2],cov1=t.var,failcode=1,cencode=0)
      out[j]=2*(g$loglik-g$loglik.null)
    }
    return(list(best=best,best.def=best.def,iidd=iidd,pvalue=c(mean(out[1]<=out[-1]),
    1-(1-mean(out[1]<=out[-1]))^ncol(x.))))
}

#' Print an CRrpart object
#'
#' Print an CRrpart object
#' @param x CRrpart object
#' @param failcode Code of fstatus corresponding to the event of interest (default = 1)
#' @param tim The time at which to report cumulative incidence of event failcode
#' @param ... arguments to be passed to methods
#' @examples
#' \dontrun{
#' library(CRrpart)
#' g=CRrpart(minbucket=NULL,support=.2,ftime=bmtcrr$ftime,fstatus=bmtcrr$Status,
#' cencode=0,failcode=1,x=bmtcrr[,c('Sex','D','Phase','Age','Source')],
#' n.splits=500, sig.level=0.05,iter=2000,p.adj=FALSE)
#' }
#' @export

print.CRrpart=function(x,tim=24,failcode=object$call$failcode,...)
{
  object=x
  if (object$def[1]=='None')
  {
    out=matrix('',1,6)
    colnames(out)=c('Mother','Daughter','Cut',paste('CumInc_',tim,sep=''),'N','Terminal')
    g=cuminc(ftime=y[,1],fstatus=y[,2])
    g1=timepoints(g,tim)$est[failcode]
    out[1,]=c('','','',sprintf(g1,fmt='%.3f'),nrow(y),'')
    out[1,6]='***'
    out[1,1]='0'
    print(noquote(out))
  }
  if (object$def[1]!='None')
  {
    y=eval(parse(text=object$call$y))
    x=eval(parse(text=object$call$x))
    n=max(nchar(object$class.id))-2
    g0=model.matrix(lm(y[,1]~factor(object$class.id)))[,-1,drop=F]
    colnames(g0)=levels(factor(object$class.id))[-1]
    g=cuminc(ftime=y[,1],fstatus=y[,2])
    g1=timepoints(g,min(max(y[,1]),tim))$est[failcode]

    out=matrix('',n*8-1,6)
    colnames(out)=c('Mother','Daughter','Cut',paste('CumInc_',tim,sep=''),'N','Terminal')
    out[1,]=c('','','',sprintf(g1,fmt='%.3f'),nrow(y),'')

    k=2
    for (i in 1:n)
    {
      if(i==1) tmp.class.m=rep('0',length(object$class.id))
      if(i>1) tmp.class.m=paste('0',substr(object$class.id,3,2+i-1),sep='')
      tmp.class.d=substr(object$class.id,2+i,2+i)
      use.lev=unique(tmp.class.m[which(tmp.class.d %in% c('0','1'))])
      for (j in 1:length(use.lev))
      {
        out[k,1]=out[k+1,1]=use.lev[j]
        out[k,2]=0
        out[k+1,2]=1
        def.id=which(gsub('\\.','',object$def)==use.lev[j])
        out[k,3]=paste('!(',object$def[def.id,3],')',sep='')
        out[k+1,3]=paste('(',object$def[def.id,3],')',sep='')
        ids=sort(unique(paste('0.',substr(tmp.class.m,2,max(nchar(tmp.class.m))),tmp.class.d,sep='')[tmp.class.d %in% c('0','1') & tmp.class.m==use.lev[j]]))
        n1=nchar(ids[1])
        id1=which(substr(object$class.id,1,n1) %in% ids[1])
        id2=which(substr(object$class.id,1,n1) %in% ids[2])

        if(sum(y[id1,2]==failcode)>0)
        {
          g=cuminc(ftime=y[,1],fstatus=y[,2],subset=id1)
          tp=timepoints(g,min(max(y[id1,1]),tim))$est
          tp.id=which(rownames(tp)==paste('1',failcode))
          out[k,4]=sprintf(timepoints(g,min(max(y[id1,1]),tim))$est[tp.id],fmt='%.3f')
        }
        if(sum(y[id1,2]==failcode)==0)
        {
          out[k,4]=sprintf(0,fmt='%.3f')
        }
        if(sum(y[id2,2]==failcode)>0)
        {
          g=cuminc(ftime=y[,1],fstatus=y[,2],subset=id2)
          tp=timepoints(g,min(max(y[id2,1]),tim))$est
          tp.id=which(rownames(tp)==paste('1',failcode))
          out[k+1,4]=sprintf(timepoints(g,min(max(y[id2,1]),tim))$est[tp.id],fmt='%.3f')
        }
        if(sum(y[id2,2]==failcode)==0)
        {
          out[k+1,4]=sprintf(0,fmt='%.3f')
        }
        out[k,5]=length(id1)
        out[k+1,5]=length(id2)
        k=k+2
      }
    }
    out=out[which(out[,4]!=''),]
    rownames(out)=rep('',nrow(out))
    out[which(paste('0.',substr(out[,1],2,max(nchar(out[,1]))),out[,2],sep='') %in% unique(object$class.id)),6]='***'

    md=matrix(NA,(nrow(out)-1)/2,3)
    colnames(md)=c('M','D0','D1')
    lev=paste(out[,1],out[,2],sep='')
    lev[1]='0'
    for (i in seq(2,nrow(out)-1,by=2))
    {
      md[i/2,1]=out[i,1]
      md[i/2,2]=lev[i]
      md[i/2,3]=lev[i+1]
    }
    colnames(md)=c('M','D0','D1')
    rd=matrix(NA,nrow(md),2)
    colnames(rd)=c('R0','R1')
    rd[1,]=c(1,2)
    if (nrow(rd)>1)
    {
      for (i in 2:nrow(rd))
      {
        tmp=md[i,1]
        ind=which(md[,-1]==tmp,arr.ind=T)
        rd[i,]=paste(rd[ind],c(1,2),sep='')
      }
    }
    for (i in 1:nrow(rd))
    {
      for (j in 1:2)
      {
        rd[i,j]=paste(substr(rd[i,j],1,1),'.',substr(rd[i,j],2,nchar(rd[i,j])),sep='')
      }
    }
    o=order(matrix(t(rd)))
    out=rbind(out[1,],out[-1,][o,])
    for (i in 2:nrow(out))
    {
      out[i,1]=paste(paste(rep(' ',nchar(out[i,1])),collapse=''),out[i,1],sep='')
    }
    print(noquote(out))
  }
}

#' Predict output classes to a new database for CRrpart object
#'
#' Predict output classes to a new database
#' @param object CRrpart object
#' @param newdata new data with same input predictors as the CRrpart object
#' @return outclass: output class for observation
#' @examples
#' \dontrun{
#' library(CRrpart)
#' g=CRrpart(minbucket=NULL,support=.2,ftime=bmtcrr$ftime,fstatus=bmtcrr$Status,
#' cencode=0,failcode=1,x=bmtcrr[,c('Sex','D','Phase','Age','Source')],
#' n.splits=500,sig.level=0.05,iter=2000,p.adj=FALSE)
#' predict(object=g)
#' }
#' @export

predict.CRrpart=function(object,newdata=NULL)
{
  if (is.null(newdata)) newdata=eval(parse(text=object$call$x))
  if (sum(colnames(newdata)%in% colnames(eval(parse(text=object$call$x))))<ncol(newdata)) stop('Colnames in test dataset are not in training data set')
  out.tf=matrix(F,nrow(newdata),nrow(object$def))
  for (i in 1:ncol(out.tf))
  {
    if (length(grep('%in%',object$def[i,3]))==1)
    {
      txt.tmp=strsplit(object$def[i,3],' %in% ')[[1]]
      tt=paste('c("',paste(strsplit(txt.tmp[2],';')[[1]],collapse='","'),'")',sep='')
      out.tf[,i]=eval(parse(text=paste('newdata$',txt.tmp[1],'%in%',tt,sep='')))
    }
    if (length(grep('%in%',object$def[i,3]))==0)
    {
      out.tf[,i]=eval(parse(text=paste('newdata$',object$def[i,3],sep='')))
    }
  }
  md=matrix('',nrow(object$def),3)
  md[,1]=object$def[,1]
  for (i in 1:nrow(object$def))
  {
    md[i,2:3]=strsplit(object$def[i,2],'_')[[1]]
  }
  out.class=rep('',nrow(newdata))
  for (i in 1:nrow(newdata))
  {
    begin='0.'
    col.sel=which(md[,1]==begin)
    while(length(col.sel)>0)
    {
      begin=md[col.sel,as.numeric(out.tf[i,col.sel])+2]
      col.sel=which(md[,1]==begin)
    }
    out.class[i]=begin
  }
  return(out.class)
}

#' Plot the cumulative incidence by output classes
#'
#' Predict output classes to a new database
#' @param x CRrpart object
#' @param failcode Code of fstatus corresponding to the event of interest (default = 1)
#' @param ... arguments to be passed to methods
#' @examples
#' \dontrun{
#' library(CRrpart)
#' g=CRrpart(minbucket=NULL,support=.2,ftime=bmtcrr$ftime,fstatus=bmtcrr$Status,
#' cencode=0,failcode=1,x=bmtcrr[,c('Sex','D','Phase','Age','Source')],
#' n.splits=500,sig.level=0.05,iter=2000,p.adj=FALSE)
#' plot(object=g)
#' }
#' @export

plot.CRrpart=function(x,failcode=1,...)
{
   object=x
   y=eval(parse(text=object$call$y))

   ti=y[,1]
   ev=y[,2]
   a=as.numeric(as.factor(object$class.id))

   h3=cuminc(ftime=ti,fstatus=ev,a,cencode=0)

   par(mar=c(8,4,4,2))
   plot(h3,ylab='',xaxt='n',yaxt='n',
   main='Cumulative incidence function',
   curvlab=rep('',3),col='white',xlab='')
   box()
   clz=c('blue','gold','purple','darkgreen','lightgreen','black','red','orange','cyan','pink')

   st=which(rownames(timepoints(h3,1)$est)==paste('1',failcode))-1
   for (i in 1:max(a))
   {
     points(h3[[st+i]]$time,h3[[st+i]]$est,ty='l',lwd=2,lty=1,col=clz[i])
   }

   m.y=max(y[,1])
   if (m.y<=24) spc=3
   if (m.y>24 & m.y<=48) spc=6
   if (m.y>48 & m.y<=84) spc=12
   if (m.y>84 & m.y<=192) spc=24
   if (m.y>192 & m.y<=288) spc=36
   if (m.y>288 & m.y<=384) spc=48
   if (m.y>384 & m.y<=480) spc=60
   if (m.y>480) stop('Rescale the event time, the maximum time is greater than 480')

   seqs=seq(0,m.y,by=spc)

   l0=matrix(0,length(seqs),max(a))
   for (j in 1:max(a))
   {
     for (i in 1:nrow(l0))
     {
       l0[i,j]=sum(ti[a==j]>=seqs[i],na.rm=T)
     }
   }
   par(cex.axis=.8)
   axis(2,at=seq(0,1,by=.2),labels=F)
   axis(2,at=seq(0,1,by=.2),tick=F,line=-.5,cex.axis=.7)
   mtext('Cumulative incidence',side=2,line=1.2,cex=.75,font=2)

   axis(1,at=seqs,line=-.75,tick=F,cex.axis=.7)
   axis(1,at=seqs,labels=F,cex.axis=.7)
   mtext('Time (months)',side=1,line=1,cex=.7,font=2)

   at.=-spc*5/12
   mtext(side=1,line=1.4,at=at.,text='N at risk',cex=.7)

   for (i in 1:max(a))
   {
     axis(1,line=1+i/2,at=seqs,l0[,i],tick=F,col.axis=clz[i],cex.axis=.6)
     mtext(side=1,line=2+i/2,at=at.,text=levels(factor(object$class.id))[i],col=clz[i],cex=.6)
   }
   legend('topright',col=clz[1:max(a)],lty=rep(c(1),each=max(a)),legend=c(paste(levels(factor(object$class.id)),' event =',failcode)),cex=.7)
}
