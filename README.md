# CRrpart

The objective of this CRrpart is the provide an R package to perform recursive partitioning with hypothesis testing at each split for a competing risks endpoint.

It can be installed using:
```
library(remotes)
install_github("sigpvalue/CRrpart")
```

An example:

```
library(CRrpart)
g=CRrpart(minbucket=NULL,support=.2,ftime=bmtcrr$ftime,fstatus=bmtcrr$Status, cencode=0,failcode=1,
x=bmtcrr[,c('Sex','D','Phase','Age','Source')], n.splits=500, sig.level=0.05,iter=2000,p.adj=FALSE)
print(g)
plot(g)
```
