# CRrpart

## Summary

Recursive partitioning (RP) is a statistical procedure used to classify observations into subgroups that have a consistent outcome based on features of the observations. Starting with the complete set of observations, RP will split the input dataset into two using a binary (T/F) split based on 1 feature (e.g., Age>=70; Blood Type  = "A" | "B") that is optimal according to a statistical heuristic. Each resultant subgroup is similarly split until the potential splits are exhausted and/or other stopping conditions are met. Software has been developed to perform RP with hypothesis testing at each split for continuous and categorical outcome measures. However, competing risks (CR) outcomes, which are a special type of survival (time-to-event) analysis where multiple event types are possible to be observed, have no software implementation nor hypothesis testing included. The R package CRrpart includes RP for CR outcomes (along with survival outcomes) with statistical hypothesis testing at each split.

## CRrpart package

The objective of CRrpart is the provide an R package to perform recursive partitioning with hypothesis testing at each split for a competing risks endpoint.

It can be installed using:
```
library(remotes)
install_github("sigpvalue/CRrpart")
```
The CRrpart package uses the CR implementation from Fine and Gray (Fine and Gray, 1999) through the cmprsk R package (Gray, 2024). The main function for performing the analysis for the package is `CRrpart`, with `plot`, `predict`, and `print` routines also included (see example below using a transplant (bmtcrr) dataset from the casebase package (Bhatnagar, 2022). The naming conventions for CR response inputs (e.g., ftime, fstatus, cencode, failcode) established by the cmprsk package are followed for this package to ease integration. There are no methods to account for missing data for `CRrpart`, so only observations with complete data can be analzyed. The additional user inputted parameters for `CRrpart` include

- <u>n.splits</u> Maximum number of splits
- <u>minbucket</u> Minimum number of observations needed in each daughter node. Either minbucket or support needs to be specified.
- <u>support</u> Minimum percent of observations in the mother node needed in each daughter node. Either minbucket or support needs to be specified.
- <u>sig.level</u> Significance level 
- <u>iter</u> Number of iterations for the p-value calculation
- <u>p.adj</u> Adjust p-value for multiple testing?

Note that as the `crr` function in the cmprsk package (Gray, 2024) can analyze standard survival models (only 1 event type), `CRrpart` also will work for typical survival endpoints.

An example:

```
library(CRrpart)
g=CRrpart(minbucket=NULL,support=.2,ftime=bmtcrr$ftime,fstatus=bmtcrr$Status, cencode=0,
failcode=1,x=bmtcrr[,c('Sex','D','Phase','Age','Source')], n.splits=500, sig.level=0.05,
iter=2000,p.adj=FALSE)
print(g)
plot(g)
```

## Citations

Bhatnagar S, Turgeon M, Islam J, Saarela O, Hanley J (2022). casebase: An Alternative Framework for Survival Analysis and Comparison of Event Rates. The R Journal, 14(3).

Fine J, Gray RJ (1999). A Proportional Hazards Model for the Subdistribution of a Competing Risk. Journal of the American Statistical Association, 94(446), 496–509.

Gray RJ (2024). cmprsk: Subdistribution Analysis of Competing Risks. R package version 2.2-12. https://CRAN.R-project.org/package=cmprsk
