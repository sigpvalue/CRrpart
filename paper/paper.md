---
title: 'Recursive partitioning with hypothesis testing at each split for a competing risks endpoint'
author:
  - name: Greg Dyson
    orcid: 0000-0001-9783-8464
  - affiliation: Department of Oncology, Wayne State University, United States
 
date: 06 June 2025

bibliography: 
  - paper.bib

tags:
  - competing risks
  - recursive partitioning
  - hypothesis testing

## Summary

Recursive partitioning (RP) is a statistical procedure used to classify observations into subgroups that have a consistent outcome based on features of the observations. Starting with the complete set of observations, RP will split the input dataset into two using a binary (T/F) split based on 1 feature (e.g., Age>=70; Blood Type  = "A" | "B") that is optimal according to a statistical heuristic. Each resultant subgroup is similarly split until the potential splits are exhausted and/or other stopping conditions are met. Software has been developed to perform RP with hypothesis testing at each split for continuous and categorical outcome measures. However, competing risks (CR) outcomes, which are a special type of survival (time-to-event) analysis where multiple event types are possible to be observed, have no software implementation nor hypothesis testing included. The R package CRrpart includes RP for CR outcomes (along with survival outcomes) with statistical hypothesis testing at each split.

## Statement of need

In medicine, progression free survival (PFS) is an endpoint used to evaluate a patient's response to treatment. PFS is a composite endpoint, comprising of disease progression events and death events. A CR analysis [@FINE1999], which can determine the impact of the treatment on progression utilizing death events as competing events, may be more valuable to evaluate the efficacy of the treatment than an analysis of PFS.

RP (also known as classification and regression trees) was developed [@BREIMAN1984] as a method to create subgroups of observations with a consistent outcome based on features of the observations. In R [@RCORE2024], RP is implemented in the rpart package [@THERNEAU2023]. RP models have easy to interpret output (akin to a flowchart) and can handle larger dimension data than standard regression models. However, they are also prone to overfitting (particularly without a hypothesis test conducted at each split) and can be difficult to generalize to a new data set. CR endpoints cannot be analyzed by RP algorithms for survival outcomes because of the multiple possible event types.

Other implementations of the utilizing RP for a CR outcome have been developed (e.g., [@XU2016, @DAUDA2019]), but the code used to define the subgroups are not freely available and hypothesis testing (which was developed for continuous and categorical outcomes [@DYSON2018]) is not included in their methods. The R package CRrpart will recursively partition a dataset with a competing risks endpoint while conducting hypothesis testing at each split. 

## CRrpart package

The CRrpart package uses the CR implementation from Fine and Gray [@FINE1999] through the cmprsk R package [@GRAY2024]. The main function for performing the analysis for the package is `CRrpart`, with `plot`, `predict`, and `print` routines also included (see example below using a transplant (bmtcrr) dataset from the casebase package [@BHATNAGAR2022]). The naming conventions for CR response inputs (e.g., ftime, fstatus, cencode, failcode) established by the cmprsk package are followed for this package to ease integration. There are no methods to account for missing data for `CRrpart`, so only observations with complete data can be analzyed. The additional user inputted parameters for `CRrpart` include

- <u>n.splits</u> Maximum number of splits
- <u>minbucket</u> Minimum number of observations needed in each daughter node. Either minbucket or support needs to be specified.
- <u>support</u> Minimum percent of observations in the mother node needed in each daughter node. Either minbucket or support needs to be specified.
- <u>sig.level</u> Significance level 
- <u>iter</u> Number of iterations for the p-value calculation
- <u>p.adj</u> Adjust p-value for multiple testing?

Note that as the `crr` function in the cmprsk package [@GRAY2024] can analyze standard survival models (only 1 event type), `CRrpart` also will work for typical survival endpoints.

## P-value calculation 

To compute the p-value for each split based upon 2 times the difference in "pseudo" log likelihoods [@FINE1999] between the full and reduced models, multiple options were evaluated. Permutation testing is feasible but it cannot explore the entire space of possible subgroups given that it is restricted by the distributions of the input predictor variables. It is also computationally the most intensive. An alternative approach computed random splits for the observed data that met the minbucket/support criteria. Computing the log likelihood ratio (LLR) test statistic over all of the iterations provides a null distribution for the given split. This quantity mostly (but not always) follows a $\chi^2_1$ distribution. Therefore, to ensure that accurate p-values are computed for all input datasets, the reported p-value will be derived from the null distribution defined by the random splits.

## Multiple testing correction

If the p.adj parameter is TRUE, a multiple testing adjustment will be applied to the nominal p-value, based upon the cumulative distribution of the maximum LLR statistic $F(x)^n$, where $n$ is the number of input variables. The number of valid data splits, which was used previously [@DYSON2018] in a similar situation, is not employed as $n$ here because it results in too conservative of an p-value adjustment with the number of potential splits existing for continuous variables.

## Example

```
library(CRrpart)
g=CRrpart(minbucket=NULL,support=.2,ftime=bmtcrr$ftime,fstatus=bmtcrr$Status, cencode=0,failcode=1,x=bmtcrr[,c('Sex','D','Phase','Age','Source')], n.splits=500, sig.level=0.05,iter=2000,p.adj=FALSE)
print(g)
plot(g)
```

## Acknowledgements
Development of this package was supported by R01CA200864 and P30CA022453.
