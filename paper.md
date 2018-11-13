---
title: 'iRF: Extracting interactions from random forests'
tags:
  - R
  - Random Forests
  - Interpretable machine learning
authors:
  - name: Sumanta Basu
    affiliation: "&#8727;, 1, 2"
  - name: Karl Kumbier
    affiliation: "&#8727;, 3"
  - name: James B. Brown
    affiliation: "3, 4, 5"
  - name: Bin Yu
    affiliation: "3, 6"
affiliations:
 - name: Department of Biological Statistics and Computational Biology, Cornell University
   index: 1
 - name: Department of Statistical Science, Cornell University
   index: 2
 - name: Statistics Department, University of California, Berkeley
   index: 3
 - name: Centre for Computational Biology, School of Biosciences, University of Birmingham
   index: 4
 - name: Molecular Ecosystems Biology Department, Biosciences Area, Lawrence Berkeley National Laboratory 
   index: 5
 - name: Department of Electrical Engineering and Computer Sciences, University of California, Berkeley 
   index: 6
 - name: Denotes equal contribution
   index: &#8727;;
date: 2 October 2018
bibliography: paper.bib
---

# Summary
Random forests [@breiman2001random] are a popular class of supervised learning
models that have demonstrated impressive empirical success across a wide variety
of problems. The predictive accuracy of random forests stems from their ability
to learn high-order, non-linear interactions in large datasets. Although
approaches exist for evaluating the importance of individual features in a
fitted random forest, identifying interactions that drive predictive accuracy
remains a challenge. This challenge is in large part due to the enourmous number
of interactions that must be considered (i.e. there are $O(p^s)$ possible
interactions of size $s$ among $p$ features) and the instability of random
forest decision paths.

The iterative Random Forest algorithm (iRF), and corresponding `iRF` R package,
take a step towards addressing these issues with a computationally tractable
approach to search for important interactions in a fitted random forest
[@basu2018iterative]. Our algorithm grows a series of feature weighted random
forests [@amaratunga2008enriched] to perform soft regularization on the model
based on predictive features. We then search for prevalent interactions in the
fitted random forest using a generalization of random intersection trees
[@shah2014random].  Finally, we assess the stability of recovered interactions
by repeating this search across random forests trained on bootstrap samples of
the data. The `iRF` R package combines these steps into a single workflow. It is
based on the source codes from the R packages `randomForest`
[@liaw2002classification] and `FSInteract` [@shah2014random]. A detailed
vignette is available on
[github](https://cdn.rawgit.com/sumbose/iRF/master/vignettes/vignette2.html).

# Acknowledgements
This research was supported in part by grants NHGRI U01HG007031, ARO
W911NF1710005, ONR N00014-16-1-2664, DOE DE-AC02-05CH11231, NHGRI R00 HG006698,
DOE (SBIR/STTR) Award DE-SC0017069, DOE DE-AC02-05CH11231, and NSF DMS-1613002.
We thank the Center for Science of Information (CSoI), a US NSF Science and
Technology Center, under grant agreement CCF-0939370. Research reported in this
publication was supported by the National Library Of Medicine of the NIH under
Award Number T32LM012417. The content is solely the responsibility of the
authors and does not necessarily represent the official views of the NIH. BY
acknowledges support from the Miller Institute for her Miller Professorship in
2016-2017. SB acknowledges the support of UC Berkeley and LBNL, where he
conducted most of his work on this paper as a postdoc. We thank P. Bickel and S.
Shrotriya for helpful discussions and comments.

# References
