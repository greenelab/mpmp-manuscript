---
title: Prediction of cancer mutation states using multiple data modalities reveals the utility and consistency of gene expression and DNA methylation
keywords:
- machine learning
- functional genomics
- gene expression
- epigenetics
- cancer
- pan-cancer
lang: en-US
date-meta: '2021-09-16'
author-meta:
- Jake Crawford
- Brock C. Christensen
- Maria Chikina
- Casey S. Greene
header-includes: |-
  <!--
  Manubot generated metadata rendered from header-includes-template.html.
  Suggest improvements at https://github.com/manubot/manubot/blob/main/manubot/process/header-includes-template.html
  -->
  <meta name="dc.format" content="text/html" />
  <meta name="dc.title" content="Prediction of cancer mutation states using multiple data modalities reveals the utility and consistency of gene expression and DNA methylation" />
  <meta name="citation_title" content="Prediction of cancer mutation states using multiple data modalities reveals the utility and consistency of gene expression and DNA methylation" />
  <meta property="og:title" content="Prediction of cancer mutation states using multiple data modalities reveals the utility and consistency of gene expression and DNA methylation" />
  <meta property="twitter:title" content="Prediction of cancer mutation states using multiple data modalities reveals the utility and consistency of gene expression and DNA methylation" />
  <meta name="dc.date" content="2021-09-16" />
  <meta name="citation_publication_date" content="2021-09-16" />
  <meta name="dc.language" content="en-US" />
  <meta name="citation_language" content="en-US" />
  <meta name="dc.relation.ispartof" content="Manubot" />
  <meta name="dc.publisher" content="Manubot" />
  <meta name="citation_journal_title" content="Manubot" />
  <meta name="citation_technical_report_institution" content="Manubot" />
  <meta name="citation_author" content="Jake Crawford" />
  <meta name="citation_author_institution" content="Genomics and Computational Biology Graduate Group, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA, USA" />
  <meta name="citation_author_orcid" content="0000-0001-6207-0782" />
  <meta name="twitter:creator" content="@jjc2718" />
  <meta name="citation_author" content="Brock C. Christensen" />
  <meta name="citation_author_institution" content="Department of Epidemiology, Geisel School of Medicine, Dartmouth College, Lebanon, NH, USA" />
  <meta name="citation_author_institution" content="Department of Molecular and Systems Biology, Geisel School of Medicine, Dartmouth College, Lebanon, NH, USA" />
  <meta name="citation_author_orcid" content="0000-0003-3022-426X" />
  <meta name="citation_author" content="Maria Chikina" />
  <meta name="citation_author_institution" content="Department of Computational and Systems Biology, School of Medicine, University of Pittsburgh, Pittsburgh, PA, USA" />
  <meta name="citation_author_orcid" content="0000-0003-2550-5403" />
  <meta name="citation_author" content="Casey S. Greene" />
  <meta name="citation_author_institution" content="Department of Biochemistry and Molecular Genetics, University of Colorado School of Medicine, Aurora, CO, USA" />
  <meta name="citation_author_institution" content="Center for Health AI, University of Colorado School of Medicine, Aurora, CO, USA" />
  <meta name="citation_author_orcid" content="0000-0001-8713-9213" />
  <meta name="twitter:creator" content="@GreeneScientist" />
  <link rel="canonical" href="https://greenelab.github.io/mpmp-manuscript/" />
  <meta property="og:url" content="https://greenelab.github.io/mpmp-manuscript/" />
  <meta property="twitter:url" content="https://greenelab.github.io/mpmp-manuscript/" />
  <meta name="citation_fulltext_html_url" content="https://greenelab.github.io/mpmp-manuscript/" />
  <meta name="citation_pdf_url" content="https://greenelab.github.io/mpmp-manuscript/manuscript.pdf" />
  <link rel="alternate" type="application/pdf" href="https://greenelab.github.io/mpmp-manuscript/manuscript.pdf" />
  <link rel="alternate" type="text/html" href="https://greenelab.github.io/mpmp-manuscript/v/1cddd2284c81c33e53f40346ee7ff1e101eb58f0/" />
  <meta name="manubot_html_url_versioned" content="https://greenelab.github.io/mpmp-manuscript/v/1cddd2284c81c33e53f40346ee7ff1e101eb58f0/" />
  <meta name="manubot_pdf_url_versioned" content="https://greenelab.github.io/mpmp-manuscript/v/1cddd2284c81c33e53f40346ee7ff1e101eb58f0/manuscript.pdf" />
  <meta property="og:type" content="article" />
  <meta property="twitter:card" content="summary_large_image" />
  <link rel="icon" type="image/png" sizes="192x192" href="https://manubot.org/favicon-192x192.png" />
  <link rel="mask-icon" href="https://manubot.org/safari-pinned-tab.svg" color="#ad1457" />
  <meta name="theme-color" content="#ad1457" />
  <!-- end Manubot generated metadata -->
bibliography:
- content/manual-references.json
manubot-output-bibliography: output/references.json
manubot-output-citekeys: output/citations.tsv
manubot-requests-cache-path: ci/cache/requests-cache
manubot-clear-requests-cache: false
...






<small><em>
This manuscript
([permalink](https://greenelab.github.io/mpmp-manuscript/v/1cddd2284c81c33e53f40346ee7ff1e101eb58f0/))
was automatically generated
from [greenelab/mpmp-manuscript@1cddd22](https://github.com/greenelab/mpmp-manuscript/tree/1cddd2284c81c33e53f40346ee7ff1e101eb58f0)
on September 16, 2021.
</em></small>

## Authors



+ **Jake Crawford**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0001-6207-0782](https://orcid.org/0000-0001-6207-0782)
    · ![GitHub icon](images/github.svg){.inline_icon}
    [jjc2718](https://github.com/jjc2718)
    · ![Twitter icon](images/twitter.svg){.inline_icon}
    [jjc2718](https://twitter.com/jjc2718)<br>
  <small>
     Genomics and Computational Biology Graduate Group, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA, USA
     · Funded by National Institutes of Health's National Cancer Institute (R01 CA237170); National Institutes of Health's National Human Genome Research Institute (R01 HG010067)
  </small>

+ **Brock C. Christensen**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0003-3022-426X](https://orcid.org/0000-0003-3022-426X)<br>
  <small>
     Department of Epidemiology, Geisel School of Medicine, Dartmouth College, Lebanon, NH, USA; Department of Molecular and Systems Biology, Geisel School of Medicine, Dartmouth College, Lebanon, NH, USA
     · Funded by National Institutes of Health's National Cancer Institute (R01 CA216265, R01 CA253976)
  </small>

+ **Maria Chikina**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0003-2550-5403](https://orcid.org/0000-0003-2550-5403)
    · ![GitHub icon](images/github.svg){.inline_icon}
    [mchikina](https://github.com/mchikina)<br>
  <small>
     Department of Computational and Systems Biology, School of Medicine, University of Pittsburgh, Pittsburgh, PA, USA
  </small>

+ **Casey S. Greene**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [0000-0001-8713-9213](https://orcid.org/0000-0001-8713-9213)
    · ![GitHub icon](images/github.svg){.inline_icon}
    [cgreene](https://github.com/cgreene)
    · ![Twitter icon](images/twitter.svg){.inline_icon}
    [GreeneScientist](https://twitter.com/GreeneScientist)<br>
  <small>
     Department of Biochemistry and Molecular Genetics, University of Colorado School of Medicine, Aurora, CO, USA; Center for Health AI, University of Colorado School of Medicine, Aurora, CO, USA
     · Funded by National Institutes of Health's National Cancer Institute (R01 CA237170); National Institutes of Health's National Human Genome Research Institute (R01 HG010067)
  </small>



## Abstract {.page_break_before}

In studies of cellular function in cancer, researchers are increasingly able to choose from many -omics assays as functional readouts.
Choosing the correct readout for a given study can be difficult, and which layer of cellular function is most suitable to capture the relevant signal may be unclear.
In this study, we consider prediction of cancer mutation status (presence or absence) from functional -omics data as a representative problem.
Since functional signatures of cancer mutation have been identified across many data types, this problem presents an opportunity to quantify and compare the ability of different -omics readouts to capture signals of dysregulation in cancer.
The TCGA Pan-Cancer Atlas contains genetic alteration data including somatic mutations and copy number variants (CNVs), as well as several -omics data types.
From TCGA, we focus on RNA sequencing, DNA methylation arrays, reverse phase protein arrays (RPPA), microRNA, and somatic mutational signatures as -omics readouts.

Across a collection of cancer-associated genetic alterations, RNA sequencing and DNA methylation were the most effective predictors of alteration state.
Surprisingly, we found that for most alterations, they were approximately equally effective predictors.
The target gene was the primary driver of performance, rather than the data type, and there was little difference between the top data types for the majority of genes.
We also found that combining data types into a single multi-omics model often provided little or no improvement in predictive ability over the best individual data type.
Based on our results, for the design of studies focused on the functional outcomes of cancer mutations, we recommend focusing on gene expression or DNA methylation as first-line readouts.



## Introduction

Although cancer can be initiated and driven by many different genetic alterations, these tend to converge on a limited number of pathways or signaling processes [@doi:10.1016/j.cell.2018.03.035].
As driver mutation status alone confers limited prognostic information, a comprehensive understanding of how diverse genetic alterations perturb central pathways is vital to precision medicine and biomarker identification efforts [@doi:10.7554/eLife.39217; @doi:10.1038/ncomms12096].
While many methods exist to distinguish driver mutations from passenger mutations based on genomic sequence characteristics [@doi:10.1073/pnas.1616440113; @doi:10.1038/s41467-019-11284-9; @doi:10.1371/journal.pcbi.1006658], until recently it has been a challenge to connect driver mutations to downstream changes in gene expression and cellular function within individual tumor samples.

The Cancer Genome Atlas (TCGA) Pan-Cancer Atlas provides uniformly processed, multi-platform -omics measurements across tens of thousands of samples from 33 cancer types [@doi:10.1038/ng.2764].
Enabled by this publicly available data, a growing body of work on linking the presence of driving genetic alterations in cancer to downstream gene expression changes has emerged.
Recent studies have considered Ras pathway alteration status in colorectal cancer [@doi:10.1158/1078-0432.CCR-13-1943], alteration status across many cancer types in Ras genes [@doi:10.1016/j.celrep.2018.03.046; @doi:10.1093/bib/bbaa258], _TP53_ [@doi:10.1016/j.celrep.2018.03.076], and _PIK3CA_ [@doi:10.1371/journal.pone.0241514], and alteration status across cancer types in frequently mutated genes [@doi:10.1186/s13059-020-02021-3].
More broadly, other groups have drawn on similar ideas to distinguish between the functional effects of different alterations in the same driver gene [@doi:10.1101/2020.06.02.128850], to link alterations with similar gene expression signatures within cancer types [@doi:10.1142/9789811215636_0031], and to identify trans-acting expression quantitative trait loci (trans-eQTLs) in germline genetic studies [@doi:10.1101/2020.05.07.083386].

These studies share a common thread: they each combine genomic (point mutation and copy number variation) data with transcriptomic (RNA sequencing) data within samples to interrogate the functional effects of genetic variation.
RNA sequencing is ubiquitous and cheap, and its experimental and computational methods are relatively mature, making it a vital tool for generating insight into cancer pathology [@doi:10.1038/nrg.2017.96].
Some driver mutations, however, are known to act indirectly on gene expression through varying mechanisms.
For example, oncogenic _IDH1_ and _IDH2_ mutations in glioma have been shown to interfere with histone demethylation, which results in increased DNA methylation and blocked cell differentiation [@doi:10.1016/j.ccr.2010.03.017; @doi:10.1093/jnci/djq497; @doi:10.1056/NEJMoa0808710; @doi:10.1038/nature10860].
Other genes implicated in aberrant DNA methylation in cancer include the TET family of genes [@doi:10.1016/j.tig.2014.07.005] and _SETD2_ [@doi:10.1101/cshperspect.a026468].
Certain driver mutations, such as those in DNA damage repair genes, may lead to detectable patterns of somatic mutation [@doi:10.1038/nrg3729].
Additionally, correlation between gene expression and protein abundance in cancer cell lines is limited, and proteomics data could correspond more directly to certain cancer phenotypes and pathway perturbations [@doi:10.1016/j.cell.2019.12.023].
In these contexts and others, integrating different data modalities or combining multiple data modalities could be more effective than relying solely on gene expression as a functional signature.

Here, we compare -omics data types profiled in the TCGA Pan-Cancer Atlas to evaluate use as a multivariate functional readout of genetic alterations in cancer.
We focus on gene expression (RNA sequencing data), DNA methylation (27K and 450K probe chips), reverse phase protein array (RPPA), microRNA expression, and mutational signatures data [@doi:10.1038/s41586-020-1943-3] as possible readouts.
Prior studies have identified univariate correlations of CpG site methylation [@doi:10.1371/journal.pcbi.1005840; @doi:10.1186/s12920-018-0425-z] and correlations of RPPA protein profiles [@doi:10.1186/s13073-018-0591-9] with the presence or absence of certain driver mutations.
Other relevant past work includes linking point mutations and copy number variants (CNVs) with changes in methylation and expression at individual genes [@doi:10.1093/bioinformatics/btr019; @doi:10.1093/bib/bbw037] and identifying functional modules that are perturbed by somatic mutations [@doi:10.1093/bioinformatics/btq182; @doi:10.1038/ncomms9554].
However, direct comparison among different data types for this application is lacking, particularly in the multivariate case where we consider changes to -omics-derived gene signatures rather than individual genes in isolation.

We select a collection of potential cancer drivers with varying functions and roles in cancer development [@doi:10.1126/science.1235122].
We use mutation status in these genes as labels to train classifiers, using each of the data types listed as training data, in a pan-cancer setting; we follow similar methods to the elastic net logistic regression approach described in Way et al. 2018 [@doi:10.1016/j.celrep.2018.03.046] and Way et al. 2020 [@doi:10.1186/s13059-020-02021-3].
We show that there is considerable predictive signal for many genes relative to a random baseline, and that gene expression and DNA methylation generally provide the best predictions of mutation state.
Surprisingly, we find that across a variety of target genes, gene expression and DNA methylation are approximately equally effective predictors; the target gene, rather than the data type, is the primary determinant of performance.
We observe similar results for pan-cancer survival prediction across the same data types, with gene expression, DNA methylation, and RPPA data providing the most predictive ability, and little separation between these top-performing data types.
In addition, we observe that combining data types into a single multi-omics model for mutation prediction provides little, if any, performance benefit over the most performant model using a single data type.
Our results will help to inform the design of future functional genomics studies in cancer, suggesting that for many strong drivers with clear functional signatures, gene expression and DNA methylation measurements provide similar information content.



## Methods

### Mutation data download and preprocessing

To generate binary mutated/non-mutated gene labels for our machine learning model, we used mutation calls for TCGA samples from MC3 [@doi:10.1016/j.cels.2018.03.002] and copy number threshold calls from GISTIC2.0 [@doi:10.1186/gb-2011-12-4-r41].
MC3 mutation calls were downloaded from the Genome Data Commons (GDC) of the National Cancer Institute, at [`https://gdc.cancer.gov/about-data/publications/pancanatlas`](https://gdc.cancer.gov/about-data/publications/pancanatlas).
Copy number threshold calls are from an older version of PanCanAtlas, and are available here: [`https://figshare.com/articles/dataset/TCGA_PanCanAtlas_Copy_Number_Data/6144122`](https://figshare.com/articles/dataset/TCGA_PanCanAtlas_Copy_Number_Data/6144122).
We removed hypermutated samples (defined as five or more standard deviations above the mean non-silent somatic mutation count) from our dataset to reduce the number of false positives (i.e., non-driver mutations).
In total, this resulted in 9,074 TCGA samples with mutation and copy number data.
Any sample with a non-silent somatic variant in the target gene was included in the positive set.
We also included copy number gains in the target gene for oncogenes, and copy number losses in the target gene for tumor suppressor genes, in the positive set; all remaining samples were considered negative for mutation in the target gene.

### Omics data download and preprocessing

RNA sequencing, 27K and 450K methylation array, and RPPA datasets for TCGA samples were all downloaded from GDC, at the same link provided above.
Mutational signatures information for TCGA samples with whole-exome sequencing data was downloaded from the International Cancer Genome Consortium (ICGC) data portal, at [`https://dcc.icgc.org/releases/PCAWG/mutational_signatures/Signatures_in_Samples/SP_Signatures_in_Samples`](https://dcc.icgc.org/releases/PCAWG/mutational_signatures/Signatures_in_Samples/SP_Signatures_in_Samples).
For our experiments, we used only the "single base signature" (SBS) mutational signatures, generated in [@doi:10.1038/s41586-020-1943-3].
We standardized (took z-scores of) each column of RNA sequencing and RPPA data; methylation data and mutational signatures data were left untransformed (beta values and mutation counts respectively), except in multi-omics experiments where all data types were standardized.
For the RNA sequencing dataset, we used only the top 8,000 gene features by mean absolute deviation as predictors in our models, except in multi-omics experiments where all 15,639 genes were used.

To remove missing values from the methylation datasets, we removed the 10 samples with the most missing values, then performed mean imputation for probes with 1 or 2 values missing.
All probes with missing values remaining after sample filtering and imputation were dropped from the analysis.
This left us with 20,040 CpG probes in the 27K methylation dataset, and 370,961 CpG probes in the 450K methylation dataset.
For experiments where "raw" methylation data was used, we used the top 100,000 probes in the 450K dataset by mean absolute deviation for computational efficiency, and we used all of the 20,040 probes in the 27K dataset.
For experiments where "compressed" methylation data was used, we used principal component analysis (PCA), as implemented in the `scikit-learn` Python library [@url:https://jmlr.csail.mit.edu/papers/v12/pedregosa11a.html], to extract the top 5,000 principal components from the methylation datasets.
We initially applied the beta-mixture quantile normalization (BMIQ) method [@doi:10.1093/bioinformatics/bts680] to correct for variability in signal intensity between type I and type II probes, but we observed that this had no effect on our results.
We report uncorrected results in the main paper for simplicity.

### Comparing data modalities

We made three main comparisons in this study, listed in the results section: one between different sets of genes using only expression data, one comparing expression and DNA methylation data types, and one comparing all data types.
This was mainly due to sample size limitations -- running a single comparison using all data types would force us to use only samples that are profiled for every data type, which would discard a large number of samples that lack profiling on only one or a few data types.
Thus, for each of the three comparisons, we used the intersection of TCGA samples having measurements for all of the datasets being compared in that experiment.
This resulted in three distinct sets of samples: 9,074 samples shared between {expression, mutation} data, 7,981 samples shared between {expression, mutation, 27K methylation, 450K methylation}, and 5,226 samples shared between {expression, mutation, 27K methylation, 450K methylation, RPPA, microRNA, mutational signatures}.
When we dropped samples between experiments as progressively more data types were added, we observed that the dropped samples had approximately the same cancer type proportions as the dataset as a whole.
In other words, samples that were profiled for one data type but not another did not tend to come exclusively from one or a few cancer types.
Exceptions included acute myeloid leukemia (LAML) which had no samples profiled in the RPPA data, and ovarian cancer (OV) which had only 8 samples with 450K methylation data.
More detailed information on cancer type proportions profiled for each data type is provided in Supplementary Figure {@fig:cancer_type_proportions} and Supplementary Table {@tbl:cancer_type_proportions}.

For each target gene, in order to ensure that the training dataset was reasonably balanced (i.e. that there would be enough mutated samples to train an effective classifier), we included only cancer types with at least 15 mutated samples and at least 5% mutated samples, which we refer to here as "valid" cancer types.
After applying these filters, the number of valid cancer types remaining for each gene varied based on the set of samples used: more data types resulted in fewer shared samples, and fewer samples generally meant fewer valid cancer types.
In some cases, this resulted in genes with no valid cancer types, which we dropped from the analysis.
Out of the 127 genes from the original cancer gene set described in Vogelstein et al. 2013 [@doi:10.1126/science.1235122], for the analysis using {expression, mutation} data we retained 85 target genes, for the {expression, mutation, 27k methylation, 450k methylation} analysis we retained 84 genes, and for the analysis using all data types we retained 75 genes.

We additionally explored mutation prediction from gene expression alone using 3 gene sets of equal size: the cancer-related genes from Vogelstein et al. 2013 described above, a set of frequently mutated genes in TCGA, and a set of random genes with mutations profiled by MC3.
To match the size of the Vogelstein et al. gene set, we took the 85 most frequently mutated genes in TCGA as quantified by MC3, all of which had at least one valid cancer type.
For the random gene set, we first filtered to the set of all genes with 2 or more valid cancer types by the above criteria, then sampled 85 of these genes uniformly at random.
Based on the results of the gene expression experiments, we used the Vogelstein et al. gene set for all subsequent experiments comparing -omics data types.

### Training classifiers to detect cancer mutations

We trained logistic regression classifiers to predict whether or not a given sample has a mutational event in a given target gene, using data from various -omics datasets as explanatory variables.
Our model is trained on -omics data ($X$) to predict mutation presence or absence ($y$) in a target gene.
To control for varying mutation burden per sample, and to adjust for potential cancer type-specific expression patterns, we included one-hot encoded cancer type and log~10~(sample mutation count) in the model as covariates.
Since our -omics datasets tend to have many dimensions and comparatively few samples, we used an elastic net penalty to prevent overfitting [@doi:10.1111/j.1467-9868.2005.00503.x], in line with the approach used in Way et al. 2018 [@doi:10.1016/j.celrep.2018.03.046] and Way et al. 2020 [@doi:10.1186/s13059-020-02021-3].
Elastic net logistic regression finds the feature weights $\hat{w} \in \mathbb{R}^{p}$ solving the following optimization problem:

$$\hat{w} = \text{argmin}_{w} \ \ell(X, y; w) + \alpha \lambda||w||_1 + \frac{1}{2}\alpha (1 - \lambda) ||w||_2$$

where $i \in \{1, \dots, n\}$ denotes a sample in the dataset, $X_i \in \mathbb{R}^{p}$ denotes features (omics measurements) from the given sample, $y_i \in \{0, 1\}$ denotes the label (mutation presence/absence) for the given sample, and $\ell(\cdot)$ denotes the negative log-likelihood of the observed data given a particular choice of feature weights, i.e.

$$\ell(X, y; w) = -\sum_{i=1}^{n} y_i \log\left(\frac{1}{1 + e^{-w^{\top}X_i}}\right) + (1 - y_i) \log\left(1 - \frac{1}{1 + e^{-w^{\top}X_i}}\right)$$

This optimization problem leaves two hyperparameters to select: $\alpha$ (controlling the tradeoff between the data log-likelihood and the penalty on large feature weight values), and $\lambda$ (controlling the tradeoff between the L1 penalty and L2 penalty on the weight values).
Although the elastic net optimization problem does not have a closed form solution, the loss function is convex, and iterative optimization algorithms are commonly used for finding reasonable solutions.
For fixed values of $\alpha$ and $\lambda$, we solved for $\hat{w}$ using stochastic gradient descent, as implemented in `scikit-learn`'s `SGDClassifier` method.

Given weight values $\hat{w}$, it is straightforward to predict the probability of a positive label (mutation in the target gene) $P(y^{*} = 1 \mid X^{*}; \hat{w})$ for a test sample $X^{*}$:

$$P(y^{*} = 1 \mid X^{*}; \hat{w}) = \frac{1}{1 + e^{-\hat{w}^{\top}X^{*}}}$$

and the probability of no mutation in the target gene, $P(y^{*} = 0 \mid X^{*}; \hat{w})$, is given by (1 - the above quantity).

For each target gene, we evaluated model performance using 2 replicates of 4-fold cross-validation, where train and test splits were stratified by cancer type and sample type.
That is, each training set/test set combination had equal proportions of each cancer type (BRCA, SKCM, COAD, etc) and each sample type (primary tumor, recurrent tumor, etc).
To choose the elastic net hyperparameters, we used 3-fold nested cross-validation, with a grid search over the following hyperparameter ranges: $\lambda$ = [0.0, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] and $\alpha$ = [0.0001, 0.001, 0.01, 0.1, 1, 10].
Using the grid search results, for each evaluation fold we selected the set of hyperparameters with the optimal area under the precision-recall curve (AUPR), averaged over the three inner folds.

### Evaluating mutation prediction classifiers

To quantify classification performance for a continuous or probabilistic output, such as that provided by logistic regression, the area under the receiver-operator curve (AUROC) [@doi:10.1016/j.patrec.2005.10.010] and the area under the precision-recall curve (AUPR) [@doi:10.1145/65943.65945] metrics are frequently used.
These metrics summarize performance across a variety of binary label thresholds, rather than requiring choice of a single threshold to determine positive or negative predictions.
In the main text, we report results using AUPR, summarized using average precision.
AUPR has been shown to distinguish between models more accurately than AUROC when there are few positively labeled samples [@doi:10.1371/journal.pone.0118432; @arxiv:2006.11278].
As an additional correction for imbalanced labels, in many of the results in the main text we report the difference in AUPR between a classifier fit to true mutation labels, and a classifier fit to data where the mutation labels are randomly shuffled.
In cases where mutation labels are highly imbalanced (very few mutated samples and many non-mutated samples), a classifier with shuffled labels may perform well simply by chance, e.g. by predicting the negative/non-mutated class for most samples.
To maintain the same label balance for the classifiers with shuffled labels as the classifiers with the true labels, we shuffled labels separately in the train and test sets for each cross-validation split.

Recall that for each target gene and each -omics dataset, we ran 2 replicates of 4-fold cross-validation, for a total of 8 performance results.
To make a statistical comparison between two models using these performance distributions, we used paired-sample _t_-tests, where performance measurements derived from the same cross-validation fold are considered paired measurements.
We used this approach to compare a model trained on true labels with a model trained on shuffled labels (addressing the question, "for the given gene using the given data type, can we predict mutation status better than random"), and to compare a model trained on data type A with a model trained on data type B (addressing the question, "for the given gene, can we make more effective mutation status predictions using data type A or data type B").

We corrected for multiple tests using a Benjamini-Hochberg false discovery rate correction.
For experiments where we chose a binary threshold for accepting/rejecting $H_0$ we set a conservative corrected threshold of $p = 0.001$; we were able to estimate the number of false positives by examining genes with better performance for shuffled mutation labels than true labels.
We chose this threshold to ensure that none of the observed false positive genes were considered significant, since we would never expect permuting labels to improve performance.
However, our results were not sensitive to the choice of this threshold, and we display cutoffs of $p = 0.05$ and $p = 0.01$ in many of our plots as well.

### Survival prediction using -omics datasets

As a complementary comparison to mutation prediction, we constructed predictors of patient survival using the clinical data available from the GDC, in the `TCGA-CDR-SupplementalTableS1.xlsx` file.
Following the methods described in [@doi:10.1101/2021.06.01.446243], as the clinical endpoint we used overall survival (OS), except in 9 cancer types with few deaths observed where we used progression-free intervals (PFI) as the clinical endpoint: BRCA, DLBC, LGG, PCPG, PRAD, READ, TGCT, THCA and THYM.
For prediction, we used elastic net Cox regression as implemented in the `scikit-survival` Python package [@url:https://jmlr.org/papers/v21/20-729.html], with patient age at diagnosis and log~10~(sample mutation count) included as covariates, as well as a one-hot encoded variable for cancer type in the pan-cancer case.
To ensure that the per-feature information content was comparable between data types, we preprocessed the -omics datasets using PCA and extracted the top $k$ principal components, for $k \in \{10, 100, 500, 1000, 5000\}$; in the case where the number of features in the original dataset was less than $k$ we used all available PCs (that is, we set $k = \min(p, k)$ where $p$ is the number of features in the unprocessed dataset).

To select hyperparameters for the elastic net Cox regression model, we performed a grid search over $\lambda$ = [0.0, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0], and used the default $\alpha$ selection procedure implemented in `scikit-survival` to select a range of $\alpha$ values based on the data.
This procedure begins by deriving the maximum $\alpha$ value as the smallest value for which all coefficients are 0 (call this $\alpha_{\text{max}}$), then it selects 100 possibilities for alpha spaced evenly on a log scale between $\alpha_{\text{max}}$ and $0.01 \cdot \alpha_{\text{max}}$.
For the survival models, we found that this data-driven procedure resulted in more consistent and stable model convergence than choosing a fixed set of alphas to search over as in the mutation prediction experiments.

We measured survival prediction performance using the censored concordance index (c-index) [@pubmed:8668867], which quantifies agreement between the order of survival time predictions and true outcomes for a held-out dataset: higher c-index values indicate more accurate survival prediction performance.
Similar to the mutation prediction experiments, we calculated c-index values on held-out subsets of the data for 2 replicates of 4-fold cross-validation, resulting in 8 performance measurements for each model.
For the results shown in the main text, we report the difference in c-index relative to a model where the survival endpoints are permuted, so values close to 0 suggest that the -omics features and other covariates are uninformative, and values greater than 0 suggest that they add predictive value.

### Multi-omics mutation prediction experiments

To predict mutation presence or absence in cancer genes using multiple data types simultaneously, we concatenated individual datasets into a large feature matrix, then used the same elastic net logistic regression method described previously.
For this task, we considered only the gene expression, 27K methylation, and 450K methylation datasets.
We used only these data types to limit the number of multi-omics combinations: the expression and methylation datasets resulted in the best overall performance across the single-omics experiments so we limited combinations to those datasets here.
In the main text, we report results using the top 5,000 principal components for each dataset.
In the supplement, we also report results using "raw" features: for gene expression we used all 15,639 genes available in our RNA sequencing dataset, for the 27K methylation dataset we used all 20,040 CpG probes, and for the 450K methylation dataset we used the top 5,000 principal components.

To construct the multi-omics models, we considered each of the pairwise combinations of the datasets listed above, as well as a combination of all 3 datasets.
When combining multiple datasets, we concatenated along the column axis, including covariates for cancer type and sample mutation burden as before.
For all multi-omics experiments, we used only the samples from TCGA with data for all three data types (i.e. the same 7,981 samples used in the single-omics experiments comparing expression and methylation data types).
We considered a limited subset of well-performing genes from the Vogelstein et al. gene set as target genes, including _EGFR_, _IDH1_, _KRAS_, _PIK3CA_, _SETD2_, and _TP53_.
We selected these genes because we have previously observed that they have good predictive performance, and they represent a combination of alterations that have strong gene expression signatures (_KRAS_, _EGFR_, _IDH1_, _TP53_) and strong DNA methylation signatures (_IDH1_, _SETD2_, _TP53_).

### Data and code availability

All analyses were implemented in the Python programming language and are available in the following GitHub repository: [`https://github.com/greenelab/mpmp`](https://github.com/greenelab/mpmp), under the open-source BSD 3-clause license.
Scripts to download large data files from GDC and other sources are located in the `00_download_data` directory.
Scripts to run experiments comparing data modalities used individually are located in the `02_classify_mutations` directory, and scripts to run multi-omics experiments are located in the `05_classify_mutations_multimodal` directory.
The Python environment was managed using `conda`, and directions for setting up the environment can be found in the `README.md` file. All analyses were run locally on a CPU.
This manuscript was written using Manubot [@doi:10.1371/journal.pcbi.1007128] and is available on GitHub at [`https://github.com/greenelab/mpmp-manuscript`](https://github.com/greenelab/mpmp-manuscript).


## Results

### Using diverse data modalities to predict cancer alterations

We collected five different data modalities from cancer samples in the TCGA Pan-Cancer Atlas, capturing five steps of cellular function that are perturbed by genetic alterations in cancer (Figure {@fig:overview}A).
These included gene expression (RNA-seq data), DNA methylation (27K and 450K Illumina BeadChip arrays), protein abundance (RPPA data), microRNA expression data, and patterns of somatic mutation (mutational signatures).
To link these diverse data modalities to changes in mutation status, we used elastic net logistic regression to predict the presence or absence of mutations in cancer genes, using these readouts as predictive features (Figure {@fig:overview}B).
We evaluated the resulting mutation status classifiers in a pan-cancer setting, preserving the proportions of each of the 33 cancer types in TCGA for 8 train/test splits (4 folds x 2 replicates) in each of 85 cancer genes (Figure {@fig:overview}C).

We sought to compare classifiers against a baseline where mutation labels are permuted (to identify genes whose mutation status correlates strongly with a functional signature in a given data type), and also to compare classifiers trained on true labels across different data types (to identify data types that are more or less predictive of mutations in a given gene).
To account for variation between dataset splits in making these comparisons, we treat classification metrics from the 8 train/test splits as performance distributions, which we compare using _t_-tests.
We summarize performance across all genes in our cancer gene set using a similar approach to a volcano plot, in which each point is a gene.
In our summary plots, the x-axis shows the magnitude of the change in the classification metric between conditions, and the y-axis shows the _p_-value for the associated _t_-test (Figure {@fig:overview}C).

![
**A.** Cancer mutations can perturb cellular function via a variety of cellular processes.
Arrows represent major potential paths of information flow from a somatic mutation in DNA to its resulting cell phenotype; circular arrow represents the ability of certain mutations (e.g. in DNA damage repair genes) to alter somatic mutation patterns.
Note that this does not reflect all possible relationships between cellular processes: for instance, changes in gene expression can lead to changes in somatic mutation rates.
**B.** Predicting presence/absence of somatic alterations in cancer from diverse data modalities.
In this study, we use functional readouts from TCGA as predictive features and the presence or absence of mutation in a given gene as labels.
This reverses the primary direction of information flow shown in Panel A.
**C.** Schematic of evaluation pipeline.
](images/figure_1.png){#fig:overview}

### Selection of cancer-related genes improves predictive signal

As a baseline, we evaluated prediction of mutation status from gene expression data across several different gene sets.
Past work has evaluated mutation prediction for the top 50 most mutated genes in TCGA [@doi:10.1186/s13059-020-02021-3], and we sought to extend this to a broader list of gene sets.
We compared a set of cancer-related genes (n=85) from Vogelstein et al. 2013 [@doi:10.1126/science.1235122] with an equal number of random genes (n=85) and an equal number of the most mutated genes in TCGA (n=85).
For all gene sets, we used only the set of TCGA samples for which both gene expression and somatic mutation data exists, resulting in a total of 9,074 samples from all 33 cancer types.
This set of samples was further filtered for each target gene to cancer types containing at least 15 mutated samples and at least 5% of samples mutated for that cancer type.
We then evaluated the performance for each target gene in each of the three gene sets.

Genes from the Vogelstein et al. set were more predictable than randomly chosen genes or those selected by total mutation count (Figure {@fig:expression_gene_sets}A).
In total, for a significance threshold of $\alpha = 0.001$, 47/85 genes (55.3%) in the Vogelstein et al. gene set are significantly predictable from gene expression data, compared to 12/85 genes (14.1%) in the random gene set and 15/85 genes (17.6%) in the most mutated gene set.
Of the 12 significantly predictable genes in the random gene set, 9 of them are also in the Vogelstein gene set (highlighted in red in Figure {@fig:expression_gene_sets}B), and of the 15 significantly predictable genes in the most mutated gene set, 9 of them are also in the Vogelstein gene set (highlighted in red in Figure {@fig:expression_gene_sets}C).
Additionally, many of the significant genes in the most mutated gene set are clustered close to the significance threshold (Figure {@fig:expression_gene_sets}C), while the significant genes in the Vogelstein et al. gene set tend to be further from the threshold (Figure {@fig:expression_gene_sets}D, higher AUPR differences and lower _p_-values).
These results suggest that selecting target genes for mutation prediction based on prior knowledge of their involvement in cancer pathways and processes, rather than randomly or based on mutation frequency alone, can improve predictive signal and identify more highly predictable mutations from gene expression data.

![
**A.** Overall distribution of performance across three gene sets, using gene expression (RNA-seq) data to predict mutations.
Each data point represents the mean cross-validated AUPR difference, compared with a baseline model trained on permuted mutation presence/absence labels, for one gene in the given gene set; notches show bootstrapped 95% confidence intervals.
"random" = 85 random genes, "most mutated" = 85 most mutated genes, "Vogelstein et al." = 85 cancer related genes from Vogelstein et al. 2013 gene set.
Significance stars indicate results of Bonferroni-corrected pairwise Wilcoxon tests: \*\*: _p_ < 0.01, \*\*\*: _p_ < 0.001, ns: not statistically significant for a cutoff of _p_ = 0.05.
**B, C, D.** Volcano-like plots showing mutation presence/absence predictive performance for each gene in each of the 3 gene sets.
The _x_-axis shows the difference in mean AUPR compared with a baseline model trained on permuted labels, and the _y_-axis shows _p_-values for a paired _t_-test comparing cross-validated AUPR values within folds.
Points (genes) highlighted in red are overlapping between the Vogelstein gene set and either the random or most mutated gene set.
](images/figure_2.png){#fig:expression_gene_sets width="90%"}

### Gene expression and DNA methylation have similar mutation prediction performance

We compared gene expression with DNA methylation as downstream readouts of the effects of cancer alterations.
In these analyses, we considered both the 27K probe and 450K probe methylation datasets generated for the TCGA Pan-Cancer Atlas.
We performed this comparison using the cancer-related gene set derived from Vogelstein et al. [@doi:10.1126/science.1235122].
We used samples that had data for each of the data types being compared, including somatic mutation data to generate mutation labels.
This process retained 7,981 samples in the intersection of the expression, 27K methylation, 450K methylation, and mutation datasets, which we used for subsequent analyses.
The most frequent missing data types were somatic mutation data (1,114 samples) and 450K methylation data (1,072 samples) (Figure {@fig:methylation}A).

For most genes, predictions are better than our baseline model where labels are permuted (most values greater than 0 in the box plots), suggesting that there is considerable predictive signal in both expression and methylation datasets across the Vogelstein et al. gene set (Figure {@fig:methylation}B).
Performance distributions are similar for expression and methylation, and aggregate performance is also similar for models using both 8,000 raw features (genes or CpG probes for expression and methylation respectively, selected using mean absolute deviation) or 5,000 principal components.
Both before and after filtering for genes that exceed the significance threshold, gene expression with raw gene features provides a significant performance improvement relative to the methylation datasets (Figure {@fig:methylation}B, C).
This provides evidence that gene expression measurements may be a slightly better predictor of mutation status than DNA methylation, although we observed no significant difference between aggregate performance using 27k methylation and 450k methylation.
Results were similar with PCA-compressed gene expression features or raw CpG probes as predictors (Supplementary Figure {@fig:me_compress_boxes}).

Considering each target gene in the Vogelstein gene set individually, we observed that 42/84 genes significantly outperformed the permuted baseline using gene expression data, as compared to 42/84 genes for 27K methylation and 39/84 genes for 450K methylation (Figure {@fig:methylation}D-F, more information in Supplementary Figure {@fig:methylation_heatmap}).
In most cases, these "well-predicted" genes that outperformed the permuted baseline tended to be similar between data types (Figure {@fig:methylation}D-F; genes in the top right of each plot).
For example, _TP53_, _BRAF_, and _PTEN_ appear in the top right of all 3 plots, suggesting that mutations in these genes have strong gene expression and DNA methylation signatures, and these signatures tend to be preserved across cancer types.

In addition to comparing mutation classifiers trained on different data types to the permuted baseline, we also compared classifiers trained on true labels directly to each other, using a similar methodology (Figure {@fig:methylation}G-H).
We observed that 6/100 genes were significantly more predictable from expression data than 27K methylation data, and 2/100 genes were significantly more predictable from expression data than 450K methylation data.
In both cases, 0/100 genes were significantly more predictable using the methylation data types.
For both comparisons (expression vs. 27K methylation and expression vs. 450K methylation), we observed that the majority of points were clustered around the origin, indicating that the data types appear to confer similar information about mutation status.
Additionally, many genes near the origin are significantly predictable vs. the shuffled baseline (labeled with an "X" in Figure {@fig:methylation}G-H), but equally predictable between data types (blue shading in Figure {@fig:methylation}G-H).
That is, in many cases, matching the gene being studied with the "correct" data modality seems to be unimportant: mutation status has a strong signature which can be extracted from both expression and DNA methylation data roughly equally.

We additionally compared pan-cancer survival prediction performance using principal components derived from each data type; results were comparable across the three data types (Figure {@fig:methylation}I).
The 450K methylation predictor appears to benefit slightly more from higher numbers of PCs than the other data modalities, with performance for gene expression and 27K methylation remaining fairly constant as more PCs were added to the model.
However, all three data types outperform the permuted-labels baseline, and confidence intervals between the best- and worst-performing data types overlap at each PC count, suggesting that similarly to mutation prediction, the three data types have comparable effectiveness for pan-cancer survival prediction.

![
**A.** Count of overlapping samples between gene expression, 27K methylation, 450K methylation, and somatic mutation data used from TCGA.
Only non-zero overlap counts are shown.
**B.** Predictive performance for genes in the Vogelstein et al. gene set, using each of the three data types as predictors.
The gene expression predictor uses the top 8000 gene features by mean absolute deviation, and the methylation predictors use the top 5000 principal components as predictive features.
Significance stars indicate results of Bonferroni-corrected pairwise Wilcoxon tests: \*\*: _p_ < 0.01, \*\*\*: _p_ < 0.001, ns: not statistically significant for a cutoff of _p_ = 0.05.
**C.** Predictive performance for genes where at least one of the considered data types predicts mutation labels significantly better than the permuted baseline.
**D-F.** Predictive performance for each gene in the Vogelstein et al. gene set, for each data type, compared with a baseline model trained on permuted labels.
**G-H.** Predictive performance for each gene in the Vogelstein et al. gene set, comparing gene expression directly to each methylation dataset (with classifiers trained on true labels).
**I.** Pan-cancer survival prediction performance, quantified using c-index relative to a label-permuted baseline on the _y_-axis, for gene expression, 27K methylation, and 450K methylation. The _x_-axis shows results for varying numbers of principal components included for each data type. Models also included covariates for patient age, sample mutation burden, and sample cancer type.

](images/figure_3.png){#fig:methylation width="90%"}

Focusing on several selected genes of interest, we observed that relative classifier performance varies by gene (Figure {@fig:methylation_genes}).
Past work has indicated that mutations in _TP53_ are highly predictable from gene expression data [@doi:10.1016/j.celrep.2018.03.076], and we observed that the methylation datasets provided similar predictive performance (Figure {@fig:methylation_genes}A).
Similarly, for _IDH1_ both expression and methylation features result in similar performance, consistent with IDH1's known role in regulating both DNA methylation and gene expression (Figure {@fig:methylation_genes}D) [@doi:10.1038/ng.3457].
Mutations in _KRAS_ and _ERBB2_ (_HER2_) were most predictable from gene expression data, and in both cases the methylation datasets significantly outperformed the baseline as well (Figure {@fig:methylation_genes}B and {@fig:methylation_genes}E).
Gene expression signatures of _ERBB2_ alterations are historically well-studied in breast cancer [@doi:10.1038/sj.onc.1207361], and samples with activating _ERBB2_ mutations have recently been shown to share sensitivities to some small-molecule inhibitors across cancer types [@doi:10.1016/j.ccell.2019.09.001].
These observations are consistent with the pan-cancer _ERBB2_ mutant-associated expression signature that we observed in this study.
_NF1_ mutations were also most predictable from gene expression data, although the gene expression-based _NF1_ mutation classifier did not significantly outperform the baseline with permuted labels at a cutoff of $\alpha = 0.001$ (Figure {@fig:methylation_genes}C).
_SETD2_ is an example of a gene that is more predictable from the methylation datasets than from gene expression, although gene expression with raw gene features significantly outperformed the permuted baseline as well (Figure {@fig:methylation_genes}F).
_SETD2_ is widely mutated across cancer types and affects H3K36 histone methylation most directly, but SETD2-mediated changes in H3K36 methylation have been linked to dysregulation of diverse cellular processes including DNA methylation and RNA splicing [@doi:10.1101/cshperspect.a026468; @doi:10.1007/s00018-017-2517-x].

![
Performance across varying PCA dimensions for specific genes of interest. Dotted lines represent results for "raw" features (8,000 gene features for gene expression data and 8,000 CpG probes for both methylation datasets, selected by largest mean absolute deviation). Error bars and shaded regions show bootstrapped 95% confidence intervals. Stars in boxes show statistical testing results compared with permuted baseline model: each box refers to the model using the number of PCA components it is over (far right box = models with raw features). \*\*: _p_ < 0.01, \*\*\*: _p_ < 0.001, no stars: not statistically significant for a cutoff of _p_ = 0.05.
](images/figure_4.png){#fig:methylation_genes width="90%"}

### Comparing six different readouts favors expression and DNA methylation

Next, we expanded our comparison to all five functional data modalities (six total readouts, since there are two DNA methylation platforms) available in the TCGA Pan-Cancer Atlas.
As with previous experiments, we limited our comparison to the set of samples profiled for each readout, resulting in 5,226 samples with data for all readouts.
The data types with the most missing samples were RPPA data (2,215 samples that were missing RPPA data) and 450K methylation (630 samples that were missing 450K methylation data) (Figure {@fig:all_data}A).
Summarized over all genes in the Vogelstein et al. dataset, we observed that gene expression and both methylation datasets tended to produce similar quality predictions, and these were significantly better than the remaining data types (Figure {@fig:all_data}B).
For the set of genes having at least one significant predictor (i.e. "well-predicted" genes), median performance using gene expression was slightly higher than for the methylation data types, although this difference was not statistically significant (Figure {@fig:all_data}C).

On the individual gene level, mutations in 26/75 genes were significantly predictable from RPPA data relative to the permuted baseline, compared to 23/75 genes from microRNA data and 15/75 genes from mutational signatures data (Figure {@fig:all_data}D-F).
For the remaining data types on this smaller set of samples, 35/75 genes outperformed the baseline for gene expression data, 34/75 for 27k methylation, and 31/75 for 450k methylation.
Compared to the methylation experiments (Figure {@fig:methylation}), we observed slightly fewer "well-predicted" genes for the expression and methylation datasets here (likely due to the considerably smaller sample size) but relative performance was comparable (Supplementary Figure {@fig:all_volcano_me}).
Direct comparisons between each added data type and gene expression data showed that RPPA, microRNA and mutational signatures data generally provide similar or worse performance than the remaining data types (Figure {@fig:all_data}G-I).

Performance using RPPA data (Figure {@fig:all_data}G) is notable because of its drastically smaller dimensionality than the other data types (190 proteins, compared to thousands of features for the expression and methylation data types).
This suggests that each protein abundance measurement provides a high information content, although this is by design as the antibody probes used for the TCGA analysis were selected to cover established cancer-related pathways [@doi:10.1038/nmeth.2650].
Mutations that are more predictable using RPPA data include _PIK3R1_ and _MAP2K1_ (Figure {@fig:all_data}G), although neither classifier significantly outperforms the permuted baseline.
Both genes are kinases involved in phosphorylation-mediated signal transduction.
The ability of RPPA technology to quantify protein phosphorylation status may thus provide an advantage in identifying mutations in these genes, relative to the other data types we used that cannot directly measure protein phosphorylation.
Similarly, the scope of the features captured by the mutational signatures data we used is limited to single-base substitution signatures; a broader spectrum of possible signatures is described in previous work [@doi:10.1038/s41586-020-1943-3; @doi:10.1038/s41586-019-1913-9] including doublet-substitution signatures, small indel signatures, and signatures of structural variation, but these were not readily available for the TCGA exome sequencing data.
The relatively poor predictive ability of mutational signatures likely stems from a combination of biological and technical factors, as there is no reason to expect that changes in somatic mutation patterns would be directly caused by most cancer driver mutations (excepting, for instance, mutations in genes involved in DNA damage repair and similar processes).

As in the expression/methylation comparison, we compared pan-cancer survival prediction performance between all six readouts, using the top principal components derived from each data type to ensure comparable information content (Figure {@fig:all_data}J).
In this case, RPPA data resulted in comparable prediction to the expression and methylation-derived feature sets, performing slightly better at low numbers of PCs possibly due to its higher information content over a smaller dimensionality.
The microRNA and mutational signatures datasets generally proved ineffective for pan-cancer survival prediction (values near zero in Figure {@fig:all_data}J), although in survival analyses for individual cancer types these data modalities were sometimes effective predictors (Supplementary Figure {@fig:survival_cancer_types}).

![
**A.** Overlap of TCGA samples between data types used in mutation prediction comparisons. Only overlaps with more than 100 samples are shown.
**B.** Overall distribution of performance per data type across 75 genes from Vogelstein et al. gene set.
Each data point represents mean cross-validated AUPR difference, compared with a baseline model trained on permuted labels, for one gene; notches show bootstrapped 95% confidence intervals.
Significance stars indicate results of Bonferroni-corrected pairwise Wilcoxon tests: \*\*: _p_ < 0.01, \*\*\*: _p_ < 0.001, ns: not statistically significant for a cutoff of _p_ = 0.05.
All pairwise tests were run (and corrected for) but only neighboring test results are shown.
**C.** Overall performance distribution per data type for genes where the permuted baseline model is significantly outperformed for one or more data types, resulting in a total of 39 genes.
**D, E, F.** Volcano-like plots showing predictive performance for each gene in the Vogelstein et al. gene set, in each of the added data types (RPPA, microRNA, mutational signatures). The _x_-axis shows the difference in mean AUPR compared with a baseline model trained on permuted labels, and the _y_-axis shows _p_-values for a paired _t_-test comparing cross-validated AUPR values within folds.
**G, H, I.** Volcano-like plots comparing predictive performance between data types for each gene in the Vogelstein et al. gene set. The _x_-axis shows the difference in mean AUPR between gene expression and another data type (positive values = better mean performance using gene expression features), and the _y_-axis shows _p_-values for a paired _t_-test comparing cross-validated AUPR values within folds.
**J.** Pan-cancer survival prediction performance, quantified using c-index relative to a label-permuted baseline on the _y_-axis, for all data types. The _x_-axis shows results for varying numbers of principal components included for each data type. Models also included covariates for patient age, sample mutation burden, and sample cancer type.
](images/figure_5.png){#fig:all_data width="90%"}

When we constructed a heatmap depicting predictive performance for each gene across data types, we found that very few genes tended to be well-predicted exclusively by one or two data types (Figure {@fig:heatmap}).
Of the 39 genes that are well-predicted using at least one data type (blue or red highlighted boxes in Figure {@fig:heatmap}), only three of them are well-predicted exclusively by a single data type, meaning that mutations in the other 37 genes can be predicted effectively using at least two different data sources.
This supports our observation that choosing the "correct" data modality is often unimportant for driver genes with strong functional signatures.
Notable exceptions included _NF1_ (only well-predicted using gene expression data), _SETD2_ (only well-predicted using the two methylation datasets), and _TSC1_ (only well-predicted using gene expression data).
Gene expression provided the best performance in 32/39 genes with at least one significant data type (red highlighted boxes in Figure {@fig:heatmap}), but only 2 of those 32 genes did not have any other significantly predictive data types (_NF1_ and _TSC1_); in the other 23 genes one or more non-expression data types also outperformed the permuted baseline.

![
Heatmap displaying predictive performance for mutations in each of the 75 genes from the Vogelstein et al. gene set, across all 6 TCGA data modalities. Each cell quantifies performance for a target gene, using predictive features derived from a particular data type. Blue highlights indicate that the given data type provides significantly better predictions than the permuted baseline for the given gene; red highlights indicate the same and additionally that the given data type provides statistically equivalent performance to the data type with the best average performance (determined by pairwise _t_-tests with FDR correction).
](images/figure_6.png){#fig:heatmap width="90%"}

### Simple multi-omics baseline provides no performance benefit

We also trained "multi-omics" classifiers to predict mutations in six well-studied and widely mutated driver genes across various cancer types: _EGFR_, _IDH1_, _KRAS_, _PIK3CA_, _SETD2_, and _TP53_.
Each of these genes is well-predicted from several data types in our earlier experiments (Figure {@fig:heatmap}), consistent with having strong pan-cancer driver effects.
For the multi-omics classifiers, we considered all pairwise combinations of the top three performing individual data types (gene expression, 27K methylation, and 450K methylation), in addition to a model using all three data types.
We trained a classifier for multiple data types by concatenating features from the individual data types, then fitting the same elastic net logistic regression model as we used for the single-omics models.
Here, we show results using the top 5000 principal components from each data type as predictive features, to ensure that feature count and scale is comparable among data types; results for raw features are shown in Supplementary Figure {@fig:multi_omics_raw_feats}.

For each of the six target genes, we observed comparable performance between the best single-omics classifier (blue boxes in Figure {@fig:multi_omics}A) and the best multi-omics classifier (orange boxes in Figure {@fig:multi_omics}A).
Across all classifiers and data types, we found varied patterns based on the target gene.
For _IDH1_ and _TP53_ performance is relatively consistent regardless of data type(s), suggesting that baseline performance is high and there is little room for improvement as data is added (Figure {@fig:multi_omics}C, G).
For _EGFR_, _KRAS_, and _PIK3CA_, combining gene expression with methylation data results in statistically equivalent performance to gene expression alone; classifiers trained only on methylation data perform worse (Figure {@fig:multi_omics}B, D, E).
The best classifiers for _SETD2_ use methylation data alone; results are comparable whether 27K methylation or 27K+450K methylation features are used (Figure {@fig:multi_omics}F).
Overall, we saw that combining data types in a relatively simple manner (i.e. by concatenating features from each individual data type) provided little or no improvement in predictive ability over the best individual data type.
This supports our earlier observations of the redundancy of gene expression and methylation data as functional readouts, since our multi-omics classifiers are not in general able to extract gains in predictive performance as more data types are added for this set of cancer drivers.

![
**A.** Comparing the best-performing model (i.e. highest mean AUPR relative to permuted baseline) trained on a single data type against the best "multi-omics" model for each target gene. None of the differences between single-omics and multi-omics models were statistically significant using paired-sample Wilcoxon tests across cross-validation folds, for a threshold of 0.05.
**B-G.** Classifier performance, relative to baseline with permuted labels, for mutation prediction models trained on various combinations of data types. Each panel shows performance for one of the six target genes; box plots show performance distribution over 8 evaluation sets (4 cross-validation folds x 2 replicates).
](images/figure_7.png){#fig:multi_omics width="90%"}


## Discussion

We carried out a large-scale comparison of data types in the TCGA Pan-Cancer Atlas as functional readouts of genetic alterations in cancer, integrating results across cancer types and across driver genes.
Overall, we found that gene expression and DNA methylation capture signatures of mutation state equally effectively in most cases, and that other data types (RPPA, microRNA, mutational signatures) were less effective at predicting mutation presence or absence.
For pan-cancer survival prediction, we found that RPPA data had comparable effectiveness to the expression and methylation data types, and microRNA and mutational signatures datasets were ineffective.
Our multi-omics modeling experiment indicated that the mutation state information captured by gene expression and DNA methylation is highly redundant, as added data types resulted in no gain or modest gains in classifier performance.

Comparing mutation status prediction using raw and PCA compressed expression and DNA methylation data, we observed that feature extraction using PCA provided no benefit compared to using raw gene or CpG probe features.
Other studies using DNA methylation array data have found that nonlinear dimension reduction methods, such as variational autoencoders and capsule networks, can be effective for extracting predictive features [@doi:10.1186/s12859-020-3443-8; @pubmed:34417465].
The latter approach is especially interesting because capsule networks and "capsule-like methods" can be constrained to extract features that align with known biology (i.e. that correspond to known disease pathways or CpG site annotations).
This can improve model interpretability as well as predictive performance.
Similar methods have been applied to extract biologically informed features from gene expression data (see, for instance, [@doi:10.1016/j.ccell.2020.09.014; @doi:10.1101/2021.05.25.445604]).
A more comprehensive study of dimension reduction methods in the context of mutation prediction, including the features selected by these methods and their biological relevance and interpretation, would be a beneficial area of future work.

In contrast to many other studies demonstrating the benefits of integrating multiple -omics data types for various cancer-related prediction problems [@doi:10.1038/nrc3721; @doi:10.1016/j.jbi.2015.05.019; @doi:10.1186/s13040-018-0184-6; @doi:10.1093/bioinformatics/btz318; @doi:10.1371/journal.pcbi.1008878], we found that combining multiple data types to predict mutation status was generally not effective for this problem.
The method we used to integrate different data types by concatenating feature sets is sometimes referred to as "early" data integration (discussed in more detail in [@doi:10.1098/rsif.2015.0571] and [@doi:10.1016/j.inffus.2018.09.012]).
It is possible that more sophisticated data integration methods, such as "intermediate" integration methods that learn a set of features jointly across datasets, would produce improved predictions.
We do not interpret our results as concrete evidence that multi-omics integration is not effective for this problem; rather, we see them as an indication that this is a challenging data integration problem for which further investigation is needed.
We also present this problem as a set of benchmark tasks on which multi-omics integration methods can be evaluated.
In addition to the methodological questions, the issue of data integration also has implications for the underlying biology: a more nuanced understanding of when different data readouts provide redundant information, and when they can contribute unique information about cancer pathology and development, could have many translational applications.

One limitation of the current study is that, for the mutation prediction problem, we only evaluated classifiers that were trained on pan-cancer data.
Considering every possible combination of target gene and TCGA cancer type (85 target genes x 33 cancer types x 6 data types) would have drastically increased the computational load and presented a large multiple testing burden.
Alternatively, choosing only a subset of gene/cancer type combinations to study would have biased our results toward known driver gene/cancer type relationships, which we aimed to avoid.
In future work it would be interesting to identify classifiers that perform well in a certain cancer type but not in the pan-cancer context, and to compare these instances across different cancer types.
As a motivating example, other studies have shown that activating mutations in Ras isoforms (_HRAS_, _KRAS_, _NRAS_) tend to have similar effects to one another in thyroid cancer, producing similar gene expression signatures [@doi:10.1142/9789811215636_0031].
In multiple myeloma, however, activating _KRAS_ and _NRAS_ mutations produce distinct expression signatures, necessitating separate classifiers [@doi:10.1182/bloodadvances.2019000303].
A high-throughput computational pipeline to identify cases where functional signatures of a particular cancer driver are either concordant or discordant between cancer types could identify opportunities for context-specific protein function prediction, improve biomarker identification, and suggest cases where drugs targeting specific alterations might produce discordant results in different cancer types.

As with any study relying on observational, cross-sectional data such as the TCGA Pan-Cancer Atlas, the conclusions that we can draw are limited by the data.
In particular, for any of our "well-predicted" genes (i.e. genes that, when mutated, have strong signatures in one or more data types), we cannot definitively distinguish correlation from causation.
To directly assess the effects of particular mutations on various data modalities, some studies use cell line data from sources such as the Cancer Cell Line Encyclopedia (CCLE) [@doi:10.1038/s41586-019-1186-3].
While this approach could help to isolate the causal effect of a given mutation on a given cell line, cell lines are sometimes an imperfect match for the cancers they are derived from [@doi:10.1038/s41467-019-11415-2].
We are also limited in that we cannot assign timing or clonal status to mutations, or fully characterize intratumor heterogeneity, with certainty from the bulk sequencing data generated by TCGA (although some features of tumor mutational processes over time can be estimated from bulk data, e.g. [@doi:10.1126/scitranslmed.aaa1408]).
As methods for generating large longitudinal datasets at single-cell resolution mature and scale, we will need to revise the way we think about cellular function and dysregulation in cancer cells, as dynamic and adaptive processes rather than a single representative snapshot of a tumor.

Based on our results, for studies focused on the functional consequences of cancer mutations, we recommend that researchers designing large-scale studies aiming to understand the regulatory state of cancers prioritize gene expression or DNA methylation as downstream readouts.
On balance, prediction of mutation status was slightly better using gene expression data.
However, the finding that DNA methylation profiles contain much of the same information will be useful for some study designs, given the varying stability of mRNA across storage methods and times.
Results using RPPA measurements as predictive features were also promising, especially considering the substantially lower dimensionality of the RPPA dataset compared to other data types.
Future technology advances, in both quality and quantity of data, are likely to improve our understanding of the full picture of functional consequences of mutations in cancer cells.


## Acknowledgements

We would like to thank Alexandra Lee, Ben Heil, Milton Pividori, and Natalie Davidson for reviewing the software associated with this work and providing insightful feedback.
Figure 1 (the schematic of the background and evaluation pipeline) was created using BioRender.com.



## Supplementary Information

A version of the main paper figures using the area under the receiver-operator curve (AUROC) metric rather than AUPR is available at [`https://doi.org/10.6084/m9.figshare.14919729`](https://doi.org/10.6084/m9.figshare.14919729).

![
Proportion of samples from each TCGA cancer type that are "dropped" as more data types are added to our analyses. We started with gene expression data, and for each added data type, we took the intersection of samples that were profiled for that data type and the previous data types, dropping all samples that were missing 1 or more data types. Overall, at each step, the proportions of "dropped" samples appear to be fairly evenly spread between cancer types, showing that in general we are not disproportionately losing one or several cancer types as more data modalities are added to our analyses.
](images/supp_figure_8.png){#fig:cancer_type_proportions width=80%}

| cancer type | mutation | me_27k | me_450k | mut_sigs | mirna | rppa | base |
|-------------|----------|--------|---------|----------|-------|------|------|
| ACC | 3 | 0 | 0 | 0 | 0 | 32 | 44 |
| BLCA | 29 | 0 | 0 | 20 | 3 | 65 | 310 |
| BRCA | 237 | 16 | 291 | 44 | 9 | 109 | 512 |
| CESC | 38 | 0 | 0 | 18 | 0 | 120 | 134 |
| CHOL | 9 | 0 | 3 | 0 | 0 | 6 | 27 |
| COAD | 154 | 0 | 77 | 15 | 0 | 33 | 216 |
| DLBC | 11 | 0 | 0 | 6 | 0 | 10 | 21 |
| ESCA | 27 | 0 | 0 | 0 | 0 | 53 | 116 |
| GBM | 48 | 0 | 73 | 12 | 39 | 0 | 0 |
| HNSC | 79 | 0 | 0 | 39 | 3 | 267 | 178 |
| KICH | 26 | 0 | 0 | 15 | 0 | 2 | 48 |
| KIRC | 254 | 0 | 97 | 16 | 0 | 18 | 221 |
| KIRP | 52 | 0 | 13 | 17 | 0 | 53 | 188 |
| LAML | 58 | 1 | 0 | 0 | 0 | 114 | 0 |
| LGG | 23 | 0 | 1 | 8 | 4 | 85 | 409 |
| LIHC | 75 | 0 | 1 | 46 | 2 | 179 | 120 |
| LUAD | 74 | 1 | 56 | 20 | 1 | 125 | 299 |
| LUSC | 89 | 1 | 113 | 21 | 4 | 94 | 231 |
| MESO | 5 | 0 | 0 | 0 | 0 | 23 | 59 |
| OV | 132 | 2 | 166 | 0 | 0 | 0 | 8 |
| PAAD | 31 | 0 | 2 | 2 | 0 | 51 | 97 |
| PCPG | 26 | 0 | 0 | 1 | 0 | 83 | 77 |
| PRAD | 71 | 0 | 0 | 10 | 4 | 134 | 331 |
| READ | 53 | 0 | 29 | 1 | 0 | 20 | 68 |
| SARC | 36 | 0 | 0 | 28 | 1 | 33 | 167 |
| SKCM | 111 | 0 | 0 | 52 | 14 | 92 | 205 |
| STAD | 67 | 0 | 40 | 25 | 0 | 67 | 251 |
| TGCT | 12 | 0 | 0 | 0 | 0 | 30 | 114 |
| THCA | 92 | 0 | 1 | 16 | 1 | 124 | 338 |
| THYM | 3 | 0 | 0 | 0 | 0 | 32 | 87 |
| UCEC | 60 | 0 | 109 | 22 | 0 | 84 | 292 |
| UCS | 1 | 0 | 0 | 0 | 1 | 9 | 46 |
| UVM | 0 | 0 | 0 | 0 | 0 | 68 | 12 |

Table: Number of samples from each TCGA cancer type that are "dropped" as more data types are added to the analysis. The "base" column indicates the number of samples that are present per cancer type in the final intersection of all data types (i.e. each sample counted in the last column has data for each of the 7 data types, including gene expression (not listed here) and somatic mutations). {#tbl:cancer_type_proportions}

![
Heatmap displaying predictive performance for mutations in each of the 84 genes from the Vogelstein et al. gene set, across gene expression and the two DNA methylation arrays. Each cell quantifies performance for a target gene, using predictive features derived from a particular data type. Blue highlights indicate that the given data type provides significantly better predictions than the permuted baseline for the given gene; red highlights indicate the same and additionally that the given data type provides statistically equivalent performance to the data type with the best average performance (determined by pairwise _t_-tests with FDR correction).
](images/supp_figure_9.png){#fig:methylation_heatmap width=90%}

![
Volcano-like plots showing predictive performance for each gene in the Vogelstein et al. gene set for expression and DNA methylation, on the sample set used for the "all data types" experiments. The first row shows performance relative to the permuted baseline, and the second row shows direct comparisons between data types. The _x_-axis shows the difference in mean AUPR compared with a baseline model trained on permuted labels, and the _y_-axis shows _p_-values for a paired _t_-test comparing cross-validated AUPR values within folds.
](images/supp_figure_10.png){#fig:all_volcano_me width=90%}

![
Predictive performance for genes in the Vogelstein et al. gene set, using each of the three data types as predictors.
The _x_-axis shows the number of PCA components used as features, "raw" = no PCA compression.
](images/supp_figure_11.png){#fig:me_compress_boxes width=90%}

![
Survival prediction performance across increasing numbers of PCA components, for individual cancer types.
The top row shows the five cancer types with the highest average coefficient of variation across data types, for expression, 27K methylation, and 450K methylation.
The bottom row shows the five cancer types with the highest average coefficient of variation across data types, for the comparison using all 6 data types.
Each point shows the mean concordance index, relative to baseline predictors with permuted labels, across 8 performance measurements (2 cross-validation replicates x 4 folds); error bars show bootstrapped 95% confidence intervals for the mean.
](images/supp_figure_12.png){#fig:survival_cancer_types width=90%}

![
Top plot: comparing the best-performing model (i.e. highest mean AUPR relative to permuted baseline) trained on a single data type against the best "multi-omics" model for each target gene, using raw (not PCA compressed) features. For feature parity between data types, the top 20,000 features by mean absolute deviation were used for each feature type. None of the differences between single-omics and multi-omics models were statistically significant using paired-sample Wilcoxon tests across cross-validation folds, for a threshold of 0.05.
Bottom plots: classifier performance, relative to baseline with permuted labels, for individual genes. Each panel shows performance for one of the six target genes; box plots show performance distribution over 8 evaluation sets (4 cross-validation folds x 2 replicates).
](images/supp_figure_13.png){#fig:multi_omics_raw_feats width=90%}


## References {.page_break_before}

<!-- Explicitly insert bibliography here -->
<div id="refs"></div>
