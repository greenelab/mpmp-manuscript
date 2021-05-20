---
title: Manuscript Title
keywords:
- markdown
- publishing
- manubot
lang: en-US
date-meta: '2021-05-20'
author-meta:
- John Doe
- Jane Roe
header-includes: |-
  <!--
  Manubot generated metadata rendered from header-includes-template.html.
  Suggest improvements at https://github.com/manubot/manubot/blob/main/manubot/process/header-includes-template.html
  -->
  <meta name="dc.format" content="text/html" />
  <meta name="dc.title" content="Manuscript Title" />
  <meta name="citation_title" content="Manuscript Title" />
  <meta property="og:title" content="Manuscript Title" />
  <meta property="twitter:title" content="Manuscript Title" />
  <meta name="dc.date" content="2021-05-20" />
  <meta name="citation_publication_date" content="2021-05-20" />
  <meta name="dc.language" content="en-US" />
  <meta name="citation_language" content="en-US" />
  <meta name="dc.relation.ispartof" content="Manubot" />
  <meta name="dc.publisher" content="Manubot" />
  <meta name="citation_journal_title" content="Manubot" />
  <meta name="citation_technical_report_institution" content="Manubot" />
  <meta name="citation_author" content="John Doe" />
  <meta name="citation_author_institution" content="Department of Something, University of Whatever" />
  <meta name="citation_author_orcid" content="XXXX-XXXX-XXXX-XXXX" />
  <meta name="twitter:creator" content="@johndoe" />
  <meta name="citation_author" content="Jane Roe" />
  <meta name="citation_author_institution" content="Department of Something, University of Whatever" />
  <meta name="citation_author_institution" content="Department of Whatever, University of Something" />
  <meta name="citation_author_orcid" content="XXXX-XXXX-XXXX-XXXX" />
  <link rel="canonical" href="https://greenelab.github.io/mpmp-manuscript/" />
  <meta property="og:url" content="https://greenelab.github.io/mpmp-manuscript/" />
  <meta property="twitter:url" content="https://greenelab.github.io/mpmp-manuscript/" />
  <meta name="citation_fulltext_html_url" content="https://greenelab.github.io/mpmp-manuscript/" />
  <meta name="citation_pdf_url" content="https://greenelab.github.io/mpmp-manuscript/manuscript.pdf" />
  <link rel="alternate" type="application/pdf" href="https://greenelab.github.io/mpmp-manuscript/manuscript.pdf" />
  <link rel="alternate" type="text/html" href="https://greenelab.github.io/mpmp-manuscript/v/0146d1530e9ee906d634dd42aadc0545aba0d77c/" />
  <meta name="manubot_html_url_versioned" content="https://greenelab.github.io/mpmp-manuscript/v/0146d1530e9ee906d634dd42aadc0545aba0d77c/" />
  <meta name="manubot_pdf_url_versioned" content="https://greenelab.github.io/mpmp-manuscript/v/0146d1530e9ee906d634dd42aadc0545aba0d77c/manuscript.pdf" />
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
([permalink](https://greenelab.github.io/mpmp-manuscript/v/0146d1530e9ee906d634dd42aadc0545aba0d77c/))
was automatically generated
from [greenelab/mpmp-manuscript@0146d15](https://github.com/greenelab/mpmp-manuscript/tree/0146d1530e9ee906d634dd42aadc0545aba0d77c)
on May 20, 2021.
</em></small>

## Authors



+ **John Doe**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [XXXX-XXXX-XXXX-XXXX](https://orcid.org/XXXX-XXXX-XXXX-XXXX)
    · ![GitHub icon](images/github.svg){.inline_icon}
    [johndoe](https://github.com/johndoe)
    · ![Twitter icon](images/twitter.svg){.inline_icon}
    [johndoe](https://twitter.com/johndoe)<br>
  <small>
     Department of Something, University of Whatever
     · Funded by Grant XXXXXXXX
  </small>

+ **Jane Roe**<br>
    ![ORCID icon](images/orcid.svg){.inline_icon}
    [XXXX-XXXX-XXXX-XXXX](https://orcid.org/XXXX-XXXX-XXXX-XXXX)
    · ![GitHub icon](images/github.svg){.inline_icon}
    [janeroe](https://github.com/janeroe)<br>
  <small>
     Department of Something, University of Whatever; Department of Whatever, University of Something
  </small>



## Abstract {.page_break_before}




## Introduction

Although cancer can be initiated and driven by many different genetic alterations, these tend to converge on a limited number of pathways or signaling processes [@doi:10.1016/j.cell.2018.03.035].
A comprehensive understanding of how diverse genetic alterations perturb these central pathways is vital to precision medicine and biomarker identification efforts, as driver mutation status alone confers limited prognostic information [@doi:10.7554/eLife.39217; @doi:10.1038/ncomms12096].
While many methods exist to distinguish driver mutations from passenger mutations based on genomic sequence characteristics [@doi:10.1073/pnas.1616440113; @doi:10.1038/s41467-019-11284-9; @doi:10.1371/journal.pcbi.1006658], until recently it has been a challenge to connect driver mutations to downstream changes in gene expression and cellular function within individual tumor samples.

The Cancer Genome Atlas (TCGA) Pan-Cancer Atlas provides uniformly processed, multi-platform -omics measurements across tens of thousands of samples from 33 cancer types [@doi:10.1038/ng.2764].
Enabled by this publicly available data, a growing body of work on linking the presence of driving genetic alterations in cancer to downstream gene expression changes has emerged.
Recent studies have considered Ras pathway alteration status in colorectal cancer [@doi:10.1158/1078-0432.CCR-13-1943], alteration status across many cancer types in Ras genes [@doi:10.1016/j.celrep.2018.03.046; @doi:10.1093/bib/bbaa258], TP53 [@doi:10.1016/j.celrep.2018.03.076], and PIK3CA [@doi:10.1371/journal.pone.0241514], and alteration status across cancer types in frequently mutated genes [@doi:10.1186/s13059-020-02021-3].
More broadly, other groups have drawn on similar ideas to distinguish between the functional effects of different alterations in the same driver gene [@doi:10.1101/2020.06.02.128850], to link alterations with similar gene expression signatures within cancer types [@doi:10.1142/9789811215636_0031], and to identify trans-acting expression quantitative trait loci (trans-eQTLs) in germline genetic studies [@doi:10.1101/2020.05.07.083386].

These studies share a common thread: they each combine genomic (point mutation and copy number variation) data with transcriptomic (RNA sequencing) data within samples to interrogate the functional effects of genetic variation.
RNA sequencing is ubiquitous and cheap, and its experimental and computational methods are relatively mature, making it a vital tool for generating insight into cancer pathology [@doi:10.1038/nrg.2017.96].
Some driver mutations, however, are known to act indirectly on gene expression through varying mechanisms.
For example, oncogenic IDH1 and IDH2 mutations in glioma have been shown to interfere with histone demethylation, which results in increased DNA methylation and blocked cell differentiation [@doi:10.1056/NEJMoa0808710; @doi:10.1038/nature10860].
Other genes implicated in aberrant DNA methylation in cancer include the TET family of genes [@doi:10.1016/j.tig.2014.07.005] and SETD2 [@doi:10.1101/cshperspect.a026468].
Certain driver mutations, such as those in DNA damage repair genes, may lead to detectable patterns of somatic mutation [@doi:10.1038/nrg3729].
Additionally, correlation between gene expression and protein abundance in cancer cell lines is limited, and proteomics data could correspond more directly to certain cancer phenotypes and pathway perturbations [@doi:10.1016/j.cell.2019.12.023].
In these contexts and others, integrating different data modalities or combining multiple data modalities could be more effective than relying solely on gene expression as a functional signature.

Here, we seek to compare -omics data types profiled in the TCGA Pan-Cancer Atlas for use as a multivariate functional readout of genetic alterations in cancer.
We focus on DNA methylation (27K and 450K probe chips), reverse phase protein array (RPPA), and mutational signatures data [@doi:10.1038/s41586-020-1943-3] as alternative readouts.
Prior studies have identified univariate correlations of CpG site methylation [@doi:10.1371/journal.pcbi.1005840; @doi:10.1186/s12920-018-0425-z] and correlations of RPPA protein profiles [@doi:10.1186/s13073-018-0591-9] with the presence or absence of certain driver mutations.
Other relevant past work includes linking point mutations and copy number variants (CNVs) with changes in methylation and expression at individual genes [@doi:10.1093/bioinformatics/btr019; @doi:10.1093/bib/bbw037] and identifying functional modules that are perturbed by somatic mutations [@doi:10.1093/bioinformatics/btq182; @doi:10.1038/ncomms9554].
However, no direct comparison has been made between different data types for this application, particularly in the multivariate case where we consider changes to -omics-derived gene signatures rather than individual genes in isolation.

We select a wide-ranging collection of potential cancer drivers with varying functions and roles in cancer development [@doi:10.1126/science.1235122].
We use mutation status in these genes as labels to train classifiers, using each of the data types listed as training data, in a pan-cancer setting; we follow similar methods to the elastic net logistic regression approach described in Way et al. 2018 [@doi:10.1016/j.celrep.2018.03.046] and Way et al. 2020 [@doi:10.1186/s13059-020-02021-3].
We show that although there is considerable predictive signal for many genes in each dataset relative to a random baseline, gene expression data tends to provide more effective predictions than the other data types in the vast majority of cases.
In addition, we observe that combining data types into a single multi-omics model provides little, if any, performance benefit over the most performant model using a single data type.
Our results will help to inform the design of future functional genomics studies in cancer, suggesting that RNA sequencing can serve as a broadly effective first-line readout for a variety of genetic alterations.



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

In order to remove missing values from the methylation datasets, we removed the 10 samples with the most missing values, then performed mean imputation for probes with 1 or 2 values missing.
All probes with missing values remaining after sample filtering and imputation were dropped from the analysis.
This left us with 20,040 CpG probes in the 27K methylation dataset, and 370,961 CpG probes in the 450K methylation dataset.
For experiments where "raw" methylation data was used, we used the top 100,000 probes in the 450K dataset by mean absolute deviation for computational efficiency, and we used all of the 20,040 probes in the 27K dataset.
For experiments where "compressed" methylation data was used, we used principal component analysis (PCA), as implemented in the `scikit-learn` Python library [@url:https://jmlr.csail.mit.edu/papers/v12/pedregosa11a.html], to extract the top 5,000 principal components from the methylation datasets.
We initially applied the beta-mixture quantile normalization (BMIQ) method [@doi:10.1093/bioinformatics/bts680] to correct for variability in signal intensity between type I and type II probes, but we observed that this had no effect on our results.
We report uncorrected results in the main paper for simplicity.

To make a fair comparison in each of the experiments displayed in the results, we used the intersection of TCGA samples having measurements for all of the datasets being compared in that experiment.
This resulted in 3 distinct sets of samples: 9,074 samples shared between {expression, mutation} data, 7,981 samples shared between {expression, mutation, 27K methylation, 450K methylation}, and 5,282 samples shared between {expression, mutation, 27K methylation, 450K methylation, RPPA, mutational signatures}.
When we dropped samples between experiments as progressively more data types were added, we observed that the dropped samples had approximately the same cancer type proportions as the dataset as a whole.
In other words, samples that were profiled for one data type but not another did not tend to come exclusively from one or a few cancer types.
Exceptions included acute myeloid leukemia (LAML) which had no samples profiled in the RPPA data, and ovarian cancer (OV) which had only 8 samples with 450K methylation data.
More detailed information on cancer type proportions profiled for each data type is provided in (the supplement).

### Training classifiers to detect cancer mutations

We trained logistic regression classifiers to predict whether or not a given sample has a mutational event in a given target gene, using data from various -omics datasets as explanatory variables.
We explored mutation prediction from gene expression alone using 3 gene sets of equal size: a cancer-related gene dataset of 124 genes described in Vogelstein et al. 2013 [@doi:10.1126/science.1235122], the most mutated genes in TCGA in descending order, and a set of random genes with mutations profiled by MC3.
For each target gene, in order to ensure that the training dataset was reasonably balanced (i.e. that there would be enough mutated samples to train a classifier), we included only cancer types with at least 15 mutated samples and at least 5% mutated samples.
After filtering for sufficient mutated samples, 17 of the genes from the Vogelstein et al. gene set had no valid cancer types remaining, leaving 107 genes with one or more valid cancer types to use in further analyses.
To match the size of this gene set, we took the 107 most frequently mutated genes in TCGA as quantified by MC3, all of which had at least one valid cancer type.
For our random gene set, we first filtered to the set of all genes with 2 or more valid cancer types by the above criteria, then sampled 107 of these genes uniformly at random.
Based on the results of the gene expression experiments, we used the Vogelstein et al. cancer gene set for all subsequent experiments comparing -omics data types.

Our model is trained on -omics data ($X$) to predict mutation presence or absence ($y$) in a target gene.
To control for varying mutation burden per sample, and to adjust for potential cancer type-specific expression patterns, we included one-hot encoded cancer type and log~10~(sample mutation count) in the model as covariates.
Since our -omics datasets tend to have many dimensions and comparatively few samples, we used logistic regression with an elastic net penalty to prevent overfitting [@doi:10.1111/j.1467-9868.2005.00503.x], in line with the approach used in Way et al. 2018 [@doi:10.1016/j.celrep.2018.03.046] and Way et al. 2020 [@doi:10.1186/s13059-020-02021-3].
Elastic net logistic regression finds the feature weights $\hat{w} \in \mathbb{R}^{p}$ solving the following optimization problem:  
$$\hat{w} = \text{argmin}_{w} \ \ell(X, y; w) + \alpha \lambda||w||_1 + \frac{1}{2}\alpha (1 - \lambda) ||w||_2$$

where $i \in \{1, \dots, n\}$ denotes a sample in the dataset, $X_i \in \mathbb{R}^{p}$ denotes features (omics measurements) from the given sample, $y_i \in \{0, 1\}$ denotes the label (mutation presence/absence) for the given sample, and $\ell(\cdot)$ denotes the negative log-likelihood of the observed data given a particular choice of feature weights, i.e.  
$$\ell(X, y; w) = -\sum_{i=1}^{n} y_i \log\left(\frac{1}{1 + e^{-w^{\top}X_i}}\right) + (1 - y_i) \log\left(1 - \frac{1}{1 + e^{-w^{\top}X_i}}\right)$$

This optimization problem leaves two hyperparameters to select: $\alpha$ (controlling the tradeoff between the data log-likelihood and the penalty on large feature weight values), and $\lambda$ (controlling the tradeoff between the L1 penalty and L2 penalty on the weight values).
Although the elastic net loss function does not have a closed form solution, it is convex, and iterative optimization algorithms are commonly used for finding reasonable solutions.
For fixed values of $\alpha$ and $\lambda$, we solved for $\hat{w}$ using stochastic gradient descent, as implemented in `scikit-learn`'s `SGDClassifier` method.

Given weight values $\hat{w}$, it is straightforward to predict the probability of a positive label (mutation in the target gene) $P(y^{*} = 1 \mid X^{*}; \hat{w})$ for a test sample $X^{*}$:  
$$P(y^{*} = 1 \mid X^{*}; \hat{w}) = \frac{1}{1 + e^{-\hat{w}^{\top}X^{*}}}$$

and the probability of no mutation in the target gene, $P(y^{*} = 0 \mid X^{*}; \hat{w})$, is given by (1 - the above quantity).

For each target gene, we evaluated model performance using 2 replicates of 4-fold cross-validation, where train and test splits were stratified by cancer type and sample type.
That is, each training set/test set combination had equal proportions of each cancer type (BRCA, SKCM, COAD, etc) and each sample type (primary tumor, recurrent tumor, etc).
To choose the elastic net hyperparameters, we used 3-fold nested cross-validation, with a grid search over the same hyperparameter ranges used in Way et al. 2020 [@doi:10.1186/s13059-020-02021-3]: $\lambda$ = [0.15, 0.16, 0.2, 0.25, 0.3, 0.4] and $\alpha$ = [0.1, 0.13, 0.15, 0.2, 0.25, 0.3].
Using the grid search results, for each evaluation fold we selected the set of hyperparameters with the optimal area under the receiver-operator curve (AUROC), averaged over the three inner folds.

### Evaluating mutation prediction classifiers

To quantify classification performance for a continuous or probabilistic output, such as that provided by logistic regression, the area under the receiver-operator curve (AUROC) [@doi:10.1016/j.patrec.2005.10.010] and the area under the precision-recall curve (AUPR) [@doi:10.1145/65943.65945] metrics are frequently used.
These metrics summarize performance across a variety of binary label thresholds, rather than requiring choice of a single threshold to determine positive or negative predictions.
In the main text, we report results using AUPR, summarized using average precision.
AUPR has been shown to distinguish between models more accurately than AUROC when there are few positively labeled samples [@doi:10.1371/journal.pone.0118432; @arxiv:2006.11278].
As an additional correction for imbalanced labels, in many of the results in the main text we report the difference in AUPR between a classifier fit to true mutation labels, and a classifier fit to data where the mutation labels are randomly shuffled.
In cases where mutation labels are highly imbalanced (very few mutated samples and many non-mutated samples), a classifier with shuffled labels may perform well simply by chance, e.g. by predicting the negative/non-mutated class for most samples.

Recall that for each target gene and each -omics dataset, we ran 2 replicates of 4-fold cross-validation, for a total of 8 performance results.
To make a statistical comparison between two models using these performance distributions, we used paired-sample _t_-tests, where performance measurements derived from the same cross-validation fold are considered paired measurements.
We used this approach to compare a model trained on true labels with a model trained on shuffled labels (addressing the question, "for the given gene using the given data type, can we predict mutation status better than random"), and to compare a model trained on data type A with a model trained on data type B (addressing the question, "for the given gene, can we make more effective mutation status predictions using data type A or data type B").
We corrected for multiple tests using a Benjamini-Hochberg false discovery rate correction.
For all of our experiments, we set a conservative corrected threshold of $p = 0.001$; we were able to estimate the number of false positives by examining genes with better performance for shuffled mutation labels than true labels.
We chose our threshold to ensure that none of these genes were considered significant, since we would never expect permuting labels to improve performance.
However, our results were not sensitive to the choice of this threshold.

### Multi-omics prediction experiments

To predict mutation presence or absence in cancer genes using multiple data types simultaneously, we concatenated individual datasets into a large feature matrix, then used the same elastic net logistic regression method described previously.
For this task, we considered only the gene expression, 27K methylation, and 450K methylation datasets.
We used only these data types to limit the number of multi-omics combinations: the expression and methylation datasets resulted in the best overall performance across the single-omics experiments so we limited combinations to those datasets here.
For gene expression we used all 15,639 genes available in our RNA sequencing dataset, for the 27K methylation dataset we used all 20,040 CpG probes, and for the 450K methylation dataset we used the top 5,000 principal components.

To construct the multi-omics models, we considered each of the pairwise combinations of the datasets listed above, as well as a combination of all 3 datasets.
When combining multiple datasets, we concatenated along the column axis, including covariates for cancer type and sample mutation burden as before.
For all multi-omics experiments, we used only the samples from TCGA with data for all three data types (i.e. the same 7,981 samples used in the single-omics experiments comparing expression and methylation data types).
Due to computational demands, we considered only a limited subset of the Vogelstein et al. genes as target genes, including EGFR, IDH1, KRAS, PIK3CA, SETD2 and TP53.
We selected these genes because we have previously observed that they have good predictive performance, and they represent a combination of alterations that have strong gene expression signatures (KRAS, EGFR, IDH1, TP53) and strong DNA methylation signatures (IDH1, SETD2, TP53 to some degree).

### Data and code availability

All analyses were implemented in the Python programming language and are available in the following GitHub repository: [`https://github.com/greenelab/mpmp`](https://github.com/greenelab/mpmp), under the open-source BSD 3-clause license. Scripts to download large data files from GDC and other sources are located in the `00_download_data` directory. Scripts to run experiments comparing data modalities used individually are located in the `02_classify_mutations` directory, and scripts to run multi-omics experiments are located in the `05_classify_mutations_multimodal` directory. The Python environment was managed using `conda`, and directions for setting up the environment can be found in the `README.md` file. All analyses were run locally on a CPU.


## Results

### Using diverse data modalities to predict cancer alterations

We collected five different data modalities from cancer samples in the TCGA Pan-Cancer Atlas, capturing five steps of cellular function that are perturbed by genetic alterations in cancer (Figure {@fig:overview}A).
These included gene expression (RNA-seq data), DNA methylation (27K and 450K Illumina BeadChip arrays), protein abundance (RPPA data), microRNA expression data, and patterns of somatic mutation (mutational signatures).
To link these diverse data modalities to changes in mutation status, we used elastic net logistic regression to predict the presence or absence of mutations in cancer genes, using these readouts as predictive features (Figure {@fig:overview}B).
We evaluated the resulting mutation status classifiers in a pan-cancer setting, preserving the proportions of each of the 33 cancer types in TCGA for 8 train/test splits (4 folds x 2 replicates) in each of 107 cancer genes (Figure {@fig:overview}C).

We sought to compare classifiers against a baseline where mutation labels are permuted (to identify genes whose mutation status correlates strongly with a functional signature in a given data type), and also to compare classifiers trained on true labels across different data types (to identify data types that are more or less predictive of mutations in a given gene).
To account for variation between dataset splits in making these comparisons, we treat classification metrics from the 8 train/test splits as performance distributions, which we compare using _t_-tests.
We summarize performance across all genes in our cancer gene set using a similar approach to a volcano plot, in which each point is a gene.
In our summary plots, the x-axis shows the magnitude of the change in the classification metric between conditions, and the y-axis shows the _p_-value for the associated _t_-test.
Figure {@fig:overview}C illustrates this evaluation pipeline.

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
We compared a set of cancer-related genes from Vogelstein et al. 2013 [@doi:10.1126/science.1235122] with a set of random genes having equal size, and a set of the most mutated genes in TCGA having equal size.
For all gene sets, we used only the set of TCGA samples for which both gene expression and somatic mutation data exists, resulting in a total of 9,074 samples from all 33 cancer types (Figure {@fig:expression_gene_sets}A).
This set of samples was further filtered for each target gene to cancer types containing at least 15 mutated samples and at least 5% of samples mutated for that cancer type.
We then evaluated the performance for each target gene in each of the three gene sets.

Figure {@fig:expression_gene_sets}B summarizes performance for all genes in each gene set.
Genes from the Vogelstein et al. set were more predictable than randomly chosen genes or those selected by total mutation count.
Figures {@fig:expression_gene_sets}C, {@fig:expression_gene_sets}D, and {@fig:expression_gene_sets}E show performance for each gene in the random gene set, most mutated gene set, and Vogelstein et al. gene set respectively.
<!-- TODO: I'll fix these numbers once the missing genes stuff gets figured out/finishes running (gene counts for all three datasets will be the same then) -->
<!-- https://github.com/greenelab/mpmp/issues/44 -->
In total, for a significance threshold of $\alpha = 0.001$, 60/98 genes (61.2%) in the Vogelstein et al. gene set are significantly predictable from gene expression data, compared to 10/107 genes (9.35%) in the random gene set and 43/107 genes (40.2%) in the most mutated gene set.
Additionally, Figure {@fig:expression_gene_sets}D shows that many of the significant genes in the most mutated gene set are clustered close to the significance threshold, while the significant genes in the Vogelstein et al. gene set tend to be further from the threshold in Figure {@fig:expression_gene_sets}E (i.e. higher AUPR differences and lower _p_-values).
These results suggest that selecting target genes for mutation prediction based on prior knowledge of their involvement in cancer pathways and processes, rather than randomly or based on mutation frequency alone, can improve predictive signal and identify more highly predictable mutations from gene expression data.

![
**A.** Overlap of TCGA samples between gene expression and MC3 somatic mutation data.
**B.** Overall distribution of performance across 3 gene sets. Each data point represents mean cross-validated AUPR difference, compared with a baseline model trained on permuted mutation presence/absence labels, for one gene in the given gene set; notches show bootstrapped 95% confidence intervals. "random" = 107 random genes, "most mutated" = 107 most mutated genes, "Vogelstein et al." = 107 cancer related genes from Vogelstein et al. 2013 gene set.
**C, D, E.** Volcano-like plots showing mutation presence/absence predictive performance for each gene in each of the 3 gene sets. The _x_-axis shows the difference in mean AUPR compared with a baseline model trained on permuted labels, and the _y_-axis shows _p_-values for a paired _t_-test comparing cross-validated AUPR values within folds.
](images/figure_2.png){#fig:expression_gene_sets width="90%"}

### Comparing six different readouts favors expression and DNA methylation

Next, we expanded our comparison to all five functional data modalities (six total readouts, since there are two DNA methylation datasets) available in the TCGA Pan-Cancer Atlas.
As with previous experiments, we limited our comparison to the set of samples profiled for each readout, resulting in 5,226 samples with data for all readouts.
The data types with the most missing samples were RPPA data (2,215 samples that were missing RPPA data) and 450K methylation (630 samples that were missing 450K methylation data) (Figure {@fig:all_data}A).
Summarized over all genes in the Vogelstein et al. dataset, we observed that gene expression and both methylation datasets tended to produce similar quality predictions, and these tended to be slightly better than the remaining data types (Figure {@fig:all_data}B).
For the set of genes having at least one significant predictor (i.e. "well-predicted" genes), median performance using gene expression was slightly higher than for the methylation data types, although this difference was not statistically significant (Figure {@fig:all_data}C).

On the individual gene level, mutations in 25/75 genes were significantly predictable from RPPA data relative to the permuted baseline, compared to 24/75 genes from microRNA data and 3/75 genes from mutational signatures data (Figure {@fig:all_data}D-F).
For the remaining data types on this smaller set of samples, 34/75 genes outperformed the baseline for gene expression data, 32/75 for 27k methylation, and 30/75 for 450k methylation.
Compared to the methylation experiments shown in Figure {@fig:methylation}, we observed slightly fewer "well-predicted" genes for the expression and methylation datasets here (likely due to the considerably smaller sample size) but relative performance was comparable.
Direct comparisons between each added data type and gene expression data showed that RPPA data produces comparable or improved predictions for a small number of genes (Figure {@fig:all_data}G), but microRNA and mutational signatures data generally provide similar or worse performance (Figure {@fig:all_data}H-I).
Performance using RPPA data is notable because of its drastically smaller dimensionality than the other data types (190 proteins, compared to thousands of features for the expression and methylation data types), suggesting that each protein abundance measurement provides a high information content.
Mutations that are more predictable using RPPA data include _PIK3R1_ and _MAP2K1_ (Figure {@fig:all_data}G).
Both genes are kinases involved in phosphorylation-mediated signal transduction.
The ability of RPPA technology to quantify protein phosphorylation status may thus provide an advantage in identifying mutations in these genes, relative to the other data types we used that cannot directly measure protein phosphorylation.

![
**A.** Overlap of TCGA samples between data types used in mutation prediction comparisons. Only overlaps with more than 100 samples are shown.
**B.** Overall distribution of performance per data type across 75 genes from Vogelstein et al. gene set. Each data point represents mean cross-validated AUPR difference, compared with a baseline model trained on permuted labels, for one gene; notches show bootstrapped 95% confidence intervals.
**C.** Overall performance distribution per data type for genes where the permuted baseline model is significantly outperformed for one or more data types, resulting in a total of 39 genes.
**D, E, F.** Volcano-like plots showing predictive performance for each gene in the Vogelstein et al. gene set, in each of the added data types (RPPA, microRNA, mutational signatures). The _x_-axis shows the difference in mean AUPR compared with a baseline model trained on permuted labels, and the _y_-axis shows _p_-values for a paired _t_-test comparing cross-validated AUPR values within folds.
**G, H, I.** Volcano-like plots comparing predictive performance between data types for each gene in the Vogelstein et al. gene set. The _x_-axis shows the difference in mean AUPR between gene expression and another data type (positive values = better mean performance using gene expression features), and the _y_-axis shows _p_-values for a paired _t_-test comparing cross-validated AUPR values within folds.
](images/figure_5.png){#fig:all_data width="90%"}


## References {.page_break_before}

<!-- Explicitly insert bibliography here -->
<div id="refs"></div>
