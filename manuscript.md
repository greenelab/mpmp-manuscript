---
title: Manuscript Title
keywords:
- markdown
- publishing
- manubot
lang: en-US
date-meta: '2021-04-07'
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
  <meta name="dc.date" content="2021-04-07" />
  <meta name="citation_publication_date" content="2021-04-07" />
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
  <link rel="alternate" type="text/html" href="https://greenelab.github.io/mpmp-manuscript/v/296df815ee993739f3e34e72f0bae27f5c20ab41/" />
  <meta name="manubot_html_url_versioned" content="https://greenelab.github.io/mpmp-manuscript/v/296df815ee993739f3e34e72f0bae27f5c20ab41/" />
  <meta name="manubot_pdf_url_versioned" content="https://greenelab.github.io/mpmp-manuscript/v/296df815ee993739f3e34e72f0bae27f5c20ab41/manuscript.pdf" />
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
([permalink](https://greenelab.github.io/mpmp-manuscript/v/296df815ee993739f3e34e72f0bae27f5c20ab41/))
was automatically generated
from [greenelab/mpmp-manuscript@296df81](https://github.com/greenelab/mpmp-manuscript/tree/296df815ee993739f3e34e72f0bae27f5c20ab41)
on April 7, 2021.
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



## References {.page_break_before}

<!-- Explicitly insert bibliography here -->
<div id="refs"></div>
