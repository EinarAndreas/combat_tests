Methods that remove batch effects while retaining group differences may induce false confidence in downstream analyses
========================================================

## Abstract
Removal of, or adjustment for, batch effects or centre differences is generally required when such effects are present in data. In particular, when preparing microarray gene expression data from multiple cohorts, array platforms, or batches for later analyses, batch effects can have confounding effects. Many methods and tools exist for this purpose. One method, ComBat which is part of the R package sva, is particularly popular due to its ability to remove batch differences even when batches are small and heterogeneous, while retaining the difference between study groups, e.g. treatments. Unfortunately, this recommended and frequently used approach to batch correction will systematically induce incorrect group differences in downstream analyses when groups are unevenly distributed between the batches. The scientific community seems to be largely unaware of this problem, which most likely has contributed to false discoveries being presented in the published literature.

## Introduction


```r
"\nKort oppsummere behovet for batch-korrigering: generelt og i microarraydata.\n\nKort beskrive metoder som gjoer dette (batch-sentrering, ANOVA, ComBat).\n\nForklare tilsiktet effekt disse paa batch+gruppe-forskjeller (balansert versus ikke-balansert).\n\nForklare kort hva faktisk effekt er: ANOVA gir redusert gruppe-forskjeller, mens for ComBat overestimeres den statistiske styrken/sikkerheten av forskjellen.\n"
```

```
## [1] "\nKort oppsummere behovet for batch-korrigering: generelt og i microarraydata.\n\nKort beskrive metoder som gjoer dette (batch-sentrering, ANOVA, ComBat).\n\nForklare tilsiktet effekt disse paa batch+gruppe-forskjeller (balansert versus ikke-balansert).\n\nForklare kort hva faktisk effekt er: ANOVA gir redusert gruppe-forskjeller, mens for ComBat overestimeres den statistiske styrken/sikkerheten av forskjellen.\n"
```


## Simple sanity check
1. Vise til sanity checken med random numbers paa eksempelet fra sva
pakka (bladderbatch) med p-verdi-histogram og mange DE-gener.

## Formler og ord
2. Forklare hva som gaar galt med formler og ord.

## Examples of Undesired consequences
As the amount of false positive results as a consequence of trying to retain group differences depends on the group/batch balance, we will show two examples with varying degree of unbalancedness. The first is a "worst-case"scenario where a few batches are completly missing groups coupled with a a priori assumption of no group difference (24421904). The second compares different cell types and is only mildly unbalanced ("data set 2" from 16632515). We tested their most highlighted group difference for differentially expressed genes using the described approach with ComBat and an alternative approach handling batch inside the test. In addition we performed one sanity check permuting the group labels whitin batches, and one sanity check using random numbers. The different p-value distributions are shown in Figure 2.

![fig. 1](../reanalysis/ComBat_original/figure/pvaluesjohnson.svg)


```r
"\n3. Vise at det er variasjon i stoerrelsen paa feilen introdusert av ComBat\nved a)worst-case scenario (MS-artikel)og b) mild-case (data set 2 fra\nComBat). Vises med p-verdi plot med/uten ComBat og eventuelt antall\nsignifikante prober.\n"
```

```
## [1] "\n3. Vise at det er variasjon i stoerrelsen paa feilen introdusert av ComBat\nved a)worst-case scenario (MS-artikel)og b) mild-case (data set 2 fra\nComBat). Vises med p-verdi plot med/uten ComBat og eventuelt antall\nsignifikante prober.\n"
```

## Discussion

## Sliding application area of ComBat
In the original ComBat article(2007) it is clear that the primary motivation behind ComBat was to employ an Emprical Bayes method for batch-effect removal. The feature of retaining group differences for unbalanced designs is optional and seems to be subordinate, only examplified in the supplementary information. However, over the years this feature became more important judging from advice given by the author to users
(
[1](https://groups.google.com/d/msg/combat-user-forum/eSVSKwGtuyE/ZIWV2juYmmAJ)
,[2](https://groups.google.com/d/msg/combat-user-forum/Vkb9p7wekd4/h5Etie7FTVQJ)
,[3](https://groups.google.com/d/msg/combat-user-forum/fpBTcgDjiR8/uo4QIZL4sZgJ)
,[4](https://groups.google.com/d/msg/combat-user-forum/26FZlgU2LFQ/W6U_Lhh_64EJ)
).
And when ComBat was incorporated in the sva package ( [22257669](http://bioinformatics.oxfordjournals.org/cgi/pmidlookup?view=long&pmid=22257669), 2012), inclusion of group labels was made almost mandatory (an undocumented option of passing a NULL value exists). Thus, this problematic use has likely been common.


## Motivation for this Warning
Our knowledge of ComBat came through a typical use case when trying to salvage unbalanced data which had batch effects. Upon realizing that the confidence on our group differences were exagerated, the literature was searched for a better understanding of correct use and potential overseen limitations of ComBat. But the authors of ComBat and the sva package recomended our usage (16632515, 22257669). In addition, other works looking into the problem of batch effects were mostly recomending ComBat without much concern (22682473, 22133085). A brief inquiry into  some of the articles citing ComBat (574 Google Scholar) revealed few problems, and their method descriptions regarding ComBat were mostly sparse, limited to one or two sentences (24391845, 18414638, 21731603,23630272, 24584070). A further indication of their carefree use of this potentially devastating procedure was the frequent neglect to state the program parameters, i.e batch-labels (18414638, 21731603) or group labels (24391845, 18414638, 21731603, 23630272, 23482648, 24584070). Often no effort was done in order to substantiate the existance of batch effects in their data, ecxept for stating the presence of batches (24391845?, 18414638, 21731603, 23630272). In one instance ComBat was even applied on data were effects due to batch were investigated, but not found (husker ikke ref). The incorporation of the method into analysis pipelines (16642009,TCGA) and other packages (23452776, 21937664) could make its usage even more trivial and parameters setting harder to perceive.
Taken together we fear that many published results from data adjusted by ComBat are completly or partially false. And knowing that scientist don't give up when one analysis fails, a method that almost ensures a result given a sufficently unbalanced design will continue to be used.




```r
"\n\nMange siteringer\nI pipelines 16642009,TCGA  pakker.   23452776, 21937664\nrecomondations, review vinnere. 22682473, 22133085\nmangelfull forstaaelse fra forfattere. 16632515, 22257669\nFaa setninger, 24391845, 18414638, 21731603,23630272, 24584070\nikke oppgit parametere, dvs batch(18414638, 21731603) eller covariates(24391845, 18414638, 21731603, 23630272, 23482648, 24584070),\nIkke vist batch problem (24391845?, 18414638, 21731603, 23630272).\nBrukt uansett.\nP-hacking. Gjoer ComBat populaer.\n"
```

```
## [1] "\n\nMange siteringer\nI pipelines 16642009,TCGA  pakker.   23452776, 21937664\nrecomondations, review vinnere. 22682473, 22133085\nmangelfull forstaaelse fra forfattere. 16632515, 22257669\nFaa setninger, 24391845, 18414638, 21731603,23630272, 24584070\nikke oppgit parametere, dvs batch(18414638, 21731603) eller covariates(24391845, 18414638, 21731603, 23630272, 23482648, 24584070),\nIkke vist batch problem (24391845?, 18414638, 21731603, 23630272).\nBrukt uansett.\nP-hacking. Gjoer ComBat populaer.\n"
```





## Guide for "worried" scientist
6. Gi en rettledning for hvordan en bekymret forsker kanskje kan
estimere feilen p.g.a ComBat-bruk i en eller flere utvalgte artikler.

likely small or no difference (24391845, 18414638?, 21731603, 23630272)

## Reproducible research

## References

Combat org
http://biostatistics.oxfordjournals.org/content/8/1/118.abstract
Adjusting batch effects in microarray expression data using empirical Bayes methods.
Johnson WE, Li C, Rabinovic A.
Biostatistics. 2007 Jan;8(1):118-27. Epub 2006 Apr 21.
PMID: 16632515 [PubMed - indexed for MEDLINE] Free Article

1
Molecular profiling of single Sca-1+/CD34+,- cells--the putative murine lung stem cells.
Hittinger M, Czyz ZT, Huesemann Y, Maneck M, Botteron C, Kaeufl S, Klein CA, Polzer B.
PLoS One. 2013 Dec 31;8(12):e83917. doi: 10.1371/journal.pone.0083917. eCollection 2013.
PMID: 24391845 [PubMed - in process] 
http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0083917
" Based on the hybridization date the data set is separated into two main batches of microarrays (early 2007/2008 and late 2009), which degrade further into small sub-batches. We used the ComBat algorithm [15] to adjust for the main batches. "


2
Drinking-water arsenic exposure modulates gene expression in human lymphocytes from a U.S. population.
Andrew AS, Jewell DA, Mason RA, Whitfield ML, Moore JH, Karagas MR.
Environ Health Perspect. 2008 Apr;116(4):524-31. doi: 10.1289/ehp.10861.
PMID: 18414638 [PubMed - indexed for MEDLINE] Free PMC Article
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2290973/
"followed by empirical Bayes adjustment procedures implemented on the R platform, as described previously (Johnson et al. 2007)"


3
Epigenetic predictor of age.
Bocklandt S, Lin W, Sehl ME, Sánchez FJ, Sinsheimer JS, Horvath S, Vilain E.
PLoS One. 2011;6(6):e14821. doi: 10.1371/journal.pone.0014821. Epub 2011 Jun 22.
PMID: 21731603 [PubMed - indexed for MEDLINE] Free PMC Article
http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0014821
"Batch effect were removed using the Combat algorithm [9], and one outlier sample was removed."

4.
Childhood maltreatment is associated with distinct genomic and epigenetic profiles in posttraumatic stress disorder.
Mehta D, Klengel T, Conneely KN, Smith AK, Altmann A, Pace TW, Rex-Haffner M, Loeschner A, Gonik M, Mercer KB, Bradley B, Müller-Myhsok B, Ressler KJ, Binder EB.
Proc Natl Acad Sci U S A. 2013 May 14;110(20):8302-7. doi: 10.1073/pnas.1217750110. Epub 2013 Apr 29.
PMID: 23630272 [PubMed - indexed for MEDLINE] Free PMC Article
http://www.pnas.org/content/110/20/8302.full
"To correct for confounding as a result of batch effects, the data were normalized using an empirical Bayes method for batch correction " (expression)
"Hybridization and chip batch effects were accounted for using an empirical Bayes method (43)." (methylation)



5.
Patterns of regulatory activity across diverse human cell types predict tissue identity, transcription factor binding, and long-range interactions.
Sheffield NC, Thurman RE, Song L, Safi A, Stamatoyannopoulos JA, Lenhard B, Crawford GE, Furey TS.
Genome Res. 2013 May;23(5):777-88. doi: 10.1101/gr.152140.112. Epub 2013 Mar 12.
PMID: 23482648 [PubMed - indexed for MEDLINE] Free PMC Article
http://genome.cshlp.org/content/23/5/777.full
"We improved on the previously published open chromatin measurements by accounting for batch affects that grouped the data by laboratory rather than by biological signal (see Methods). We used ComBat (Johnson et al. 2007) to remove these batch effects, after which both the DNase I and expression data clustered according to expected biological relationships (Supplemental Fig. S1)."
"Counts were quantile-normalized and scaled, and protocol batch effects were corrected using ComBat (Supplemental Fig. S1; Johnson et al. 2007)."
"Artifactual differences between protocols (Figure S1) were corrected using ComBat (Johnson et al. 2007). Because the data was too large for the ComBat model, we divided it into subsets of about 300,000 DHSs and corrected the batch affect individually on each subset. "
"To make the arrays comparable, we used ComBat to correct for this batch affect (Johnson et al. 2007). After correction with ComBat, the arrays grouped according to expected biological similarity (Figure S1). "


6.
A recurrent inactivating mutation in RHOA GTPase in angioimmunoblastic T cell lymphoma.
Yoo HY, Sung MK, Lee SH, Kim S, Lee H, Park S, Kim SC, Lee B, Rho K, Lee JE, Cho KH, Kim W, Ju H, Kim J, Kim SJ, Kim WS, Lee S, Ko YH.
Nat Genet. 2014 Mar 2. doi: 10.1038/ng.2916. [Epub ahead of print]
PMID: 24584070 [PubMed - as supplied by publisher]
http://www.nature.com/ng/journal/vaop/ncurrent/pdf/ng.2916.pdf
http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.2916.html
http://www.nature.com/ng/journal/vaop/ncurrent/extref/ng.2916-S1.pdf
"We combined the RNA-seq and microarray data using the Combat algorithm36 in the Bioconductor sva package to remove batch biases due to platform differences"


## anbefaler ComBat

10
Batch correction of microarray data substantially improves the identification of genes differentially expressed in rheumatoid arthritis and osteoarthritis.
Kupfer P, Guthke R, Pohlers D, Huber R, Koczan D, Kinne RW.
BMC Med Genomics. 2012 Jun 8;5:23. doi: 10.1186/1755-8794-5-23.
PMID: 22682473 [PubMed - indexed for MEDLINE] Free PMC Article
Related citations
http://www.biomedcentral.com/1755-8794/5/23

11
Relative impact of key sources of systematic noise in Affymetrix and Illumina gene-expression microarray experiments.
Kitchen RR, Sabine VS, Simen AA, Dixon JM, Bartlett JM, Sims AH.
BMC Genomics. 2011 Dec 1;12:589. doi: 10.1186/1471-2164-12-589.
PMID: 22133085 [PubMed - indexed for MEDLINE] Free PMC Article
http://www.biomedcentral.com/1471-2164/12/589/


## skeptisk til ComBat

20 - HAr kjoert permutering av covariatelabler og finner at combat finner masse gener uansett. Fokus paa artikkel er ikke combat spesifikkt.
Expanding the understanding of biases in development of clinical-grade molecular signatures: a case study in acute respiratory viral infections.
Lytkin NI, McVoy L, Weitkamp JH, Aliferis CF, Statnikov A.
PLoS One. 2011;6(6):e20662. doi: 10.1371/journal.pone.0020662. Epub 2011 Jun 1.
PMID: 21673802 [PubMed - indexed for MEDLINE] Free PMC Article
http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0020662


## ComBat i pipelines

TCGA Batch Effects Tool
http://bioinformatics.mdanderson.org/main/TCGABatchEffects:Overview

GenePattern 2.0.
Reich M, Liefeld T, Gould J, Lerner J, Tamayo P, Mesirov JP.
Nat Genet. 2006 May;38(5):500-1. No abstract available.
PMID: 16642009 [PubMed - indexed for MEDLINE]
http://www.nature.com/ng/journal/v38/n5/full/ng0506-500.html

Uverifisert.
insilicodb? bioconductopakke eller hva
http://bib.oxfordjournals.org/content/early/2012/07/31/bib.bbs037.short


The sva package for removing batch effects and other unwanted variation in high-throughput experiments.
Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD.
Bioinformatics. 2012 Mar 15;28(6):882-3. doi: 10.1093/bioinformatics/bts034. Epub 2012 Jan 17.
PMID: 22257669 [PubMed - indexed for MEDLINE] Free PMC Article
http://bioinformatics.oxfordjournals.org/cgi/pmidlookup?view=long&pmid=22257669

virtualArray: a R/bioconductor package to merge raw data from different microarray platforms.
Heider A, Alt R.
BMC Bioinformatics. 2013 Mar 2;14:75. doi: 10.1186/1471-2105-14-75.
PMID: 23452776 [PubMed - indexed for MEDLINE] Free PMC Article
http://bioconductor.org/packages/release/bioc/html/virtualArray.html


inSilicoDb: an R/Bioconductor package for accessing human Affymetrix expert-curated datasets from GEO.
Taminau J, Steenhoff D, Coletta A, Meganck S, Lazar C, de Schaetzen V, Duque R, Molter C, Bersini H, Nowé A, Weiss Solís DY.
Bioinformatics. 2011 Nov 15;27(22):3204-5. doi: 10.1093/bioinformatics/btr529. Epub 2011 Sep 21.
PMID: 21937664 [PubMed - indexed for MEDLINE] Free Article


## eksempler med reporduserte p-value

Comparing the biological impact of glatiramer acetate with the biological impact of a generic.
Towfic F, Funt JM, Fowler KD, Bakshi S, Blaugrund E, Artyomov MN, Hayden MR, Ladkani D, Schwartz R, Zeskind B.
PLoS One. 2014 Jan 8;9(1):e83757. doi: 10.1371/journal.pone.0083757. eCollection 2014.
PMID: 24421904 [PubMed - in process] Free PMC Article
http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0083757

### user group recomandations
https://groups.google.com/d/msg/combat-user-forum/eSVSKwGtuyE/ZIWV2juYmmAJ
https://groups.google.com/d/msg/combat-user-forum/Vkb9p7wekd4/h5Etie7FTVQJ
https://groups.google.com/d/msg/combat-user-forum/fpBTcgDjiR8/uo4QIZL4sZgJ
https://groups.google.com/d/msg/combat-user-forum/26FZlgU2LFQ/W6U_Lhh_64EJ

Alternativer til tittel
Methods that remove batch effects while retaining group differences may lead to incorrect conclusions in downstream analyses
Methods that remove batch effects while retaining group differences may lead to incorrect conclusions






