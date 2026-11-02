# Study bias inflates cross-phenotype gene overlaps: a macrocephaly case study

**Hafsa Fathima Khajamainudeen**<br>
M.Sc. Bioinformatics & Data Science, Sathyabama Institute of Science and Technology, Chennai, India<br>
Correspondence: hafsafathima2004@gmail.com<br>
ORCID: https://orcid.org/0009-0007-5168-4050

---

## Abstract

Gene sets drawn from phenotype ontologies are widely used to ask how much two conditions share their genetic basis. These comparisons rest on an assumption that is rarely tested: that a gene is equally likely to be annotated to any phenotype. In practice a small number of genes have been studied so intensively that they are annotated to hundreds of phenotypes, while most genes have only a few annotations. Shared membership of these heavily studied genes can create the appearance of a shared genetic basis where none exists. Using macrocephaly as a case study, I built gene sets for macrocephaly, microcephaly, autism and seizure from the Human Phenotype Ontology and measured their pairwise overlaps with the standard hypergeometric test. All three overlaps with macrocephaly were strongly significant (1.5–2.2 fold, p < 1e-4). I then retested each overlap against a permutation null matched to each gene's annotation load, so that the comparison asks whether the overlap exceeds what an equally well-studied random gene set would give. Under this test the macrocephaly–microcephaly overlap fell to 1.00 fold (p = 0.51) and the macrocephaly–autism overlap to 1.07 fold (p = 0.29); only the macrocephaly–seizure overlap remained significant, and modestly so (1.14 fold, p < 1e-4). The 37 genes shared by all four phenotypes, which include several of the most heavily annotated genes in the ontology, were also consistent with study bias alone (1.03 fold, p = 0.45). These results show that ontology-based overlap statistics can be dominated by how much each gene has been studied, and that a study-bias-matched null is needed before an overlap is interpreted as biological.

**Keywords:** Human Phenotype Ontology, gene set overlap, ascertainment bias, macrocephaly, neurodevelopmental phenotypes, permutation testing

---

## 1. Introduction

Comparing the genes behind two conditions is a common first step in understanding whether they share a mechanism. Phenotype ontologies such as the Human Phenotype Ontology (HPO) make this easy: each phenotype maps to a set of genes, and the overlap between two sets can be tested for significance with a hypergeometric test. Overlaps found this way are often read as evidence of shared biology.

This reading depends on an assumption that the hypergeometric test makes but does not check — that every gene has an equal chance of being annotated to a phenotype. Phenotype annotations come from published case reports and curated literature, and that literature is uneven. A gene such as MECP2 or SHANK3 has been studied in thousands of papers and is annotated to a very large number of phenotypes, while most genes are annotated to only a handful. When two phenotype sets both contain these heavily studied genes, they overlap for a reason that has nothing to do with shared biology: the genes are shared because they are well studied, not because the two conditions converge on them.

Macrocephaly is a useful case to examine this problem. It is genetically heterogeneous, it has reported links to several neurodevelopmental phenotypes, and some of those links — for example to autism through the PI3K/AKT/mTOR pathway — are well established (Yeung et al., 2017). This makes it a setting where real overlap and study-bias-driven overlap are both plausible, and where the two need to be told apart. Here I measure the overlap of macrocephaly with microcephaly, autism and seizure using the standard test, then retest each overlap against a null that controls for how heavily each gene has been studied, and show that most of the apparent overlap does not survive.

## 2. Methods

**Gene sets.** Gene sets were built from the Human Phenotype Ontology (Gargano et al., 2024) for macrocephaly (HP:0000256), microcephaly (HP:0000252), autism (HP:0000717) and seizure (HP:0001250), each propagated to include all descendant terms. The background was the set of all HPO-annotated genes. Using HPO-annotated genes as the background, rather than the whole genome, removes the coarsest form of bias, because it asks what is specific to a phenotype among genes known to cause disease.

**Overlap tests.** Pairwise overlaps between macrocephaly and each other phenotype were tested with the hypergeometric distribution, with Benjamini–Hochberg correction and Jaccard similarity reported alongside.

**Study-bias-matched null.** The study load of each gene was taken as the number of HPO terms it is annotated to, obtained from PyHPO. The background pool was binned into deciles of study load. For each permutation, a pseudo-macrocephaly gene set was drawn from the background matched to the real macrocephaly set's study-load distribution — the same number of genes from each decile — and its overlap with the fixed comparison set was recorded. A size-only null, drawing genes at random without study-load matching, was run in parallel as a control; it reproduces the hypergeometric result. Each null used 10,000 permutations with a fixed seed. The empirical p-value is the proportion of permutations with an overlap at least as large as observed. Fold enrichment is the observed overlap divided by the mean of the null. All sets were restricted to genes with study-load data so that the observed and permuted values are measured on the same gene universe.

**Descriptive network.** For context, STRING interactions among the macrocephaly genes were retrieved at a combined score ≥ 0.7, tested for PPI enrichment, and partitioned into modules with the Louvain method (seed 42). Pathway enrichment used g:Profiler with the HPO-annotated genes as a custom background.

**Code and reproducibility.** The analysis is in five ordered notebooks. The databases update over time, so counts and p-values shift slightly between runs; the analysis reported here was run in October 2026.

## 3. Results

**The standard test finds strong overlaps.** By the hypergeometric test, macrocephaly overlaps significantly with all three phenotypes: microcephaly (1.54 fold), autism (2.17 fold) and seizure (1.63 fold), each at p < 1e-4. These values match the analytic hypergeometric result and are consistent with the overlaps reported for these phenotypes elsewhere.

**Most overlaps do not survive study-bias control.** Against the study-bias-matched null, the picture changes (Table 1, Figure 1). The macrocephaly–microcephaly overlap falls to 1.00 fold (p = 0.51) — exactly what study bias alone predicts. The macrocephaly–autism overlap, the strongest under the standard test, falls to 1.07 fold (p = 0.29). Only the macrocephaly–seizure overlap remains significant, and the real excess is modest at 1.14 fold (p < 1e-4), far smaller than the 1.63 fold the standard test reported.

**Table 1.** Overlap of macrocephaly with each phenotype under the two nulls.

| Comparison | Observed | Size-only fold (p) | Bias-matched fold (p) | Survives |
|---|---|---|---|---|
| vs microcephaly | 192 | 1.54 (<1e-4) | 1.00 (0.51) | No |
| vs autism | 55 | 2.17 (<1e-4) | 1.07 (0.29) | No |
| vs seizure | 301 | 1.63 (<1e-4) | 1.14 (<1e-4) | Yes |

**The shared "core" is a study-bias artifact.** Thirty-seven genes are annotated to all four phenotypes (31 have study-load data and could be tested). Under the size-only null this four-way overlap is 2.16 fold and highly significant, but under the bias-matched null it is 1.03 fold (p = 0.45). These genes carry a median of 113 HPO annotations each, against a background median of 37 (Figure 1, right). They are shared across the four phenotypes because they are among the most studied genes in the ontology, not because the four conditions specifically converge on them.

**Figure 1.** Left: fold enrichment of each macrocephaly overlap under the size-only null (≈ hypergeometric) and the study-bias-matched null, with bias-matched p-values. Right: distribution of HPO annotation load for all genes versus the shared-core genes.

**Macrocephaly's own gene set has a two-route structure (descriptive).** Independent of the overlap question, the macrocephaly genes separate into two groups: a growth-signalling group (PI3K/AKT/mTOR, Ras/MAPK, Hedgehog, Wnt) and a metabolic/structural group (peroxisome biogenesis, O-mannose glycosylation, lysosomal storage). The STRING network is denser than a random set of the same size (1,268 vs 565 expected edges) and splits into modules matching these pathways. This structure reproduces the known clinical distinction between megalencephaly and other causes of macrocephaly and is reported here as context rather than as a new finding; edge density and hub prominence are themselves partly driven by study, so they are not treated as independent evidence.

## 4. Discussion

The central result is a cautionary one. Three overlaps that look convincing under the standard test — one of them more than twofold — largely disappear when the test controls for how much each gene has been studied. The overlap that the field would most expect to be real, macrocephaly with autism, is among those that do not survive. The one overlap that does survive, with seizure, is real but much smaller than the raw number implied.

This does not mean the phenotypes are genetically unrelated. It means that ontology-based overlap counts are not, on their own, evidence of a shared genetic basis, because they are confounded by the uneven study of genes. The confounder is strong: the genes shared by all four phenotypes are annotated several times more heavily than average, and a null that simply matches set size cannot separate that from biology. A study-bias-matched null can, and it is inexpensive to run.

The practical recommendation is that any claim of shared genetics from phenotype-ontology overlap should be checked against a study-bias-matched null before it is interpreted, and that the surviving fold enrichment, rather than the raw count or the hypergeometric p-value, should be reported.

## 5. Limitations

- The bias-matched null controls for each gene's annotation load, but annotation load is a proxy for study intensity and does not capture every source of ascertainment bias.
- The null resamples the macrocephaly set against a fixed comparison set. The comparison sets are themselves study-biased, so a two-sided matched resampling would be a stricter test.
- STRING edge density and hub status are partly study-driven, so the descriptive network is not independent evidence of biological connectivity.
- Overlaps describe shared reported genetics and point to candidate shared mechanisms; they do not establish them.

## Data and code availability

All gene sets, result tables, the network file and the five analysis notebooks are available at the project repository: https://github.com/HafsaFathimakhaja/macrocephaly-gene-landscape

## Funding and conflicts

No specific funding was received for this work. The author declares no competing interests.

## References

1. Gargano MA, et al. The Human Phenotype Ontology in 2024: phenotypes around the world. *Nucleic Acids Research* 52, D1333–D1346 (2024).
2. Marcello J. PyHPO: a Python library to work with, analyze, filter and inspect the Human Phenotype Ontology. Software, https://pypi.org/project/pyhpo/ (annotation load used as the study-bias covariate).
3. Kolberg L, et al. g:Profiler: interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update). *Nucleic Acids Research* 51, W207–W212 (2023).
4. Szklarczyk D, et al. The STRING database in 2023. *Nucleic Acids Research* 51, D638–D646 (2023).
5. Blondel VD, et al. Fast unfolding of communities in large networks. *Journal of Statistical Mechanics* P10008 (2008).
6. Yeung KS, Tso WWY, Ip JJK, et al. Identification of mutations in the PI3K-AKT-mTOR signalling pathway in patients with macrocephaly and developmental delay and/or autism. *Molecular Autism* 8, 66 (2017). doi:10.1186/s13229-017-0182-4
