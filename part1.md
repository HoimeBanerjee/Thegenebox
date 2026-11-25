## Section 1

### 1.1

* **GGB9712 (Genus-level Genome Bin #9712):** A GGB groups related SGBs at a 15% whole-genome genetic distance threshold approximating the genus level.
* **SGB15212 (Species-level Genome Bin #15212):** Corresponds to the standard prokaryotic species boundary.
* **1.2% (Relative Abundance):** During sample profiling, metagenomic reads align specifically to the unique core marker genes defining SGB15212. Because these markers are exclusive to that specific MAG cluster, the classifier quantifies the organism at 1.2% relative abundance despite the lack of a formal species name.
* **Conclusion:** According to the classification, this indicates the presence of a distinct, reproducible prokaryotic species-level lineage (SGB15212) within an uncharacterized genus (GGB9712) that has not yet been isolated, cultured, and formally named in taxonomic repositories.

> **Reference:** Edwin, N. R., Fitzpatrick, A. H., Brennan, F., Abram, F., & O’Sullivan, O. (2024). An in-depth evaluation of metagenomic classifiers for soil microbiomes. *Environmental Microbiome*, 19(1), 19.

---

### 1.2

Since microbial kingdoms differ dramatically in genome size, gene density, and ploidy, these metrics are mathematically and biologically non-equivalent.

| Metrics | Marker-Gene Profiling | K-mer Classification |
| :--- | :--- | :--- |
| **Numerator** | Length-normalized read depth (RPK) summed across unique clade-specific marker genes for a given kingdom. | Total number of sequencing reads (or read pairs) assigned to a given kingdom via lowest common ancestor (LCA) k-mer matching. |
| **Denominator** | Sum of length-normalized read depths across all detected clades (all profiled kingdoms/SGBs) in the sample. | Either total sequenced reads in the library (including unclassified and host reads) or total classified microbial reads (excluding unclassified reads, depending on reporting mode). |

> **Reference:** Bokulich, N. A. (2025). Integrating sequence composition information into microbial diversity analyses with k-mer frequency counting. *mSystems*, 10(3), e01550-24.

---

### 1.3

This happens because the summation double-counts every functional feature by adding unstratified community totals directly to their own taxon-stratified breakdowns.

In some functional metagenomic profilers like HUMAnN, each gene family, pathway, or reaction is reported across two hierarchical tiers within the exact same table as **Unstratified (Community-Level) Rows** followed by **Stratified (Taxon-Specific) Rows**. Therefore, when a column normalized to relative abundance ($\text{sum} = 1.0$ or $100\%$) is summed incorrectly across all rows without filtering by stratification level, it yields $\approx 200\%$.

> **Reference:** Franzosa, E. A., McIver, L. J., Rahnavard, G., Thompson, L. R., Schirmer, M., Weingart, G., ... & Huttenhower, C. (2018). Species-level functional profiling of metagenomes and metatranscriptomes. *Nature Methods*, 15(11), 962-968.

---

### 1.4

| Analytical Aspect | Mechanism / Finding | Summary Explanation |
| :--- | :--- | :--- |
| **Hybridization & Calling Artifact** | BAF Clustering Shift | The single intact allele hybridizes only to the reference (A) probe ($I_B = 0$), yielding $\text{BAF} = 0.0$. Standard callers map this directly to the diploid `AA` cluster while ignoring the $\sim 50\%$ drop in total signal intensity ($\text{LRR} \approx -0.5$). |
| **True Genotype Representation** | Hemizygous State | Written as `A/-` (or `*1/*5` in PGx star-allele notation, such as *CYP2D6*), reflecting one intact reference allele and one deleted allele rather than two functional copies (`*1/*1`). |
| **Clinical Risks in PGx Reporting** | Clinical Risks in PGx Reporting | Inflates the predicted metabolizer status (e.g., calling an Intermediate Metabolizer a Normal Metabolizer), risking toxicity from standard doses of active drugs or treatment failure from under-activated prodrugs. |

> **Reference:** Cantsilieris, S., Baird, P. N., & White, S. J. (2013). Molecular methods for genotyping complex copy number polymorphisms. *Genomics*, 101(2), 86-93.

---

### 1.5

**Answering using the example of Kraken 2:**

| Aspect | Issue / Action | Details |
| :--- | :--- | :--- |
| **Root Cause 1: Schema Mismatch** | 8-Column vs. 6-Column Layout | `--report-minimizer-data` inserts two additional columns (total minimizers and distinct minimizers) after column 3 (direct reads). This shifts rank code, NCBI taxID, and scientific name to columns 6–8, breaking Bracken's fixed 6-column parser. |
| **Root Cause 2: Model Violation** | Minimizer Counts vs. Read Counts | Bracken's Bayesian redistribution (`databaseXmers.kmer_distrib`) requires discrete read counts for a known read length. Minimizer sketch counts represent sub-sampled sequence diversity and violate its probabilistic model. |
| **Remediation Step 1** | Restore Standard 6-Column Report | Run Kraken 2 without `--report-minimizer-data`, or remove columns 4–5 from an existing report:<br>`cut -f1-3,6-8 kraken2_minimizer.report > kraken2_standard.report` |
| **Remediation Step 2** | Run Bracken on Standard Report | Execute Bracken using the cleaned 6-column report:<br>`bracken -d /path/to/db -i kraken2_standard.report -o sample.bracken -r 150 -l S` |
| **Remediation Step 3** | Post-Hoc Minimizer Filtering | Keep the 8-column minimizer report separate and join it to `sample.bracken` using NCBI taxID after re-estimation. Use distinct minimizer thresholds to filter out low-coverage false positives. |

> **Reference:** Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20(1), 257.