# CYP2D6 Clinical Findings Summary

## 1. Highest Clinical Risk Samples
* **Samples:** **S10** and **S11** (tied with identical results: `*2;*2;*41 | *2;*41;*41 | *2;*41;*68 | *41;*41;*68 | *41;*68;*68`).
* **Why It Is High Risk:** 
  * Both samples have an extra gene copy (`total_cn = 3`) and an unassayed deletion/duplication marker (`4125_4133dup`), resulting in five possible diplotypes.
  * These five possibilities span **three different metabolizer types**: **Intermediate** (AS `0.5`), **Normal** (AS `1.5`–`2.0`), and **Ultrarapid** (AS `2.5`).
  * Guessing a single result could cause severe dosing errors—such as life-threatening opioid toxicity (e.g., with codeine or tramadol) if an Ultrarapid metabolizer is treated as Normal, or treatment failure/drug toxicity if an Intermediate metabolizer is misclassified.
  * *(Note: Among standard 2-copy samples, **S13** [`*2;*4|*4;*68`] and **S18** [`*2;*5|*5;*68`] also carry high risk because they span **Intermediate** [AS `1.0`] and **Poor** [AS `0.0`] metabolizer status.)*

### Required Clinical Report Statement
> **Clinical Reporting Notice (S10 / S11):** Definitive *CYP2D6* diplotype and metabolizer phenotype assignment cannot be resolved from current array data due to an unphased copy-number gain (`total_cn = 3`) and missing coverage for the non-functional `*68` variant (`4125_4133dup`). Consistent diplotypes span an Activity Score of **0.5 to 2.5** (**Intermediate**, **Normal**, and **Ultrarapid** metabolizer phenotypes). Do not assume Normal Metabolizer status; avoid CYP2D6-dependent prodrugs with narrow therapeutic indices (e.g., codeine, tramadol) or confirm allele-specific copy number and `*68` status via targeted follow-up testing before prescribing.

---

## 2. Identical Platform Results from Different Genotypes
* **Samples:** **S17** and **S24** (both output `*1;*2 | *1;*68`, Activity Score `1`–`2`, `Normal|Intermediate`). *(Note: **S10** and **S11** also share identical platform outputs.)*
* **Why They Look Identical:**
  * Both samples have 2 gene copies (`total_cn = 2`) and heterozygous (`0/1`) calls at `2850C>T` and `4180G>C`.
  * Because the platform does not test for `4125_4133dup`—the only marker separating normal-function `*2` (AS `1.0`) from non-functional `*68` (AS `0.0`)—a true **`*1;*2`** sample (Normal Metabolizer) and a true **`*1;*68`** sample (Intermediate Metabolizer) produce the exact same test data.

### Most Cost-Effective Follow-Up Test
* **Recommended Assay:** Targeted PCR with Sanger sequencing or fragment-length analysis (capillary electrophoresis) specifically targeting the **`4125_4133dup`** locus in exon 9.
* **Capabilities:** Directly confirms whether `4125_4133dup` is present or absent at low cost, cleanly separating `*68` from `*2` in standard 2-copy samples like S17 and S24.
* **Limitations:**
  * **No chromosome phasing:** Short-amplicon testing cannot tell which chromosome carries the duplication in complex compound heterozygotes without long-range allele-specific PCR.
  * **No copy-number quantification:** Cannot measure exact allele dosage in multi-copy samples (`total_cn > 2`, such as S10/S11) or detect whole-gene deletions (`*5`) and hybrid tandem genes (`*36+*10`, `*68+*4`), which require ddPCR, MLPA, or long-read sequencing.