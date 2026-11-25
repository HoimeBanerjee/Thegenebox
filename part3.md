## Part 3

### 3.1

#### (a) Positive Predictive Value (PPV) Derivation

1. **Given Parameters:**
   * Prevalence ($P$) $= \frac{1}{10,000} = 0.0001$
   * Probability of true negative population ($1 - P$) $= 1 - 0.0001 = 0.9999$
   * Sensitivity ($\text{Se}$) $= 0.99$
   * Specificity ($\text{Sp}$) $= 0.999 \implies \text{False Positive Rate (FPR)} = 1 - \text{Sp} = 1 - 0.999 = 0.001$

2. **Bayes' Theorem Formulation:**
   $$\text{PPV} = \frac{\text{Se} \cdot P}{(\text{Se} \cdot P) + (1 - \text{Sp})(1 - P)}$$

3. **Step-by-Step Substitution:**
   * **True Positive Fraction ($\text{TP}$):**
     $$\text{TP} = 0.99 \times 0.0001 = 0.000099$$
   * **False Positive Fraction ($\text{FP}$):**
     $$\text{FP} = 0.001 \times 0.9999 = 0.0009999$$
   * **Total Positive Call Fraction ($\text{TP} + \text{FP}$):**
     $$\text{TP} + \text{FP} = 0.000099 + 0.0009999 = 0.0010989$$
   * **PPV Calculation:**
     $$\text{PPV} = \frac{0.000099}{0.0010989} = \frac{99}{1098.9} = \frac{1}{11.1} \approx 0.09009009 \quad (9.01\%)$$

---

#### (b) Target Specificity Derivation for 95% PPV

1. **Target Equation:** Set $\text{PPV}_{\text{target}} = 0.95$, keeping $\text{Se} = 0.99$ and $P = 0.0001$ constant:
   $$0.95 = \frac{0.000099}{0.000099 + (1 - \text{Sp}_{\text{req}})(0.9999)}$$

2. **Algebraic Rearrangement for Required False Positive Rate ($1 - \text{Sp}_{\text{req}}$):**
   $$0.95 \times [0.000099 + (1 - \text{Sp}_{\text{req}})(0.9999)] = 0.000099$$
   $$0.00009405 + 0.95(0.9999)(1 - \text{Sp}_{\text{req}}) = 0.000099$$
   $$0.949905(1 - \text{Sp}_{\text{req}}) = 0.000099 - 0.00009405$$
   $$0.949905(1 - \text{Sp}_{\text{req}}) = 0.00000495$$
   $$1 - \text{Sp}_{\text{req}} = \frac{0.00000495}{0.949905} = \frac{495}{94,990,500} = \frac{1}{191,900} \approx 5.211047 \times 10^{-6}$$

3. **Required Target Specificity ($\text{Sp}_{\text{req}}$):**
   $$\text{Sp}_{\text{req}} = 1 - \frac{1}{191,900} = \frac{191,899}{191,900} \approx 0.99999479 \quad (99.999479\%)$$

4. **Final Threshold Expression:** Since $\text{FPR}_{\text{req}} = 1 - \text{Sp}_{\text{req}} = \frac{1}{191,900}$, the required threshold is **one false positive per 191,900 negatives**.

---

### 3.2

#### (a) Expected Total Unsequenced Bases (Zero Coverage)

1. **Given Parameters:**
   * Mean sequencing coverage ($\lambda$) $= 5\times$
   * Genome length ($G$) $= 3.1\text{ Gb} = 3.1 \times 10^9\text{ bp}$

2. **Poisson Probability of Zero Coverage:**
   According to the Lander–Waterman Poisson model, the probability that any genomic position has read depth $K = k$ is given by:
   $$P(K = k) = \frac{\lambda^k e^{-\lambda}}{k!}$$
   Putting $k = 0$ and $\lambda = 5$:
   $$P(K = 0) = \frac{5^0 e^{-5}}{0!} = e^{-5} \approx 0.00674$$

3. **Expected Unsequenced Bases ($E[B_0]$):**
   By linearity of expectation across all $G$ independent genomic positions:
   $$E[B_0] = G \cdot P(K = 0) = (3.1 \times 10^9\text{ bp}) \cdot e^{-5}$$
   $$E[B_0] = 3.1 \times 10^9 \times 0.00674 = 20,887,600\text{ bp} \quad \text{or } (20.89\text{ Mb})$$

---

#### (b) Exact Probability of Observing Both Alleles at a Heterozygous Locus

1. **Conditional Probability Given Read Depth $K = k$:**
   At a diploid heterozygous site (alleles $A$ and $B$), each read independently samples allele $A$ with $p = 0.5$ or allele $B$ with $q = 1 - p = 0.5$.
   * For $k = 0$ or $k = 1$, it is impossible to observe both alleles: $P(\text{Both} \mid K = 0) = P(\text{Both} \mid K = 1) = 0$.
   * For $k \ge 2$, the only failure modes are observing exclusively allele $A$ (probability $0.5^k$) or exclusively allele $B$ (probability $0.5^k$). Thus:
     $$P(\text{Both} \mid K = k) = 1 - [(0.5)^k + (0.5)^k] = 1 - 2(0.5)^k$$

2. **Integration (Summation) Across the Full Poisson Depth Distribution:**
   Applying the Law of Total Probability over all $k \ge 1$:
   $$P(\text{Both}) = \sum_{k=1}^{\infty} P(K = k) \cdot P(\text{Both} \mid K = k) = \sum_{k=1}^{\infty} \frac{\lambda^k e^{-\lambda}}{k!} [1 - 2(0.5)^k]$$
   Split the infinite series into two standard Taylor exponential series:
   $$P(\text{Both}) = e^{-\lambda} \sum_{k=1}^{\infty} \frac{\lambda^k}{k!} - 2e^{-\lambda} \sum_{k=1}^{\infty} \frac{(0.5\lambda)^k}{k!}$$
   Using the identity $\sum_{k=1}^{\infty} \frac{x^k}{k!} = e^x - 1$:
   $$P(\text{Both}) = e^{-\lambda}(e^\lambda - 1) - 2e^{-\lambda}(e^{0.5\lambda} - 1)$$
   $$P(\text{Both}) = (1 - e^{-\lambda}) - (2e^{-0.5\lambda} - 2e^{-\lambda})$$
   $$P(\text{Both}) = 1 - 2e^{-0.5\lambda} + e^{-\lambda}$$

3. **Numerical Evaluation at $\lambda = 5$:**
   $$P(\text{Both}) = 1 - 2e^{-0.5(5)} + e^{-5} = 1 - 2e^{-2.5} + e^{-5}$$
   Substitute the exponential terms:
   * $e^{-2.5} = 0.0820849986 \implies 2e^{-2.5} \approx 0.1641$
   * $e^{-5} \approx 0.006738$
   $$P(\text{Both}) = 1 - 0.1641 + 0.0068 = 0.842 \quad (84.2\%)$$

---

### 3.3

$$\text{Absolute abundance} = \text{Relative abundance} \times \text{Total microbial load}$$

* **Sample A:**
  $$0.04 \times 1.0 \times 10^{11} = 4.0 \times 10^9\text{ cells/g}$$
* **Sample B:**
  $$0.02 \times 2.5 \times 10^{11} = 5.0 \times 10^9\text{ cells/g}$$
* **Fold Change:**
  $$\frac{5.0 \times 10^9}{4.0 \times 10^9} = 1.25$$
* **$\log_2$ Fold Change ($\log_2\text{FC}$):**
  $$\log_2\text{FC} = \log_2(1.25) = 0.322$$

---

### 3.4

#### (a) Expected Baseline Concordance Across a Uniform MAF Distribution

1. **Hardy–Weinberg Genotype Concordance Function:**
   For a biallelic locus with minor allele frequency (MAF) $q \in (0, 0.5]$ and major (reference) allele frequency $p = 1 - q$, the genotype probabilities under Hardy–Weinberg equilibrium (HWE) are:
   * Homozygous reference ($p^2$): $(1 - q)^2$
   * Heterozygous ($2pq$): $2q(1 - q)$
   * Homozygous alternate ($q^2$): $q^2$

   Because the baseline caller unconditionally predicts homozygous reference at every locus, its concordance $C(q)$ at a locus with MAF $q$ equals the true frequency of homozygous reference individuals:
   $$C(q) = (1 - q)^2 = 1 - 2q + q^2$$

2. **Integration Over $q \sim \text{Uniform}(0, 0.5]$:**
   The probability density function (PDF) of $q$ on the interval $(0, 0.5]$ is:
   $$f(q) = \frac{1}{0.5 - 0} = 2 \quad \text{for } 0 < q \le 0.5$$
   Taking the expected value $E[C(q)]$:
   $$E[C] = \int_{0}^{0.5} C(q)\,f(q)\,dq = \int_{0}^{0.5} (1 - 2q + q^2) \cdot 2\,dq$$
   Evaluate the definite integral step-by-step:
   $$E[C] = 2\left[q - q^2 + \frac{q^3}{3}\right]_{0}^{0.5}$$
   $$E[C] = 2\left[\left(0.5 - (0.5)^2 + \frac{(0.5)^3}{3}\right) - (0 - 0 + 0)\right]$$
   $$E[C] = 2\left[\frac{1}{2} - \frac{1}{4} + \frac{1}{24}\right]$$
   Convert to a common denominator of 24:
   $$E[C] = 2\left[\frac{12 - 6 + 1}{24}\right] = 2\left(\frac{7}{24}\right) = \frac{7}{12} \approx 0.5833 \quad (58.33\%)$$

---

#### (b) Fixed MAF Required for a Baseline Concordance of 0.90

1. **Target Concordance Equation:** Set the HWE homozygous reference frequency equal to the target concordance $C(q) = 0.90$:
   $$(1 - q)^2 = 0.90$$

2. **Step-by-Step Algebraic Solution for $q$:**
   Take the principal square root of both sides (since $q \in (0, 0.5] \implies 1 - q > 0$):
   $$1 - q = \sqrt{0.90}$$
   $$1 - q = \sqrt{\frac{9}{10}} = \frac{3}{\sqrt{10}} = \frac{3\sqrt{10}}{10} \approx 0.949$$
   Solve for the minor allele frequency $q$:
   $$q = 1 - \sqrt{0.90} = 1 - 0.949 = 0.051 \quad (= 5.10\%)$$

---

### 3.5

#### (a) Total Sequencing Yield Required for $1\times$ Mean Coverage

1. **Parameters:**
   * Target organism genome size ($G$) $= 3\text{ Mb} = 3 \times 10^6\text{ bp}$
   * Relative DNA abundance ($r$) $= 0.01\% = \frac{0.01}{100} = 10^{-4}$
   * Desired mean sequencing coverage ($\lambda_{\text{target}}$) $= 1\times$

2. **Required Target-Mapped Sequencing Bases ($B_{\text{target}}$):**
   To achieve $1\times$ mean coverage across a $3\text{ Mb}$ genome, the total sequenced bases originating from the target organism must equal its genome size:
   $$B_{\text{target}} = \lambda_{\text{target}} \times G = 1 \times (3 \times 10^6\text{ bp}) = 3 \times 10^6\text{ bp}$$

3. **Total Metagenomic Sequencing Yield ($Y_{\text{total}}$):**
   Assuming sequencing reads sample DNA proportionally to relative abundance ($B_{\text{target}} = r \times Y_{\text{total}}$):
   $$Y_{\text{total}} = \frac{B_{\text{target}}}{r} = \frac{3 \times 10^6\text{ bp}}{10^{-4}} = 3 \times 10^{10}\text{ bp} = 30\text{ Gb}$$

---

#### (b) Expected Fraction of Target Genome Covered by $\ge 1$ Read at 5 Gb Yield

1. **Target-Mapped Bases at 5 Gb Total Yield:**
   Given total metagenomic yield $Y_{\text{given}} = 5\text{ Gb} = 5 \times 10^9\text{ bp}$:
   $$B_{\text{target}} = r \times Y_{\text{given}} = 10^{-4} \times (5 \times 10^9\text{ bp}) = 5 \times 10^5\text{ bp} \quad (0.5\text{ Mb})$$

2. **Mean Coverage Depth ($\lambda$) of the Target Genome:**
   $$\lambda = \frac{B_{\text{target}}}{G} = \frac{5 \times 10^5\text{ bp}}{3 \times 10^6\text{ bp}} = \frac{1}{6}\times \approx 0.167\times$$

3. **Poisson Breadth of Coverage ($P(K \ge 1)$):**
   Under the Lander–Waterman Poisson read distribution model, the probability that any base receives zero reads ($K = 0$) is $P(K = 0) = e^{-\lambda}$. Therefore, the expected fraction of the genome covered by at least one read ($K \ge 1$) is:
   $$F_{\ge 1} = 1 - P(K = 0) = 1 - e^{-\lambda}$$
   Substitute $\lambda = \frac{1}{6}$:
   $$F_{\ge 1} = 1 - e^{-1/6} \approx 1 - 0.846 = 0.154 \quad (15.40\%)$$

---

### 3.6

#### (a) Detection Power at $60\times$ Sequencing Depth

1. **Binomial Sampling Model (Fixed Depth $N = 60$):**
   At a fixed sequencing depth $N = 60$, each read independently samples the variant allele with probability equal to the variant allele frequency $p = \text{VAF} = 0.05$ and reference allele with $1 - p = 0.95$.
   Thus, the number of alternate reads $X$ follows $X \sim \text{Binomial}(N = 60, p = 0.05)$:
   $$P(X = k) = \binom{60}{k} (0.05)^k (0.95)^{60-k}$$
   Detection requires $X \ge 3$. Using the complement rule:
   $$\text{Power} = P(X \ge 3) = 1 - [P(X = 0) + P(X = 1) + P(X = 2)]$$

   Calculating each failure term step-by-step:
   * **$P(X = 0)$:**
     $$P(X = 0) = \binom{60}{0} (0.05)^0 (0.95)^{60} = (0.95)^{60} \approx 0.046$$
   * **$P(X = 1)$:**
     $$P(X = 1) = \binom{60}{1} (0.05)^1 (0.95)^{59} = 60 \times 0.05 \times (0.95)^{59} \approx 0.138$$
   * **$P(X = 2)$:**
     $$P(X = 2) = \frac{60 \times 59}{2} \times 0.0025 \times (0.95)^{58} = 1770 \times 0.0025 \times 0.051 \approx 0.226$$

   Sum the non-detection probabilities ($\beta$, Type II error rate):
   $$\beta = P(X \le 2) = 0.046 + 0.138 + 0.226 = 0.41$$

   Calculate detection power ($1 - \beta$):
   $$\text{Power} = 1 - 0.41 = 0.59 \quad (59\%)$$

---

#### (b) Expected False-Negative Calls Across 20 Subclonal Variants

1. **Given Parameters:**
   * Total true subclonal variants ($M$) $= 20$
   * Per-variant false-negative probability ($\beta = 1 - \text{Power}$) $= 0.41$

2. **Expectation Calculation:**
   By linearity of expectation across $M = 20$ independent variant loci:
   $$E[\text{FN}] = M \times \beta = 20 \times 0.41 = 8.2\text{ missed variants}$$

---

#### (c) Minimum Sequencing Depth Required for 95% Detection Power

1. **Exact Binomial Threshold Derivation:**
   $$(0.95)^N + N(0.05)(0.95)^{N-1} + \frac{N(N-1)}{2}(0.05)^2(0.95)^{N-2} \le 0.05$$
   Factoring $(0.95)^{N-2}$:
   $$(0.95)^{N-2} \left[(0.95)^2 + 0.05(0.95)N + \frac{0.0025N(N-1)}{2}\right] \le 0.05$$
   $$(0.95)^{N-2} [0.9025 + 0.0475N + 0.00125N^2 - 0.00125N] \le 0.05$$
   $$(0.95)^{N-2} [0.00125N^2 + 0.04625N + 0.9025] \le 0.05$$

2. **Boundary Verification at $N = 123$ vs. $N = 124$:**
   * **Test $N = 123$:**
     $$(0.95)^{121} [0.00125(123)^2 + 0.046(123) + 0.9025] = 0.002 \times [18.91 + 5.68 + 0.9025]$$
     $$\beta(123) = 0.002 \times 25.5 = 0.050 \implies \text{Power} = 1 - 0.050 = 0.95$$
   * **Test $N = 124$:**
     $$(0.95)^{122} [0.00125(124)^2 + 0.046(124) + 0.90] = 0.0019 \times [19.22 + 5.70 + 0.90]$$
     $$\beta(124) = 0.0019 \times 25.82 = 0.049 \implies \text{Power} = 1 - 0.049 = 0.951 \ge 0.95$$
     The condition is strictly satisfied at $N = 124$.

3. **Minimum Depth Result:**
   Using the binomial distribution with $p = 0.05$:
   $$P(X \le 2) = (0.95)^n + n(0.05)(0.95)^{n-1} + \binom{n}{2}(0.05)^2(0.95)^{n-2}$$
   Solving numerically for the smallest integer $n$ satisfying $P(X \le 2) \le 0.05$ gives:
   $$\boxed{n = 124}$$

---

### 3.7

1. **Yield Derivations ($2 \times 150\text{ bp}$ Paired-End, $1\text{ Gb} = 10^9\text{ bp}$):**
   * **Definition A — Read Pairs ($2 \times 150\text{ bp} = 300\text{ bp/unit}$):**
     * Case 1 ($N_1 = 5 \times 10^6$): $Y_{\text{pairs}, 1} = \frac{(5 \times 10^6) \times 300\text{ bp}}{10^9} = 1.5\text{ Gb}$
     * Case 2 ($N_2 = 18 \times 10^6$): $Y_{\text{pairs}, 2} = \frac{(18 \times 10^6) \times 300\text{ bp}}{10^9} = 5.4\text{ Gb}$
   * **Definition B — Total Single Reads ($150\text{ bp/unit}$):**
     * Case 1 ($N_1 = 5 \times 10^6$): $Y_{\text{single}, 1} = \frac{(5 \times 10^6) \times 150\text{ bp}}{10^9} = 0.75\text{ Gb}$
     * Case 2 ($N_2 = 18 \times 10^6$): $Y_{\text{single}, 2} = \frac{(18 \times 10^6) \times 150\text{ bp}}{10^9} = 2.7\text{ Gb}$

2. **Comparison:**
   Both definitions yield gigabase metrics, but interpreting labels as total single reads results in a negative underestimation by an exact factor of $2\times$:
   $$\frac{Y_{\text{single}}}{Y_{\text{advertised}}} = 0.5\times \quad (\text{or } -50\%)$$

---

### 3.8

#### (a) Distinct Possible Diplotypes at Copy Number 2 and Copy Number 3

1. **Combinatorial Model (Combinations with Repetition / Stars and Bars):**
   $$\left(\!\binom{k}{n}\!\right) = \binom{k + n - 1}{n} = \frac{(k + n - 1)!}{n!(k - 1)!}$$

2. **Calculation for Copy Number $n = 2$ ($k = 9$):**
   $$\binom{9 + 2 - 1}{2} = \binom{10}{2} = \frac{10!}{2! \cdot 8!} = \frac{10 \times 9}{2 \times 1} = 45\text{ distinct diplotypes}$$

3. **Calculation for Copy Number $n = 3$ ($k = 9$):**
   $$\binom{9 + 3 - 1}{3} = \binom{11}{3} = \frac{11!}{3! \cdot 8!} = \frac{11 \times 10 \times 9}{3 \times 2 \times 1} = \frac{990}{6} = 165\text{ distinct diplotypes}$$

---

#### (b) Why Total Combinations Do Not Equal $9^n$

* **Order Does Not Matter in Clinical Diplotypes:** Using $9^n$ assumes that each gene copy sits in a specific, labeled slot—for instance, distinguishing the maternal copy from the paternal copy or a tandem duplication. In standard unphased genotyping, however, we only care about which alleles are present, not which chromosome carries them. Because the order makes no biological difference, `(*1, *2)` and `(*2, *1)` represent the exact same diplotype (`*1/*2`), causing $9^n$ to overcount.
* **Uneven Overcounting Prevents Simple Division:** Dividing $9^n$ by $n!$ (the number of ways to arrange $n$ items) fails because different genotypes are overcounted by different amounts. A homozygous genotype like `(*1, *1)` can only be written one way, whereas a heterozygous genotype like `(*1, *2)` can be arranged in $n!$ ways. Dividing by $n!$ does not even yield a whole number (e.g., $\frac{9^2}{2!} = 40.5$ and $\frac{9^3}{3!} = 121.5$).
* **Exact Mathematical Relationship:** Each distinct unordered diplotype with allele counts $(m_1, m_2, \dots, m_k)$ (where $\sum m_i = n$) maps to $\frac{n!}{m_1!\,m_2!\cdots m_k!}$ ordered arrangements. Summing those arrangements across every valid diplotype reconstructs $9^n$:
  $$9^n = \sum_{\text{all diplotypes}} \frac{n!}{m_1!\,m_2!\cdots m_k!} > \binom{9 + n - 1}{n} \quad (\text{for all } n \ge 2)$$

---

#### (c) Distinct Ordered Allele-to-Copy Assignments for `{*2, *2, *41}`

1. **Multinomial Coefficient Derivation:**
   For an unordered multiset of size $n = 3$ containing allele `*2` with multiplicity $m_1 = 2$ and allele `*41` with multiplicity $m_2 = 1$:
   $$N_{\text{ordered}} = \frac{n!}{m_1!\,m_2!} = \frac{3!}{2! \cdot 1!} = \frac{6}{2 \times 1} = 3\text{ distinct ordered assignments}$$

2. **Explicit Enumeration of the 3 Ordered Assignments:**
   * `(*2, *2, *41)` — `*41` on copy 3
   * `(*2, *41, *2)` — `*41` on copy 2
   * `(*41, *2, *2)` — `*41` on copy 1

---

### 3.9

Under HWE, true genotype frequencies for alternate allele frequency $p$ (with reference allele frequency $q = 1 - p$) are:
* **Homozygous reference ($G = 0$):** $P(G = 0) = q^2$
* **Heterozygous ($G = 1$):** $P(G = 1) = 2pq$
* **Homozygous alternate ($G = 2$):** $P(G = 2) = p^2$

Because the uninformative caller assigns a constant dosage of $2p$ to every sample, nearest-integer rounding produces a single identical hard call across the entire cohort. Hard-call concordance equals the true HWE frequency of that predicted genotype.

1. **Alternate Allele Frequency $p = 0.01$**
   * **Assigned dosage:** $2p = 2(0.01) = 0.02$.
   * **Derived hard call:** Rounding $0.02$ to the nearest integer yields a hard call of $0$ (homozygous reference) for every sample.
   * **Hard-call concordance:** The proportion of samples matching the predicted genotype $G = 0$ is:
     $$q^2 = (1 - 0.01)^2 = (0.99)^2 = 0.9801 \quad (98.01\%)$$

2. **Alternate Allele Frequency $p = 0.30$**
   * **Assigned dosage:** $2p = 2(0.30) = 0.60$.
   * **Derived hard call:** Rounding $0.60$ to the nearest integer yields a hard call of $1$ (heterozygous) for every sample.
   * **Hard-call concordance:** The proportion of samples matching the predicted genotype $G = 1$ is:
     $$2pq = 2(0.30)(0.70) = 0.4200 \quad (42.00\%)$$

---

### 3.10

Doctors and researchers use epigenetic clocks to calculate a person's biological age versus their chronological age. If biological age is significantly higher than chronological age—referred to as **age acceleration**—it signals an elevated risk for age-related conditions like cardiovascular disease, neurodegeneration, and all-cause mortality.

An epigenetic clock estimates biological age as a weighted linear combination of methylation beta ($\beta$) values across its target CpG probes:
$$\text{Predicted Age} = b_0 + \sum_{i=1}^{68} w_i \beta_i$$

This is the foundational linear model introduced by Steve Horvath (2013) and Gregory Hannum et al. (2013), where every first-generation epigenetic clock calculates age as an intercept ($b_0$) plus a weighted sum of methylation beta values ($\beta_i$) multiplied by their learned regression weights ($w_i$). When 14 probes are absent from a newer array, their contribution to the linear predictor is altered, introducing a systematic shift ($\text{Bias} = \text{Observed} - \text{True}$) across the cohort.

#### (a) Systematic Prediction Bias from Probe Exclusion

Excluding the 14 missing probes effectively sets their methylation values to zero ($\beta_i = 0$) in the linear model. The true average contribution of these 14 probes in the cohort is given by the product of their summed weights ($\sum w_{\text{missing}} = -18.6$) and their mean methylation level ($\bar{\beta} = 0.62$):
$$\text{True Contribution} = (-18.6) \times 0.62 = -11.53$$

Because dropping the probes omits this $-11.53\text{ year}$ contribution from the final sum, the predicted ages are systematically shifted upward:
$$\text{Bias} = \text{Predicted}_{\text{excluded}} - \text{Predicted}_{\text{true}} = 0 - (-11.53) = +11.53\text{ years}$$

Simply excluding the missing probes causes the clock to overestimate biological age by an average of **11.53 years**.

---

#### (b) Prediction Bias Using Constant Beta Imputation ($\beta = 0.50$)

Instead of treating the missing probes as zero, imputing a constant $\beta = 0.50$ replaces their true average contribution with an imputed baseline:
$$\text{Imputed Contribution} = (-18.6) \times 0.50 = -9.30$$

Comparing this imputed contribution against the true cohort contribution ($-11.53$) gives the remaining systematic bias:
$$\text{Bias} = \text{Imputed Contribution} - \text{True Contribution} = -9.30 - (-11.53) = +2.23\text{ years}$$

Imputing a constant beta of $0.50$ reduces the distortion substantially, though it still overestimates biological age by an average of **2.23 years**.

---

#### (c) Preferred Strategy and Required Reporting Disclosures

* **Methodologically Preferred Strategy:**
  Rather than dropping probes or using a flat $0.50$ constant, the preferred strategy is **reference-based or population-specific mean imputation** (such as using the original training dataset's probe-specific mean $\beta$ values, $k$-nearest neighbors [$k\text{NN}$] imputation, or specialized methylation imputation tools like *SeSAMe*). Because CpG sites are typically either heavily methylated or unmethylated rather than sitting at $0.50$, probe-specific reference imputation preserves the baseline intercept and biological scaling of the original model without distorting age acceleration metrics.

* **Required Supplementary Reporting Disclosures:**
  To ensure reproducibility and allow accurate cross-platform interpretation, the following details must be reported in the methods or supplementary materials:
  1. **Probe Overlap & Missingness:** The exact number and proportion of model probes present versus missing on the new array (e.g., 54 of 68 probes retained; 20.59% missing), along with the specific CpG identifiers.
  2. **Summed Weight Impact:** The cumulative weight of the missing probes (e.g., $-18.6$) to clarify how heavily the missing features influence the model's baseline offset.
  3. **Imputation Methodology:** The exact algorithm, reference dataset, or constant values used to handle missing probes prior to clock calculation.
  4. **Sensitivity & Benchmarking Analysis:** A comparison of predicted ages and epigenetic age acceleration residuals on a subset of samples measured across both array versions (or simulated probe dropout on the original array) to quantify any residual offset or loss of variance.

> **References:**
> * Dhingra, R., Kwee, L. C., Diaz-Sanchez, D., Devlin, R. B., Cascio, W., Hauser, E. R., ... & Ward-Caviness, C. K. (2019). Evaluating DNA methylation age on the Illumina MethylationEPIC bead chip. *PLOS ONE*, 14(4), e0207834.
> * Hannum, G., Guinney, J., Zhao, L., Zhang, L. I., Hughes, G., Sadda, S., ... & Zhang, K. (2013). Genome-wide methylation profiles reveal quantitative views of human aging rates. *Molecular Cell*, 49(2), 359–367.
> * Horvath, S. (2013). DNA methylation age of human tissues and cell types. *Genome Biology*, 14(10), R115.
> * Zhou, W., Triche, T. J., Laird, P. W., & Shen, H. (2018). SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions. *Nucleic Acids Research*, 46(20), e123.