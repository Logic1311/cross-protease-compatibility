# Censor-aware compatibility budgets for paired protein-folding assays

**Research manuscript, version 1.0 | 26 November 2026 | AI-assisted research; human authorship pending; not peer reviewed**

## Abstract

Thermodynamic coupling in double-mutant landscapes is well established, but inference from incomplete, range-limited measurements requires distinguishing nonadditivity from assay disagreement. We develop an exact deterministic diagnostic for two interval-valued landscapes: the minimum combined additional tolerance needed for a common unrestricted landscape, allowing additive row and column offsets between assays, is the minimax additive defect of their interval difference. The common additive-landscape budget instead equals the sum of the two individual additive defects. A negative mean cycle in a bipartite difference-constraint graph certifies the individual defects, including one-sided bounds and missing cells. Applying these established graph principles to paired trypsin and chymotrypsin estimates in the published MegaScale data, we independently reconstruct 210,118 records in 559 panels, exclude 434 duplicated-coordinate records, and retain 209,684 unique records. In 549 of 559 panels, a nonempty equal-tolerance interval permits a shared latent landscape while each assay separately excludes every additive landscape. Its median width is 0.505 kcal mol^-1. At the unexpanded interval setting, a cycle chosen exclusively from trypsin remains positive in chymotrypsin in 428 panels, with 488 assessable out of 554 positive source witnesses; the reverse direction gives 429 positive transfers among 460 assessable of 556. A direct joint optimization checks 116 panel/mask combinations, and integer certificates check 3,354 individual or difference optimizations. These are model-conditional compatibility statements, not simultaneous confidence tests or claims of newly discovered molecular mechanisms.

## 1. Introduction

Protein stability depends on both individual amino-acid effects and interactions. Double-mutant cycles are standard tools for quantifying coupling [1]. The MegaScale study [2,3] measured selected site pairs exhaustively using cDNA display proteolysis, including 559 panels over 190 protein domains. It already reported thermodynamic couplings, additive-model residuals and comparisons between trypsin and chymotrypsin. Subsequent studies evaluated predictive epistasis models [4], characterized energetic architectures [5], and developed reference-free decomposition with nonlinear phenotypes [6]. Minimum-epistasis interpolation [7] addresses unmeasured genotypes through a different optimization objective.

Here the question is narrower: given incomplete, interval-valued observations from two assays, what additional tolerance is needed to make the assays compatible with one another, and at what tolerance do they no longer exclude additivity? An explicit incompatibility witness can also be transferred between assays without reselecting its amino-acid states on the target. The complete-data uniform-approximation problem is classical [8], and directed minimum mean cycles are classical [9]. We do not claim a new graph algorithm or the discovery of epistasis. The contribution is the paired-assay budget identity, its one-sided interval formulation, and an auditable empirical reanalysis of the published paired-protease data.

## 2. Mathematical formulation

### 2.1 Individual interval-valued additive defect

Let I be a finite rectangular array of nonempty closed intervals I_ij = [L_ij,U_ij], allowing L_ij = -infinity or U_ij = +infinity. An unobserved entry is the whole real line. The additive subspace A consists of arrays a_i + b_j. For x in the real line, dist(x,I_ij) denotes ordinary distance to the interval. Define the additive defect

$$
D(I)=\min_{a,b}\max_{i,j}\operatorname{dist}(a_i+b_j,I_{ij}).
$$

This is a minimum extra uniform tolerance, not a standard deviation. For a fixed epsilon >= 0, introduce vertex potentials p(A_i) = a_i and p(B_j) = -b_j. A finite lower bound gives the edge A_i to B_j with weight -L_ij; a finite upper bound gives B_j to A_i with weight U_ij. Add epsilon to every existing edge. The inequalities p(v)-p(u) <= w(u,v)+epsilon are exactly the expanded-interval constraints.

Proposition 1. If the directed graph has a cycle, then

$$
D(I)=\max\left(0,-\min_C\frac{\sum_{e\in C}w(e)}{|C|}\right).
$$

For an acyclic graph D(I)=0. Here C ranges over directed simple cycles. Proof: summing the potential inequalities around C makes their left side zero, so its total expanded weight must be nonnegative. Conversely, if all directed cycles have nonnegative weight, append a super-source with zero-weight edges to all vertices. Shortest-path distances are finite and satisfy each potential inequality. A negative closed walk contains a negative simple cycle, so checking simple cycles suffices. The finite collection of inequalities is closed at the minimum feasible epsilon; the graph argument constructs an attained solution. The result recovers classical uniform additive approximation on exact arrays [8,9] and handles missing and one-sided bounds by omitting the relevant edges.

For any finite directed cycle define h(C;I) as minus its mean weight. A positive value is a certificate: no additive landscape fits every used bound at a tolerance smaller than h. A four-edge cycle recovers the interval form of a double-mutant contrast divided by four. The best four-edge certificate is a lower bound on D, not necessarily equal to D. For the exact 3-by-3 array with rows (0.5,1,0), (0,0.5,1), (1,0,0.5), the full optimum is 0.5 whereas the largest quartet bound is 0.375. Longer cycles refer to amino-acid alternatives at two sites, not interactions among more than two residues.

### 2.2 Exact paired-assay tolerance theorem

Let T and C be two interval arrays on the same sites. Allow each assay an unknown additive offset and its own nonnegative uniform extra tolerance e_T and e_C. The common latent landscape G is initially unrestricted. After absorbing the first additive offset into G, the two constraints can be written

$$
G_{ij}\in T_{ij}^{e_T},\quad G_{ij}+a_i+b_j\in C_{ij}^{e_C},
$$

where the superscript denotes interval expansion by the specified tolerance. Define the cellwise difference J=C-T, i.e., [L_C-U_T,U_C-L_T]. Infinite lower and upper endpoints are interpreted in the extended-real sense; unobserved intervals are unconstrained. Define D_T=D(T), D_C=D(C), and D_J=D(J).

Theorem 2. A common unrestricted G and additive offset exist if and only if

$$
e_T+e_C\geq D_J.
$$

A common additive G exists if and only if e_T >= D_T and e_C >= D_C.

Proof. For a given additive offset A_ij=a_i+b_j, a feasible G_ij is precisely an element of the intersection T_ij expanded by e_T and (C_ij expanded by e_C) shifted by -A_ij. This intersection is nonempty exactly when L_C-U_T-(e_C+e_T) <= A_ij <= U_C-L_T+(e_C+e_T), with absent endpoints imposing no inequality. Thus selecting a common G is equivalent to fitting an additive array to J with total tolerance e_T+e_C. Cellwise selection gives sufficiency even with infinite endpoints. If G is additive, both predicted arrays are additive. Conversely, choose any feasible additive prediction for each assay separately; their difference is additive, so they admit a shared additive G in the chosen gauge. Proposition 1 ensures attainment. This proves both statements.

Consequently the minimal total unrestricted budget is D_J, while the minimal total additive budget is D_T+D_C. In particular, D_J <= D_T+D_C. With equal per-assay tolerance eta, the interval

$$
\left[D_J/2,\ \min(D_T,D_C)\right)
$$

is precisely the set of eta for which a common unrestricted landscape is feasible and neither assay admits any additive landscape. The width W is max(0,min(D_T,D_C)-D_J/2). The weaker condition that a common additive landscape is impossible uses max rather than min. The width is an identifiable diagnostic of the recorded intervals, not a measurement of the true experimental error. Sign reversal of J leaves its additive defect invariant.

If a common G is feasible at (e_T,e_C), its distance to the additive subspace in the full-array sup norm is at least max(0,D_T-e_T,D_C-e_C). Indeed, adding an additive offset does not change this distance; the triangle inequality and cellwise tolerance bounds imply each individual inequality. No actual biochemical interaction energy is identified without the observation model.

### 2.3 Fixed-witness transfer

Select a positive optimal cycle using one assay only. Hold its ordered vertices and direction fixed, then compute h on the other assay. This contrast is invariant to additive shifts in the target assay. If every required target endpoint is finite, a positive transferred score provides a target-side incompatibility certificate. A missing required endpoint makes the witness unassessable; it is not a negative result and is never silently imputed. Selection is not repeated on the target. Compare the transferred full cycle with the transferred best source quartet selected by the same policy. If a source optimum is nonunique, the solver's dual-support traversal chooses a deterministic witness in the tested software environment. This is cross-assay validation on the same sequences, not an independent test on other proteins.

## 3. Data and analysis

Data source: Tsuboyama et al. [2], Zenodo record 7992926 version v2_230420 [3]. The 1,013,678,473-byte Processed_K50_dG_datasets.zip matches the published MD5 f7e8c553efee734cf161ee6f2b0a09cf. The processed table contains assay-specific deltaG_t and deltaG_c estimates and their lower and upper 95% interval fields. We use these published model estimates, not a new re-fit of raw sequencing counts.

We select rows with a nonempty pair_name. Parse mut_type into the two one-indexed positions and their substituted amino acids. Check against each aa_seq, require a constant remaining sequence background within each panel and a constant ordered position pair. Exclude all rows at duplicated panel/state coordinates rather than choosing a replicate on its outcome. This yields 209,684 unique observations from 210,118 records across 559 panels. The independent reconstruction of the earlier combined-quality analysis reproduces 207,381 combined observations. All panels pass the sequence/background audit.

For each assay separately, keep finite estimates with valid finite ordered interval endpoints. The source paper describes a main reliable range of -1 to 5 kcal mol^-1 and clips estimates for its figures [2]. Our deliberately conservative construction widens each source interval to contain the estimate clipped to this range. For an estimate at or below -1 discard its lower bound; at or above 5 discard its upper bound. This differs from the source's combined dG_ML categorical rule, which is a quality annotation; it does not imply all retained measurements are of equal reliability. Missing or invalid records impose no bounds. For sensitivity we additionally expand every finite endpoint outward by 0.25 or 0.5 kcal mol^-1, and repeat with the published combined-quality mask (exclude dG_ML='-'). This latter mask uses both assays and is therefore not independent target-blind selection. We also report a high-coverage subset (at least 350 retained states) and a subset with at least 250 jointly two-sided observations. Reported interval widths reflect sequencing-count/model uncertainty and not all systematic errors [2].

Linear programming on 40 potentials and one tolerance calculates each optimum. The negative dual-support cycle provides a witness, and direct enumeration evaluates every 4-cycle. The paired criterion uses the same operation on the interval difference. A separate direct LP in 400 latent-cell variables, 40 additive offset parameters and two nonnegative tolerances checks the combined-budget identity on 58 panels selected by a SHA-256 panel-name rule, for both masks (116 tests). Independent integer verification rounds every finite constraint weight outward at 10^-9 kcal mol^-1, recovers an exact rational mean-cycle value and verifies a shortest-path feasible potential using Python integers. The outward rounding weakens constraints by approximately one grid unit; the exact rounded and original floating-point optima are compared. The accompanying Lean 4.28 proof verifies the interval intersection and difference criterion for integer endpoints (arbitrary scaled integer units), including unbounded endpoints, and its cellwise landscape extension. The full real-variable proof is given above. No Lean verification of the empirical data, LP implementation, or complete graph theorem is claimed.

## 4. Results

### 4.1 Shared compatibility and separate nonadditivity

| Subset | Panels | Positive W | Median W |
|---|---:|---:|---:|
| All unique observations | 559 | 549 | 0.5047 |
| At least 350 states | 518 | 516 | 0.5145 |
| At least 350 states and 250 jointly two-sided | 396 | 396 | 0.4894 |
| Combined-quality mask | 559 | 549 | 0.5022 |
| Combined mask, at least 350 states | 510 | 508 | 0.5118 |
| Combined mask, at least 350 states and 250 jointly two-sided | 392 | 392 | 0.4894 |

Table 1. Positive equal-tolerance compatibility windows and median widths (kcal mol^-1). A positive width is defined numerically as greater than 10^-8; these counts describe the measured panels, not prevalence estimates over all proteins. Subset selection is based on data coverage, and the combined mask is assay-conditioned.

In the full unique-record analysis, medians D_T, D_C and D_J are 0.7691, 0.7145 and 0.2913 kcal mol^-1, respectively. The median required common per-assay budget D_J/2 is 0.1456. Window widths exceed 0.25 in 487 panels and 0.5 in 287 panels. This separation is preserved by the stricter coverage/quality subsets above. The median of the window is computed per panel; it is not the difference between independently reported medians.

![Figure 1. Paired-assay budget endpoints for all 559 panels. Points above the diagonal have a positive window; histogram includes zero-width panels. Units: kcal/mol.](../figures/agreement.png)

### 4.2 Transfer of a fixed cycle

| Source -> target | Extra margin | Positive source | Assessable | Positive transfer |
|---|---:|---:|---:|---:|
| Trypsin -> chymotrypsin | 0 | 554 | 488 | 428 |
| Chymotrypsin -> trypsin | 0 | 556 | 460 | 429 |
| Trypsin -> chymotrypsin | 0.25 | 551 | 485 | 368 |
| Chymotrypsin -> trypsin | 0.25 | 539 | 448 | 380 |
| Trypsin -> chymotrypsin | 0.5 | 487 | 422 | 244 |
| Chymotrypsin -> trypsin | 0.5 | 433 | 366 | 265 |

Table 2. All 559 panels, all unique observations. Extra margins widen both assays' finite bounds before discovery and evaluation. Each row selects only positive source witnesses, so denominators vary. Unassessable cycles are not counted as confirmations.

Without extra margins, the median transferred score among assessable positive-source witnesses is 0.5115 for trypsin to chymotrypsin and 0.5720 for the reverse. Confirmation fractions among assessable cases are 87.7% and 93.3%, respectively. Relative to all 559 panels they are 76.6% and 76.7%. Among 518 panels with at least 350 states, 416 and 419 confirm in the two directions, with 460 and 439 assessable out of 516 and 518 positive source cycles. At a source-to-target margin of 0.5, the all-panel counts fall to 244 and 265, so the finding depends quantitatively on the tolerated error.

A longer-cycle advantage over the best four-cycle in the source (>0.01 kcal mol^-1) occurs in 105 trypsin-source panels and 80 chymotrypsin-source panels. Of these, 69 and 53 transfer positively; 45 and 28 retain a >0.01 advantage over the same source-selected quartet on the target. The latter comparison requires both witnesses to be assessable. This is not a demonstration that every longer witness beats every target-selected quartet.

For an illustrative real panel, 1A0N.pdb_hnet0, the trypsin-selected 6-cycle certifies 1.1201 kcal mol^-1 versus the source best quartet 1.0454. The fixed 6-cycle transfers to chymotrypsin at 1.2689 versus the transferred source quartet at 1.1034. This is an example, not a preselected representative sample.

![Figure 2. Fixed source-witness scores versus target scores at zero added margin, with missing required target endpoints omitted from the scatter. Negative scores are retained.](../figures/transfer.png)

### 4.3 Verification

All 116 sampled direct joint-LP values agree with the reduced interval-difference optimum to floating-point precision (maximum difference 3.33 x 10^-16 kcal mol^-1). The 3,354 integer certificate checks cover three optima per panel for each of the two masks. The maximum difference from the floating-point reported values is 9.80 x 10^-10 kcal mol^-1. The numerical checks test the theorem's application to these inputs; they are not a substitute for the proof. Unit tests cover exact additivity, a strict longer-cycle advantage, invariance to additive shifts, identical and antagonistic paired landscapes, and exact rounded-constraint feasibility.

## 5. Scope, limitations, and relation to prior work

The source work [2] established nonadditive thermodynamic coupling in these panels and used the two protease conditions for reliability analysis. Dieckhaus and Kuhlman [4] showed that current models fail to capture double-mutant epistasis reliably; their predictive benchmark and our deterministic tolerance budget answer different questions. Faure et al. [5], Park et al. [6] and Zhou and McCandlish [7] address genetic architecture, nonlinear phenotype mapping or interpolation. Classical approximation and graph algorithms [8,9] underlie Proposition 1; novelty is not claimed for these facts. Targeted literature searches and public code search did not identify the specific paired interval-difference compatibility budget and fixed-witness census in this dataset. This does not establish priority over all publications or inaccessible material.

The data are selected site pairs, many chosen using prior structural/interaction information, rather than a random sample of all pairs. Several panels share domains and backgrounds. Counts and medians are descriptive; no independent-panel confidence interval, p-value, causal inference or population-level prevalence claim is made. The reported intervals are individual model uncertainty summaries, not simultaneous coverage of an entire landscape. A positive cycle excludes additivity within the stated bounds, conditional on the observation model. It does not prove that all observed bounds simultaneously contain true thermodynamic values. Two protease measurements share experimental material, sequence-dependent modeling assumptions and possible systematic errors. Shared nonadditive bias can remain even when the difference landscape is additive. Thus a wide compatibility window is not evidence that actual errors necessarily lie in it, and cross-assay agreement is not independent biochemical confirmation of a mechanism. Extra margins are transparent sensitivity parameters, not calibrated error estimates.

The paired criterion remains useful precisely because it separates two questions: whether the intervals admit one landscape up to additive offsets, and whether they admit an additive one. The fixed witnesses are short lists of observed states whose signed bounds can be checked directly and prioritized for follow-up measurements. Independent purified-protein experiments or different measurement modalities would be required to claim mechanistic or broadly transferable physical conclusions.

## 6. Data, code and disclosures

Original measurements and model outputs belong to the source authors [2,3]; cite them when reusing this work. The package includes executable analysis, complete derived panel summaries and witness tables, exact-certificate results, source audit, Lean source, figures, manuscript and file hashes. The original 1 GB archive is not redistributed. See README for pinned input and reproduction commands. The manuscript and analysis were generated and reviewed with AI assistance; no human author or affiliation is asserted here. This is an unreviewed research preprint, not a journal publication.

## References

[1] Pagano L, Toto A, Malagrino F, Visconti L, Jemth P, Gianni S. Double mutant cycles as a tool to address folding, binding, and allostery. Int J Mol Sci 22, 828 (2021). doi:10.3390/ijms22020828.

[2] Tsuboyama K et al. Mega-scale experimental analysis of protein folding stability in biology and design. Nature 620, 434-444 (2023). doi:10.1038/s41586-023-06328-6.

[3] Tsuboyama K et al. Mega-scale experimental analysis of protein folding stability in biology and design. Zenodo data record (2023). doi:10.5281/zenodo.7992926.

[4] Dieckhaus H, Kuhlman B. Protein stability models fail to capture epistatic interactions of double point mutations. Protein Science 34, e70003 (2025). doi:10.1002/pro.70003.

[5] Faure AJ et al. The genetic architecture of protein stability. Nature 634, 995-1003 (2024). doi:10.1038/s41586-024-07966-0.

[6] Park Y, Metzger BPH, Thornton JW. The simplicity of protein sequence-function relationships. Nature Communications 15, 7953 (2024). doi:10.1038/s41467-024-51895-5.

[7] Zhou J, McCandlish DM. Minimum epistasis interpolation for sequence-function relationships. Nature Communications 11, 1782 (2020). doi:10.1038/s41467-020-15512-5.

[8] Diliberto SP, Straus EG. On the approximation of a function of several variables by the sum of functions of fewer variables. Pacific J Math 1, 195-210 (1951). doi:10.2140/pjm.1951.1.195.

[9] Karp RM. A characterization of the minimum cycle mean in a digraph. Discrete Mathematics 23, 309-311 (1978). doi:10.1016/0012-365X(78)90011-0.
