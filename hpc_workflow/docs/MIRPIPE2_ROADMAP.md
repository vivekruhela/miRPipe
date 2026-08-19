# miRPipe2: proposed scientific and engineering roadmap

## Executive direction

miRPipe2 should not be only a faster wrapper around miRDeep\*. Its core scientific contribution should be an **evidence-calibrated small-RNA identity engine** that separates four questions:

1. What sequence and genomic locus generated the reads?
2. Does the locus show canonical miRNA or piRNA biogenesis?
3. Is it known, a paralogue/isomiR/editing event, or a truly novel candidate?
4. Is it reproducible, disease-relevant, and functionally supported?

The current HPC workflow solves the execution and provenance layer. The following phases address the biological gaps.

## Current gaps and proposed changes

| Gap in miRPipe v1 | Consequence | miRPipe2 direction |
|---|---|---|
| One primary miRNA caller | Caller-specific false positives and false negatives | Optional ensemble of miRDeep\*, miRDeep2, miRge3.0, and sRNAbench; retain caller-level evidence |
| Selenium/DASHR lookup | Fragile, slow, non-versioned, difficult to reproduce | Local indexed snapshots of miRBase, MirGeneDB, RNAcentral/DASHR-compatible exports, with checksums and release metadata |
| Seed clustering dominates reannotation | Same seed does not prove the same locus, precursor, or function | Hierarchical identity model using sequence, strand, locus, precursor, arm, seed, and processing evidence |
| Sequential names for novel candidates | Names may be mistaken for official miRBase accessions | Stable provisional IDs; official nomenclature only after database acceptance |
| Limited contamination model | rRNA/tRNA/snoRNA/Y-RNA fragments can appear novel | Ordered competitive annotation against rRNA, tRNA/tRF, snoRNA, snRNA, Y-RNA, repeats, and microbial contaminants |
| No UMI/isomiR/editing model | PCR duplicates and biologically meaningful variants are conflated | Protocol-aware UMI extraction, mirGFF3 output, 5′/3′ isomiRs, non-templated addition, and A-to-I editing |
| piRNA = length + database overlap | Somatic degradation fragments may be overcalled | Evidence for clusters, 1U/10A, 10-nt ping-pong overlap, phasing, strand, transposons, and PIWI expression |
| Two-group DE only | Batch, pairing, ancestry, sex, RIN, and clinical covariates are ignored | Arbitrary designs, contrasts, paired/repeated designs, independent filtering, effect-size shrinkage, and QC |
| Accuracy evaluated mainly by synthetic reads | Simulations may not reproduce all library and mapping artifacts | Spike-ins, curated positives, hard negatives, cross-study validation, and independent simulated seeds |

## A biologically driven miRNA evidence matrix

Use one row per **candidate locus × mature arm**, not only one row per observed sequence.

| Feature block | Example features |
|---|---|
| Read support | total and unique reads, CPM, number of samples, condition recurrence, library complexity |
| Processing precision | 5′ homogeneity, 3′ dispersion, dominant mature fraction, loop-read depletion |
| Duplex evidence | mature/star support, 2-nt 3′ overhang, arm assignment, duplex pairing |
| Structure | precursor MFE, MFEI, paired bases in mature region, bulges, loop length, ensemble diversity |
| Mapping | number of loci, strand consistency, repeat overlap, mappability, soft-clipping/mismatch pattern |
| Caller evidence | each caller's score and pass/fail, consensus count, caller disagreement |
| Annotation | exact/reverse-complement catalogue match, MirGeneDB status, Rfam class, other ncRNA overlap |
| Evolution | phyloP/phastCons, synteny, homologous hairpin, lineage specificity |
| Functional | AGO CLIP support, target-site evidence, miRNA–mRNA anti-correlation, perturbation evidence |
| Technical | adapter confidence, UMI deduplication rate, batch, read-length peak, contamination fraction |

MirGeneDB's curation criteria are especially valuable: expression from both hairpin arms, homogeneous 5′ ends, imperfect precursor complementarity, and the expected Drosha/Dicer offsets should be encoded explicitly rather than approximated by seed similarity alone.

## Candidate ranking model

Start with an interpretable rule-based model, then compare it with calibrated machine learning. A useful conceptual score is:

\[
S = f(E_{processing}, E_{duplex}, E_{structure}, E_{mapping}, E_{reproducibility}, E_{conservation}, E_{functional}) - P_{contamination}
\]

Recommended labels:

- **Positive**: high-confidence MirGeneDB loci, held out by miRNA family and genomic locus.
- **Hard negative**: expressed tRNA/rRNA/snoRNA/Y-RNA/repeat fragments matched for length, GC, abundance, and mappability.
- **Unlabeled**: unvalidated novel predictions; do not automatically treat them as negatives.

Avoid random read-level train/test splits because reads from the same precursor leak nearly identical sequence and structure into both sets. Prefer family-held-out, chromosome-held-out, study-held-out, and species-held-out evaluation. Report calibrated precision-recall, not accuracy alone, because true novel miRNAs are rare.

## Better novel-miRNA identification

1. Collapse identical reads after optional UMI correction to reduce compute and PCR bias.
2. Perform ordered competitive annotation to remove known small-RNA fragments before discovery.
3. Discover candidate loci jointly across all samples, then quantify the fixed catalogue per sample. This avoids giving deep samples more candidate opportunities.
4. Record caller evidence separately and create consensus tiers rather than requiring blind intersection.
5. Require canonical biogenesis evidence: precise 5′ processing, plausible hairpin, mature/star duplex and overhang, and cross-sample support.
6. Distinguish:
   - known mature miRNA;
   - isomiR or edited miRNA;
   - known-family paralogue at a new locus;
   - conserved novel locus;
   - lineage-specific novel candidate;
   - other-small-RNA fragment/rejected candidate.
7. Validate top candidates with AGO-IP/CLIP where available and orthogonal RT-qPCR or Northern blot for a focused set.

## Better piRNA identification

The existing overlap-based branch should be renamed **piRNA annotation**. A high-confidence piRNA branch should add:

- de novo piRNA-cluster discovery and cluster-level quantification;
- 1U and 10A nucleotide-bias summaries;
- enrichment of 10-nt 5′ overlaps between opposite-strand read pairs (ping-pong Z-score);
- phased primary piRNA signatures;
- transposable-element family and orientation assignment;
- locus mappability and explicit multi-mapping treatment;
- tissue-aware PIWI-pathway gene expression evidence;
- replication of signatures across samples, not only pooled reads.

Return separate labels such as `database_overlap`, `cluster_supported`, `ping_pong_supported`, and `high_confidence_biogenesis` rather than one binary piRNA call.

## Biological prioritization after identification

For each robust differential miRNA:

1. Combine conserved sequence prediction with experimentally supported AGO/target databases.
2. If matched mRNA-seq is available, prioritize negatively correlated miRNA–mRNA pairs while adjusting for phenotype and batch.
3. Test pathway enrichment on the supported target set, not all predicted targets.
4. Add disease/tissue context, survival or treatment response, and external-cohort replication.
5. Where genotype data exist, integrate miRNA-eQTL, target-eQTL, colocalization, and mediation as downstream evidence—not as proof of miRNA identity.

## Speed and scale

- Keep the current sample-level SLURM scatter/gather and Nextflow cache.
- Collapse reads early and quantify unique sequences once per cohort.
- Build each reference/index once, checksum it, and share it read-only across runs.
- Discover loci on pooled/collapsed evidence, then quantify the locked catalogue in parallel.
- Use adaptive resources from prior trace files; retry only memory/time failures.
- Avoid FASTQ splitting unless a single file exceeds the efficient range of a node; scheduler-level sample parallelism is simpler and safer.

## Reproducibility and software quality

1. Pin Nextflow, containers by digest, reference releases, and reference SHA-256 values.
2. Add `nf-test` unit and integration tests plus a tiny permissible FASTQ truth fixture.
3. Run continuous integration for schema validation, Python/R tests, linting, container build, and a stub workflow.
4. Emit RO-Crate or equivalent machine-readable provenance and a methods paragraph generated from run parameters.
5. Version the output schema and use mirGFF3 for interoperable isomiR reporting.
6. Never modify or delete user inputs; retain a manifest of every input checksum.

## Suggested releases

| Release | Deliverable | Exit criterion |
|---|---|---|
| `2.0-alpha` | Current Nextflow/SLURM conversion | Real-data run completes and reproduces v1 count matrices within explained differences |
| `2.0-beta` | Local reference bundle, UMI/isomiR, contamination hierarchy, multi-caller evidence | Locked synthetic + real truth benchmarks pass in CI |
| `2.0-rc` | Biogenesis scoring, piRNA signatures, report/dashboard | Independent external cohorts preserve calibrated precision |
| `2.0` | Versioned container/reference release and documentation | Reproducible archived run with DOI and complete provenance |

## Primary references and implementations

- Ruhela V. et al. miRPipe. *Frontiers in Bioinformatics* (2022). https://doi.org/10.3389/fbinf.2022.842051
- Nextflow cache and resume documentation. https://www.nextflow.io/docs/latest/cache-and-resume.html
- nf-core/smrnaseq 2.4.1. https://nf-co.re/smrnaseq/2.4.1/
- MirGeneDB structural criteria. https://mirgenedb.org/information
- Patil A.H. et al. miRge3.0. *Nucleic Acids Research* (2021). https://pmc.ncbi.nlm.nih.gov/articles/PMC8294687/
- Pawlina-Tyszko K. et al. Benchmarking small-RNA tools (2023). https://pmc.ncbi.nlm.nih.gov/articles/PMC10687144/
