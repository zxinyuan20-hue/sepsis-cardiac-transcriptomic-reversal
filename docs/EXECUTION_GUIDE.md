# Execution stages and dependency boundaries

Read each module docstring and its associated configuration before running. Historical release modules verify stage-specific artifacts; they are not prerequisites for installing or reading the code. Keep protocol seeds and family definitions unchanged.

| Stage | Modules | Input and output role | Private-data dependence |
|---|---|---|---|
| Public preparation | 00–07 | GEO/LINCS/reference acquisition, dictionaries and input QC | Public; large downloads |
| H9c2 import and preprocessing | 08, 10_rat_qc.R, 12 | Original archives and explicit C/L sample contract; one-to-one orthology | Requires unavailable private input |
| Human discovery | 09, 10, 11, 13 | Raw GSE79962 RMA; discovery/comparator models; 100-UP/90-DOWN query | Human analysis public; some stage audits also inspect rat files |
| LINCS screening | 14–19 | Landmark profiles; disease bootstraps; null family; all 482 entries | Public after human preprocessing |
| Cross-model context | 20–23, 27 | Separate mouse/H9c2 effects and ortholog models | Mixed; rat branch requires private input |
| Cardiac perturbation | 24–26, 28–35 | GSE217421 identity, donors and fixed 21-drug family | Public |
| Disease specificity | 36–38 | Sepsis versus IHD/DCM disease-weighted drug effects | Public |
| Literature context | 39–42, 57–58 | Public bibliographic evidence and manual-review records | No raw private input; historical context |
| External human cohorts | 43–45 | GSE141864 and initial GSE237861 acquisition/analysis | Public; initial partial GSE237861 stage is superseded |
| Complete human reconstruction | 46–51 | All 14 GSE237861 heart FASTQ; GENCODE v47, Subjunc/featureCounts | Public; compute/storage intensive; authoritative complete stage |
| External drug reassessment | 52–56 | Fixed 21 drugs with two disease backgrounds and donor/patient omissions | Public |
| Inferential resolution | 69 | Exact sign-flip/BH attainability | Public / mathematical audit |
| Robustness extensions | 77–82 | Structural identity, alternative scores and toxicity/magnitude | Public; post-result extensions |
| Discovery technical sensitivity | 86–89 | Clinical annotation and CEL scan-date sensitivity | Public |
| Hallmark concordance | 92 | Six effect tables and licensed Hallmark GMT; common/full backgrounds | Requires H9c2 effect table; cannot reproduce complete six-model matrix publicly |
| Adult human heart localization | 93–97, 100 | Public adult heart reference, donor/cell annotations, frozen human query | Public; large reference or range-based reads |
| Plotting archive | 60, 64–65, 82, 96 | Frozen upstream outputs and figure-specific source tables | Mixed where H9c2 shown; requires upstream artifacts |

## Important ordering exceptions

1. R environment setup modules 03 and 05 accompany acquisition rather than alphabetical execution. Numbered .R files are wrapper dependencies.
2. Module 51 orchestrates the complete raw-heart pipeline. Avoid simultaneously running both the orchestrator and its child modules. Inspect its stage checks before restarting.
3. Module 92 is the post-result Hallmark extension; adult localization requires 93–97 before 100. Figure assembly 96 also consumes previously generated publication assets and is historical, not a standalone graphics demo.
4. Modules ending release/audit may require original manifests or prior output directories. This archive preserves those checks. It does not contain their historical result files or assert a turnkey end-to-end rerun.
5. GSE190856 in the initial study configuration is optional historical context, not the adult human heart reference used in the final manuscript.

## Required private schema

Module 08 reads Summary.tar.gz members listed in config/local_rat.json, including gene_count_matrix.txt (gene rows; six C1–C6 and six L1–L6 sample columns), gene annotation, mapping metrics and read-QC tables; Report.tar.gz supplies sample_info.txt. Preserve original identifiers and biological replicate mapping. Counts can be fractional estimates; the pipeline flags this and does not silently round values. Do not generate mock data to pass gates.

## Manually obtained resources

Hallmark: MSigDB release 2025.1.Hs, Entrez GMT, path specified by biology_extension_v1.json. Adult heart: public source endpoints and exact row-identity verification are in modules 93–97. Sources used to seed curated literature records are supplied as URLs/identifiers in the relevant config, not as copied copyrighted full texts.
