# Context dependence and inferential limits of transcriptomic drug reversal

Code archive for the sepsis-associated myocardial transcriptomic reversal study, prepared for Computational and Structural Biotechnology Journal. Version: v0.1.0-prepublication (2026-10-06).

The study distinguishes global reversal ranking, shortlist stability, and formal inferential support. Expression reversal is hypothesis-generating evidence and is not evidence of cardioprotection.

## What this archive contains

Numbered Python/R acquisition, preprocessing, analysis and audit modules; frozen scientific configurations; Python dependency snapshots and an R/Bioconductor lock; input-contract tests; a module inventory and stage guide. Historical modules are retained for provenance. Their numeric names reflect development stages, not a command to execute every file in lexical order. Manuscript builders, machine-specific reports, private data and generated results are excluded.

## Install and check

Validated analysis platform: Windows, Python 3.11, R 4.4.3 and Bioconductor 3.20. Raw-read quantification includes Windows-specific Subread 2.1.1 installation; other operating systems require corresponding binaries and an explicit configuration review.

```powershell
python -m venv .venv
.venv/Scripts/python.exe -m pip install -r environment/python_final.lock
.venv/Scripts/python.exe setup_runtime.py
.venv/Scripts/python.exe preflight.py
.venv/Scripts/python.exe -m unittest discover -s tests -v
```

Set RSCRIPT or edit the ignored config/runtime.local.json to your Rscript executable. The original R package versions are recorded in environment/renv.lock. With R 4.4.3 and renv installed, use `renv::restore(lockfile="environment/renv.lock", library="environment/R-library", prompt=FALSE)` from the repository root. Do not install the project R library into a shared global library.

## Run the study

Start with docs/EXECUTION_GUIDE.md and docs/MODULE_INVENTORY.tsv. Study settings, seeds, paths, contrasts and family definitions are in config/. Each numbered Python module uses the common logger and writes outputs/; its corresponding R file is normally launched by that Python module. Do not launch R files independently without checking the wrapper arguments.

The archive is an analysis-source release, not a one-command data bundle. Full execution requires large public downloads, the correct external resources and the private H9c2 inputs. Some audit/release modules require original stage manifests and completed upstream outputs. Those archival assertions are not disabled to manufacture a successful rerun.

## Data access and reproducibility boundary

No H9c2 raw counts, normalized expression, gene-level results or private source reports are included. The H9c2 data remain unavailable publicly while the investigators' research is ongoing. Cross-model modules and figures using H9c2 cannot be reproduced from public inputs alone. This restriction is not a journal-approved data-sharing exemption. See data/README.md and docs/EXECUTION_GUIDE.md.

Local validation covers syntax, package integrity, portability of the shared runtime configuration and six input-contract regression checks. The full study was not rerun in a clean environment for this archive. Source code availability and complete computational reproducibility are distinct.

## Citation and licensing

Manuscript title: Context dependence and inferential limits of transcriptomic drug reversal in sepsis-associated myocardial alterations. Authors, journal publication details and DOI are pending. Do not cite an invented DOI. A code licence has not yet been assigned by the rights holders; do not describe this archive as licensed open-source software. Third-party dependencies and data retain their original terms.

AI assistance was used during code development, source checking, drafting and figure preparation. Final author review and responsibility must be confirmed before journal submission.
