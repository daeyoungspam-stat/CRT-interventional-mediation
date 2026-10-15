# Verification

Local verification completed on 2026-10-14.

| Check | Result |
|---|---|
| Archived source comparison | All seven engine-file MD5 values matched |
| Final bundle comparison | Estimator signature and full default configuration matched |
| Data generation | Deterministic by replication identifier; cluster sizes within 10 through 20 |
| Target values | Agreement with the recorded values to absolute tolerance 1e-12 |
| Workspace protection | Configuration mismatches and occupied unrelated directories rejected |
| Implementation assertions | Passed |
| Complete-method smoke study | K = 50, one replication, one worker |
| Canonical smoke output | 288 method-scenario-target rows |
| Low-dimensional smoke output | 396 method-scenario-target rows |
| Report generation | CSV summaries, LaTeX tables, and 32 PDF figures generated |
| Resume behavior | Completed replication-file hashes remained identical |
| R package check | R 4.6.1, Status: OK |

The fitting smoke study ran with R 4.5.1, SuperLearner 2.0.29, glmnet 4.1.10,
and ranger 0.17.0. The study exercised the full fitting configuration, including
200 g-computation bootstrap resamples, 10 reference draws, both reference laws,
all three classifiers, and all scenarios. The integration test is available through
`CRTMEDSIM_INTEGRATION=true` when running `R CMD check` with fitting dependencies
installed.

These checks establish source provenance and executable reproduction. Monte
Carlo performance is assessed from the repeated-replication study outputs.

