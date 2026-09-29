# Curated intermediate results

These are small outputs of the workflow, enough to regenerate every figure and table with `make figures`. File
names follow the step identifiers of the Methods.

| Folder | Contents |
|---|---|
| `A3/`, `A3_h16a/`, `A3_h16b/` | 100-kb window tables per branch: W/S substitution counts, parent W/S sites, species telomere distance. Input to the model. |
| `A4*/` | pooled model: fits, profile likelihoods, observations, calibration on real ends, two-flank LR test (`A4b_sides_*.tsv`) |
| `A5/` | deCODE recombination at real ends, interior and fusion flanks |
| `A6/`, `A12c/` | present-day gBGC strength *B* by region and by paternal-DSB quantile (site-frequency spectra) |
| `A7/` | time strata of human-branch substitutions, composition-corrected (`A7c_corrected.tsv`) |
| `A8/`, `A16*/` | ILS profiles; flank diagnostics (profile, callability, orthologue end strength, junction ILS) |
| `A10/`, `A10b/` | exact Wright–Fisher bias (constant size, bottleneck), msprime lineage-sorting predictions, SLiM estimates |
| `A12/` | recombination calibration (models, per-end) and the window table joined to deCODE |
| `A13/` | LD recombination profiles of chimpanzee, bonobo and gorilla ends and flank orthologues |
| `A14/`, `A15*/` | assembly/haplotype cells; single-reference estimators applied to real ends |
| `figures/stats.json` | statistics shown in the figures and quoted in the text |
