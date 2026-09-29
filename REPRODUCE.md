# Reproducing the analysis

Two levels:
1. **Figures and tables from the shipped results.** Run `make figures`. It needs only this repository.
2. **Everything from public data.** Follow the stages below.

Stage numbers match `workflow/`, and step identifiers (A1, A4, …) match the Methods.

## Setup

```bash
conda env create -f environment.yml && conda activate hsa2
export HSA2_ROOT=/scratch/hsa2        # working directory; needs about 150 GB (alignments, VCF slices)
make install                          # copies the workflow into $HSA2_ROOT/scripts
cd $HSA2_ROOT
```

All scripts run from `$HSA2_ROOT` and read and write relative paths (`data/`, `results/`, `logs/`). For long steps,
`scripts/run_job.sh <name> <command>` records `logs/jobs/<name>.{done,fail}`. `scripts/overnight_monitor.sh`
appends a health snapshot every 10 minutes. The reference run used an 80-core Linux server.

## Stage 0: data (a few hours, download-bound)

```bash
bash scripts/00_download_core.sh      # 8-way T2T ape HAL, CHM13, chains, Cactus v2.9.3 binaries
bash scripts/00_download_archaic.sh   # Altai, Vindija, Chagyrskaya, Denisova (region-streamed later)
```

- The diploid 16-way HAL for the haplotype cells comes from the same Cactus release (`data/cactus16/`).
- The deCODE maps (Palsson *et al.* 2025 supplement) go to `data/decode/`.
- The great-ape LD maps (github.com/lstevison/great-ape-recombination) go to `data/pan_recomb/`.
- NCBI assembly reports for species placement go to `data/asm_reports/`.

## Stage 1: branch-specific substitutions (A1; about 8 h on 16 cores)

```bash
awk -v L=5000000 '$1 ~ /^chr([0-9]+|X)$/ {for (s = 0; s < $2; s += L) print $1, s, L}' data/ref/hs1.chrom.sizes > results/A1/chunks.txt
cat results/A1/chunks.txt | xargs -P 16 -L 1 sh -c 'scripts/A1_run_chunk.sh $0 $1 $2'
```

Output: 624 chunks in `results/A1/chunks/`, each with `*.subs.tsv.gz`, `*.callable.tsv` and `*.ils.tsv`.

- **Other references and haplotypes.** GRCh38 as reference: `scripts/A2_hg38.sh`. Haplotype cells h16a/h16b on
  the 16-way HAL: `A1_run_chunk_ref.sh` with `config/h16{a,b}.json`, as in `utils/launch_parallel2.sh`.

## Stage 2: window tables (A3; minutes)

```bash
python scripts/A3_aggregate_chunks.py results/A1/chunks data/asm_reports results/A3/windows.tsv.gz
```

`results/A3/windows.tsv.gz` is shipped. It is the input to the model, calibration and diagnostics.

## Stage 3: present-day status (A5, A6, A12, A13; about 3 h)

```bash
python scripts/A5_recombination_today.py data results/A5
python scripts/A6_regions.py && cat results/A6/regions.txt | xargs -P 5 -L 1 scripts/A6_run_region.sh
bash scripts/A6c_fit_regions.sh
python scripts/A12_recomb_calibration.py results/A3/windows.tsv.gz data results/A12
python scripts/A12b_recomb_calibration_v2.py results/A12
bash scripts/A12c_B_by_recomb.sh 40
bash scripts/A13a_pan_recomb_to_hs1.sh
for sp in chimp bonobo gorilla; do python scripts/A13b_pan_recomb_profiles.py $sp results/A13 results/A3/windows.tsv.gz; done
```

## Stage 4: time strata (A7; about 2 h)

```bash
bash scripts/A7_run.sh
python scripts/A7b_classify.py results/A7 data/ref/hs1.chrom.sizes
python scripts/A7c_composition_correct.py results/A7/A7_summary.tsv results/A12/A12_windows.tsv.gz results/A7/A7c_corrected.tsv
```

## Stage 5: single-reference estimators and robustness (A14, A15; minutes)

```bash
python scripts/A15_pseudo_fusion.py results/A3/windows.tsv.gz results/A15 2000 6.0
python scripts/A14_logodds_cells.py data/ref/hs1.chrom.sizes results/A14/A14_logodds_cells.tsv \
    chm13_8way=results/A1/chunks h16a=results/A1_h16a/chunks h16b=results/A1_h16b/chunks
```

## Stage 6: pooled model (A4, A4b; about 1 h on 16 cores)

```bash
python scripts/A4_pooled_model.py results/A3/windows.tsv.gz results/A4 chimp            # main fit + profile
python scripts/A4_pooled_model.py results/A3/windows.tsv.gz results/A4 chimp --pseudo   # calibration on real ends
python scripts/A4_pooled_model.py results/A3/windows.tsv.gz results/A4 bonobo           # robustness (also --no-pan, --no-gorilla)
cd scripts && python A4b_sides.py ../results/A3/windows.tsv.gz ../results/A4 chimp 60 16 && cd ..   # two-flank LR test
```

## Stage 7: theory and simulations (A10; about 1 h)

```bash
python scripts/A10a_nonlinear_transient.py results/A10 400 8      # also 200 8 and 200 16
python scripts/A10d_bottleneck.py results/A10 300
python scripts/A10c_ils_structured.py results/A10 50000 16
bash scripts/A10b_run.sh 48 40                                    # SLiM, 600 runs
```

## Stage 8: lineage sorting and diagnostics (A8, A16, B1; minutes)

```bash
python scripts/A8_ils_profiles.py results/A1/chunks data/ref/hs1.chrom.sizes results/A8
python scripts/A16_2b_investigation.py results/A3/windows.tsv.gz results/A1/chunks results/A12/A12_windows.tsv.gz results/A16
bash scripts/B1_cen_feasibility.sh
```

## Stage 9: figures

From the repository root, `make figures`; or, in `$HSA2_ROOT`, run the scripts of `workflow/09_figures` as in the
`Makefile`.

## Expected key numbers

- Main fit: `T_off = 2.8` Mya, profile 95% CI 1.5–3.8.
- Calibration: 37/39 single ends and 38/40 pairs cover `T = 0`.
- SLiM: 15/15 scenarios cover the true value.
- Two-flank LR = 5.9, with empirical p = 0.049.

Small differences can arise from random seeds in the bootstrap and simulation steps.
