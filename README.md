# GuaCAMOLE Diagnostics Viewer

An interactive **Streamlit app** for visualizing per-iteration diagnostic output from [GuaCAMOLE](https://github.com/Cibiv/GuaCAMOLE) — a GC-aware species abundance estimation tool for metagenomic data ([Nature Communications, 2025](https://doi.org/10.1038/s41467-025-65530-4)).

## What it does

GuaCAMOLE's `--verbose` flag writes one JSON file per iteration of its false-positive removal algorithm. This app reads those files and lets you scrub through each step to explore how:

- **Sequencing efficiencies** (GC-bin curves) evolve as taxa are removed
- **Residuals** per taxon reveal which species fit the model poorly
- **Abundance estimates** shift between iterations
- **Taxa removal** unfolds across cycles

## Demo data

The `guac_output/` directory contains example output from a real GuaCAMOLE run on sample **SRR12996245** — upload these files in the app to explore it immediately.

## Running the app

Requires [pixi](https://pixi.sh):

```bash
pixi run streamlit run app.py
```

The app opens at `http://localhost:8501`.

## Usage

1. **Upload files** in the sidebar:
   - All `iteration_*.json` files (and optionally `iteration_manifest.json`) from a GuaCAMOLE `--verbose` run
   - Optionally the `.guac` output file for species names instead of taxids

2. **Scrub through iterations** with the slider at the top

3. **Explore four tabs**:
   - **Sequencing Efficiency** — GC-bin efficiency curve, with previous step overlay
   - **Residuals** — per-taxon residual lines; taxa exceeding the removal threshold shown in red
   - **Abundances** — horizontal bar chart with log2 fold-change vs. previous step
   - **Taxa Removal Timeline** — bar chart of removals per cycle and a detailed table

## Tech stack

- [Streamlit](https://streamlit.io) — UI
- [Plotly](https://plotly.com/python/) — interactive charts
- [pandas](https://pandas.pydata.org) / [NumPy](https://numpy.org) — data handling
- [pixi](https://pixi.sh) — environment management

## Data format

Each `iteration_{cycle}_rep_{rep}.json` contains:

| Field | Description |
|---|---|
| `iteration` | False-positive removal cycle number |
| `rep` | Repetition within cycle (threshold halves each cycle) |
| `threshold` | Current residual-range threshold for taxon removal |
| `taxids` | NCBI taxids of taxa still in the analysis |
| `abundances` | Estimated relative abundances (parallel to `taxids`) |
| `efficiencies` | GC-dependent sequencing efficiency (101 values, GC bins 0–100) |
| `residuals` | 2D array (101 GC bins × n\_taxa); NaN where no reads expected |
| `res_range` | Per-taxon residual range (max − min of non-NaN residuals) |
| `removed_taxids` | Taxids removed at this step |

## Related

- [GuaCAMOLE source](https://github.com/Cibiv/GuaCAMOLE)
- [GuaCAMOLE paper](https://doi.org/10.1038/s41467-025-65530-4)
