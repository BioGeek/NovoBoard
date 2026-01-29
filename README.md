# NovoBoard

A comprehensive framework for evaluating the false discovery rate and accuracy of de novo peptide sequencing.

## Features

- Calculate fragment ion, amino acid, and peptide accuracies
- Generate decoy spectra for FDR estimation
- Validate FDR estimation against known database matches
- Visualization of FDR validation results
- **Software-agnostic**: Works with any de novo sequencing tool

## Installation

### Using uv (recommended)

```bash
# Clone the repository
git clone https://github.com/BioGeek/NovoBoard.git
cd NovoBoard

# Create virtual environment and install dependencies
uv sync

# Activate the virtual environment
source .venv/bin/activate
```

### Using pip

```bash
pip install -e .
```

## Usage

NovoBoard provides a command-line interface with three main commands:

### 1. Calculate Accuracy

Compare de novo predictions against database search results:

```bash
novoboard accuracy \
    --db-file db_results.csv \
    --denovo-file denovo_predictions.csv \
    --spectrum-file spectra.mgf \
    --score-column "ALC (%)" \
    --aa-score-column "local confidence (%)"
```

**Options:**
- `--db-file`: Database search results CSV (required)
- `--denovo-file`: De novo sequencing results CSV (required)
- `--spectrum-file`: MGF spectrum file (required)
- `--score-column`: Column name for peptide score (default: "ALC (%)")
- `--aa-score-column`: Column name for AA-level scores (default: "local confidence (%)")

### 2. Generate Decoy Spectra

Create decoy MGF files for FDR estimation:

```bash
novoboard decoy \
    --spectrum-file spectra.mgf \
    --sampling-strategy random \
    --sampling-rate 0.5 \
    --seed 99
```

**Options:**
- `--spectrum-file`: Input MGF file(s) (required, accepts multiple)
- `--sampling-strategy`: Strategy for peak sampling (default: random)
  - `random`: Random peak sampling
  - `intensity`: Sample by peak intensity
  - `intensity_mass`: Intensity + peptide mass
  - `permutation`: Shuffle intensities
  - `500Da`: Remove peaks < 500 Da
  - `distance`: Remove by AA mass distances
- `--sampling-rate`: Fraction of peaks to keep (default: 0.5)
- `--seed`: Random seed for reproducibility (default: 99)

### 3. Validate FDR

Validate FDR estimation using target-decoy approach:

```bash
novoboard fdr \
    --target-file target_denovo.csv \
    --decoy-files decoy1.csv decoy2.csv decoy3.csv \
    --db-file db_results.csv \
    --spectrum-file spectra.mgf \
    --output-file fdr_validation.png \
    --ion-threshold 0.90
```

**Options:**
- `--target-file`: Target de novo results CSV (required)
- `--decoy-files`: Decoy de novo results CSV(s) (required, accepts multiple)
- `--db-file`: Database search results CSV (required)
- `--spectrum-file`: MGF spectrum file (required)
- `--output-file`: Output plot path (default: fdr_validation.png)
- `--score-column`: Column name for peptide score (default: "ALC (%)")
- `--aa-score-column`: Column name for AA-level scores (default: "local confidence (%)")
- `--ion-threshold`: Ion matching threshold (default: 0.90)

### 4. Download Example Data

Download the ABRF example dataset:

```bash
novoboard download --output-dir data
```

## Example Workflow

```bash
# 1. Download example data
novoboard download --output-dir data

# 2. Generate decoy spectra
novoboard decoy --spectrum-file data/spectra.mgf

# 3. Run de novo sequencing on target and decoy spectra
# (using your preferred tool: PEAKS, Casanovo, InstaNovo, etc.)

# 4. Calculate accuracy
novoboard accuracy \
    --db-file data/db_results.csv \
    --denovo-file results/denovo.csv \
    --spectrum-file data/spectra.mgf

# 5. Validate FDR
novoboard fdr \
    --target-file results/target.csv \
    --decoy-files results/decoy*.csv \
    --db-file data/db_results.csv \
    --spectrum-file data/spectra.mgf \
    --output-file figures/fdr.png
```

## Jupyter Notebook

The original notebook `aa.fdr_github.ipynb` is still available for interactive analysis:

```bash
source .venv/bin/activate
jupyter notebook aa.fdr_github.ipynb
```

## Development

### Running Tests

```bash
source .venv/bin/activate
pytest tests/ -v
```

### Project Structure

```
novoboard/
├── src/novoboard/
│   ├── __init__.py
│   ├── cli.py           # Command-line interface
│   ├── config.py        # Vocabulary and mass definitions
│   ├── accuracy.py      # Accuracy calculation (WorkerTest)
│   ├── decoy.py         # Decoy MGF generation
│   ├── fdr.py           # FDR calculation and validation
│   ├── mgf.py           # MGF file parser
│   └── plotting.py      # Visualization functions
├── tests/               # Unit tests
├── data/                # Data directory (not in git)
├── aa.fdr_github.ipynb  # Original Jupyter notebook
├── config.py            # Compatibility shim for notebook
├── download_data.py     # Data download script
└── pyproject.toml       # Project configuration
```

## Input File Formats

### De novo Results CSV

Required columns:
- `Source File`: Source spectrum file name
- `Scan`: Scan number
- `Peptide`: Peptide sequence with modifications
- Score column (configurable, e.g., `ALC (%)`)
- AA score column (configurable, e.g., `local confidence (%)`)

### Database Search Results CSV

Required columns:
- `Source File`: Source spectrum file name
- `Scan`: Scan number
- `Peptide`: Peptide sequence with modifications

### MGF Spectrum File

Standard MGF format with:
- `BEGIN IONS` / `END IONS` markers
- `TITLE`, `PEPMASS`, `CHARGE`, `SCANS` headers
- Peak list as `m/z intensity` pairs

## Citation

If you use NovoBoard in your research, please cite:

> "NovoBoard: a comprehensive framework for evaluating the false discovery rate and accuracy of de novo peptide sequencing"
> https://doi.org/10.1101/2024.04.16.589668
