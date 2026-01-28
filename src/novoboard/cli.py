"""Command-line interface for NovoBoard."""

from __future__ import annotations

import argparse
import logging
import os
import sys

from novoboard.accuracy import WorkerTest
from novoboard.decoy import generate_decoy_mgf
from novoboard.fdr import validate_FDR
from novoboard.plotting import plot_fdr_validation

logger = logging.getLogger(__name__)


def download_data(data_dir: str) -> None:
    """Download data from Google Drive using gdown.
    
    Args:
        data_dir: Directory to download data into
    """
    import gdown
    
    folder_url = "https://drive.google.com/drive/folders/1_6azR4-YjTUfRYdsXbFhZL9lFvjdrIDh"
    
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    
    logger.info(f"Downloading data to {data_dir}...")
    gdown.download_folder(folder_url, output=data_dir, quiet=False)
    logger.info("Download complete.")


def run_accuracy(data_dir: str, col_score: str, col_aa_score: str) -> None:
    """Run accuracy calculation (Cell 2 equivalent).
    
    Args:
        data_dir: Path to data directory
        col_score: Column name for score values
        col_aa_score: Column name for AA score values
    """
    folder = f"{data_dir}/"
    target_file = f"{folder}pd_merged.csv.db.psms.csv"
    spectrum_file = f"{folder}2017-12-4_ABRF_200_DDA1.mgf"
    
    for x in range(10, 10 + 1):
        predicted_file = f"{folder}PEAKS/Sample {x}.denovo.csv"
        worker_test = WorkerTest(target_file, predicted_file, spectrum_file, col_score, col_aa_score)
        worker_test.test_accuracy()


def run_decoy_generation(
    data_dir: str,
    peak_sampling: str = 'random',
    sampling_rate: float = 0.5,
) -> None:
    """Run decoy MGF generation (Cell 4 equivalent).
    
    Args:
        data_dir: Path to data directory
        peak_sampling: Peak sampling strategy
        sampling_rate: Fraction of peaks to sample
    """
    folder = f"{data_dir}/"
    input_mgf_list = [
        '2017-12-4_ABRF_200_DDA1.mgf',
    ]
    input_mgf_list = [f"{folder}{x}" for x in input_mgf_list]
    
    generate_decoy_mgf(input_mgf_list, peak_sampling, sampling_rate)


def run_fdr_validation(
    data_dir: str,
    output_dir: str,
    col_score: str,
    col_aa_score: str,
) -> None:
    """Run FDR validation (Cell 8 equivalent).
    
    Args:
        data_dir: Path to data directory
        output_dir: Path for output files
        col_score: Column name for score values
        col_aa_score: Column name for AA score values
    """
    folder = f"{data_dir}/"
    db_csv = f"{folder}pd_merged.csv.db.psms.csv"
    spectrum_file = f"{folder}2017-12-4_ABRF_200_DDA1.mgf"

    p_decoy = [x / 1000. for x in range(0, 50, 1)]
    T_pct = 0.90

    samples = range(3, 7 + 1)
    target_csv = f"{folder}PEAKS/Sample 10.denovo.csv"
    decoy_csv_list = [f"{folder}PEAKS/Sample {x}.denovo.csv" for x in samples]
    engine_score = col_score
    
    results_list = [
        validate_FDR(target_csv, decoy_csv, engine_score, db_csv, spectrum_file, p_decoy, T_pct, col_score, col_aa_score) 
        for decoy_csv in decoy_csv_list
    ]

    output_path = os.path.join(output_dir, 'fig.decoy_fdr_valid_X_random_abrf_peaks.png')
    plot_fdr_validation(results_list, samples, output_path)


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for CLI usage.
    
    Args:
        verbose: If True, set DEBUG level; otherwise INFO level
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def main() -> None:
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        description='NovoBoard - Framework for evaluating de novo peptide sequencing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  novoboard --data-dir data --output-dir fig
  novoboard --download --data-dir data --output-dir fig
        """
    )
    parser.add_argument(
        '--data-dir', 
        default='data',
        help='Path to data directory (default: data/)'
    )
    parser.add_argument(
        '--output-dir', 
        default='fig',
        help='Path for output files/plots (default: fig/)'
    )
    parser.add_argument(
        '--download',
        action='store_true',
        help='Download data using gdown before running'
    )
    parser.add_argument(
        '--skip-accuracy',
        action='store_true',
        help='Skip accuracy calculation step'
    )
    parser.add_argument(
        '--skip-decoy',
        action='store_true',
        help='Skip decoy MGF generation step'
    )
    parser.add_argument(
        '--skip-fdr',
        action='store_true',
        help='Skip FDR validation step'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose (debug) logging'
    )
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(verbose=args.verbose)
    
    # Column names for score values
    col_score = "ALC (%)"
    col_aa_score = "local confidence (%)"
    
    # Download data if requested
    if args.download:
        download_data(args.data_dir)
    
    # Check if data directory exists
    if not os.path.exists(args.data_dir):
        print(f"Error: Data directory '{args.data_dir}' does not exist.")
        print("Use --download to download the data first, or specify a valid --data-dir.")
        sys.exit(1)
    
    # Create output directory if needed
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Run all steps sequentially (like the notebook)
    if not args.skip_accuracy:
        print("\n" + "=" * 80)
        print("Step 1: Running accuracy calculation...")
        print("=" * 80)
        run_accuracy(args.data_dir, col_score, col_aa_score)
    
    if not args.skip_decoy:
        print("\n" + "=" * 80)
        print("Step 2: Generating decoy MGF...")
        print("=" * 80)
        run_decoy_generation(args.data_dir)
    
    if not args.skip_fdr:
        print("\n" + "=" * 80)
        print("Step 3: Running FDR validation...")
        print("=" * 80)
        run_fdr_validation(args.data_dir, args.output_dir, col_score, col_aa_score)
    
    print("\n" + "=" * 80)
    print("NovoBoard analysis complete!")
    print("=" * 80)


if __name__ == '__main__':
    main()
