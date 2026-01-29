"""Command-line interface for NovoBoard."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys
from typing import Sequence

from novoboard import config
from novoboard.accuracy import WorkerTest
from novoboard.decoy import generate_decoy_mgf
from novoboard.fdr import validate_FDR
from novoboard.plotting import plot_fdr_validation

logger = logging.getLogger(__name__)


def download_data(data_dir: Path) -> None:
    """Download example ABRF data from Google Drive using gdown.
    
    Args:
        data_dir: Directory to download data into
    """
    import gdown
    
    folder_url = "https://drive.google.com/drive/folders/1_6azR4-YjTUfRYdsXbFhZL9lFvjdrIDh"
    
    data_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Downloading data to {data_dir}...")
    gdown.download_folder(folder_url, output=str(data_dir), quiet=False)
    logger.info("Download complete.")


def run_accuracy(
    db_file: Path,
    denovo_file: Path,
    spectrum_file: Path,
    col_score: str,
    col_aa_score: str,
) -> None:
    """Calculate accuracy of de novo predictions against database search results.
    
    Args:
        db_file: Path to database search results CSV
        denovo_file: Path to de novo sequencing results CSV
        spectrum_file: Path to MGF spectrum file
        col_score: Column name for peptide score
        col_aa_score: Column name for AA-level scores
    """
    logger.info(f"Database file: {db_file}")
    logger.info(f"De novo file: {denovo_file}")
    logger.info(f"Spectrum file: {spectrum_file}")
    
    worker_test = WorkerTest(
        str(db_file),
        str(denovo_file),
        str(spectrum_file),
        col_score,
        col_aa_score,
    )
    worker_test.test_accuracy()


def run_decoy_generation(
    spectrum_files: Sequence[Path],
    peak_sampling: str = 'random',
    sampling_rate: float = config.DEFAULT_SAMPLING_RATE,
    seed: int = 99,
) -> None:
    """Generate decoy MGF files for FDR estimation.
    
    Args:
        spectrum_files: List of input MGF spectrum files
        peak_sampling: Peak sampling strategy
        sampling_rate: Fraction of peaks to sample
        seed: Random seed for reproducibility
    """
    logger.info(f"Generating decoys for {len(spectrum_files)} spectrum file(s)")
    input_mgf_str_list = [str(p) for p in spectrum_files]
    generate_decoy_mgf(input_mgf_str_list, peak_sampling, sampling_rate, seed)


def run_fdr_validation(
    target_file: Path,
    decoy_files: Sequence[Path],
    db_file: Path,
    spectrum_file: Path,
    output_file: Path,
    col_score: str,
    col_aa_score: str,
    ion_threshold: float = 0.90,
) -> None:
    """Validate FDR estimation using target-decoy approach.
    
    Args:
        target_file: Path to target de novo results CSV
        decoy_files: Paths to decoy de novo results CSVs
        db_file: Path to database search results CSV (ground truth)
        spectrum_file: Path to MGF spectrum file
        output_file: Path for output plot
        col_score: Column name for peptide score
        col_aa_score: Column name for AA-level scores
        ion_threshold: Ion matching threshold percentage (default: 0.90)
    """
    logger.info(f"Target file: {target_file}")
    logger.info(f"Decoy files: {len(decoy_files)}")
    logger.info(f"Database file: {db_file}")
    
    p_decoy = [x / 1000. for x in range(0, 50, 1)]
    
    results_list = [
        validate_FDR(
            str(target_file),
            str(decoy_file),
            col_score,  # engine_score
            str(db_file),
            str(spectrum_file),
            p_decoy,
            ion_threshold,
            col_score,
            col_aa_score,
        ) 
        for decoy_file in decoy_files
    ]
    
    # Use file indices as sample labels
    samples = range(1, len(decoy_files) + 1)
    plot_fdr_validation(results_list, samples, str(output_file))


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
  # Calculate accuracy of de novo predictions
  novoboard accuracy --db-file db_results.csv --denovo-file denovo.csv --spectrum-file spectra.mgf

  # Generate decoy spectra
  novoboard decoy --spectrum-file spectra.mgf --sampling-rate 0.5

  # Validate FDR estimation
  novoboard fdr --target-file target.csv --decoy-files decoy1.csv decoy2.csv \\
                --db-file db_results.csv --spectrum-file spectra.mgf

  # Download example ABRF dataset
  novoboard download --output-dir data
        """
    )
    
    # Global arguments
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose (debug) logging'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # =========================================================================
    # DOWNLOAD command
    # =========================================================================
    download_parser = subparsers.add_parser(
        'download',
        help='Download example ABRF dataset from Google Drive'
    )
    download_parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('data'),
        help='Directory to download data into (default: data/)'
    )
    
    # =========================================================================
    # ACCURACY command
    # =========================================================================
    accuracy_parser = subparsers.add_parser(
        'accuracy',
        help='Calculate accuracy of de novo predictions against database search'
    )
    accuracy_parser.add_argument(
        '--db-file',
        type=Path,
        required=True,
        help='Path to database search results CSV'
    )
    accuracy_parser.add_argument(
        '--denovo-file',
        type=Path,
        required=True,
        help='Path to de novo sequencing results CSV'
    )
    accuracy_parser.add_argument(
        '--spectrum-file',
        type=Path,
        required=True,
        help='Path to MGF spectrum file'
    )
    accuracy_parser.add_argument(
        '--score-column',
        type=str,
        default='ALC (%)',
        help='Column name for peptide score (default: "ALC (%%)")'
    )
    accuracy_parser.add_argument(
        '--aa-score-column',
        type=str,
        default='local confidence (%)',
        help='Column name for AA-level scores (default: "local confidence (%%)")'
    )
    
    # =========================================================================
    # DECOY command
    # =========================================================================
    decoy_parser = subparsers.add_parser(
        'decoy',
        help='Generate decoy MGF files for FDR estimation'
    )
    decoy_parser.add_argument(
        '--spectrum-file',
        type=Path,
        required=True,
        nargs='+',
        help='Path(s) to input MGF spectrum file(s)'
    )
    decoy_parser.add_argument(
        '--sampling-strategy',
        type=str,
        default='random',
        choices=['random', 'intensity', 'intensity_mass', 'permutation', '500Da', 'distance'],
        help='Peak sampling strategy (default: random)'
    )
    decoy_parser.add_argument(
        '--sampling-rate',
        type=float,
        default=config.DEFAULT_SAMPLING_RATE,
        help=f'Fraction of peaks to keep (default: {config.DEFAULT_SAMPLING_RATE})'
    )
    decoy_parser.add_argument(
        '--seed',
        type=int,
        default=99,
        help='Random seed for reproducibility (default: 99)'
    )
    
    # =========================================================================
    # FDR command
    # =========================================================================
    fdr_parser = subparsers.add_parser(
        'fdr',
        help='Validate FDR estimation using target-decoy approach'
    )
    fdr_parser.add_argument(
        '--target-file',
        type=Path,
        required=True,
        help='Path to target de novo results CSV'
    )
    fdr_parser.add_argument(
        '--decoy-files',
        type=Path,
        required=True,
        nargs='+',
        help='Path(s) to decoy de novo results CSV(s)'
    )
    fdr_parser.add_argument(
        '--db-file',
        type=Path,
        required=True,
        help='Path to database search results CSV (ground truth)'
    )
    fdr_parser.add_argument(
        '--spectrum-file',
        type=Path,
        required=True,
        help='Path to MGF spectrum file'
    )
    fdr_parser.add_argument(
        '--output-file',
        type=Path,
        default=Path('fdr_validation.png'),
        help='Path for output plot (default: fdr_validation.png)'
    )
    fdr_parser.add_argument(
        '--score-column',
        type=str,
        default='ALC (%)',
        help='Column name for peptide score (default: "ALC (%%)")'
    )
    fdr_parser.add_argument(
        '--aa-score-column',
        type=str,
        default='local confidence (%)',
        help='Column name for AA-level scores (default: "local confidence (%%)")'
    )
    fdr_parser.add_argument(
        '--ion-threshold',
        type=float,
        default=0.90,
        help='Ion matching threshold percentage (default: 0.90)'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(verbose=args.verbose)
    
    # Handle no command
    if args.command is None:
        parser.print_help()
        sys.exit(0)
    
    # Execute command
    if args.command == 'download':
        download_data(args.output_dir)
        
    elif args.command == 'accuracy':
        # Validate files exist
        for f in [args.db_file, args.denovo_file, args.spectrum_file]:
            if not f.exists():
                logger.error(f"File not found: {f}")
                sys.exit(1)
        
        run_accuracy(
            args.db_file,
            args.denovo_file,
            args.spectrum_file,
            args.score_column,
            args.aa_score_column,
        )
        
    elif args.command == 'decoy':
        # Validate files exist
        for f in args.spectrum_file:
            if not f.exists():
                logger.error(f"File not found: {f}")
                sys.exit(1)
        
        run_decoy_generation(
            args.spectrum_file,
            args.sampling_strategy,
            args.sampling_rate,
            args.seed,
        )
        
    elif args.command == 'fdr':
        # Validate files exist
        files_to_check = [args.target_file, args.db_file, args.spectrum_file] + list(args.decoy_files)
        for f in files_to_check:
            if not f.exists():
                logger.error(f"File not found: {f}")
                sys.exit(1)
        
        # Create output directory if needed
        args.output_file.parent.mkdir(parents=True, exist_ok=True)
        
        run_fdr_validation(
            args.target_file,
            args.decoy_files,
            args.db_file,
            args.spectrum_file,
            args.output_file,
            args.score_column,
            args.aa_score_column,
            args.ion_threshold,
        )
    
    logger.info("NovoBoard analysis complete!")


if __name__ == '__main__':
    main()
