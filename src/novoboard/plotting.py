"""Plotting functions for FDR validation results."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from matplotlib import pyplot

if TYPE_CHECKING:
    from novoboard.fdr import FDRValidationResult

logger = logging.getLogger(__name__)


def plot_fdr_validation(
    results_list: list[FDRValidationResult],
    samples: range,
    output_path: str,
) -> None:
    """Plot FDR validation results.
    
    Args:
        results_list: List of FDRValidationResult from validate_FDR
        samples: Range of sample numbers used
        output_path: Path to save the output figure
    """
    # Create output directory if needed
    output_dir = Path(output_path).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    fig, ax = pyplot.subplots(1, 2, figsize=(9, 4))
    labels = [f'decoy {10-x}0%' for x in samples]
    colors9 = ['b', 'c', 'g', 'k', 'm', 'r', 'y', 'orange', 'pink']
    colors = [colors9[x-1] for x in samples]
    
    for results, c, label in zip(results_list, colors, labels):
        ax[0].plot(results.estimated_fdr, results.true_fdr_T, color=c, label=label)
        ax[1].plot(results.cumsum, results.estimated_fdr, color=c, label=label)
    
    ax[0].plot([0, 0.05], [0, 0.05], color='black', linestyle='--', label='True FDR')
    ax[0].set_xlim(0, 0.05)
    ax[0].set_xlabel('Estimated FDR')
    ax[0].set_ylabel('True FDR')
    ax[0].legend()
    
    ax[1].plot(results_list[2].cumsum, results_list[2].true_fdr_T, color='black', label='True FDR', linestyle='--')
    ax[1].set_ylim(0, 0.05)
    ax[1].set_xlabel('Number of PSMs')
    ax[1].set_ylabel('FDR')
    ax[1].legend()
    
    fig.tight_layout()
    fig.savefig(output_path)
    logger.info(f"Figure saved to {output_path}")
    pyplot.close(fig)

