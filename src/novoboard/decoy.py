"""Generate decoy MGF files for FDR estimation."""

from __future__ import annotations

import re
import random
import numpy as np
from novoboard import config


def generate_decoy_mgf(
    input_mgf_list: list[str],
    peak_sampling: str = 'random',
    sampling_rate: float = 0.5,
) -> None:
    """Generate decoy MGF files with specified peak sampling strategy.
    
    Creates decoy spectra by sampling/shuffling peaks from target spectra,
    used for estimating false discovery rate in de novo sequencing.
    
    Args:
        input_mgf_list: List of input MGF file paths
        peak_sampling: Sampling strategy. Options:
            - 'random': Random peak sampling
            - 'intensity': Sample based on peak intensity
            - 'intensity_mass': Sample based on intensity and peptide mass
            - 'permutation': Shuffle peak intensities
            - '500Da': Remove peaks under 500 Da
            - 'distance': Remove peaks matching AA mass differences
        sampling_rate: Fraction of peaks to keep (0.0-1.0)
    """
    print(f"peak_sampling = {peak_sampling}")
    print(f"sampling_rate = {sampling_rate}")
    
    for input_mgf in input_mgf_list:

        if peak_sampling == 'permutation':
            output_mgf = f"{input_mgf}.permutation.mgf"
        elif peak_sampling == '500Da':
            output_mgf = f"{input_mgf}.500Da.mgf"
        else:
            output_mgf = f"{input_mgf}.decoy_{sampling_rate:.2f}.mgf"

        # collect peak_distr and removed_peaks_distr for noise sampling
        peaks_distr: list[list[float]] = []
        removed_peaks_distr: list[list[float]] = []
        
        if peak_sampling in ('random', 'permutation', 'distance'):
            with open(input_mgf, 'r') as f:
                peaks_distr = [[float(x) for x in re.split(r' |\r|\n', line)[:2]] for line in f if line[0].isdigit()]
            print(f"len(peaks_distr) = {len(peaks_distr)}")
        elif peak_sampling == '500Da':
            with open(input_mgf, 'r') as f:
                peaks_distr = [[float(x) for x in re.split(r' |\r|\n', line)[:2]] for line in f if line[0].isdigit()]
            print(f"len(peaks_distr) = {len(peaks_distr)}")
            removed_peaks_distr = [[x, y] for x, y in peaks_distr if x < 500]
            print(f"len(removed_peaks_distr) = {len(removed_peaks_distr)}")
        elif peak_sampling in ('intensity', 'intensity_mass'):
            with open(input_mgf, 'r') as f_in:
                while True:
                    line = f_in.readline()
                    if not line:  # end of file
                        break
                    if line == '\n':  # empty line
                        continue
                    peak_list: list[list[float]] = []
                    peptide_mass = 0.0
                    while "END IONS" not in line:
                        # parse header lines
                        if 'BEGIN IONS' in line or '=' in line:
                            line = f_in.readline()
                            if 'PEPMASS' in line:
                                mz = float(re.split(r'=| |\r|\n', line)[1])
                            if 'CHARGE' in line:
                                z = float(re.split(r'=|\+|\r|\n', line)[1])
                                peptide_mass = mz * z - z * 1.0078
                            continue
                        # parse ions
                        mz_str, intensity_str = re.split(r' |\r|\n', line)[:2]
                        peak_list.append([float(mz_str), float(intensity_str)])
                        line = f_in.readline()

                    # peak removal and noise sampling by intensity
                    num_peaks = len(peak_list)
                    if peak_sampling == 'intensity':
                        num_sampling = int(num_peaks * sampling_rate)
                    elif peak_sampling == 'intensity_mass':
                        est_len = int(peptide_mass / 122.8652943)
                        num_noise = int(min((est_len - 1) * 2, num_peaks) * (1 - sampling_rate))
                        num_sampling = num_peaks - num_noise
                    peak_list_sorted = sorted(peak_list, key=lambda x: x[1])
                    removed_peaks = peak_list_sorted[num_sampling:]
                    removed_peaks_distr += removed_peaks
                    peaks_distr += peak_list
            print(f"len(peaks_distr) = {len(peaks_distr)}")
            print(f"len(removed_peaks_distr) = {len(removed_peaks_distr)}")

        sampling_peaks_distr: list[list[float]] = []
        noise_peaks_distr: list[list[float]] = []
        decoy_peaks_distr: list[list[float]] = []
        
        with open(input_mgf, 'r') as f_in:
            with open(output_mgf, 'w') as f_out:
                while True:
                    line = f_in.readline()
                    if not line:  # end of file
                        break
                    if line == '\n':  # empty line
                        continue
                    peak_list = []
                    peptide_mass = 0.0
                    while "END IONS" not in line:
                        # parse header lines
                        if 'BEGIN IONS' in line or '=' in line:
                            f_out.write(line)
                            line = f_in.readline()
                            if 'PEPMASS' in line:
                                mz = float(re.split(r'=| |\r|\n', line)[1])
                            if 'CHARGE' in line:
                                z = float(re.split(r'=|\+|\r|\n', line)[1])
                                peptide_mass = mz * z - z * 1.0078
                            if 'SCANS=' in line:
                                scan = re.split(r'=|\r|\n', line)[1]
                            continue
                        # parse ions
                        mz_str, intensity_str = re.split(r' |\r|\n', line)[:2]
                        peak_list.append([float(mz_str), float(intensity_str)])
                        line = f_in.readline()
                    
                    num_peaks = len(peak_list)
                    num_sampling = int(num_peaks * sampling_rate)
                    num_noise = num_peaks - num_sampling
                    random.seed(99)
                    np.random.seed(99)
                    
                    sampling_peaks: list = []
                    noise_peaks: list = []
                    
                    # random peak sampling
                    if peak_sampling == 'random':
                        sampling_peaks = random.sample(peak_list, num_sampling)
                        noise_peaks = random.sample(peaks_distr, num_noise)
                    # peak removal and noise sampling by intensity
                    elif peak_sampling == 'intensity':
                        peak_list_sorted = sorted(peak_list, key=lambda x: x[1])
                        sampling_peaks = peak_list_sorted[:num_sampling]
                        noise_peaks = random.sample(removed_peaks_distr, num_noise)
                    # peak removal and noise sampling by intensity and peptide mass
                    elif peak_sampling == 'intensity_mass':
                        est_len = int(peptide_mass / 122.8652943)
                        num_noise = int(min((est_len - 1) * 2, num_peaks) * (1 - sampling_rate))
                        num_sampling = num_peaks - num_noise
                        peak_list_sorted = sorted(peak_list, key=lambda x: x[1])
                        sampling_peaks = peak_list_sorted[:num_sampling]
                        noise_peaks = random.sample(removed_peaks_distr, num_noise)
                    # peak permutation
                    elif peak_sampling == 'permutation':
                        mz_list = [x[0] for x in peak_list]
                        intensity_list = [x[1] for x in peak_list]
                        sampling_peaks = []
                        random.shuffle(intensity_list)
                        noise_peaks = list(zip(mz_list, intensity_list))
                    # remove peaks under 500 Da and replace by noise peaks
                    elif peak_sampling == '500Da':
                        sampling_peaks = [[x, y] for x, y in peak_list if x >= 500]
                        num_sampling = len(sampling_peaks)
                        num_noise = num_peaks - num_sampling
                        noise_peaks = random.sample(removed_peaks_distr, num_noise)
                    # peak removal by distance
                    elif peak_sampling == 'distance':
                        mz_array = np.array([peak[0] for peak in peak_list])
                        pair_distance = np.absolute(np.reshape(mz_array, (num_peaks, 1)) - np.reshape(mz_array, (1, num_peaks)))
                        aa_masses = config.mass_ID_np[3:].reshape(1, -1)
                        pair_aa_match = np.absolute(np.expand_dims(pair_distance, axis=2) - aa_masses)
                        pair_aa_match = np.any(pair_aa_match <= 0.02, axis=2)
                        match_peak_indices = list(np.flatnonzero(np.any(pair_aa_match, axis=1)))
                        num_noise = int(len(match_peak_indices) * (1 - sampling_rate))
                        removed_indices = random.sample(match_peak_indices, num_noise)
                        sampling_peaks = [peak for index, peak in enumerate(peak_list) if index not in removed_indices]
                        noise_peaks = random.sample(peaks_distr, num_noise)
                    else:
                        sampling_peaks = peak_list
                        noise_peaks = []
                    
                    sampling_peaks_distr += sampling_peaks
                    noise_peaks_distr += noise_peaks
                    decoy_peaks_distr += sampling_peaks + noise_peaks

                    # write ion lines
                    sorted_peaks = sorted(sampling_peaks + noise_peaks, key=lambda x: x[0])
                    for x, y in sorted_peaks:
                        f_out.write(f"{x:.5f} {y:.5f}\n")
                    f_out.write(line)  # END IONS line
                    f_out.write(f_in.readline())  # empty line between spectra
        
        print(f"len(sampling_peaks_distr) = {len(sampling_peaks_distr)}")
        print(f"len(noise_peaks_distr) = {len(noise_peaks_distr)}")
        print(f"len(decoy_peaks_distr) = {len(decoy_peaks_distr)}")
        print()
