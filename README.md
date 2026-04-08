# Implementation and Calibration of a 3GPP-Compliant ISAC Channel Model Simulator

This repository presents an implementation of the Integrated Sensing and Communication (ISAC) channel model based on 3GPP TR 38.901 (v19).

The full source code will be released after the acceptance of the corresponding paper.

<p align="center">
  <img src="https://github.com/user-attachments/assets/ab331751-f0b5-403d-abff-9277f0e265f4" width="800">
</p>

## Calibration Results

The monostatic ISAC calibration results for each scenario defined in 3GPP TR 38.901 are provided, including:

- Table 7.9.6.1-1 / 7.9.6.1-2 / 7.9.6.1-3 / 7.9.6.1-4 (Large-scale metrics)
- Table 7.9.6.2-1 / 7.9.6.2-2 / 7.9.6.2-3 / 7.9.6.2-4 (Full-scale metrics)

All results are stored in the `Result` directory.

## Data Structure

For each calibration scenario, we provide:
- The number of elements for each calibration matrix
- Corresponding data organized in separate folders

## Usage

To reproduce the calibration plots:

1. Modify the file path in `Calibration.m`
2. Run the script to generate the corresponding calibration results for each scenario



