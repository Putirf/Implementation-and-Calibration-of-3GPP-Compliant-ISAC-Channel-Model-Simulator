# Implementation and Calibration of a 3GPP-Compliant ISAC Channel Model Simulator

This repository presents an implementation of the Integrated Sensing and Communication (ISAC) channel model based on 3GPP TR 38.901 (v19).

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
- TThe required parameters for each calibration matrix
- Corresponding data organized in separate folders

## Usage

To reproduce the calibration plots:

1. Modify the file path in `Calibration.m`
2. Run the script to generate the corresponding calibration results for each scenario


## Publication

This project is associated with the following published paper:

**Implementation and Calibration of 3GPP-Compliant ISAC Channel Simulator**

If you find this repository, calibration results, or implementation useful for your research, please cite the corresponding paper.


## Usage Notice

This repository is provided for academic and research reference purposes.

Unless otherwise stated in a separate license file, all rights to the source code, calibration data, figures, and related materials are reserved by the authors. Users may view the repository contents for research reference, but no permission is granted for redistribution, modification, or commercial use without prior written permission from the authors.

If you use this repository, calibration results, or implementation details in academic work, please cite the corresponding paper:

**Implementation and Calibration of 3GPP-Compliant ISAC Channel Simulator**




