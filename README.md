# Implementation and Calibration of a 3GPP-Compliant ISAC Channel Model Simulator

This repository presents an implementation of the Integrated Sensing and Communication (ISAC) channel model based on 3GPP TR 38.901 (v19).

<p align="center">
  <img src="https://github.com/user-attachments/assets/ab331751-f0b5-403d-abff-9277f0e265f4" width="800">
</p>

## Calibration Results

The ISAC calibration workflows cover TRP-monostatic (`TRP_mono`), UE-monostatic (`UE_mono`), and bistatic (`TRP_TRP`, `TRP_UE`, and `UE_UE`) sensing, including:

- Table 7.9.6.1-1 / 7.9.6.1-2 / 7.9.6.1-3 / 7.9.6.1-4 (Large-scale metrics)
- Table 7.9.6.2-1 / 7.9.6.2-2 / 7.9.6.2-3 / 7.9.6.2-4 (Full-scale metrics)
- Table 7.9.6.3-1 (Spatial-consistency calibration for Urban Grid Vehicle and InH Human)

Snapshot calibration includes large-scale coupling loss and, with full-channel generation enabled, full coupling loss, delay spread (DS), and angular spreads (ASD, ZSD, ASA, and ZSA). The ISAC frequency presets are `ISAC_FR1` (6 GHz, 100 MHz bandwidth) and `ISAC_FR2` (30 GHz, 400 MHz bandwidth).

Calibration figures are saved under `Results`. Small-sample workflow checks do not establish agreement of every full Monte Carlo CDF with the company reference curves.

## Data Structure

For each calibration scenario, we provide:

- The required parameters for each calibration matrix
- Corresponding data organized in separate folders

The current directory layout is:

```text
ISAC_channel/                         Channel models, scenarios, sensing types,
                                     network layouts, antennas, and example driver
Calibrations/
  ISAC_channel_calibrations/          ISAC calibration drivers and plotting
  ISAC_benchmark/                     ISAC company reference curves
  Comm_channel_calibrations/          Communication LS/Full/spatial calibration
  Comm_benchmark/                     Communication company reference curves
Results/                             Calibration figures and existing result sets
3GPP規格書/                          Local specification documents
```

## Usage

### Calibration

The calibration workflows use MATLAB, Signal Processing Toolbox, Statistics and Machine Learning Toolbox, and Communications Toolbox. The workflows were checked with MATLAB R2025b.

The external `multiprod.m` function is required.

To reproduce the calibration plots, open the appropriate driver under `Calibrations/ISAC_channel_calibrations/`, edit its controls, and run it in MATLAB:

| Driver | Sensing modes / purpose |
|---|---|
| `calibration_ISAC_channel_monostatic.m` | `TRP_mono`, `UE_mono` |
| `calibration_ISAC_channel_bistatic.m` | `TRP_TRP`, `TRP_UE`, `UE_UE` |
| `calibration_ISAC_spatial_consistency.m` | Procedure B spatial-consistency calibration |

The snapshot drivers share the following case IDs:

| ID | Target and deployment |
|---|---|
| 1 | Small UAV / UMa-AV |
| 2 | Human / UMa |
| 3 | Human / UMi |
| 4 | Human / InH |
| 5 | Human / InF-SH |
| 6 | AGV / InF-SH |
| 7 | Vehicle / Urban Grid |

`UE_mono` supports Human, InF AGV, and Urban Grid Vehicle; its all-case run skips UAV case 1, which is available in `TRP_mono` and the three bistatic modes. Shared paths are resolved relative to the driver by `get_ISAC_calibration_paths.m`; no machine-specific path edit is required for the ISAC drivers.

Set `full_cali = false` for large-scale-only calibration or `true` to include full-channel metrics. `do_background` controls additional background-channel evaluation. `use_traditional_rms` selects the angular-spread calculation. Both snapshot drivers select links using large-scale target-channel power and retain that selection for Full, DS, and angular-spread statistics.

For example, to run UMa Human in `TRP_UE` at both frequencies, edit the bistatic controls:

```matlab
sensing_mode = 'TRP_UE';
case_id = 2;
use_isac_frequency_preset = true;
isac_frequency_preset = {'ISAC_FR1', 'ISAC_FR2'};
run_all_calibration_cases = false;
```

The bistatic driver also accepts multiple case IDs, such as `case_id = [2, 3]`. `run_all_calibration_cases = true` runs all seven cases at `calibration_frequency_presets`. `calibration_simulation_times` can override the scenario's drop count. TRP-UE uses directed Best-4 pairs; UE-UE uses distinct unordered Best-4 pairs with balanced STX/SRX roles. Physical-leg caches preserve shared propagation states across selected pairs. See [Bistatic calibration controls](Calibrations/ISAC_channel_calibrations/calibration_ISAC_channel_bistatic.md) for implementation and validation details.

For the monostatic driver, use `run_all_frequency_presets = true` to run both entries in `calibration_frequency_presets`. Urban Grid `UE_mono` additionally supports `urban_ut_type = 1` (Pedestrian), `2` (RSU), or `3` (Vehicle), one type per run. The bistatic Urban Grid calibration currently uses pedestrian UTs. Monostatic RSU/Vehicle layout generation and candidate LS Best-4 selection occur before the channel-calculation waitbar is created.

Current outdoor Human calibration trials use 100% outdoor UTs in `UE_mono` and `TRP_UE`. The effective UE height used for the breakpoint is 0.8 m for UMi Human in both modes and for UMa Human in `TRP_UE`; UMa Human `UE_mono` retains 0.5 m. In UMa/UMi Human `TRP_UE`, the Human-UE leg uses shadow-fading standard deviation 7 dB for both LOS and NLOS. These are explicit trial settings in the current implementation, not general specification defaults.

Spatial consistency is evaluated separately from snapshot LS/Full calibration, which disables spatial consistency. In `calibration_ISAC_spatial_consistency.m`, set `selected_case_ids` to `1` for Urban Grid Vehicle, `2` for InH Human, or `[]` for both. Configure `frequency_presets`, `num_ue_per_drop`, `num_drops`, and `full_channel` there. Communication-channel calibration drivers remain under `Calibrations/Comm_channel_calibrations/`.

The shared snapshot plotter reads the company benchmark workbooks and exports PNG/MATLAB FIG files to `Results`. Snapshot drivers pass results directly to plotting rather than automatically saving every result structure as a MAT file.

### ISAC Channel

For detailed ISAC channel configuration and usage, see `ISAC_channel/ISAC_channel.m`.

The script supports custom BS, UE, and ST positions. Additional parameters can be manually adjusted in the corresponding scenario, sensing type, network layout, channel, and antenna files.

This example driver exposes `TRP_mono` and `UE_mono`. Use the dedicated bistatic calibration driver for `TRP_TRP`, `TRP_UE`, and `UE_UE` calibration.

## Publication

This project is associated with the following published paper:

**Implementation and Calibration of 3GPP-Compliant ISAC Channel Simulator**

https://arxiv.org/abs/2606.07328

If you find this repository, calibration results, or implementation useful for your research, please cite the corresponding paper.

## Usage Notice

This repository .is provided, and the use of source for academic purpose is allowed only if the correct citation to the associated paper is conducted.

Unless otherwise stated in a separate license file, all rights to the source code, calibration data, figures, and related materials are reserved by the authors. Users may view the repository contents for research reference, but no permission is granted for redistribution, modification, or commercial use without prior written permission from the authors.

If you use this repository, calibration results, or implementation details in academic work, please cite the corresponding paper:

**Implementation and Calibration of 3GPP-Compliant ISAC Channel Simulator**
