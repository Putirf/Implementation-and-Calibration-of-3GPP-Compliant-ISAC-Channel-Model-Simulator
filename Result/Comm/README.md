# Communication Channel Calibration

MATLAB implementations for communication-only large-scale and full-scale calibration experiments in the ISAC channel simulator.

The calibration workflow has two stages:

* **Simulation:** generate BS-UE communication-channel statistics and save them as `.mat` files.
* **Plotting:** compare the generated results with benchmark data and save the CDF figures as `.png` files.

The outer entry scripts run both stages automatically.

## Project Structure

```text
Comm_channel_calibration/
|-- README.md
|-- Calibration_largescale_comm/
|   |-- Calibration_largescale_comm.m
|   |-- Calibration_largescale_comm_run.m
|   |-- Plot_Comm_channel_calibration_largeScale.m
|   |-- Phase1Calibration_v42_CMCC.xlsx
|   `-- results/comm_large_scale/
|-- Calibration_fullscale_comm/
|   |-- Calibration_full_config1/
|   |   |-- Calibration_fullscale_config1.m
|   |   |-- Calibration_fullscale_config1_run.m
|   |   |-- Plot_Comm_channel_calibration_fullScale_config1.m
|   |   |-- Phase2Config1Calibration_v35_CMCC.xlsx
|   |   `-- results/comm_full_scale_config1/
|   `-- Calibration_full_config2/
|       |-- Calibration_fullscale_config2.m
|       |-- Calibration_fullscale_config2_run.m
|       |-- Plot_Comm_channel_calibration_fullScale_config2.m
|       |-- Phase2Config2Calibration_v28_CMCC.xlsx
|       `-- results/comm_full_scale_config2/
```

## Requirements

* MATLAB with support for package folders such as `+channel` and `+comm_scenario`.
* The project root added to the MATLAB path.
* The benchmark `.xlsx` file located in the corresponding calibration folder.
* Permission to create files under each `results` directory.

## Calibration Sets

| Set | Scenarios | Frequencies | Main purpose |
|---|---|---|---|
| Large-scale | UMa, UMi, InH | 6, 30, 70 GHz | Large-scale communication calibration |
| Full-scale Config 1 | UMa, UMi, InH | 6, 30, 60, 70 GHz | Full calibration, configuration 1 |
| Full-scale Config 2 | UMa, UMi, InH | 6, 30, 60, 70 GHz | Full calibration, configuration 2 |

## Running a Calibration

Use the outer entry script when you want the `.mat` results and `.png` figures to be generated automatically.

### Large-scale calibration

```matlab
cd('...\\ISAC_channel\\Comm_channel_calibration\\Calibration_largescale_comm');
addpath(pwd);
Calibration_largescale_comm
```

### Full-scale Config 1

```matlab
cd('...\\ISAC_channel\\Comm_channel_calibration\\Calibration_fullscale_comm\\Calibration_full_config1');
addpath(pwd);
Calibration_fullscale_config1
```

### Full-scale Config 2

```matlab
cd('...\\ISAC_channel\\Comm_channel_calibration\\Calibration_fullscale_comm\\Calibration_full_config2');
addpath(pwd);
Calibration_fullscale_config2
```

## Selecting a Scenario

The supported scenario names are:

```matlab
'UMa'
'UMi'
'InH'
```

For a selected scenario, the runner and plotter can also be called separately. For example:

```matlab
Calibration_fullscale_config1_run('UMi');
Plot_Comm_channel_calibration_fullScale_config1('UMi');
```

The outer entry scripts use the following default scenarios:

| Entry script | Default scenario |
|---|---|
| `Calibration_largescale_comm.m` | UMi |
| `Calibration_fullscale_config1.m` | UMa |
| `Calibration_fullscale_config2.m` | InH |

An environment variable named `CALIB_SELECTED_SCENARIO` can override the default selected by the entry script.

## Runner and Plotter Roles

The `_run.m` files only generate numerical results:

```text
<result folder>/<scenario>_<frequency>.mat
```

The outer entry scripts then call the plotter:

```matlab
Calibration_*_run(selected_scenario);
Plot_Comm_channel_calibration_*(selected_scenario);
```

The plotter reads the `.mat` files, loads the benchmark workbook, creates CDF comparisons, and saves the figures into the same `results` folder.

Running only a `_run.m` file will therefore create `.mat` files but will not create `.png` figures.

## Output

Typical output files include:

```text
UMa_6GHz.mat
UMa_30GHz.mat
UMa_60GHz.mat
UMa_70GHz.mat
```

The exact files depend on the calibration set. The output folders are:

```text
Calibration_largescale_comm/results/comm_large_scale/
Calibration_fullscale_comm/Calibration_full_config1/results/comm_full_scale_config1/
Calibration_fullscale_comm/Calibration_full_config2/results/comm_full_scale_config2/
```

## Calibration Data Flow

```text
scenario configuration
        |
BS / UE placement
        |
Comm_channel generation
        |
Coupling Loss, SIR/SINR, DS and angular statistics
        |
MAT-file result
        |
benchmark comparison and CDF plots
        |
PNG figures in results/
```

These are communication-only calibration tools. They use `Comm_channel` for BS-UE links and do not generate the sensing `Target_channel` links used by `ISAC_channel_calibration.m`.

## Troubleshooting

### Only `.mat` files are created

Run the outer entry script instead of the `_run.m` file, or call the plotter manually:

```matlab
Plot_Comm_channel_calibration_fullScale_config1('UMa');
```

### A benchmark workbook cannot be found

Check that the expected `.xlsx` file exists in the calibration folder and that MATLAB is running from the correct folder.

### Figures are visible but not saved

Check the MATLAB Command Window for missing `.mat` files, missing benchmark sheets, or file-permission errors. The plotter must reach its `saveas` command for the PNG file to be written.
