# Implement Communication Channel Calibration

MATLAB implementation for communication-only large-scale and full-scale calibration experiments in the ISAC channel simulator.

The folder has three main calibration workflows:

- **Large-scale calibration**: evaluate communication path loss, coupling loss, SIR, and SINR statistics.
- **Full-scale Config 1**: evaluate full communication-channel statistics using the Config 1 calibration assumptions.
- **Full-scale Config 2**: evaluate full communication-channel statistics using the Config 2 calibration assumptions.

All three workflows generate BS-UE communication links through `channel.Comm_channel`. They do not generate the sensing `Target_channel` links used by `ISAC_channel_calibration.m`.

## Project Structure

```text
Calibration_largescale_comm/
    Calibration_largescale_comm.m                 Main entry script
    Calibration_largescale_comm_run.m             Simulation runner
    Plot_Comm_channel_calibration_largeScale.m    Result plotting utility
    Phase1Calibration_v42_CMCC.xlsx               Benchmark workbook
    results/comm_large_scale/                     MAT and PNG outputs

Calibration_fullscale_comm/
    Calibration_full_config1/
        Calibration_fullscale_config1.m           Main entry script
        Calibration_fullscale_config1_run.m       Simulation runner
        Plot_Comm_channel_calibration_fullScale_config1.m
                                                     Result plotting utility
        Phase2Config1Calibration_v35_CMCC.xlsx    Benchmark workbook
        results/comm_full_scale_config1/           MAT and PNG outputs

    Calibration_full_config2/
        Calibration_fullscale_config2.m           Main entry script
        Calibration_fullscale_config2_run.m       Simulation runner
        Plot_Comm_channel_calibration_fullScale_config2.m
                                                     Result plotting utility
        Phase2Config2Calibration_v28_CMCC.xlsx    Benchmark workbook
        results/comm_full_scale_config2/           MAT and PNG outputs
```

Shared model components are located in the parent `ISAC_channel` folder:

```text
+comm_scenario/                 Communication scenario definitions
+network_layout/                BS and UE placement functions
+channel/                       Communication and target channel models
+antennas/                      Antenna element, panel, and array models
+elements/                      Equipment and sector objects
```

## Requirements

- MATLAB with object-oriented package folder support.
- The project root added to the MATLAB path.
- The benchmark Excel workbook in the corresponding calibration folder.
- MATLAB functions used by the channel model, including `randsrc`, `pow2db`, and `db2pow`.
- Permission to create files under the calibration `results` folders.

## Calibration

Use the outer entry script for each calibration set. The entry script runs the simulation, saves the numerical results, loads the benchmark data, and saves the calibration figures.

| Calibration set | Entry script | Scenarios | Frequencies |
| --- | --- | --- | --- |
| Large-scale | `Calibration_largescale_comm.m` | UMa, UMi, InH | 6, 30, 70 GHz |
| Full-scale Config 1 | `Calibration_fullscale_config1.m` | UMa, UMi, InH | 6, 30, 60, 70 GHz |
| Full-scale Config 2 | `Calibration_fullscale_config2.m` | UMa, UMi, InH | 6, 30, 60, 70 GHz |

### Large-scale Calibration

```matlab
cd('...\\ISAC_channel\\Comm_channel_calibration\\Calibration_largescale_comm');
addpath(pwd);
Calibration_largescale_comm
```

The default scenario is `UMi`.

### Full-scale Config 1

```matlab
cd('...\\ISAC_channel\\Comm_channel_calibration\\Calibration_fullscale_comm\\Calibration_full_config1');
addpath(pwd);
Calibration_fullscale_config1
```

The default scenario is `UMa`.

### Full-scale Config 2

```matlab
cd('...\\ISAC_channel\\Comm_channel_calibration\\Calibration_fullscale_comm\\Calibration_full_config2');
addpath(pwd);
Calibration_fullscale_config2
```

The default scenario is `InH`.

### Runner and Plotter Separation

Each calibration set has a runner and a plotter. The `_run.m` file only generates and saves `.mat` results:

```matlab
Calibration_fullscale_config1_run('UMi');
```

To plot an existing result manually:

```matlab
Plot_Comm_channel_calibration_fullScale_config1('UMi');
```

The outer entry script calls both functions automatically:

```matlab
Calibration_fullscale_config1_run(selected_scenario);
Plot_Comm_channel_calibration_fullScale_config1(selected_scenario);
```

Use the outer entry script when both `.mat` and `.png` files are required.

## Scenario Selection

The communication scenarios used by these calibration runners are:

- `UMa`
- `UMi`
- `InH`

Scenario objects are created through the `scenario_configs` list in each runner. For example:

```matlab
scenario_configs = {
    localScenarioConfig('UMa', @() comm_scenario.UMa(), 'UMa')
    localScenarioConfig('UMi', @() comm_scenario.UMi(), 'UMi')
    localScenarioConfig('InH', @() localCreateInHScenario(), 'InH')
};
```

To run one selected scenario directly:

```matlab
Calibration_largescale_comm_run('UMa');
Plot_Comm_channel_calibration_largeScale('UMa');
```

The outer entry scripts also recognize the environment variable:

```text
CALIB_SELECTED_SCENARIO
```

## Frequency Configuration

The calibration frequency sets are defined in each runner.

Large-scale calibration uses:

```text
6 GHz, 30 GHz, 70 GHz
```

Full-scale Config 1 and Config 2 use:

```text
6 GHz, 30 GHz, 60 GHz, 70 GHz
```

Each frequency configuration includes the carrier frequency, bandwidth, and calibration tag. For example:

```matlab
localFrequencyConfig(6e9, 20e6, '6GHz')
localFrequencyConfig(30e9, 100e6, '30GHz')
```

Config 2 also supports optional frequency filtering when calling its runner:

```matlab
Calibration_fullscale_config2_run('InH', {'6GHz', '30GHz'});
```

## Parameters That Can Be Manually Adjusted

Communication scenario parameters are defined in `+comm_scenario/`.

Common parameters include:

- `frequency`
- `BW`
- `ISD`
- `x_range`
- `y_range`
- `BS_height`
- `BS_Tx_power`
- `BS_noise_figure`
- `BS_UE_min_d`
- `UE_TX_power`
- `UE_height`
- `UE_noise_figure`
- `UE_per_sec`
- `BS_sec_num`
- `UE_sec_num`

Shared defaults and frequency preset behavior are defined in:

```text
ISAC_channel/+comm_scenario/Comm_Scenario.m
```

Scenario-specific defaults are defined in files such as:

```text
ISAC_channel/+comm_scenario/UMa.m
ISAC_channel/+comm_scenario/UMi.m
ISAC_channel/+comm_scenario/InH.m
```

Communication-channel behavior is implemented in:

```text
ISAC_channel/+channel/Comm_channel.m
```

Placement behavior is implemented in:

```text
ISAC_channel/+network_layout/Drop_BaseStation.m
ISAC_channel/+network_layout/Drop_UE_ISAC.m
```

## Output

The runner saves one `.mat` file per scenario and frequency. The plotter saves the corresponding `.png` figures in the same result directory.

Expected output directories are:

```text
Calibration_largescale_comm/results/comm_large_scale/
Calibration_fullscale_comm/Calibration_full_config1/results/comm_full_scale_config1/
Calibration_fullscale_comm/Calibration_full_config2/results/comm_full_scale_config2/
```

Typical output names include:

```text
UMa_6GHz.mat
UMa_30GHz.mat
UMa_60GHz.mat
UMa_70GHz.mat
```

## Notes

- Use the outer entry script when the simulation results and figures are both required.
- Use the `_run.m` file when only numerical `.mat` results are required.
- Use the plotter directly to regenerate figures from existing `.mat` files.
- These scripts are communication-only calibration tools.
- `ISAC_channel_calibration.m` is the separate entry point for target and sensing calibration.
- Changes to `Comm_channel.m` affect the communication calibration results generated by these scripts.
