# Implement ISAC Channel

MATLAB implementation for Integrated Sensing and Communication (ISAC) channel generation and calibration experiments.

The project has two main workflows:

- **Calibration**: run predefined calibration cases and compare generated results with benchmark data.
- **Using ISAC Channel**: build an editable ISAC channel scenario, including custom BS, UE, and ST positions.

## Project Structure

```text
ISAC_channel.m                  Main script for building an ISAC channel scenario
ISAC_channel_calibration.m      Main script for calibration experiments
plot_TRP_mono_result.m          Calibration result plotting utility
+comm_scenario/                 Communication scenario definitions
+sensing_types/                 Sensing target type definitions
+network_layout/                BS, UE, ST, and RP placement functions
+channel/                       Communication and target channel models
+antennas/                      Antenna element, panel, and array models
+elements/                      Equipment, sector, and target objects
```

## Requirements

- MATLAB with object-oriented package folder support.
- MATLAB functions used by the channel model, including `randsrc`, `pow2db`, and `db2pow`.
- Calibration benchmark Excel files under the `ISAC calibration` directory when running calibration, unless a custom `calibration_root` is supplied.

## Calibration

Use `ISAC_channel_calibration.m` for calibration runs.

The script selects one predefined calibration case through `case_id`:

| `case_id` | Calibration case |
| --- | --- |
| `1` | `UAV-UMa-AV` |
| `2` | `Human Outdoor-UMa` |
| `3` | `Human Outdoor-UMi` |
| `4` | `Human Indoor-InH` |
| `5` | `Human Indoor-InF-SH` |
| `6` | `AGV-InF-SH` |
| `7` | `Auto-Urban Grid` |

The main calibration controls are near the top of `ISAC_channel_calibration.m`:

```matlab
case_id = 7;
full_cali = false;
do_rp = false;
plot_controller = true;
calibration_root = [];

use_isac_frequency_preset = true;
isac_frequency_preset = 'ISAC_FR1';
custom_frequency_config = struct();

run_all_calibration_cases = false;
calibration_frequency_presets = {'ISAC_FR1', 'ISAC_FR2'};
```

### Frequency Presets

The calibration script supports the same ISAC frequency presets as the normal channel script:

- `ISAC_FR1`: 6 GHz carrier frequency and 100 MHz bandwidth.
- `ISAC_FR2`: 30 GHz carrier frequency and 400 MHz bandwidth.
- `custom`: user-defined frequency-related settings through `custom_frequency_config`.

For a single calibration case, set:

```matlab
case_id = 2;
isac_frequency_preset = 'ISAC_FR1';
```

To run every calibration case at multiple frequency presets, set:

```matlab
run_all_calibration_cases = true;
calibration_frequency_presets = {'ISAC_FR1', 'ISAC_FR2'};
```

For custom settings:

```matlab
isac_frequency_preset = 'custom';
custom_frequency_config.frequency = 28e9;
custom_frequency_config.BW = 2e8;
custom_frequency_config.BS_Tx_power = 40;
custom_frequency_config.UE_TX_power = 21;
custom_frequency_config.UE_height = 1.7;
custom_frequency_config.BS_height = 24;
```

Only fields that exist as properties of the scenario object can be supplied in `custom_frequency_config`.

### Calibration Notes

- Calibration uses generated layouts from the project placement functions.
- Custom BS, UE, and ST positions are intended for `ISAC_channel.m`, not for calibration runs.
- `full_cali = true` enables additional full-channel calibration outputs such as full coupling loss, delay spread, and angular spread.
- `do_rp = true` additionally generates and plots RP-related calibration results.
- `calibration_root` can be used to point `plot_TRP_mono_result.m` to a custom benchmark data directory.

## Using ISAC Channel

Use `ISAC_channel.m` to build a general ISAC channel workspace.

This script creates:

- BS objects through `network_layout.Drop_BaseStation`.
- UE objects through `network_layout.Drop_UE_ISAC`.
- ST objects through `network_layout.Drop_ST`.
- UE communication links through `channel.Comm_channel`.
- ST target links through `channel.Target_channel`.
- Optional RP links through `network_layout.Drop_RP` and `channel.Comm_channel`.

After execution, the main objects remain in the MATLAB workspace, including:

- `scenario`
- `sensing_type`
- `BS_list`
- `UE_list`
- `ST_list`
- `link_list_UE`
- `link_list_ST`
- `link_list_RP`

### Scenario Selection

Select the communication scenario near the top of `ISAC_channel.m`:

```matlab
scenario = comm_scenario.UMa;
```

Available scenario classes are defined in `+comm_scenario/`, including:

- `UMa`
- `UMi`
- `RMa`
- `InH`
- `InF`
- `UrbanGrid`
- `HighWay`

Some scenarios support subtypes. For example:

```matlab
scenario = comm_scenario.InF('SH');
```

### Sensing Type Selection

Select the sensing target type in `ISAC_channel.m`:

```matlab
sensing_type = sensing_types.Human;
```

Available sensing types are defined in `+sensing_types/`, including:

- `Human`
- `UAV`
- `Vehicle`
- `AGV`

Some sensing types can also modify scenario properties. For example, `sensing_types.UAV(scenario)` sets UAV-related scenario options and target height.

### Custom BS, UE, and ST Positions

`ISAC_channel.m` supports custom BS, UE, and ST positions through:

```matlab
custom_bs_position = [];
custom_ue_position = [];
custom_st_position = [];
```

Each custom position matrix can be:

- `N x 2`: `[x y]`
- `N x 3`: `[x y z]`

If the height column is omitted, the script uses the default height from the scenario or sensing type:

- BS height: `scenario.BS_height`
- UE height: `scenario.UE_height`
- ST height: `sensing_type.height`

These custom positions are passed into:

- `+network_layout/Drop_BaseStation.m`
- `+network_layout/Drop_UE_ISAC.m`
- `+network_layout/Drop_ST.m`

### Parameters That Can Be Manually Adjusted

Most experiment parameters are controlled by class properties or placement functions. Edit the corresponding files directly when a run needs different assumptions.

#### Communication Scenario Parameters

Scenario-level parameters are defined in `+comm_scenario/`.

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
- `BS_ST_min_d`
- `UE_TX_power`
- `UE_height`
- `UE_noise_figure`
- `UE_per_sec`
- `ST_per_cell`
- `BS_sec_num`
- `UE_sec_num`
- `best_N`
- `RP_per_equipment`

Shared defaults and frequency preset logic are in:

```text
+comm_scenario/Comm_Scenario.m
```

Scenario-specific defaults are in files such as:

```text
+comm_scenario/UMa.m
+comm_scenario/UMi.m
+comm_scenario/InH.m
+comm_scenario/InF.m
+comm_scenario/UrbanGrid.m
```

#### Sensing Target Parameters

Sensing target parameters are defined in `+sensing_types/`.

Common adjustable parameters include:

- `height`
- `RCS_model`
- `RCS`
- `XPR`
- `sigle_STSP`
- `multi_STSP`
- `k1`
- `k2`

Target-specific files include:

```text
+sensing_types/Human.m
+sensing_types/UAV.m
+sensing_types/Vehicle.m
+sensing_types/AGV.m
```

#### Layout and Drop Parameters

Placement behavior is defined in `+network_layout/`.

Important files include:

```text
+network_layout/Drop_BaseStation.m
+network_layout/Drop_UE_ISAC.m
+network_layout/Drop_ST.m
+network_layout/Drop_RP.m
```

For indoor ST placement, `+network_layout/Drop_ST.m` currently uses:

```matlab
network_layout.drop_in_Indoor_convex(...)
```

This places indoor targets inside the convex hull of the TRP deployment. To use uniform indoor placement across the whole indoor area instead, switch the call to:

```matlab
network_layout.drop_in_Indoor_uniform(...)
```

The uniform call is already present as a commented alternative in `Drop_ST.m`.

#### Channel Parameters

Communication and target channel behavior is implemented in:

```text
+channel/Comm_channel.m
+channel/Target_channel.m
```

The ZOA filter threshold for RP generation is controlled in `+channel/Comm_channel.m`, inside `large_scale_para(obj)`:

```matlab
obj.drop_rp_angle = 0;
switch obj.scenario.name
    case {'UMa','UrbanGrid'}
        obj.drop_rp_angle = 80;
    ...
end
```

Adjust `obj.drop_rp_angle` in each scenario case to change the ZOA filter angle. Use `0` when no angle filtering is intended.

#### Antenna Parameters

Antenna behavior is implemented in:

```text
+antennas/antenna_element.m
+antennas/antenna_panel.m
+antennas/antenna_array.m
```

Equipment objects create antennas through the antenna parameter structures used by `elements.Equipment` and `elements.Sector`.

## Notes

- Use `ISAC_channel_calibration.m` when the goal is benchmark comparison.
- Use `ISAC_channel.m` when the goal is custom scenario construction or manual channel experimentation.
- Calibration and normal ISAC channel generation share scenario, sensing type, channel, antenna, and layout classes, but they use different entry scripts and control flows.
