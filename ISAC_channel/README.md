# Implement ISAC Channel

MATLAB implementation scaffold for ISAC channel generation and calibration experiments.

## Entry Points

- `ISAC_channel.m` builds the baseline ISAC scenario and leaves BS, UE, ST, and channel link objects in the workspace for extension.
- `ISAC_channel_calibration.m` runs one calibration case at a time and plots the generated results against calibration benchmark sheets through `plot_TRP_mono_result.m`.

## Frequency Presets

Frequency presets are independent of the entry point. Use either `ISAC_FR1` or `ISAC_FR2` with either script, depending on whether you want to build a baseline channel workspace or run a calibration case.

For `ISAC_channel.m`, edit the control variables near the top of the script:

```matlab
isac_frequency_preset = 'ISAC_FR1';  % 6 GHz, 100 MHz bandwidth
```

or:

```matlab
isac_frequency_preset = 'ISAC_FR2';  % 30 GHz, 400 MHz bandwidth
```

For `ISAC_channel_calibration.m`, edit the calibration control block:

```matlab
case_id = 2;
isac_frequency_preset = 'ISAC_FR1';
```

or:

```matlab
case_id = 2;
isac_frequency_preset = 'ISAC_FR2';  % 30 GHz, 400 MHz bandwidth
```

To run all calibration cases at both FR1 and FR2, set:

```matlab
run_all_calibration_cases = true;
calibration_frequency_presets = {'ISAC_FR1', 'ISAC_FR2'};
```

Custom parameters can be supplied in the same control block:

```matlab
isac_frequency_preset = 'custom';
custom_frequency_config.frequency = 28e9;
custom_frequency_config.BW = 2e8;
custom_frequency_config.BS_Tx_power = 40;
custom_frequency_config.UE_TX_power = 21;
custom_frequency_config.UE_height = 1.7;
custom_frequency_config.BS_height = 24;
```

## Custom Layout Positions

`ISAC_channel.m` accepts optional custom BS, UE, and ST positions. Use an `N x 2` matrix for `[x y]` or an `N x 3` matrix for `[x y z]`. If height is omitted, the scenario or sensing-type default height is used.

```matlab
custom_bs_position = [0 0 25];
custom_ue_position = [100 20 1.5; 120 -30 1.5];
custom_st_position = [80 10 1.7; 90 15 1.7];
```

`ISAC_channel_calibration.m` intentionally uses generated layouts for calibration runs.

For indoor scenarios, ST placement currently uses `drop_in_Indoor_convex` in `+network_layout/Drop_ST.m` so that targets are generated inside the convex hull of the TRP deployment. This is the default setting used for specification-aligned runs. To align with other validation curves that assume uniform placement over the whole indoor area, switch the indoor ST placement call from `drop_in_Indoor_convex` to `drop_in_Indoor_uniform`.

## Calibration Cases

`ISAC_channel_calibration.m` uses `case_id` to select one case:

1. `UAV-UMa-AV`
2. `Human Outdoor-UMa`
3. `Human Outdoor-UMi`
4. `Human Indoor-InH`
5. `Human Indoor-InF-SH`
6. `AGV-InF-SH`
7. `Auto-Urban Grid`

Example:

```matlab
case_id = 2;
isac_frequency_preset = 'ISAC_FR2';
```

## ZOA Filter

The ZOA filter threshold for RP generation is controlled in `+channel/Comm_channel.m`, inside `large_scale_para(obj)`:

```matlab
function large_scale_para(obj)
    fc_GHz = obj.fc/1e9; %#ok<*PROPLC>
    obj.drop_rp_angle = 0;
    switch obj.scenario.name
        case {'UMa','UrbanGrid'}
            obj.drop_rp_angle = 80;
        ...
    end
```

To adjust the filter, edit `obj.drop_rp_angle` for each scenario case in that switch block. For example, changing the value under `{'UMa','UrbanGrid'}` changes the ZOA filter angle used by those scenarios. Use `0` when no angle filtering is intended for a case.

## Requirements

- MATLAB with object-oriented package folder support.
- Toolbox functions used by the channel model, including `randsrc`, `pow2db`, and `db2pow`.
- Calibration benchmark Excel files under the `ISAC calibration` directory expected by `plot_TRP_mono_result.m`, unless a custom `calibration_root` is supplied.

## Notes

- `ISAC_channel.m` is intended as a clean starting point for new experiments.
- `ISAC_channel_calibration.m` is intended for benchmark comparison and plotting.
