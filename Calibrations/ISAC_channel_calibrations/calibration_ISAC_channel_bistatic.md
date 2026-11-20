# Bistatic calibration controls

Edit the controls near the top of `calibration_ISAC_channel_bistatic.m`, then run that script in MATLAB.

```matlab
sensing_mode = 'UE_UE'; % 'TRP_TRP', 'TRP_UE', or 'UE_UE'
case_id = 7;
```

All three modes support all seven cases. A vector such as `case_id = [1, 5, 6]` runs the selected cases. `run_all_calibration_cases = true` runs all seven at both configured calibration frequencies.

| ID | Target and deployment |
|---|---|
| 1 | Small UAV / UMa-AV |
| 2 | Human / UMa |
| 3 | Human / UMi |
| 4 | Human / InH |
| 5 | Human / InF-SH |
| 6 | AGV / InF-SH |
| 7 | Vehicle / Urban Grid |

- `full_cali = false`: LS only. `true`: LS plus Full coupling loss, DS, ASD/ZSD/ASA/ZSA.
- `isac_frequency_preset = 'ISAC_FR1'` or `'ISAC_FR2'`: 6 or 30 GHz. A cell array runs multiple frequencies.
- `calibration_simulation_times = []`: use the scenario default. Set a positive integer for a diagnostic run.
- `do_background = true`: also evaluate direct STX-SRX background channels for the same selected pairs.

The driver deliberately disables spatial consistency. This is the general LS/Full calibration, not Table 7.9.6.3-1 spatial-consistency calibration.

## TRP_UE implementation

Each drop contains one target. UAV uses 30 outdoor terrestrial UTs in the centre cell at 1.5 m and one UAV at 200 m. Indoor Human uses 20 UTs at 1 m. AGV uses 30 UTs at 1.5 m, with a 1 m UT-TRP minimum 2D distance. InF Human/AGV targets are uniform over the TRP convex hull. For TRP_UE UAV, target placement enforces at least 10 m in 3D from every candidate TRP and UT. Urban Grid uses 16 pedestrian UTs and one vehicle scattering point with the existing frequency-dependent layout.

The existing UMa/UMi Human trials retain 100% outdoor UTs; UMi retains the requested 0.8 m effective UT height. These are calibration trial settings, not changes to the specification. The existing InH target-drop choice is unchanged.

All TRP-UE candidates are ranked by LS target-channel power. Best-4 are retained; Full channels are generated for those four only. Physical-leg LS and directed-leg small-scale caches keep a shared leg consistent across pairs. TRP_TRP retains its original reciprocal physical-pair selection.

`Target_channel.m` resolves references per physical leg:

- UAV: TRP-ST uses UMa-AV; terrestrial UT-ST uses UMi-AV, with UT height 1.5 m. UMi-AV uses TR 36.777 Tables B-1/2/3 and Alternative 3 (K=15 dB). FR2 reuses the model with the corresponding frequency.
- InF Human/AGV: the high TRP leg retains InF-SH; the low UT leg uses the low-height InF reference. The UT leg adopts ASA/ZSA statistics for ASD/ZSD. The AGV UT-ST leg follows the UE-monostatic forced-LOS trial and retains the SH NLOS pathloss branch; the TRP-ST leg is not forced LOS.
- Urban Grid: TRP-Vehicle uses V2B; pedestrian UT-Vehicle uses V2P.
- InH: the UT-Human leg uses the UE-UE spread-statistics update.

UAV/Human retain sigma_D=0 dB; AGV uses the single-SPST angular RCS model path by path. Full calibration uses a -40 dB path-drop threshold.

References: TR 38.901 V19.4.0 Tables 7.9.3-1/2, 7.9.6.1-1/2/4 and 7.9.6.2-1/2/4; TR 38.858 Annex A Table A.3-3; [TR 36.777 V15.0.0 Annex B](https://www.3gpp.org/ftp/Specs/archive/36_series/36.777/36777-f00.zip).

Calibration results are passed directly to plotting without automatically saving a `.mat` file. Plotting reads the `TRP-UE-6GHz` or `TRP-UE-30GHz` benchmark sheet and saves figures under the implementation's `Results` directory.

`result.endpoint_pair_idx` lists candidate [TRP index, UE index] pairs. `result.selected_link_id` indexes those candidates per drop. `result.selected_endpoint_pairs` records the selected pairs directly. The CDF samples retain the existing sorting convention, so pairing metadata describes selection before CDF sorting.

## UE_UE implementation

`sensing_mode = 'UE_UE'` models UT-ST-UT with two different terrestrial UTs. It uses the same UT populations, antenna setup and frequency presets as the corresponding TRP-UE cases. Indoor targets are uniform over the TRP convex hull. UMa/UMi Human retain the existing outdoor-only trial and UMi's 0.8 m effective-UT-height trial.

### Outdoor Human minimum distances

All distances below are horizontal (2D). UMa uses 35 m for TRP-UT and TRP-ST; UMi uses 10 m.

| Mode | TRP-UT | TRP-ST | UT-ST | UT-UT |
|---|---|---|---|---|
| TRP_UE | UT dropping | ST dropping, relative to the centre TRP | 1 m during ST dropping, against all UTs | No additional restriction |
| UE_UE | UT dropping | No exclusion | 1 m during ST dropping, against all UTs | 1 m during pair selection |

Outdoor Human TRP_UE and UE_UE do not filter UT-ST distances again during link ranking. A rejected ST position is redrawn; the TRP-ST drop exclusion is relative to the centre TRP, not every candidate TRP. The UT-UT pair-separation check remains because ST dropping does not enforce separation between UTs. These rules do not change the UAV or indoor/AGV distance handling.

Best-4 selection ranks distinct unordered UT pairs by LS target-channel power. Reverse directions do not occupy two ranks; STX/SRX roles alternate by rank and drop, as for the TRP-TRP balanced-pair selection. The driver obtains each physical UT-ST large-scale state once, scores the pairs by the sum of their leg contributions, then generates Full channels for the selected four only. For UAV, UT placement enforces at least 1 m of horizontal spacing between all UTs, and target placement enforces at least 10 m in 3D from every UT. The other cases retain the 1 m UT-pair filter and their existing target-to-UT limits.

The two UT-ST legs use Table 7.9.3-1 independently: Case 7 (UMi-AV, Alternative 3) for UAV, Case 5 (UMi or indoor reference) for Human/AGV, and Case 6 (Urban V2P) for Vehicle. The low InF UT reference is SL and uses the UE-UE spread-statistics update. The forced-LOS AGV trials in UE_mono/TRP-UE remain limited to those modes; UE_UE samples LOS/NLOS using the low-height InF reference. UAV calibration uses the UMa absolute-delay model from Table 7.9.6.2-1. AGV/Vehicle angular RCS remains path-dependent.

`do_background = true` generates a direct UT-UT channel without target RCS. Outdoor UMa/UMi and UAV terrestrial-UT backgrounds use UMi UE-UE; indoor backgrounds use the indoor low-height reference; Urban Grid uses P2P with the V2V propagation and fast-fading equations reused by TR 37.885 Clause 6.2. Company target curves are read from `UE-UE-6GHz` / `UE-UE-30GHz`. Missing background benchmark sheets produce simulation-only background plots, as in the existing bistatic modes.

UE_UE uses the shared progress formatter and per-link updates. It also uses the common reference-column filter, so Mean/Average/Avg are excluded. UMa/UMi Human DS uses the same four-decade reference-tail cutoff as TRP-UE, preserving the original CDF percentages and all simulation samples. `result.endpoint_pair_idx` contains [UT index, UT index] candidates; `result.selected_endpoint_pairs` contains the actual directed Best-4 pairs. UTs never transmit to themselves.

## Validation limits

UE_UE was checked at 6/30 GHz for all seven cases with two drops per combination, including target/background coupling loss, DS, all four angular spreads, distinct Best-4 selection, pair metadata and benchmark-sheet loading. An eight-UT exhaustive-pair comparison checked the optimized LS ranking, reciprocal LS scores, LS/Full selection consistency and progress completion. The existing TRP-UE cases were checked at both frequencies, with additional monostatic/TRP-TRP checks for UAV, UMi Human and InF AGV. These checks do not establish agreement of full-size CDFs with the company references.

A successful run or finite output is not a calibration acceptance criterion. Compare the full CDFs with the company curves, including DS and all four angular spreads. The initial comparison report is under `Results/TRP_UE_initial_validation`.
