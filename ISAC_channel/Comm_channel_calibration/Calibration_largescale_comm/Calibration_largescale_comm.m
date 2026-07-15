clear;
clc;
close all;

% Entry script: choose a scenario, run large-scale calibration, then plot results.
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
project_dir = fileparts(fileparts(script_dir));
addpath(project_dir);
original_dir = pwd;
cleanup_dir = onCleanup(@() cd(original_dir)); %#ok<NASGU>
cd(script_dir);

selected_scenario = 'UMi';
scenario_override = strtrim(getenv('CALIB_SELECTED_SCENARIO'));
if ~isempty(scenario_override)
    selected_scenario = scenario_override;
end

% Generate simulation data first, then load the saved MAT files for plotting.
Calibration_largescale_comm_run(selected_scenario);
Plot_Comm_channel_calibration_largeScale(selected_scenario);
