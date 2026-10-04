function paths = get_Comm_calibration_paths()
%GET_COMM_CALIBRATION_PATHS Resolve shared communication calibration paths.

paths.calibrationRoot = fileparts(mfilename('fullpath'));
paths.implementationRoot = fileparts(fileparts(paths.calibrationRoot));
paths.isacRoot = fullfile(paths.implementationRoot, 'ISAC_channel');
paths.referenceRoot = fullfile( ...
    paths.calibrationRoot, 'Comm_benchmark');
paths.resultsRoot = fullfile(paths.implementationRoot, 'Results');
paths.spatialConsistencyRoot = fullfile( ...
    paths.calibrationRoot, 'Calibration__comm_spatial_consistency');
paths.spatialBenchmark = fullfile( ...
    paths.spatialConsistencyRoot, ...
    'Phase3SpatialConsistency_v15_Ericsson.xlsx');

if ~isfolder(paths.isacRoot)
    error('ISAC channel source folder was not found: %s', paths.isacRoot);
end
if ~isfolder(paths.referenceRoot)
    error('Communication calibration reference folder was not found: %s', ...
        paths.referenceRoot);
end
end
