function paths = get_ISAC_calibration_paths()
%GET_ISAC_CALIBRATION_PATHS Resolve shared ISAC calibration locations.

paths.calibrationDir = fileparts(mfilename('fullpath'));
paths.implementationRoot = fileparts(fileparts(paths.calibrationDir));
paths.isacRoot = fullfile(paths.implementationRoot, 'ISAC_channel');
paths.referenceRoot = fullfile( ...
    paths.calibrationDir, 'ISAC_benchmark');
paths.resultsRoot = fullfile(paths.implementationRoot, 'Results');

if ~isfolder(paths.isacRoot)
    error('ISAC channel source folder was not found: %s', paths.isacRoot);
end
if ~isfolder(paths.referenceRoot)
    error('ISAC calibration reference folder was not found: %s', ...
        paths.referenceRoot);
end
end
