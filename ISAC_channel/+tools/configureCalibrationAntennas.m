function equipment_list = configureCalibrationAntennas( ...
        equipment_list, polarization_angles, equipment_label)
%CONFIGURECALIBRATIONANTENNAS Apply the ISAC calibration antenna setup.
% Build one isotropic antenna location with two orthogonal polarizations.
% TRPs use +/-45 degrees and UTs use 0/90 degrees according to the
% calibration polarization convention in Clause 7.8.3.

validateattributes(polarization_angles, {'numeric'}, ...
    {'real', 'finite', 'vector', 'numel', 2}, mfilename, ...
    'polarization_angles');

for equipment_idx = 1:numel(equipment_list)
    equipment = equipment_list(equipment_idx);
    cfg = equipment.antenna_params;
    cfg.array.Mg = 1;
    cfg.array.Ng = 1;
    cfg.panel.M = 1;
    cfg.panel.N = 1;
    cfg.panel.Kv = 1;
    cfg.panel.Kh = 1;
    cfg.panel.P = 2;
    cfg.panel.X_pol = polarization_angles;
    cfg.antenna_model = 'isotropic';
    cfg.pol_model = 'model-2';
    equipment.antenna_params = cfg;

    for sector_idx = 1:numel(equipment.sector)
        old_antenna = equipment.sector(sector_idx).antenna;
        ang.alpha = old_antenna.alpha;
        ang.beta = old_antenna.beta;
        ang.gamma = old_antenna.gamma;
        attached_device = old_antenna.attachedDevice;
        attached_type = old_antenna.attachedType;

        new_antenna = antennas.antenna_array(cfg, ang);
        new_antenna.attachedDevice = attached_device;
        new_antenna.attachedType = attached_type;
        equipment.sector(sector_idx).antenna = new_antenna;
    end

    panel = equipment.sector(1).antenna.panel{1, 1};
    if panel.P ~= 2 || panel.num_element ~= 2
        error('ISACCalibration:InvalidDualPolarizedAntenna', ...
            'Calibration %s antenna must contain two polarizations.', ...
            equipment_label);
    end
end
end
