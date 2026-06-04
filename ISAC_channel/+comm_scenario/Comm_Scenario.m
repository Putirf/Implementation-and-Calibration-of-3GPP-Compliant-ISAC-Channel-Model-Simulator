classdef Comm_Scenario < handle

    properties
        % setup at each scenario
        name   
        subname = 'general'
        BW                  % bandwidth
        frequency           % frequency
        x_range             % range in x axis
        y_range             % range in y axis
        ISD                 % Intersite Distance
        % large scale parameter cross corelation
        Cross_correlation_LOS
        Cross_correlation_NLOS
        Cross_correlation_O2I

        alternative         % 1, 2, 3 TR 36.777  only for UAV


        BS_height  
        BS_Tx_power  
        BS_noise_figure
        BS_UE_min_d            % min distance between BS and ST    
        BS_ST_min_d            % min distance between BS and ST     
        total_BS_sector_num

        UE_TX_power
        UE_height
        UE_noise_figure
        UE_max_d_2D_indoor      % TR 38.901 P40
        UE_ST_min_d            % min distance between UE and ST    
        total_UE_sector_num

        % general setup
        BS_sec_num
        UE_sec_num
        UE_per_sec
        ST_per_cell        

        N0_dBm              % noise spectral density

        option              % ISAC channel step 9

        best_N              % pick up the best N-th STX-SRX pairs with the smallest power scaling factor of the target channel. 

        % RP 
        RP_per_equipment
        a_d
        b_d
        c_d
        a_h
        b_h
        c_h



    end
    
    methods   
        function obj = Comm_Scenario()
            obj.N0_dBm = -174;
            obj.UE_ST_min_d = 1;           % except for any -AV case
            obj.BS_sec_num = 1;
            obj.UE_sec_num = 1;
            obj.BS_noise_figure = 5;    % dB
            obj.UE_per_sec = 10;
            obj.ST_per_cell = 200;
            obj.option = 2;             % 1:  Each NLOS ray in the STX-SPST link is coupled with each NLOS ray in the SPST-SRX link. 
                                        % 2: The NLOS rays in the STX-SPST link are 1-by-1 randomly coupled with the NLOS rays in the STSRX link.
            obj.best_N = 4;             % except for Highway

            obj.UE_TX_power = 23;
            obj.UE_height   = 1.5;
            obj.UE_noise_figure = 9;
            obj.UE_ST_min_d = 1;
            obj.BW = 1e8;
            obj.frequency = 6e9;

            obj.RP_per_equipment = 3;
           
        end

        function applied_config = applyIsacFrequencyPreset(obj, preset_name, custom_config)
            if nargin < 2 || isempty(preset_name)
                preset_name = 'ISAC_FR1';
            end
            if nargin < 3 || isempty(custom_config)
                custom_config = struct();
            end

            preset_key = upper(string(preset_name));
            switch preset_key
                case "ISAC_FR1"
                    applied_config = obj.isacFrequencyConfig(6e9);
                case "ISAC_FR2"
                    applied_config = obj.isacFrequencyConfig(30e9);
                case "CUSTOM"
                    applied_config = struct();
                otherwise
                    error('Unknown ISAC frequency preset: %s', string(preset_name));
            end

            applied_config = obj.mergeIsacConfig(applied_config, custom_config);
            obj.applyIsacConfigFields(applied_config);
            obj.refreshIsacLayoutRange();
        end
    end

    methods(Access = private)
        function config = isacFrequencyConfig(obj, frequency)
            config = struct();
            config.frequency = frequency;

            if frequency == 6e9
                config.BW = 1e8;
                config.BS_noise_figure = 5;
                config.UE_noise_figure = 9;
                if ismember(obj.name, {'InH', 'InF'})
                    config.BS_Tx_power = 24;
                else
                    config.BS_Tx_power = 56;
                end
                if strcmp(obj.name, 'UrbanGrid')
                    config.ISD = 500;
                end
            elseif frequency == 30e9
                config.BW = 4e8;
                config.BS_noise_figure = 7;
                config.UE_noise_figure = 10;
                if ismember(obj.name, {'InH', 'InF'})
                    config.BS_Tx_power = 23;
                else
                    config.BS_Tx_power = 41;
                end
                if strcmp(obj.name, 'UrbanGrid')
                    config.ISD = 250;
                end
            else
                error('Unsupported ISAC preset frequency: %.0f Hz', frequency);
            end
        end

        function config = mergeIsacConfig(obj, config, custom_config) %#ok<INUSL>
            custom_fields = fieldnames(custom_config);
            for field_idx = 1:numel(custom_fields)
                field_name = custom_fields{field_idx};
                config.(field_name) = custom_config.(field_name);
            end
        end

        function applyIsacConfigFields(obj, config)
            config_fields = fieldnames(config);
            for field_idx = 1:numel(config_fields)
                field_name = config_fields{field_idx};
                if ~isprop(obj, field_name)
                    error('Unsupported ISAC frequency config field: %s', field_name);
                end
                obj.(field_name) = config.(field_name);
            end
        end

        function refreshIsacLayoutRange(obj)
            if ismember(obj.name, {'UMa', 'UMi'})
                layer_count = obj.('layer_num');
                cell_radius = obj.ISD / sqrt(3);
                obj.x_range = [-ceil(cell_radius * (1.5 * layer_count - 0.5)), ...
                    ceil(cell_radius * (1.5 * layer_count - 0.5))];
                obj.y_range = [-ceil(obj.ISD * (layer_count - 0.5)), ...
                    ceil(obj.ISD * (layer_count - 0.5))];
            elseif strcmp(obj.name, 'UrbanGrid')
                layer_count = obj.('layer_num');
                obj.x_range = [-obj.ISD * (2 * layer_count - 1) / 2, ...
                    obj.ISD * (2 * layer_count - 1) / 2];
                obj.y_range = [-(obj.ISD * 433 / 250) * (2 * layer_count - 1) / 2, ...
                    (obj.ISD * 433 / 250) * (2 * layer_count - 1) / 2];
            end
        end
    end
end
