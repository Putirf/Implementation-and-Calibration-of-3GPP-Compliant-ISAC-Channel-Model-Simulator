classdef Equipment < handle

    properties
        type
        ID
        scenario
        frequency
        BW
        Power
        sector_num
        inital_Position
        Position
        height
        sector
        antenna_params

        % BS
        cluster_wrapped

        % UE
        rand_LoS
        Indoor
        d_2D_in
        floor_num
        fcin
        n_fl
        O2IPL

        velocity
        theta_v
        phi_v

        % RP
        RP

        O2Isigma
        carPL

        PHI_n_m

        % Spatial consistency raw LSP Gaussian variables, indexed by peer ID.
        LSP_raw_LOS
        LSP_raw_NLOS
        LSP_raw_O2I
        SC_procB_raw
        SC_procA_comm_Xn
        SC_procA_target_Xn
        SC_procA_target_tx_Xn
        SC_procA_target_rx_Xn
        SC_procA_state
        SC_procA_comm_state
        SC_procA_target_state
        SC_time_nodes


    end
    methods
        function obj = Equipment(scenario,Type)
            obj.type = Type;
            obj.scenario = scenario.name;
            obj.frequency = scenario.frequency;
            obj.BW = scenario.BW;
            if strcmp(Type,'BS')
                obj.Power = scenario.BS_Tx_power;
                obj.sector_num = scenario.BS_sec_num;
                obj.height    = scenario.BS_height;


                % wrapping
                if ismember(scenario.name,{'RMa','UMa','UMi'})
                    if scenario.layer_num == 1
                        obj.cluster_wrapped = [0, 0];
                        return;
                    end
                    L = scenario.layer_num;
                    i = L;
                    j = L - 1;

                    e1 = [3/2, sqrt(3)/2];
                    e2 = [0,   sqrt(3)];

                    v1 = i*e1 + j*e2;
                    v2 = (i+j)*e1 - i*e2;
                    v3 = j*e1 - (i+j)*e2;

                    obj.cluster_wrapped = [
                        0,   0;
                        v1;     v2;     v3;
                        -v1;    -v2;    -v3] * scenario.R;
                end
                obj.antenna_params.array.Mg           = 1;                 %  the number of antenna panels with the same polarization in each column.
                obj.antenna_params.array.Ng           = 1;
                obj.antenna_params.array.dg_H         = 2.5;
                obj.antenna_params.array.dg_V         = 2.5;
                obj.antenna_params.panel.M            = 1;                 %  the number of antenna elements with the same polarization in each column.
                obj.antenna_params.panel.N            = 1;
                obj.antenna_params.panel.Kv           = 1;                  %  K = M or 1
                obj.antenna_params.panel.Kh           = 1;
                obj.antenna_params.panel.d_H          = 0.5;
                obj.antenna_params.panel.d_V          = 0.5;
                obj.antenna_params.panel.P            = 2;
                obj.antenna_params.panel.X_pol        = [45, -45];
                obj.antenna_params.panel.ele_downtilt = 90;
                obj.antenna_params.panel.ele_panning  = 0;
                obj.antenna_params.pol_model = 'model-2';
                if obj.sector_num == 1
                    obj.antenna_params.antenna_model = 'isotropic';
                else
                    obj.antenna_params.antenna_model = 'dipole';
                end
            elseif strcmp(Type,'UE')
                obj.Power = scenario.UE_TX_power;
                obj.sector_num = scenario.UE_sec_num;
                obj.velocity= 3/3.6;   % 3 km/h
                obj.theta_v = 90;
                obj.phi_v   = rand*360-180;
                obj.Indoor  = false;
                obj.n_fl    = 1;
                obj.height    = scenario.UE_height;

                obj.antenna_params.array.Mg           = 1;                 %  the number of antenna panels with the same polarization in each column.
                obj.antenna_params.array.Ng           = 1;
                obj.antenna_params.array.dg_H         = 2.5;
                obj.antenna_params.array.dg_V         = 2.5;
                obj.antenna_params.panel.M            = 1;                 %  the number of antenna elements with the same polarization in each column.
                obj.antenna_params.panel.N            = 1;
                obj.antenna_params.panel.Kv           = 1;                  %  K = M or 1
                obj.antenna_params.panel.Kh           = 1;
                obj.antenna_params.panel.d_H          = 0.5;
                obj.antenna_params.panel.d_V          = 0.5;
                obj.antenna_params.panel.P            = 2;
                obj.antenna_params.panel.X_pol        = [90, 0];
                obj.antenna_params.panel.ele_downtilt = 90;
                obj.antenna_params.panel.ele_panning  = 0;
                obj.antenna_params.beta               = 90;
            end
        end
    end
end
