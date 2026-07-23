classdef Comm_channel < handle
    % LINK
    properties
        fastfading_enable = false
        TX
        TX_pos_wrap
        RX
        O2I = false
        scenario
        fc
        d_2D
        d_2D_in
        d_2D_out
        d_3D
        t

        phi_LOS_AOD
        theta_LOS_ZOD
        phi_LOS_AOA
        theta_LOS_ZOA
        LOS

        PL
        sigma_SF
        CouplingLoss

        SF
        K
        DS
        ASD
        ASA
        mu_lgZSD
        sigma_lgZSD
        mu_offset_ZOD
        ZSD
        ZSA
        r_tau
        mu_XPR
        sigma_XPR
        N        % Number of clusters
        M        % Number of rays per cluster
        zeta     % Per cluster shadowing std
        c_DS
        c_ASA
        c_ZSA
        c_ASD
        c_ZSD
        drop_rp_angle
        tau_absolute
        delta_tau

        tau_order
        tau_n
        tau_n_LOS
        Pn
        Pn_LOS
        cluster_shadow_dB
        keep
        N_new
        strong_cluster_id
        map_delay

        phi_n_m_AOA
        phi_n_m_AOD
        theta_n_m_ZOA
        theta_n_m_ZOD

        XPR_n_m
        PHI_n_m
        PHI_LOS

        loss_blockage
        tau_prime
        phi_prime_AOA
        phi_prime_AOD
        theta_prime_ZOA
        theta_prime_ZOD

        S
        U
        H_cali
        H_fastfading
        H_full
        tau_channel

        isRP = false
    end

    methods
        function obj = Comm_channel(TX, RX, scenario, fastfading_enable,t,varargin)
            obj.t = t;
            obj.TX = TX;
            obj.RX = RX;
            obj.scenario = scenario;
            obj.fc = scenario.frequency;
            obj.fastfading_enable = fastfading_enable;
            if RX.Indoor
                obj.O2I = true;
            end
            if strcmp(obj.RX.type,'RP')   % UE??
                obj.TX_pos_wrap = TX.Position;
                obj.isRP  = true;
            else
                if ismember(scenario.name,{'RMa','UMa','UMi'})
                    wrapped_positions = repmat(TX.Position(1:2), size(TX.cluster_wrapped,1), 1) + TX.cluster_wrapped;
                    [~, wrapped_site_idx] = min(sum(((repmat(RX.Position(:,1:2), size(wrapped_positions,1), 1)-wrapped_positions).^2),2));
                    obj.TX_pos_wrap = wrapped_positions(wrapped_site_idx,:);
                elseif ismember(scenario.name,{'InH','InF','UrbanGrid'})
                    obj.TX_pos_wrap = TX.Position;
                end
            end


            %% step 1 calculate GLS distance and angle
            x = RX.Position(1)-obj.TX_pos_wrap(1);
            y = RX.Position(2)-obj.TX_pos_wrap(2);
            z = RX.height-TX.height;
            obj.d_2D = sqrt(x.^2+y.^2);
            obj.d_3D = sqrt(x.^2+y.^2+z.^2);

            if ismember(scenario.name,{'RMa','UMa','UMi','UrbanGrid'})
                obj.d_2D_in  = RX.d_2D_in;
                obj.d_2D_out = obj.d_2D - obj.d_2D_in;
            elseif ismember(scenario.name,{'InH','InF'})
                obj.d_2D_in  = obj.d_2D;
                obj.d_2D_out = 0;
            end

            [obj.phi_LOS_AOD, obj.theta_LOS_ZOD] = cart2sph_(x, y, z);
            [obj.phi_LOS_AOA, obj.theta_LOS_ZOA] = cart2sph_(-x, -y, -z);

            %% step 2
            obj.los_probability();

            %% step 3
            obj.pathloss();

            %% step 4
            obj.large_scale_para();

            % the following steps is needed only when fast fading is enable
            if fastfading_enable
                if obj.isSpatialConsistencyProcedureA()
                    obj.assertProcedureASupported();
                    if obj.t <= 0 || ~obj.hasProcedureAState()
                        obj.cluster_delay_procedureA_initial();
                        obj.cluster_power();
                        obj.AOA_calc();
                        obj.AOD_calc();
                        obj.ZOA_calc();
                        obj.ZOD_calc();
                    else
                        obj.cluster_delay_procedureA();
                        obj.cluster_angles_procedureA();
                        obj.cluster_power_procedureA();
                        obj.apply_angle_offsets_procedureB();
                    end
                    obj.saveProcedureAState();
                elseif obj.isSpatialConsistencyProcedureB()
                    obj.assertProcedureBSupported();
                    %% step 5-7: Procedure B for spatial consistency
                    obj.cluster_delay_procedureB();
                    obj.cluster_angles_procedureB();
                    obj.cluster_power_procedureB();
                    obj.apply_angle_offsets_procedureB();
                else
                    %% step 5: Generate cluster delays
                    obj.cluster_delay();
                    %% step 6: Generate cluster powers
                    obj.cluster_power()
                    %% step 7: Generate arrival angles and departure angles for both azimuth and elevation
                    obj.AOA_calc();
                    obj.AOD_calc();
                    obj.ZOA_calc();
                    obj.ZOD_calc();
                end
                if strcmp(obj.RX.type,'RP')
                    obj.phi_n_m_AOA = obj.phi_n_m_AOD;
                    obj.theta_n_m_ZOA = obj.theta_n_m_ZOD;
                    drop_rp = (obj.theta_n_m_ZOA < obj.drop_rp_angle);
                    obj.phi_n_m_AOA(drop_rp) = nan;
                    obj.phi_n_m_AOD(drop_rp) = nan;
                    obj.theta_n_m_ZOA(drop_rp) = nan;
                    obj.theta_n_m_ZOD(drop_rp) = nan;                    
                end

                %% step 8: Coupling of rays within a cluster for both azimuth and elevation
                obj.RandomCouplingRays();

                %% step 9: Generate XPRs
                obj.generate_XPRs();

                % The outcome of Steps 1-9 shall be identical for all the links from co-sited sectors to a UT.

                %% step 10: Draw initial random phases
                obj.initial_random_phases();

                %% step 11 generate_channel
                obj.generate_channel();
            end
        end


        % caculator the LOS probability
        % The LOS probability is derived with assuming antenna heights of
        % 3m for indoor, 10m for UMi, and 25m for UMa
        function los_probability(obj)
            d_2Dout = obj.d_2D_out;
            height = obj.RX.height;
            switch obj.scenario.name
                case {'UMa','UrbanGrid'}
                    if height>22.5 && ~strcmp(obj.RX.type,'RP')
                        if strcmp(obj.scenario.subname,'AV')
                            if height<=100
                                d1 = max(460*log10(height)-700,18);
                                p1 = 4300*log10(height)-3800;
                                if d_2Dout <= d1
                                    Pr_LOS = 1;
                                else
                                    Pr_LOS = d1/d_2Dout + exp(-d_2Dout/p1)*(1-d1/d_2Dout);
                                end
                            elseif height<=300
                                Pr_LOS = 1;
                            else
                                error('The height is over limit of TR 36.777');
                            end
                        else
                            error('The height is not support in this scenario, change to case"AV"');
                        end
                    elseif d_2Dout <= 18
                        Pr_LOS = 1;
                    else
                        C_tmp = (max((height-13)/10, 0)).^1.5;
                        Pr_LOS = (18/d_2Dout+exp(-d_2Dout/63)*(1-18/d_2Dout))*...
                            (1+C_tmp*(5/4)*((d_2Dout/100).^3)*exp(-d_2Dout/150));
                    end

                case 'UMi'
                    Pr_LOS = min(18./d_2Dout,1).*(1-exp(-d_2Dout./36))+exp(-d_2Dout./36);

                case'RMa'
                    Pr_LOS = min(exp(-(d_2Dout-10)./1000),1);

                case'InH'
                    d_2Din = obj.d_2D_in;
                    switch obj.scenario.subname
                        case {'B','open_office'}
                            Pr_LOS = min(max(exp(-(d_2Din-5)./70.8),0.54*exp(-(d_2Din-49)./211.7)),1);
                        case 'mixed_office'
                            Pr_LOS = min(max(exp(-(d_2Din-1.2)/4.7),0.32*exp(-(d_2Din-6.5)/32.6)),1);
                    end

                case'InF'
                    switch obj.scenario.subname
                        case {'SL','DL'}
                            k_subsce = -obj.scenario.d_cluster/log(1-obj.scenario.r);
                        case {'SH','DH'}
                            k_subsce = -obj.scenario.d_cluster/log(1-obj.scenario.r)*((obj.TX.height-height)/(obj.scenario.hc - height));
                        case 'HH'
                            k_subsce = inf;
                    end
                    if obj.d_2D_in < obj.scenario.d_subsce
                        Pr_LOS = obj.scenario.p_subsce;
                    else
                        Pr_LOS = obj.scenario.p_subsce*exp(-obj.d_2D_in/k_subsce);
                    end
            end
            if strcmp(obj.RX.type,'RP')
                obj.LOS = (obj.RX.rand_LoS < Pr_LOS);
            else
                obj.LOS = (obj.RX.rand_LoS(obj.TX.ID) < Pr_LOS);
            end
        end

        % caculator the path loss distance depend
        function pathloss(obj)
            tx_height = obj.TX.height;
            rx_height = obj.RX.height;
            fc_GHz = obj.fc/1e9;    % to GHz

            switch obj.scenario.name
                case {'UMa','UrbanGrid'}
                    C = 0;
                    if rx_height <= 22.5
                        rx_breakpoint_height = rx_height;
                    else
                        rx_breakpoint_height = 1.5;
                    end
                    if rx_breakpoint_height <= 13
                        C = 0;
                    elseif (rx_breakpoint_height > 13) && (rx_breakpoint_height <= 22.5)
                        if obj.d_2D <= 18
                            g = 0;
                        else
                            g = 5/4*((obj.d_2D/100).^3)*exp(-obj.d_2D/150);
                        end
                        C = (((rx_breakpoint_height-13)/10).^1.5)*g;
                    end

                    if rand < 1/(1+C) || obj.isRP
                        environment_height = 1;
                    else
                        environment_height = 3*randi([4, obj.RX.n_fl],1);
                        %                     environment_height = randsrc(1,1,[12,15,18,21]);
                    end

                    effective_tx_height = tx_height-environment_height;
                    effective_rx_height = rx_breakpoint_height-environment_height;
                    breakpoint_distance = 4*effective_tx_height*effective_rx_height*obj.fc/3e8;

                    if rx_height >22.5 && ~strcmp(obj.RX.type,'RP')
                        if strcmp (obj.scenario.subname,'AV') && obj.d_2D <= 10*1e3 
                            if obj.LOS
                                pathloss_dB = 28 + 22*log10(obj.d_3D) + 20 * log10(fc_GHz);
                                shadow_sigma_dB = 4.64 * exp(-0.0066*rx_height);
                            else
                                pathloss_dB = -17.5 + (46 - 7*log10(rx_height)) * log10(obj.d_3D) + 20*log10(40*pi*fc_GHz/3);
                                shadow_sigma_dB = 6;
                            end
                        else
                            error('distance out of range');
                        end
                    else
                        if (obj.d_2D >= 10) && (obj.d_2D <= breakpoint_distance)
                            los_pathloss_dB = 22*log10(obj.d_3D)+28+20*log10(fc_GHz);
                        elseif (obj.d_2D > breakpoint_distance) && (obj.d_2D < 5000)
                            los_pathloss_dB = 40*log10(obj.d_3D)+28+20*log10(fc_GHz)-9*log10(breakpoint_distance.^2+(tx_height-rx_height).^2);
                        elseif strcmp(obj.scenario.name,'UrbanGrid')
                            los_pathloss_dB = 22*log10(obj.d_3D)+28+20*log10(fc_GHz);
                        end
                        if obj.LOS
                            pathloss_dB = los_pathloss_dB;
                            shadow_sigma_dB = 4;
                        else
                            nlos_pathloss_dB = 13.54 + 39.08*log10(obj.d_3D) + 20*log10(fc_GHz)-0.6*(rx_height-1.5);
                            pathloss_dB = max(los_pathloss_dB, nlos_pathloss_dB);
                            shadow_sigma_dB = 6;
                            % optional
                            % nlos_pathloss_dB = 32.4+20*log10(fc_GHz)+30*log10(obj.d_3D);
                            % shadow_sigma_dB = 7.8;
                        end
                    end

                    if obj.O2I
                        [penetration_loss_dB,penetration_sigma_dB] = obj.O2IpenetrationLoss(fc_GHz,obj.RX.O2IPL);
                        shadow_sigma_dB = 7;
                        pathloss_dB = pathloss_dB + penetration_loss_dB + 0.5*obj.d_2D_in + penetration_sigma_dB*obj.RX.O2Isigma;
                    end

                case 'UMi'
                    effective_tx_height = tx_height-1;
                    effective_rx_height = rx_height-1;
                    breakpoint_distance = 4*effective_tx_height*effective_rx_height*obj.fc/3e8;
                    if (obj.d_2D <= breakpoint_distance)% && (obj.d_2D >= 10)
                        los_pathloss_dB = 32.4 + 21*log10(obj.d_3D) + 20*log10(fc_GHz);
                    elseif  (obj.d_2D > breakpoint_distance) && (obj.d_2D < 5000)
                        los_pathloss_dB = 32.4 + 40*log10(obj.d_3D) + 20*log10(fc_GHz) - 9.5*log10(breakpoint_distance.^2+(tx_height-rx_height).^2);
                    end
                    if obj.LOS
                        pathloss_dB = los_pathloss_dB;
                        shadow_sigma_dB = 4;
                    else
                        nlos_pathloss_dB = 35.3*log10(obj.d_3D)+22.4+21.3*log10(fc_GHz)-0.3*(rx_height-1.5);
                        pathloss_dB = max(los_pathloss_dB, nlos_pathloss_dB);
                        shadow_sigma_dB = 7.82;
                    end
                    if obj.O2I
                        [penetration_loss_dB,penetration_sigma_dB] = obj.O2IpenetrationLoss(fc_GHz,obj.RX.O2IPL);
                        shadow_sigma_dB = 7;
                        pathloss_dB = pathloss_dB + penetration_loss_dB + 0.5*obj.d_2D_in + penetration_sigma_dB*obj.RX.O2Isigma;
                    end
                case 'RMa'
                    W = obj.scenario.W;
                    h = obj.scenario.h;
                    breakpoint_distance_2d = (2*pi*tx_height*rx_height*obj.fc)/3e8;
                    los_pathloss_fn = @(d_3D)20*log10(40*pi*d_3D*fc_GHz/3) + min(0.03*(h^1.72),10)*log10(d_3D)...
                        -min(0.044*(h^1.72), 14.77) + 0.002*log10(h)*d_3D;
                    if (obj.d_2D > 10) && (obj.d_2D <= breakpoint_distance_2d)
                        los_pathloss_dB = los_pathloss_fn(obj.d_3D);
                        shadow_sigma_dB = 4;
                    elseif  (obj.d_2D > breakpoint_distance_2d) && (obj.d_2D <= 10000)
                        los_pathloss_dB = los_pathloss_fn(breakpoint_distance_2d) + 40*log10(obj.d_3D/breakpoint_distance_2d);
                        shadow_sigma_dB = 6;
                    end
                    if obj.LOS
                        pathloss_dB = los_pathloss_dB;
                    else
                        pathloss_dB = 161.04-7.1*log10(W)+7.5*log10(h)-(24.37-3.7*(h/tx_height)^2)*log10(tx_height)+ ...
                            (43.42-3.1*log10(tx_height))*(log10(obj.d_3D)-3)+20*log10(fc_GHz)-(3.2*(log10(11.75*rx_height))^2-4.97);
                        shadow_sigma_dB = 8;
                    end
                    if obj.O2I
                        [penetration_loss_dB,penetration_sigma_dB] = obj.O2IpenetrationLoss(fc_GHz,obj.RX.O2IPL);
                        shadow_sigma_dB = 8;
                        pathloss_dB = pathloss_dB + penetration_loss_dB + 0.5*obj.d_2D_in + penetration_sigma_dB*obj.RX.O2Isigma;
                    end
                case 'InH'
                    pathloss_dB = 32.4 + 17.3*log10(obj.d_3D) + 20*log10(fc_GHz);
                    shadow_sigma_dB = 3;
                    if ~obj.LOS
                        nlos_pathloss_dB = 17.3 + 38.3*log10(obj.d_3D) + 24.9*log10(fc_GHz);
                        pathloss_dB = max(pathloss_dB,nlos_pathloss_dB);
                        shadow_sigma_dB = 8.03;
                    end

                case 'InF'
                    pathloss_dB = 31.84+21.5*log10(obj.d_3D)+19*log10(fc_GHz);
                    shadow_sigma_dB = 4.32;
                    if ~obj.LOS
                        switch obj.scenario.subname
                            case 'SL'
                                pathloss_dB = max(pathloss_dB,33+25.5*log10(obj.d_3D)+20*log10(fc_GHz));
                                shadow_sigma_dB = 5.7;
                            case 'DL'
                                PL_dbsl = max(pathloss_dB,33+25.5*log10(obj.d_3D)+20*log10(fc_GHz));
                                pathloss_dB = max(pathloss_dB,PL_dbsl);
                                pathloss_dB = max(pathloss_dB,18.6+35.7*log10(obj.d_3D)+20*log10(fc_GHz));
                                shadow_sigma_dB = 7.2;
                            case 'SH'
                                pathloss_dB = max(pathloss_dB,32.4+23*log10(obj.d_3D)+20*log10(fc_GHz));
                                shadow_sigma_dB = 5.9;
                            case 'DH'
                                pathloss_dB = max(pathloss_dB,33.63+21.9*log10(obj.d_3D)+20*log10(fc_GHz));
                                shadow_sigma_dB = 4;
                        end
                    end
            end

            obj.PL = pathloss_dB;
            obj.sigma_SF = shadow_sigma_dB;
        end
        function [penetration_loss_dB,penetration_sigma_dB] = O2IpenetrationLoss(obj,fc,penetration_level)
            if fc<6
                penetration_loss_dB = 20; penetration_sigma_dB = 0;
            else
                if strcmp(obj.scenario.name,'RMa')
                    penetration_level = 'low';
                elseif strcmp(obj.scenario.name,'InF')
                    penetration_level = 'high';
                end
                material_loss_dB = [2+0.2*fc;23+0.3*fc;5+4*fc;4.85+0.12*fc];
                if strcmp(penetration_level,'low')
                    penetration_loss_dB = 5-10*log10(0.3*10.^(-material_loss_dB(1,:)/10) + 0.7*10.^(-material_loss_dB(3,:)/10));
                    penetration_sigma_dB = 4.4;
                elseif strcmp(penetration_level,'high')
                    penetration_loss_dB = 5-10*log10(0.7*10.^(-material_loss_dB(2,:)/10) + 0.3*10.^(-material_loss_dB(3,:)/10));
                    penetration_sigma_dB = 6.5;
                end
            end
        end

        % caculator the large scale parameters
        function large_scale_para(obj)
            fc_GHz = obj.fc/1e9; %#ok<*PROPLC>
            obj.drop_rp_angle = 0;
            switch obj.scenario.name
                case {'UMa','UrbanGrid'}
                    obj.drop_rp_angle = 80;
                    fc_GHz(fc_GHz<=6) = 6;
                    obj.delta_tau = 10^(rand*0.2-7.5); % Table 7.6.9-1
                    is_AV = strcmp(obj.scenario.subname,'AV');
                    if is_AV
                        av_alternative = obj.scenario.alternative;
                    else
                        av_alternative = [];
                    end

                    if ~is_AV || av_alternative == 3
                        if obj.O2I
                            if obj.LOS
                                obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D_out/1000)-0.01*(obj.RX.height-1.5)+0.75);
                                obj.sigma_lgZSD = 0.4;
                                obj.mu_offset_ZOD = 0;
                            else
                                obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D_out/1000)-0.01*(obj.RX.height-1.5)+0.9);
                                obj.sigma_lgZSD = 0.49;
                                obj.mu_offset_ZOD = (7.66*log10(fc_GHz)-5.96) -10^((0.208*log10(fc_GHz)-0.782)*log10(max((25), obj.d_2D_out))+(-0.13*log10(fc_GHz)+2.03)-0.07*(obj.RX.height-1.5));
                            end
                            corr_rand_vec = obj.scenario.Cross_correlation_O2I * randn(6,1);
                            % SF/DS/ASD/ASA/ZSD/ZSA
                            lsp_std_vec = [obj.sigma_SF, 0.32, 0.42, 0.16, obj.sigma_lgZSD, 0.43]';
                            lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                            obj.SF = lsp_offset_vec(1);
                            obj.DS = 10^(lsp_offset_vec(2)-6.62);
                            obj.ASD = 10^(lsp_offset_vec(3)+1.25);
                            obj.ASA = 10^(lsp_offset_vec(4)+1.76);
                            obj.ZSD = 10^(lsp_offset_vec(5)+obj.mu_lgZSD);
                            obj.ZSA = 10^(lsp_offset_vec(6)+1.01);
                            obj.r_tau = 2.2;
                            obj.mu_XPR = 9;
                            obj.sigma_XPR = 5;  % 5
                            obj.N = 12;
                            obj.M = 20;
                            obj.zeta = 4;
                            obj.c_DS  = 11;
                            obj.c_ASA = 8;
                            obj.c_ZSA = 3;
                            obj.c_ASD = 5;
                        elseif obj.LOS % 0 1
                            obj.mu_offset_ZOD = 0;
                            corr_rand_vec = obj.scenario.Cross_correlation_LOS * randn(7,1);
                            obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)-0.01*(obj.RX.height-1.5)+0.75);
                            obj.sigma_lgZSD = 0.4;
                            % SF/K/DS/ASD/ASA/ZSD/ZSA
                            lsp_std_vec = [obj.sigma_SF, 3.5, 0.66, 0.28, 0.20, obj.sigma_lgZSD, 0.16]';
                            lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                            obj.SF = lsp_offset_vec(1);
                            obj.K = lsp_offset_vec(2)+9;
                            if is_AV && av_alternative == 3
                                obj.K = 15;
                            end
                            obj.DS = 10^(lsp_offset_vec(3)+(-6.955-0.0963*log10(fc_GHz)));
                            obj.ASD = 10^(lsp_offset_vec(4)+(1.06+0.1114*log10(fc_GHz)));
                            obj.ASA = 10^(lsp_offset_vec(5)+1.81);
                            obj.ZSD = 10^(lsp_offset_vec(6)+obj.mu_lgZSD);
                            obj.ZSA = 10^(lsp_offset_vec(7)+0.95);
                            obj.r_tau = 2.5;
                            obj.mu_XPR = 8;
                            obj.sigma_XPR = 4;
                            obj.N = 12;
                            obj.M = 20;
                            obj.zeta = 3;
                            obj.c_DS  = max(0.25,6.5622-3.4084*log10(fc_GHz));
                            obj.c_ASA = 11;
                            obj.c_ZSA = 7;
                            obj.c_ASD = 5;
                        else % 1 0
                            obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)-0.01*(obj.RX.height-1.5)+0.9);
                            obj.sigma_lgZSD = 0.49;
                            obj.mu_offset_ZOD = (7.66*log10(fc_GHz)-5.96) -10^((0.208*log10(fc_GHz)-0.782)*log10(max((25), obj.d_2D))+(-0.13*log10(fc_GHz)+2.03)-0.07*(obj.RX.height-1.5));
                            corr_rand_vec = obj.scenario.Cross_correlation_NLOS * randn(6,1);
                            % SF/DS/ASD/ASA/ZSD/ZSA
                            lsp_std_vec = [obj.sigma_SF, 0.39, 0.28, 0.11, obj.sigma_lgZSD, 0.16]';
                            lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                            obj.SF = lsp_offset_vec(1);   % dB
                            obj.DS = 10^(lsp_offset_vec(2)+(-6.28-0.204*log10(fc_GHz)));
                            obj.ASD = 10^(lsp_offset_vec(3)+(1.5-0.1144*log10(fc_GHz)));
                            obj.ASA = 10^(lsp_offset_vec(4)+(2.08-0.27*log10(fc_GHz)));
                            obj.ZSD = 10^(lsp_offset_vec(5)+obj.mu_lgZSD);
                            obj.ZSA = 10^(lsp_offset_vec(6)+(-0.3236*log10(fc_GHz)+1.512));
                            obj.r_tau = 2.3;
                            obj.mu_XPR = 7;
                            obj.sigma_XPR = 3;
                            obj.N = 20;
                            obj.M = 20;
                            obj.zeta = 3;
                            obj.c_DS  = max(0.25,6.5622-3.4084*log10(fc_GHz));
                            obj.c_ASA = 15;
                            obj.c_ZSA = 7;
                            obj.c_ASD = 2;
                        end
                    else
                        if av_alternative == 1
                            error('AV alternative 1 is not implemented yet.');
                        elseif av_alternative ~= 2
                            error('Unsupported AV alternative: %g.', av_alternative);
                        end
                        if obj.LOS
                            obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)-0.01*(obj.RX.height-1.5)+0.75);
                            obj.sigma_lgZSD = 0.4;
                            obj.mu_offset_ZOD = 0;
                            corr_rand_vec = obj.scenario.Cross_correlation_LOS * randn(7,1);
                            K_std = 8.158 * exp(0.0046 * obj.RX.height);
                            DS_std = 0.7294 * exp(0.0014 * obj.RX.height);
                            ASD_std = 1.0188 * exp(-0.0001 * obj.RX.height);
                            ASA_std = 1.0389 * exp(0.0085 * obj.RX.height);
                            ZSD_std = 1.0757 * exp(0.0059 * obj.RX.height);
                            ZSA_std = 0.9576 * exp(-0.0018 * obj.RX.height);
                            lsp_std_vec = [obj.sigma_SF, K_std, DS_std, ASD_std, ASA_std, ZSD_std, ZSA_std]';
                            lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                            obj.SF = lsp_offset_vec(1);
                            K_mu = 4.217*log10(obj.RX.height)+5.787;
                            obj.K = lsp_offset_vec(2)+K_mu;
                            DS_mu = -0.31*log10(obj.RX.height)-6.845;
                            obj.DS = 10^(lsp_offset_vec(3)+DS_mu);
                            ASD_mu = -0.0135*log10(obj.RX.height)+1.345;
                            obj.ASD = 10^(lsp_offset_vec(4)+ASD_mu);
                            ASA_mu = -2.4985*log10(obj.RX.height)-1.602;
                            obj.ASA = 10^(lsp_offset_vec(5)+ASA_mu);
                            ZSD_mu = -0.2975*log10(obj.RX.height)-0.5798;
                            obj.ZSD = 10^(lsp_offset_vec(6)+ZSD_mu);
                            ZSA_mu = -0.2895*log10(obj.RX.height)+0.225;
                            obj.ZSA = 10^(lsp_offset_vec(7)+ZSA_mu);
                            obj.r_tau = 2.5;
                            obj.mu_XPR = 8;
                            obj.sigma_XPR = 4;
                            obj.N = 12;
                            obj.M = 20;
                            obj.zeta = 3;
                            obj.c_DS  = max(0.25,6.5622-3.4084*log10(fc_GHz));
                            obj.c_ASA = 11;
                            obj.c_ZSA = 7;
                            obj.c_ASD = 3.58;
                        else
                            obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)-0.01*(obj.RX.height-1.5)+0.9);
                            obj.sigma_lgZSD = 0.49;
                            obj.mu_offset_ZOD = (7.66*log10(fc_GHz) - 5.96) ...
                                - 10^((0.208*log10(fc_GHz) - 0.782) * log10(max(25, obj.d_2D))) ...
                                + (-0.13*log10(fc_GHz) + 2.03) ...
                                - 0.07*(obj.RX.height - 1.5);
                            corr_rand_vec = obj.scenario.Cross_correlation_NLOS * randn(6,1);
                            DS_std = 0.9745 * exp(-0.0045 * obj.RX.height);
                            ASD_std = 1.2387 * exp(-0.0046 * obj.RX.height);
                            ASA_std = 1.022 * exp(0.009944 * obj.RX.height);
                            ZSD_std = 1.6421 * exp(-0.0092 * obj.RX.height);
                            ZSA_std = 1.6237 * exp(-0.0076 * obj.RX.height);
                            lsp_std_vec = [obj.sigma_SF, DS_std, ASD_std, ASA_std, ZSD_std, ZSA_std]';
                            lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                            obj.SF = lsp_offset_vec(1);
                            DS_mu = 0.0965*log10(obj.RX.height)-7.503;
                            obj.DS = 10^(lsp_offset_vec(2)+DS_mu);
                            ASD_mu = 1.17*log10(obj.RX.height)-0.665;
                            obj.ASD = 10^(lsp_offset_vec(3)+ASD_mu);
                            ASA_mu = -2.266*log10(obj.RX.height)-2.666;
                            obj.ASA = 10^(lsp_offset_vec(4)+ASA_mu);
                            ZSD_mu = 0.925*log10(obj.RX.height)-2.725;
                            obj.ZSD = 10^(lsp_offset_vec(5)+ZSD_mu);
                            ZSA_mu = -0.0005*log10(obj.RX.height)-0.4695;
                            obj.ZSA = 10^(lsp_offset_vec(6)+ZSA_mu);
                            obj.r_tau = 2.3;
                            obj.mu_XPR = 7;
                            obj.sigma_XPR = 3;
                            obj.N = 20;
                            obj.M = 20;
                            obj.zeta = 3;
                            obj.c_DS  = max(0.25,6.5622-3.4084*log10(fc_GHz));
                            obj.c_ASA = 15;
                            obj.c_ZSA = 7;
                            obj.c_ASD = 1.8;
                        end
                    end

                case 'UMi'
                    fc_GHz(fc_GHz<=2) = 2;
                    obj.drop_rp_angle = 50;
                    obj.delta_tau = 10^(rand*0.5-7.5);

                    if obj.O2I %
                        if obj.LOS
                            % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D_out/1000)+0.01*abs(obj.UE.height-obj.sector.attached_BS.height)+0.75);
                            obj.mu_lgZSD = max(-0.21, -14.8*(obj.d_2D_out/1000)+0.01*abs(obj.RX.height-obj.TX.height)+0.83);
                            obj.sigma_lgZSD = 0.35; % 0.4
                            obj.mu_offset_ZOD = 0;
                        else
                            % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D_out/1000)+0.01*max(obj.UE.height-obj.sector.attached_BS.height, 0)+0.9);
                            obj.mu_lgZSD = max(-0.5, -3.1*(obj.d_2D_out/1000)+0.01*max(obj.RX.height-obj.TX.height, 0)+0.2);
                            obj.sigma_lgZSD = 0.35; %0.6
                            % obj.mu_offset_ZOD = -10^(-0.55*log10(max(10, obj.d_2D_out))+1.6);
                            obj.mu_offset_ZOD = -10^(-1.5*log10(max(10, obj.d_2D_out))+3.3);
                        end
                        raw_rand_vec = randn(6,1);
                        % corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                        corr_rand_vec =  obj.scenario.Cross_correlation_O2I*raw_rand_vec;
                        lsp_std_vec = [obj.sigma_SF, 0.32, 0.42, 0.16, obj.sigma_lgZSD, 0.43]';
                        lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                        obj.SF = lsp_offset_vec(1);
                        obj.DS = 10^(lsp_offset_vec(2)-6.62);
                        obj.ASD = 10^(lsp_offset_vec(3)+1.25);
                        obj.ASA = 10^(lsp_offset_vec(4)+1.76);
                        obj.ZSD = 10^(lsp_offset_vec(5)+obj.mu_lgZSD);
                        obj.ZSA = 10^(lsp_offset_vec(6)+1.01);
                        obj.r_tau = 2.2;
                        obj.mu_XPR = 9;
                        obj.sigma_XPR = 5;  % 5
                        obj.N = 12;
                        obj.M = 20;
                        obj.zeta = 4;
                        obj.c_DS  = 11;
                        obj.c_ASA = 8;
                        obj.c_ZSA = 3;
                        obj.c_ASD = 5;
                    elseif obj.LOS % 0 1
                        % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)+0.01*abs(obj.UE.height-obj.sector.attached_BS.height)+0.75);
                        obj.mu_lgZSD = max(-0.21, -14.8*(obj.d_2D/1000)+0.01*abs(obj.RX.height-obj.TX.height)+0.83);
                        obj.sigma_lgZSD = 0.35; % 0.4
                        obj.mu_offset_ZOD = 0;
                        raw_rand_vec = randn(7,1);
                        % corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                        corr_rand_vec = obj.scenario.Cross_correlation_LOS*raw_rand_vec;
                        % SF/K/DS/ASD/ASA/ZSD/ZSA
                        % lsp_std_vec = [obj.sigma_SF, 5, 0.4, 0.43, 0.19, obj.sigma_lgZSD, 0.16]';
                        lsp_std_vec = [obj.sigma_SF, 5, 0.38, 0.41, (0.014*log10(1+fc_GHz)+0.28), obj.sigma_lgZSD, (-0.04*log10(1+fc_GHz)+0.34)]';
                        lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                        obj.SF = lsp_offset_vec(1);
                        obj.K = lsp_offset_vec(2)+9;
                        obj.DS = 10^(lsp_offset_vec(3)+(-0.24*log10(1+fc_GHz)-7.14));
                        obj.ASD = 10^(lsp_offset_vec(4)+(-0.05*log10(1+fc_GHz)+1.21));
                        obj.ASA = 10^(lsp_offset_vec(5)+(-0.08*log10(1+fc_GHz)+1.73));
                        obj.ZSD = 10^(lsp_offset_vec(6)+obj.mu_lgZSD);
                        obj.ZSA = 10^(lsp_offset_vec(7)+(-0.1*log10(1+fc_GHz)+0.73));
                        obj.r_tau = 3;
                        obj.mu_XPR = 9;
                        obj.sigma_XPR = 3;
                        obj.N = 12;
                        obj.M = 20;
                        obj.zeta = 3;
                        obj.c_DS  = 5;
                        obj.c_ASA = 17;
                        obj.c_ZSA = 7;
                        obj.c_ASD = 3;
                    else % NLOS
                        % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)+0.01*max(obj.UE.height-obj.sector.attached_BS.height, 0)+0.9);
                        obj.mu_lgZSD = max(-0.5, -3.1*(obj.d_2D/1000)+0.01*max(obj.RX.height-obj.TX.height, 0)+0.2);
                        obj.sigma_lgZSD = 0.35;
                        % obj.mu_offset_ZOD = -10^(-0.55*log10(max(10, obj.d_2D))+1.6);
                        obj.mu_offset_ZOD = -10^(-1.5*log10(max(10, obj.d_2D))+3.3);
                        raw_rand_vec = randn(6,1);
                        % corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                        corr_rand_vec =  obj.scenario.Cross_correlation_NLOS*raw_rand_vec;
                        % SF/DS/ASD/ASA/ZSD/ZSA
                        lsp_std_vec = [obj.sigma_SF, (0.16*log10(1+fc_GHz)+0.28), (0.11*log10(1+fc_GHz)+0.33), (0.05*log10(1+fc_GHz)+0.3), obj.sigma_lgZSD, (-0.07*log10(1+fc_GHz)+0.41)]';
                        lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                        obj.SF = lsp_offset_vec(1);
                        obj.DS = 10^(lsp_offset_vec(2)+(-0.24*log10(1+fc_GHz)-6.83));
                        obj.ASD = 10^(lsp_offset_vec(3)+(-0.23*log10(1+fc_GHz)+1.53));
                        obj.ASA = 10^(lsp_offset_vec(4)+(-0.08*log10(1+fc_GHz)+1.81));
                        obj.ZSD = 10^(lsp_offset_vec(5)+obj.mu_lgZSD);
                        obj.ZSA = 10^(lsp_offset_vec(6)+(-0.04*log10(1+fc_GHz)+0.92));
                        obj.r_tau = 2.1;
                        obj.mu_XPR = 8;
                        obj.sigma_XPR = 3;
                        obj.N = 19;
                        obj.M = 20;
                        obj.zeta = 3;
                        obj.c_DS  = 11;
                        obj.c_ASA = 22;
                        obj.c_ZSA = 7;
                        obj.c_ASD = 10;
                    end
                case 'RMa'
                    if obj.O2I
                        obj.mu_lgZSD = max(-1, -0.19*(obj.d_2D_out/1000)-0.01*(obj.RX.height-1.5)+0.28);
                        obj.sigma_lgZSD = 0.30;
                        obj.mu_offset_ZOD = atand((35-3.5)/obj.d_2D_out)-atand((35-1.5)/obj.d_2D_out);
                        raw_rand_vec = randn(6,1);
                        corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                        corr_rand_vec =  obj.scenario.Cross_correlation_O2I*corr_rand_vec;
                        obj.SF = obj.sigma_SF*corr_rand_vec(1);
                        obj.DS = 10^(0.24*corr_rand_vec(2)-7.47);
                        obj.ASD = 10^(0.18*corr_rand_vec(3)+0.67);
                        obj.ASA = 10^(0.21*corr_rand_vec(4)+1.66);
                        obj.ZSD = 10^(obj.sigma_lgZSD*corr_rand_vec(5)+obj.mu_lgZSD);
                        obj.ZSA = 10^(0.22*corr_rand_vec(6)+0.93);
                        obj.r_tau = 1.7;
                        obj.mu_XPR = 7;
                        obj.sigma_XPR = 3;
                        obj.N = 10;
                        obj.M = 20;
                        obj.zeta = 3;
                        obj.c_ASA = 3;
                        obj.c_ZSA = 3;
                        obj.c_ASD = 2;
                    elseif obj.LOS
                        obj.mu_lgZSD = max(-1, -0.17*(obj.d_2D/1000)-0.01*(obj.RX.height-1.5)+0.22);
                        obj.sigma_lgZSD = 0.34;
                        obj.mu_offset_ZOD = 0;
                        raw_rand_vec = randn(7,1);
                        corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                        corr_rand_vec =  obj.scenario.Cross_correlation_LOS*corr_rand_vec;
                        obj.SF = obj.sigma_SF*corr_rand_vec(1);
                        obj.K = 4*corr_rand_vec(2)+7;
                        obj.DS = 10^(0.55*corr_rand_vec(3)-7.49);
                        obj.ASD = 10^(0.38*corr_rand_vec(4)+0.90);
                        obj.ASA = 10^(0.24*corr_rand_vec(5)+1.52);
                        obj.ZSD = 10^(obj.sigma_lgZSD*corr_rand_vec(6)+obj.mu_lgZSD);
                        obj.ZSA = 10^(0.40*corr_rand_vec(7)+0.47);
                        obj.r_tau = 3.8;
                        obj.mu_XPR = 12;
                        obj.sigma_XPR = 4;
                        obj.N = 11;
                        obj.M = 20;
                        obj.zeta = 3;
                        obj.c_ASA = 3;
                        obj.c_ZSA = 2;
                        obj.c_ASD = 2;
                    else
                        obj.mu_lgZSD = max(-1, -0.19*(obj.d_2D/1000)-0.01*(obj.RX.height-1.5)+0.28);
                        obj.sigma_lgZSD = 0.30;
                        obj.mu_offset_ZOD = atand((35-3.5)/obj.d_2D)-atand((35-1.5)/obj.d_2D);
                        raw_rand_vec = randn(6,1);
                        corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                        corr_rand_vec =  obj.scenario.Cross_correlation_NLOS*corr_rand_vec;
                        obj.SF = obj.sigma_SF*corr_rand_vec(1);
                        obj.DS = 10^(0.48*corr_rand_vec(2)-7.43);
                        obj.ASD = 10^(0.45*corr_rand_vec(3)+0.95);
                        obj.ASA = 10^(0.13*corr_rand_vec(4)+1.52);
                        obj.ZSD = 10^(obj.sigma_lgZSD*corr_rand_vec(5)+obj.mu_lgZSD);
                        obj.ZSA = 10^(0.37*corr_rand_vec(6)+0.58);
                        obj.r_tau = 1.7;
                        obj.mu_XPR = 7;
                        obj.sigma_XPR = 3;
                        obj.N = 10;
                        obj.M = 20;
                        obj.zeta = 3;
                        obj.c_ASA = 3;
                        obj.c_ZSA = 2.5;
                        obj.c_ASD = 2;
                    end
                case 'InH'
                    fc_GHz(fc_GHz<=6) = 6;
                    obj.delta_tau = 10^(rand*0.1-8.6);
                    if obj.LOS
                        % obj.mu_lgZSD = 1.02;
                        % obj.sigma_lgZSD = 0.41;
                        obj.mu_lgZSD = -1.43*log10(1+fc_GHz)+2.228;
                        obj.sigma_lgZSD = 0.13*log10(1+fc_GHz)+0.3;
                        obj.mu_offset_ZOD = 0;
                        raw_rand_vec = randn(7,1);
                        corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                        corr_rand_vec =  obj.scenario.Cross_correlation_LOS*corr_rand_vec;
                        obj.SF = obj.sigma_SF*corr_rand_vec(1);
                        obj.K = 4*corr_rand_vec(2)+7;
                        obj.DS = 10^(0.18*corr_rand_vec(3)+ (-0.01*log10(1+fc_GHz)-7.692));
                        obj.ASD = 10^(0.18*corr_rand_vec(4)+1.60);
                        obj.ASA = 10^((0.12*log10(1+fc_GHz)+0.119)*corr_rand_vec(5)+ (-0.19*log10(1+fc_GHz)+1.781));
                        obj.ZSD = 10^(obj.sigma_lgZSD*corr_rand_vec(6)+obj.mu_lgZSD);
                        obj.ZSA = 10^((-0.04*log10(1+fc_GHz)+0.264)*corr_rand_vec(7)+ (-0.26*log10(1+fc_GHz)+1.44));
                        obj.r_tau = 3.6;
                        obj.mu_XPR = 11;
                        obj.sigma_XPR = 4; % 3
                        obj.N = 15;
                        obj.M = 20;
                        obj.zeta = 6;
                        obj.c_ASA = 8;
                        obj.c_ZSA = 9;
                        obj.c_ASD = 5;
                    else
                        obj.mu_lgZSD = 1.08;
                        obj.sigma_lgZSD = 0.36;
                        obj.mu_offset_ZOD = 0;
                        raw_rand_vec = randn(6,1);
                        corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                        corr_rand_vec =  obj.scenario.Cross_correlation_NLOS*corr_rand_vec;
                        obj.SF = obj.sigma_SF*corr_rand_vec(1);
                        obj.DS = 10^((0.1*log10(1+fc_GHz)+0.055)*corr_rand_vec(2)+ (-0.28*log10(1+fc_GHz)-7.173));
                        obj.ASD = 10^(0.25*corr_rand_vec(3)+1.62);
                        obj.ASA = 10^((0.12*log10(1+fc_GHz)+0.059)*corr_rand_vec(4)+ (-0.11*log10(1+fc_GHz)+1.863));
                        obj.ZSD = 10^(obj.sigma_lgZSD*corr_rand_vec(5)+obj.mu_lgZSD);
                        obj.ZSA = 10^((-0.09*log10(1+fc_GHz)+0.746)*corr_rand_vec(6)+ (-0.15*log10(1+fc_GHz)+1.387));
                        obj.r_tau = 3;
                        obj.mu_XPR = 10;
                        obj.sigma_XPR = 4; % 3
                        obj.N = 19;
                        obj.M = 20;
                        obj.zeta = 3;
                        obj.c_ASA = 11;
                        obj.c_ZSA = 9;
                        obj.c_ASD = 5;
                    end
                case 'InF'
                    obj.delta_tau = 10^(rand*0.4-7.5);
                    % V = hall volume in m^3, S = total surface area of hall in m^2 (walls+floor+ceiling)
                    if obj.LOS
                        obj.mu_lgZSD  = 1.35;
                        obj.sigma_lgZSD  = 0.35;
                        obj.mu_offset_ZOD  = 0;
                        raw_rand_vec = randn(7,1);
                        corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                        corr_rand_vec =  obj.scenario.Cross_correlation_LOS*corr_rand_vec;
                        obj.SF  = obj.sigma_SF *corr_rand_vec(1);
                        obj.K  = 8*corr_rand_vec(2)+7;
                        obj.DS  = 10^(0.15*corr_rand_vec(3)+ (log10(26*(obj.scenario.V/obj.scenario.S)+14)-9.35));
                        obj.ASD  = 10^(0.25*corr_rand_vec(4)+1.56);
                        obj.ASA  = 10^((0.12*log10(1+fc_GHz)+0.2)*corr_rand_vec(5)+ (-0.18*log10(1+fc_GHz)+1.78));
                        obj.ZSD  = 10^(obj.sigma_lgZSD*corr_rand_vec(6)+obj.mu_lgZSD);
                        obj.ZSA  = 10^(0.35*corr_rand_vec(7)+ (-0.2*log10(1+fc_GHz)+1.5));
                        obj.r_tau  = 2.7;
                        obj.mu_XPR  = 12;
                        obj.sigma_XPR  = 6;
                        obj.N  = 25;
                        obj.M  = 20;
                        obj.zeta  = 4;
                        obj.c_ASA  = 8;
                        obj.c_ZSA  = 9;
                        obj.c_ASD  = 5;
                    else
                        obj.mu_lgZSD  = 1.2;
                        obj.sigma_lgZSD  = 0.55;
                        obj.mu_offset_ZOD  = 0;
                        raw_rand_vec = randn(6,1);
                        corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                        corr_rand_vec =  obj.scenario.Cross_correlation_NLOS*corr_rand_vec;
                        obj.SF  = obj.sigma_SF *corr_rand_vec(1);
                        obj.DS  = 10^(0.19*corr_rand_vec(2)+ (log10(30*(obj.scenario.V/obj.scenario.S)+32)-9.44));
                        obj.ASD  = 10^(0.2*corr_rand_vec(3)+1.57);
                        obj.ASA  = 10^(0.3*corr_rand_vec(4)+1.72);
                        obj.ZSD  = 10^(obj.sigma_lgZSD *corr_rand_vec(5)+obj.mu_lgZSD );
                        obj.ZSA  = 10^(0.45*corr_rand_vec(6)+ (-0.13*log10(1+fc_GHz)+1.45));
                        obj.r_tau  = 3;
                        obj.mu_XPR  = 11;
                        obj.sigma_XPR  = 6;
                        obj.N  = 25;
                        obj.M  = 20;
                        obj.zeta  = 3;
                        obj.c_ASA  = 8;
                        obj.c_ZSA  = 9;
                        obj.c_ASD  = 5;
                    end
            end
            obj.ASA = min(obj.ASA, 104);
            obj.ASD = min(obj.ASD, 104);
            obj.ZSA = min(obj.ZSA, 52);
            obj.ZSD = min(obj.ZSD, 52);
        end

        % step 5: Generate cluster delays
        function cluster_delay(obj)
            random_draw = rand(obj.N,1);
            % end
            tau_n_tmp = -obj.r_tau*obj.DS*log(random_draw.');
            [obj.tau_n, obj.tau_order] = sort(tau_n_tmp - min(tau_n_tmp));
            obj.tau_n_LOS = obj.tau_n;
            if obj.LOS && ~obj.O2I
                C_tau = 0.7705-0.0433*obj.K+0.0002*obj.K^2+0.000017*obj.K^3;
                obj.tau_n_LOS = obj.tau_n/C_tau;  % not to be used in cluster power generation
            end

        end

        % step 6: Generate cluster powers
        function cluster_power(obj)
                shadow_dB = randn(obj.N,1)*obj.zeta;
                obj.cluster_shadow_dB = shadow_dB;
            % end
            cluster_power_raw = exp(-obj.tau_n.*(obj.r_tau-1)./obj.r_tau./obj.DS).*(10.^(-shadow_dB.'/10));
            obj.Pn = cluster_power_raw/sum(cluster_power_raw);
            obj.Pn_LOS = obj.Pn;
            if obj.LOS && ~obj.O2I
                rice_linear = 10^(obj.K/10);
                P1_LOS = rice_linear/(rice_linear+1);
                obj.Pn_LOS = (1/(rice_linear+1)).*(cluster_power_raw./sum(cluster_power_raw));
                obj.Pn_LOS(1) = obj.Pn_LOS(1) + P1_LOS;  % not be used in equation step 11(7.3-22)
            end
                keep_tmp       = ((10*log10(obj.Pn_LOS/max(obj.Pn_LOS))) >= -25);
            % end
            % keep_tmp       = true(1,numel(obj.Pn_LOS));
            if sum(keep_tmp) == 1
                keep_tmp(obj.Pn_LOS == max(obj.Pn_LOS(~keep_tmp))) = true;
            end
            obj.keep  = keep_tmp;
            if ~isempty(obj.tau_absolute)
                obj.tau_absolute = obj.tau_absolute(keep_tmp);
            end
            obj.tau_n  = obj.tau_n(keep_tmp);
            obj.tau_n_LOS  = obj.tau_n_LOS(keep_tmp);
            obj.cluster_shadow_dB = obj.cluster_shadow_dB(keep_tmp);
            obj.Pn     = obj.Pn(keep_tmp);
            obj.Pn_LOS = obj.Pn_LOS(keep_tmp);
            obj.N_new  = numel(obj.Pn);
            [~,sort_idx]    = sort((1./obj.Pn));
            obj.strong_cluster_id = sort_idx(1:min(2,numel(sort_idx)));  % cluster number is 2
            obj.map_delay = zeros(1,obj.M);
            if ~isempty(obj.c_DS)
                obj.map_delay([9,10,11,12,17,18]) = obj.c_DS*1.28e-9;
                obj.map_delay(13:16) = obj.c_DS*2.56e-9;
            else
                obj.map_delay(9:18) = 3.91e-9;
            end
        end

        % caculator the AOA for every ray
        function cluster_delay_procedureB(obj)
            limits = obj.procedureBClusterLimits();
            obj.tau_prime = limits.delay*obj.drawProcedureBUniform('tau', obj.N);
            [obj.tau_n, obj.tau_order] = sort(obj.tau_prime - min(obj.tau_prime));
            obj.tau_prime = obj.tau_n;
            if obj.LOS && ~obj.O2I
                obj.tau_prime(1) = 0;
                obj.tau_n(1) = 0;
            end
            obj.tau_n_LOS = obj.tau_n;
        end

        function cluster_angles_procedureB(obj)
            limits = obj.procedureBClusterLimits();
            obj.phi_prime_AOA = limits.AOA*(2*obj.drawProcedureBUniform('AOA', obj.N)-1);
            obj.phi_prime_AOD = limits.AOD*(2*obj.drawProcedureBUniform('AOD', obj.N)-1);
            obj.theta_prime_ZOA = limits.ZOA*(2*obj.drawProcedureBUniform('ZOA', obj.N)-1);
            obj.theta_prime_ZOD = limits.ZOD*(2*obj.drawProcedureBUniform('ZOD', obj.N)-1);

            if ~isempty(obj.tau_order)
                obj.phi_prime_AOA = obj.phi_prime_AOA(obj.tau_order);
                obj.phi_prime_AOD = obj.phi_prime_AOD(obj.tau_order);
                obj.theta_prime_ZOA = obj.theta_prime_ZOA(obj.tau_order);
                obj.theta_prime_ZOD = obj.theta_prime_ZOD(obj.tau_order);
            end

            if obj.LOS && ~obj.O2I
                obj.phi_prime_AOA(1) = 0;
                obj.phi_prime_AOD(1) = 0;
                obj.theta_prime_ZOA(1) = 0;
                obj.theta_prime_ZOD(1) = 0;
            end
        end

        function cluster_delay_procedureA(obj)
            state = obj.getProcedureAState();
            if isempty(state) || ~isfield(state, 'tau_absolute') || isempty(state.tau_absolute)
                obj.cluster_delay_procedureA_initial();
                return;
            end

            obj.N = numel(state.tau_absolute);
            delta_t = max(obj.t - state.time, 0);
            tx_delta = obj.velocityVector(obj.TX) * delta_t;
            rx_delta = obj.velocityVector(obj.RX) * delta_t;
            tau_update = zeros(1, obj.N);
            for cluster_idx = 1:obj.N
                aoa_unit = obj.unitVectorFromAngles( ...
                    obj.phi_LOS_AOA + state.phi_prime_AOA(cluster_idx), ...
                    obj.theta_LOS_ZOA + state.theta_prime_ZOA(cluster_idx));
                aod_unit = obj.unitVectorFromAngles( ...
                    obj.phi_LOS_AOD + state.phi_prime_AOD(cluster_idx), ...
                    obj.theta_LOS_ZOD + state.theta_prime_ZOD(cluster_idx));
                tau_update(cluster_idx) = (-dot(rx_delta, aoa_unit) + dot(tx_delta, aod_unit)) / 3e8;
            end

            tau_abs = max(state.tau_absolute + tau_update, 0);
            tau_tmp = tau_abs - min(tau_abs);
            if obj.LOS && ~obj.O2I
                tau_tmp(1) = 0;
            end
            obj.tau_absolute = tau_abs;
            obj.tau_n = tau_tmp;
            obj.tau_prime = tau_tmp;
            obj.tau_order = state.tau_order;
            obj.tau_n_LOS = tau_tmp;
        end

        function cluster_delay_procedureA_initial(obj)
            tau_raw = -obj.r_tau*obj.DS*log(rand(1, obj.N));
            obj.tau_absolute = obj.d_3D/3e8 + tau_raw;
            [obj.tau_n, obj.tau_order] = sort(obj.tau_absolute - min(obj.tau_absolute));
            obj.tau_absolute = obj.tau_absolute(obj.tau_order);
            obj.tau_prime = obj.tau_n;
            obj.tau_n_LOS = obj.tau_n;
            if obj.LOS && ~obj.O2I
                C_tau = 0.7705-0.0433*obj.K+0.0002*obj.K^2+0.000017*obj.K^3;
                obj.tau_n_LOS = obj.tau_n/C_tau;
                obj.tau_n(1) = 0;
                obj.tau_prime(1) = 0;
            end
        end

        function cluster_angles_procedureA(obj)
            state = obj.getProcedureAState();
            if isempty(state) || ~isfield(state, 'phi_prime_AOA') || ...
                    numel(state.phi_prime_AOA) ~= obj.N
                obj.cluster_angles_procedureB();
                return;
            end

            delta_t = max(obj.t - state.time, 0);
            tx_velocity = obj.velocityVector(obj.TX);
            rx_velocity = obj.velocityVector(obj.RX);

            prev_phi_AOA = state.phi_LOS_AOA + state.phi_prime_AOA;
            prev_phi_AOD = state.phi_LOS_AOD + state.phi_prime_AOD;
            prev_theta_ZOA = state.theta_LOS_ZOA + state.theta_prime_ZOA;
            prev_theta_ZOD = state.theta_LOS_ZOD + state.theta_prime_ZOD;
            x_n = obj.procedureAXnFromState(state, numel(prev_phi_AOA));
            [rx_velocity_eff, tx_velocity_eff] = obj.procedureAEffectiveVelocities( ...
                prev_phi_AOA, prev_theta_ZOA, prev_phi_AOD, prev_theta_ZOD, ...
                tx_velocity, rx_velocity, x_n, obj.LOS && ~obj.O2I);

            [phi_AOA, theta_ZOA] = obj.updateProcedureAAnglePair( ...
                prev_phi_AOA, prev_theta_ZOA, rx_velocity_eff, obj.tau_n, delta_t);
            [phi_AOD, theta_ZOD] = obj.updateProcedureAAnglePair( ...
                prev_phi_AOD, prev_theta_ZOD, tx_velocity_eff, obj.tau_n, delta_t);

            obj.phi_prime_AOA = obj.wrapAzimuthAngles(phi_AOA - obj.phi_LOS_AOA);
            obj.phi_prime_AOD = obj.wrapAzimuthAngles(phi_AOD - obj.phi_LOS_AOD);
            obj.theta_prime_ZOA = theta_ZOA - obj.theta_LOS_ZOA;
            obj.theta_prime_ZOD = theta_ZOD - obj.theta_LOS_ZOD;

            if obj.LOS && ~obj.O2I
                obj.phi_prime_AOA(1) = 0;
                obj.phi_prime_AOD(1) = 0;
                obj.theta_prime_ZOA(1) = 0;
                obj.theta_prime_ZOD(1) = 0;
            end
        end

        function cluster_power_procedureB(obj)
            shadow_dB = obj.drawProcedureBGaussian('shadow', obj.N).' * obj.zeta;
            obj.cluster_shadow_dB = shadow_dB.';
            ds_power = obj.DS;
            asa_power = obj.ASA;
            asd_power = obj.ASD;
            zsa_power = obj.ZSA;
            zsd_power = obj.ZSD;
            rice_linear = 0;

            if obj.LOS && ~obj.O2I
                rice_linear = 10^(obj.K/10);
                ds_power = ds_power*sqrt((1 + rice_linear)/2);
                asa_power = asa_power*sqrt(1 + rice_linear);
                asd_power = asd_power*sqrt(1 + rice_linear);
                zsa_power = zsa_power*sqrt(1 + rice_linear);
                zsd_power = zsd_power*sqrt(1 + rice_linear);
            end

            cluster_power_raw = exp(-obj.tau_prime./ds_power) ...
                .* exp(-sqrt(2)*abs(obj.phi_prime_AOA)./asa_power) ...
                .* exp(-sqrt(2)*abs(obj.phi_prime_AOD)./asd_power) ...
                .* exp(-sqrt(2)*abs(obj.theta_prime_ZOA)./zsa_power) ...
                .* exp(-sqrt(2)*abs(obj.theta_prime_ZOD)./zsd_power) ...
                .* 10.^(-shadow_dB.'/10);

            if obj.LOS && ~obj.O2I
                obj.Pn = (cluster_power_raw./sum(cluster_power_raw))/(1 + rice_linear);
                obj.Pn(1) = obj.Pn(1) + rice_linear/(1 + rice_linear);
            else
                obj.Pn = cluster_power_raw./sum(cluster_power_raw);
            end

            obj.Pn_LOS = obj.Pn;
            keep_tmp = true(size(obj.Pn));
            obj.keep = keep_tmp;
            obj.tau_n = obj.tau_n(keep_tmp);
            obj.tau_n_LOS = obj.tau_n_LOS(keep_tmp);
            obj.tau_prime = obj.tau_prime(keep_tmp);
            obj.phi_prime_AOA = obj.phi_prime_AOA(keep_tmp);
            obj.phi_prime_AOD = obj.phi_prime_AOD(keep_tmp);
            obj.theta_prime_ZOA = obj.theta_prime_ZOA(keep_tmp);
            obj.theta_prime_ZOD = obj.theta_prime_ZOD(keep_tmp);
            obj.Pn = obj.Pn(keep_tmp);
            obj.Pn_LOS = obj.Pn_LOS(keep_tmp);
            obj.N_new = numel(obj.Pn);
            [~, sort_idx] = sort((1./obj.Pn));
            obj.strong_cluster_id = sort_idx(1:min(2, numel(sort_idx)));
            obj.map_delay = zeros(1, obj.M);
            if ~isempty(obj.c_DS)
                obj.map_delay([9,10,11,12,17,18]) = obj.c_DS*1.28e-9;
                obj.map_delay(13:16) = obj.c_DS*2.56e-9;
            else
                obj.map_delay(9:18) = 3.91e-9;
            end
        end

        function cluster_power_procedureA(obj)
            state = obj.getProcedureAState();
            if ~isempty(state) && isfield(state, 'cluster_shadow_dB') && ...
                    numel(state.cluster_shadow_dB) == obj.N
                obj.cluster_shadow_dB = state.cluster_shadow_dB(:);
            else
                obj.cluster_shadow_dB = randn(obj.N, 1) * obj.zeta;
            end

            ds_power = obj.DS;
            rice_linear = 0;
            if obj.LOS && ~obj.O2I
                rice_linear = 10^(obj.K/10);
            end

            cluster_power_raw = exp(-obj.tau_n.*(obj.r_tau - 1)./obj.r_tau./ds_power) ...
                .* 10.^(-obj.cluster_shadow_dB.'/10);

            if obj.LOS && ~obj.O2I
                obj.Pn = (cluster_power_raw./sum(cluster_power_raw))/(1 + rice_linear);
                obj.Pn(1) = obj.Pn(1) + rice_linear/(1 + rice_linear);
            else
                obj.Pn = cluster_power_raw./sum(cluster_power_raw);
            end

            obj.Pn_LOS = obj.Pn;
            obj.keep = true(size(obj.Pn));
            obj.tau_n = obj.tau_n(obj.keep);
            obj.tau_n_LOS = obj.tau_n_LOS(obj.keep);
            obj.tau_prime = obj.tau_prime(obj.keep);
            obj.phi_prime_AOA = obj.phi_prime_AOA(obj.keep);
            obj.phi_prime_AOD = obj.phi_prime_AOD(obj.keep);
            obj.theta_prime_ZOA = obj.theta_prime_ZOA(obj.keep);
            obj.theta_prime_ZOD = obj.theta_prime_ZOD(obj.keep);
            obj.Pn = obj.Pn(obj.keep);
            obj.Pn_LOS = obj.Pn_LOS(obj.keep);
            obj.cluster_shadow_dB = obj.cluster_shadow_dB(obj.keep);
            obj.N_new = numel(obj.Pn);
            [~, sort_idx] = sort((1./obj.Pn));
            obj.strong_cluster_id = sort_idx(1:min(2, numel(sort_idx)));
            obj.map_delay = zeros(1, obj.M);
            if ~isempty(obj.c_DS)
                obj.map_delay([9,10,11,12,17,18]) = obj.c_DS*1.28e-9;
                obj.map_delay(13:16) = obj.c_DS*2.56e-9;
            else
                obj.map_delay(9:18) = 3.91e-9;
            end
        end

        function apply_angle_offsets_procedureB(obj)
            ray_angle_offset = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];

            phi_n_AOA = obj.phi_LOS_AOA + obj.phi_prime_AOA;
            phi_n_AOD = obj.phi_LOS_AOD + obj.phi_prime_AOD;
            ZOA = obj.theta_LOS_ZOA;
            ZOA(obj.O2I) = 90;
            theta_n_ZOA = ZOA + obj.theta_prime_ZOA;
            theta_n_ZOD = obj.theta_LOS_ZOD + obj.mu_offset_ZOD + obj.theta_prime_ZOD;

            obj.phi_n_m_AOA = repmat(phi_n_AOA.', 1, obj.M) + obj.c_ASA(ones(obj.N_new,1))*ray_angle_offset;
            obj.phi_n_m_AOD = repmat(phi_n_AOD.', 1, obj.M) + obj.c_ASD(ones(obj.N_new,1))*ray_angle_offset;
            obj.theta_n_m_ZOA = repmat(theta_n_ZOA.', 1, obj.M) + obj.c_ZSA(ones(obj.N_new,1))*ray_angle_offset;

            zsd_offset_scale = (3/8)*(10^obj.mu_lgZSD);
            obj.theta_n_m_ZOD = repmat(theta_n_ZOD.', 1, obj.M) + zsd_offset_scale(ones(obj.N_new,1))*ray_angle_offset;

            obj.phi_n_m_AOA = obj.wrapAzimuthAngles(obj.phi_n_m_AOA);
            obj.phi_n_m_AOD = obj.wrapAzimuthAngles(obj.phi_n_m_AOD);
            obj.theta_n_m_ZOA = obj.wrapZenithAngles(obj.theta_n_m_ZOA);
            obj.theta_n_m_ZOD = obj.wrapZenithAngles(obj.theta_n_m_ZOD);
        end

        % caculator the AOA for every ray

        function AOA_calc(obj)
            %             if strcmp(obj.scenario.name,'InH')
            %                 nlos_az_scale_table = [1.434, 1,1,1, 1.501];
            %                 nlos_az_scale = nlos_az_scale_table(obj.N - 14);
            %                 if obj.bLOS && ~obj.O2I
            %                     az_scale = nlos_az_scale*(0.9275-0.0439*obj.K-0.0071*obj.K^2+0.0002*obj.K^3);
            %                 else
            %                     az_scale = nlos_az_scale;
            %                 end
            %                 phi_n_AOA_tmp = -obj.ASA*log(obj.Pn_LOS/max(obj.Pn_LOS))/az_scale;
            %             else  % for UMa, RMa, UMi
            %                 %                                    4       5                8           10     11     12           14     15     16                    19     20
            %                 nlos_az_scale_table = [0.5, 0.58, 0.69, 0.779, 0.860,0.92,0.975, 1.018,1.054, 1.09, 1.123, 1.146, 1.168, 1.19, 1.211, 1.226, 1.242, 1.257, 1.273,  1.289];
            % %                 nlos_az_scale_table = [0.779, 0.779, 0.779, 0.779, 0.860,0.92,0.975, 1.018,1.054, 1.09, 1.123, 1.146, 1.168, 1.19, 1.211, 1.226, 1.242, 1.257, 1.273,  1.289];
            %                 nlos_az_scale = nlos_az_scale_table(obj.N);
            %                 if obj.bLOS && ~obj.O2I
            %                     az_scale = nlos_az_scale*(1.1035-0.028*obj.K-0.002*obj.K^2+0.0001*obj.K^3);
            %                 else
            %                     az_scale = nlos_az_scale;
            %                 end
            %                 phi_n_AOA_tmp = 2*(obj.ASA/1.4)*sqrt(-log(obj.Pn_LOS/max(obj.Pn_LOS)))/az_scale;
            %             end
            n_clusters     = [4 5 8 10 11 12 14 15 16 19 20 25];
            nlos_az_scale_table = [0.779, 0.860, 1.018, 1.090, 1.123, 1.146, 1.190, 1.211, 1.226, 1.273, 1.289, 1.289];% , 1.358
            nlos_az_scale     = nlos_az_scale_table(n_clusters == obj.N);
            az_scale          = sum([nlos_az_scale ,nlos_az_scale*(0.1035-0.028*obj.K-0.002*obj.K^2+0.0001*obj.K^3)]); % for NLOS and O2I, K = [].
            phi_n_AOA_tmp  = 2*(obj.ASA/1.4)*sqrt(-log(obj.Pn_LOS/max(obj.Pn_LOS)))/az_scale;

                random_draw = randsrc(1,obj.N_new,[-1, 1]);
                angle_jitter = (obj.ASA/7)*randn(1,obj.N_new);
            % end
            phi_n_AOA = random_draw.*phi_n_AOA_tmp+angle_jitter+obj.phi_LOS_AOA;
            if obj.LOS % && ~obj.O2I
                phi_n_AOA = (random_draw.*phi_n_AOA_tmp+angle_jitter)-(random_draw(1)*phi_n_AOA_tmp(1)+angle_jitter(1)-obj.phi_LOS_AOA);
            end
            ray_angle_offset = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];
            phi_n_m_AOA_ = repmat(phi_n_AOA.',1,obj.M)+obj.c_ASA(ones(obj.N_new,1))*ray_angle_offset;
            phi_n_m_AOA_ = mod(phi_n_m_AOA_,360);
            phi_n_m_AOA_ = phi_n_m_AOA_-360*floor(phi_n_m_AOA_/180);
            obj.phi_n_m_AOA = phi_n_m_AOA_;
        end

        % caculator the AOD for every ray
        function AOD_calc(obj)
            %             if strcmp(obj.scenario.name,'InH')
            %                 nlos_az_scale_table = [1.434, 1,1,1, 1.501];
            %                 nlos_az_scale = nlos_az_scale_table(obj.N - 14);
            %                 if obj.bLOS && ~obj.O2I
            %                     az_scale = nlos_az_scale*(0.9275-0.0439*obj.K-0.0071*obj.K^2+0.0002*obj.K^3);
            %                 else
            %                     az_scale = nlos_az_scale;
            %                 end
            %                 phi_n_AOD_tmp = -obj.ASD*log(obj.Pn_LOS/max(obj.Pn_LOS))/az_scale;
            %             else  % for UMa, RMa, UMi
            %                 %                                    4       5                8           10     11     12           14     15     16                    19     20
            %                 nlos_az_scale_table = [0.5, 0.58, 0.69, 0.779, 0.860,0.92,0.975, 1.018,1.054, 1.09, 1.123, 1.146, 1.168, 1.19, 1.211, 1.226, 1.242, 1.257, 1.273,  1.289];
            % %                 nlos_az_scale_table = [0.779, 0.779, 0.779, 0.779, 0.860,0.92,0.975, 1.018,1.054, 1.09, 1.123, 1.146, 1.168, 1.19, 1.211, 1.226, 1.242, 1.257, 1.273,  1.289];
            %                 nlos_az_scale = nlos_az_scale_table(obj.N);
            %                 if obj.bLOS && ~obj.O2I
            %                     az_scale = nlos_az_scale*(1.1035-0.028*obj.K-0.002*obj.K^2+0.0001*obj.K^3);
            %                 else
            %                     az_scale = nlos_az_scale;
            %                 end
            %                 phi_n_AOD_tmp = 2*(obj.ASD/1.4)*sqrt(-log(obj.Pn_LOS/max(obj.Pn_LOS)))/az_scale;
            %             end
            n_clusters     = [4 5 8 10 11 12 14 15 16 19 20 25];
            nlos_az_scale_table = [0.779, 0.860, 1.018, 1.090, 1.123, 1.146, 1.190, 1.211, 1.226, 1.273, 1.289, 1.289];% , 1.358
            nlos_az_scale     = nlos_az_scale_table(n_clusters == obj.N);
            az_scale          = sum([nlos_az_scale ,nlos_az_scale*(0.1035-0.028*obj.K-0.002*obj.K^2+0.0001*obj.K^3)]); % for NLOS and O2I, K = [].
            phi_n_AOD_tmp  = 2*(obj.ASD/1.4)*sqrt(-log(obj.Pn_LOS/max(obj.Pn_LOS)))/az_scale;
            %             phi_n_AOD_tmp  = -obj.ASD*log(obj.Pn_LOS/max(obj.Pn_LOS))/az_scale;

                random_draw = randsrc(1,obj.N_new,[-1, 1]);
                angle_jitter = (obj.ASD/7)*randn(1,obj.N_new);
            % end
            phi_n_AOD = random_draw.*phi_n_AOD_tmp+angle_jitter+obj.phi_LOS_AOD;
            if obj.LOS % && ~obj.O2I
                phi_n_AOD = (random_draw.*phi_n_AOD_tmp+angle_jitter)-(random_draw(1)*phi_n_AOD_tmp(1)+angle_jitter(1)-obj.phi_LOS_AOD);
            end
            ray_angle_offset = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];
            phi_n_m_AOD_ = repmat(phi_n_AOD.',1,obj.M)+obj.c_ASD(ones(obj.N_new,1))*ray_angle_offset;
            phi_n_m_AOD_ = mod(phi_n_m_AOD_,360);
            phi_n_m_AOD_ = phi_n_m_AOD_-360*floor(phi_n_m_AOD_/180);
            obj.phi_n_m_AOD = phi_n_m_AOD_;
        end

        % caculator the ZOA for every ray
        function ZOA_calc(obj)
            %             %                                                                         10      11    12                    15                            19     20
            %             nlos_az_scale_table = [0.582,0.642,0.698,0.750,0.798,0.842,0.882,0.92,0.954, 0.9854, 1.013, 1.04, 1.063, 1.086, 1.1088,1.1257,1.1426,1.1595, 1.1764, 1.1918];
            % %             nlos_az_scale_table = [0.9854,0.9854,0.9854,0.9854,0.9854,0.9854,0.9854,0.9854,0.9854, 0.9854, 1.013, 1.04, 1.063, 1.086, 1.1088,1.1257,1.1426,1.1595, 1.1764, 1.1918];
            %             nlos_az_scale = nlos_az_scale_table(obj.N);
            n_clusters       = [8 10 11 12 15 19 20 25];
            nlos_zenith_scale_table = [0.889, 0.957, 1.031, 1.104, 1.1088, 1.184, 1.178, 1.178]; % 1.282
            nlos_zenith_scale     = nlos_zenith_scale_table(n_clusters == obj.N);
            zenith_scale          = sum([nlos_zenith_scale ,nlos_zenith_scale*(0.3086+0.0339*obj.K-0.0077*obj.K^2+0.0002*obj.K^3)]); % for NLOS and O2I, K = [].
            theta_n_ZOA_tmp  = -obj.ZSA*log(obj.Pn_LOS/max(obj.Pn_LOS))/zenith_scale;

                random_draw  = randsrc(1,obj.N_new,[-1, 1]);
                angle_jitter  = (obj.ZSA/7)*randn(1,obj.N_new);
            % end
            ZOA = obj.theta_LOS_ZOA;
            ZOA(obj.O2I) = 90;
            theta_n_ZOA = random_draw.*theta_n_ZOA_tmp+angle_jitter+ZOA;
            if obj.LOS  && ~obj.O2I
                theta_n_ZOA = (random_draw.*theta_n_ZOA_tmp+angle_jitter)-(random_draw(1)*theta_n_ZOA_tmp(1)+angle_jitter(1)-obj.theta_LOS_ZOA);
            end
            ray_angle_offset = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];
            theta_n_m_ZOA_ = repmat(theta_n_ZOA.',1,obj.M) + obj.c_ZSA(ones(obj.N_new,1))*ray_angle_offset;
            wrapped_zoa_angle = mod((theta_n_m_ZOA_+360),360);  % 360
            over_half_circle = ((wrapped_zoa_angle <= 360)&(wrapped_zoa_angle >= 180));
            theta_n_m_ZOA_ = (~over_half_circle).*wrapped_zoa_angle + over_half_circle.*(360 - wrapped_zoa_angle);
            obj.theta_n_m_ZOA = theta_n_m_ZOA_;
        end

        % caculator the ZOD for every ray
        function ZOD_calc(obj)
            %             %                                                                         10      11    12                    15                            19     20
            %             nlos_az_scale_table = [0.582,0.642,0.698,0.750,0.798,0.842,0.882,0.92,0.954, 0.9854, 1.013, 1.04, 1.063, 1.086, 1.1088,1.1257,1.1426,1.1595, 1.1764, 1.1918];
            % %             nlos_az_scale_table = [0.9854,0.9854,0.9854,0.9854,0.9854,0.9854,0.9854,0.9854,0.9854, 0.9854, 1.013, 1.04, 1.063, 1.086, 1.1088,1.1257,1.1426,1.1595, 1.1764, 1.1918];
            %             nlos_az_scale = nlos_az_scale_table(obj.N);
            n_clusters       = [8 10 11 12 15 19 20 25];
            nlos_zenith_scale_table = [0.889, 0.957, 1.031, 1.104, 1.1088, 1.184, 1.178, 1.178]; % 1.282
            nlos_zenith_scale     = nlos_zenith_scale_table(n_clusters == obj.N);
            zenith_scale          = sum([nlos_zenith_scale ,nlos_zenith_scale*(0.3086+0.0339*obj.K-0.0077*obj.K^2+0.0002*obj.K^3)]); % for NLOS and O2I, K = [].
            theta_n_ZOD_tmp  = -obj.ZSD*log(obj.Pn_LOS/max(obj.Pn_LOS))/zenith_scale;

                random_draw = randsrc(1,obj.N_new,[-1, 1]);
                angle_jitter = (obj.ZSD/7)*randn(1,obj.N_new);
            % end
            theta_n_ZOD = random_draw.*theta_n_ZOD_tmp+angle_jitter+obj.theta_LOS_ZOD+obj.mu_offset_ZOD;
            if obj.LOS % && ~obj.O2I
                theta_n_ZOD = (random_draw.*theta_n_ZOD_tmp+angle_jitter)-(random_draw(1)*theta_n_ZOD_tmp(1)+angle_jitter(1)-obj.theta_LOS_ZOD);
            end
            ray_angle_offset = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];
            zsd_offset_scale = (3/8)*(10^obj.mu_lgZSD);
            % if strcmp(obj.scenario.name,'UrbanGrid')
            %     c_ZSD = obj.c_ZSD;
            % end
            theta_n_m_ZOD_  = repmat(theta_n_ZOD.',1,obj.M) + zsd_offset_scale(ones(obj.N_new,1))*ray_angle_offset;
            wrapped_zod_angle = mod((theta_n_m_ZOD_+360),360);   % 360
            over_half_circle = ((wrapped_zod_angle <= 360)&(wrapped_zod_angle >= 180));
            theta_n_m_ZOD_ = (~over_half_circle).*wrapped_zod_angle + over_half_circle.*(360 - wrapped_zod_angle);
            obj.theta_n_m_ZOD = theta_n_m_ZOD_;
        end

        function RandomCouplingRays(obj)
            [rn1, rn2, rn3] = obj.drawCouplingSeeds();
            rayGroups = obj.strongClusterRayGroups();

            if obj.isRP
                obj.coupleReflectionPointRays(rayGroups);
            else
                obj.coupleCommunicationRays(rn1, rn2, rn3, rayGroups);
            end
        end

        % step 9: Generate XPRs
        function generate_XPRs(obj)
                xpr_randn = randn(obj.N_new, obj.M);
            % end
            X_n_m = obj.sigma_XPR*xpr_randn+obj.mu_XPR;
            %                 X_n_m(X_n_m<0) = 0;
            obj.XPR_n_m = 10.^(X_n_m/10);

        end
        % step 10: Draw initial random phases
        function initial_random_phases(obj)
                obj.PHI_n_m = (2*rand([2, 2, obj.N_new, obj.M])-1)*pi;
            % end
        end
        
        %% Step 11: Generate channel coefficients for each cluster n and each receiver and transmitter element pair u, s.
        function generate_channel(obj)
            lambda      = 3e8/obj.fc;
            if obj.fastfading_enable
                t0 = 0;
                if isempty(obj.TX.velocity)
                    v_tx = [0, 0, 0];
                else
                    v_tx = obj.TX.velocity * [cosd(obj.TX.phi_v),sind(obj.TX.phi_v),cosd(obj.TX.theta_v)];
                end
                if isempty(obj.RX.velocity)
                    v_rx = [0, 0, 0];
                else
                    v_rx = obj.RX.velocity * [cosd(obj.RX.phi_v),sind(obj.RX.phi_v),cosd(obj.RX.theta_v)];
                end

                tx_panel_num = obj.TX.sector.antenna.num_panel;
                Ptx = obj.TX.antenna_params.panel.P;
                E_tx = obj.TX.sector.antenna.panel{1, 1}.num_element/Ptx;
                S_ = tx_panel_num*E_tx;
                obj.S = S_;

                if strcmp(obj.RX.type,'RP')
                    Prx = Ptx;
                    E_rx = E_tx;
                    U_ = S_;
                    obj.U = U_;
                else
                    rx_panel_num = obj.RX.sector.antenna.num_panel;
                    Prx = obj.RX.antenna_params.panel.P;
                    E_rx = obj.RX.sector.antenna.panel{1, 1}.num_element/Prx;
                    U_ = rx_panel_num*E_rx;
                end

                %% BS params
                pos_panelTX = reshape(permute(obj.TX.sector.antenna.pos_panel_LCS,[3,1,2]),3,[]);
                

                dTX = zeros(3, S_);
                for s = 1:S_
                    panS = floor((s-1)/E_tx)+1;
                    idxS = mod((s-1),E_tx)+1;
                    tx_pos = reshape(permute(obj.TX.sector.antenna.panel{panS}.pos_element_LCS,[3,1,2]),3,[]);
                    dTX(:,s) = obj.TX.sector.antenna.R * (pos_panelTX(:,panS) + tx_pos(:,idxS)) * lambda;
                end

                if strcmp(obj.RX.type,'RP')
                    dRX = dTX;
                else
                    pos_panelRX = reshape(permute(obj.RX.sector.antenna.pos_panel_LCS,[3,1,2]),3,[]);
                    dRX = zeros(3, U_);
                    for u = 1:U_
                        panU = floor((u-1)/E_rx)+1;
                        idxU = mod((u-1),E_rx)+1;
                        rx_pos = reshape(permute(obj.RX.sector.antenna.panel{panU}.pos_element_LCS,[3,1,2]),3,[]);
                        dRX(:,u) = obj.RX.sector.antenna.R * (pos_panelRX(:,panU) + rx_pos(:,idxU)) * lambda;
                    end
                end

                % phase_tx(l,s) = exp(j*2*pi/lambda * r_tx(:,l)^T dTX(:,s))
                r_tx_n_m(1,:,:) = sind(obj.theta_n_m_ZOD).*cosd(obj.phi_n_m_AOD);
                r_tx_n_m(2,:,:) = sind(obj.theta_n_m_ZOD).*sind(obj.phi_n_m_AOD);
                r_tx_n_m(3,:,:) = cosd(obj.theta_n_m_ZOD);
                r_rx_n_m(1,:,:) = sind(obj.theta_n_m_ZOA).*cosd(obj.phi_n_m_AOA); % the spherical unit vector
                r_rx_n_m(2,:,:) = sind(obj.theta_n_m_ZOA).*sind(obj.phi_n_m_AOA);
                r_rx_n_m(3,:,:) = cosd(obj.theta_n_m_ZOA);

                phase_tx = exp(1j*2*pi/lambda * squeeze(sum(dTX.*r_tx_n_m,1)));   % (L x S)
                phase_rx = exp(1j*2*pi/lambda * squeeze(sum(dRX.*r_rx_n_m,1)));    % (L x U)

                doppler_tx = squeeze(sum(v_tx(:).*r_tx_n_m,1))/lambda;
                doppler_rx = squeeze(sum(v_rx(:).*r_rx_n_m,1))/lambda;

                doppler_total_ = doppler_tx + doppler_rx;
                phase_doppler_ = exp(1j*2*pi*doppler_total_*(obj.t-t0)); % leak integrted

                path_ray_mask = zeros([U_,S_,obj.N_new,obj.M,obj.N_new+4]);
                path_idx       = 0;

                if obj.LOS && ~obj.O2I
                    active_cluster_delay = obj.tau_n_LOS;
                    rice_linear = 10.^(obj.K/10);
                else
                    active_cluster_delay = obj.tau_n;
                    rice_linear = 0;
                end
                path_delay            = zeros(1,obj.N_new+4);
                for n = 1:obj.N_new
                    if ~isempty(find(obj.strong_cluster_id == n, 1))
                        subpath_count = 3;
                    else
                        subpath_count = 1;
                    end
                    for i = 1:subpath_count
                        % for next path
                        path_idx = path_idx + 1;

                        if i == 1
                            path_delay(path_idx) = active_cluster_delay(n);
                            if ~isempty(find(obj.strong_cluster_id == n, 1))
                                ray_idx_list = [1, 2, 3, 4, 5, 6, 7, 8, 19, 20];
                            else
                                ray_idx_list = 1:obj.M;
                            end
                        elseif i == 2
                            ray_idx_list = [9, 10, 11, 12, 17, 18];
                            path_delay(path_idx) = active_cluster_delay(n)+5e-9;
                        else
                            ray_idx_list = [13, 14, 15, 16];
                            path_delay(path_idx) = active_cluster_delay(n)+10e-9;
                        end
                        path_ray_mask(:,:,n,ray_idx_list,path_idx) = 1;
                    end
                end
                obj.tau_channel = path_delay(1:path_idx);



                XPR_nm(1,1,:,:)  = squeeze(obj.XPR_n_m(:,:,1));
                %             XPR_nm2(1,1,:,:) = squeeze(XPR_n_m_(:,:,2));
                pol_coupling_mat(1,1,:,:) = exp(1j*obj.PHI_n_m(1,1,:,:));                        pol_coupling_mat(1,2,:,:) = sqrt(1./XPR_nm).*exp(1j*obj.PHI_n_m(1, 2,:,:));
                pol_coupling_mat(2,1,:,:) = sqrt(1./XPR_nm).*exp(1j*obj.PHI_n_m(2, 1,:,:));      pol_coupling_mat(2,2,:,:) = exp(1j*obj.PHI_n_m(2, 2,:,:));
                cluster_power_weight              = sqrt(obj.Pn/obj.M/(rice_linear+1));
                ray_power_sqrt           = repmat(cluster_power_weight',1,obj.M);
                
                

                tx_field_pattern = complex(zeros(2, S_, obj.N_new,obj.M, Ptx));
                
                for ptx = 1:Ptx
                    % field_pattern should accept vector angles and return 2 x L (or 2 x 1 x L)
                    element_field_pattern = obj.TX.sector.antenna.panel{1}.element_list(ptx).field_pattern(obj.phi_n_m_AOD, obj.theta_n_m_ZOD); % expect 2 x L
                    element_field_pattern = reshape(element_field_pattern, 2, 1, obj.N_new,obj.M);  % 2 x 1 x L
                    tx_field_pattern(:,:,:,:,ptx) = element_field_pattern;
                end
                if strcmp(obj.RX.type,'RP')
                    rx_field_pattern  = tx_field_pattern;
                else
                    rx_field_pattern = complex(zeros(2, U_, obj.N_new,obj.M, Prx));
                    for prx = 1:Prx
                        element_field_pattern = obj.RX.sector.antenna.panel{1}.element_list(prx).field_pattern(obj.phi_n_m_AOA, obj.theta_n_m_ZOA); % 2 x L
                        element_field_pattern = reshape(element_field_pattern, 2, 1, obj.N_new,obj.M);
                        rx_field_pattern(:,:,:,:,prx) = element_field_pattern;
                    end
                end
                polarized_channel = complex(zeros(U_, S_, obj.N_new,obj.M, Prx, Ptx));

                for prx = 1:Prx
                    Frx_u2L = permute(rx_field_pattern(:,:,:,:,prx), [2 1 3 4]);   % 2 x U x L
                    for ptx = 1:Ptx
                        Ftx_2sL = tx_field_pattern(:,:,:,:,ptx);                 % 2 x S x L
                        polarized_channel(:,:,:,:,prx,ptx) =  pagemtimes(Frx_u2L, pagemtimes(pol_coupling_mat, Ftx_2sL));  % U x S x L
                    end
                end
                ray_phase_weight = reshape(ray_power_sqrt.*phase_doppler_.*phase_tx.*phase_rx, 1,1,obj.N_new,obj.M);
                ray_power_channel = sum(abs(polarized_channel.* ray_phase_weight).^2, [5 6]);   % U x S
                path_power_map = zeros(U_,S_,path_idx);
                for pa = 1:path_idx
                    path_power_map(:,:,pa) = squeeze(sum(sum(ray_power_channel.*path_ray_mask(:,:,:,:,pa),4),3));
                end
                 

                if obj.LOS && ~obj.O2I
                    r_tx_LOS = [sind(obj.theta_LOS_ZOD)*cosd(obj.phi_LOS_AOD);  % the spherical unit vector
                        sind(obj.theta_LOS_ZOD)*sind(obj.phi_LOS_AOD);
                        cosd(obj.theta_LOS_ZOD)];
                    r_rx_LOS = [sind(obj.theta_LOS_ZOA)*cosd(obj.phi_LOS_AOA);  % the spherical unit vector
                    sind(obj.theta_LOS_ZOA)*sind(obj.phi_LOS_AOA);
                    cosd(obj.theta_LOS_ZOA)];

                    phase_tx_LOS = exp(1j*2*pi/lambda * squeeze(sum(dTX.*r_tx_LOS,1)));   % (L x S)
                    phase_rx_LOS = exp(1j*2*pi/lambda * squeeze(sum(dRX.*r_rx_LOS,1)));    % (L x U)

                    doppler_tx_LOS = squeeze(sum(v_tx(:).*r_tx_LOS,1))/lambda;
                    doppler_rx_LOS = squeeze(sum(v_rx(:).*r_rx_LOS,1))/lambda;

                    doppler_total_LOS = doppler_tx_LOS + doppler_rx_LOS;
                    phase_doppler_LOS = exp(1j*2*pi*doppler_total_LOS*(obj.t-t0)); % leak integrted

                    los_pol_coupling_mat(1,1,:,:) = ones(size(obj.PHI_n_m(1,1,:,:)));       los_pol_coupling_mat(1,2,:,:) = zeros(size(obj.PHI_n_m(1,1,:,:)));
                    los_pol_coupling_mat(2,1,:,:) = zeros(size(obj.PHI_n_m(1,1,:,:)));      los_pol_coupling_mat(2,2,:,:) = ones(size(obj.PHI_n_m(1,1,:,:)));

                    tx_los_field_pattern = complex(zeros(2, S_, obj.N_new,obj.M, Ptx));
                    
                    for ptx = 1:Ptx
                        % field_pattern should accept vector angles and return 2 x L (or 2 x 1 x L)
                        element_field_pattern = obj.TX.sector.antenna.panel{1}.element_list(ptx).field_pattern(obj.phi_n_m_AOD, obj.theta_n_m_ZOD); % expect 2 x L
                        element_field_pattern = reshape(element_field_pattern, 2, 1, obj.N_new,obj.M);  % 2 x 1 x L
                        tx_los_field_pattern(:,:,:,:,ptx) = element_field_pattern;
                    end
                    
                    if strcmp(obj.RX.type,'RP')
                        rx_los_field_pattern  = tx_los_field_pattern;
                    else
                        rx_los_field_pattern = complex(zeros(2, S_, obj.N_new,obj.M, Ptx));
                        for prx = 1:Prx
                            element_field_pattern = obj.RX.sector.antenna.panel{1}.element_list(prx).field_pattern(obj.phi_n_m_AOA, obj.theta_n_m_ZOA); % 2 x L
                            element_field_pattern = reshape(element_field_pattern, 2, 1, obj.N_new,obj.M);
                            rx_los_field_pattern(:,:,:,:,prx) = element_field_pattern;
                        end
                    end
                    los_polarized_channel = complex(zeros(U_, S_,  obj.N_new,obj.M, Prx, Ptx));

                    for prx = 1:Prx
                        Frx_u2L = permute(rx_los_field_pattern(:,:,:,:,prx), [2 1 3 4]);   % 2 x U x L
                        for ptx = 1:Ptx
                            Ftx_2sL = tx_los_field_pattern(:,:,:,:,ptx);                 % 2 x S x L
                            los_polarized_channel(:,:,:,:,prx,ptx) =  pagemtimes(Frx_u2L, pagemtimes(los_pol_coupling_mat, Ftx_2sL));  % U x S x L
                            polarized_channel(:,:,:,:,prx,ptx) = polarized_channel(:,:,:,:,prx,ptx) + los_polarized_channel(:,:,:,:,prx,ptx);  % U x S x L
                        end
                    end
                    los_phase_weight = phase_doppler_LOS.*reshape(phase_tx_LOS,1,1,[]).*reshape(phase_rx_LOS,1,1,[]);
                    los_power_channel = sum(abs(los_polarized_channel.* sqrt(rice_linear/(rice_linear+1)) .*los_phase_weight).^2, [5 6]);
                    path_power_map(:,:,1) = los_power_channel(:,:,1);   % U x S
                end       

                % ===== "P part same as original": noncoherent power sum over pols =====
                calibration_power = sum(abs(polarized_channel(:,:,:,:,1,1).* ray_power_sqrt).^2,"all");   % U x S
                obj.H_fastfading = path_power_map;   % U x S                
                obj.H_cali = pow2db(calibration_power);
                obj.H_full = obj.H_fastfading * db2pow(-obj.PL+obj.SF);
            else
                obj.H_full = db2pow(-obj.PL+obj.SF);
            end
        end

        function [couplingloss_ls,couplingloss_full, pathloss] = calc_coupling_loss(obj)
             couplingloss_ls = -(obj.PL - obj.SF);
             pathloss = obj.PL - obj.SF;
             if obj.fastfading_enable
                 couplingloss_full =-obj.PL + obj.SF - pow2db(sum(sum(obj.H_fastfading))) + 10*log10(obj.U*obj.TX.antenna_params.panel.P) ;
             else
                 couplingloss_full = nan;
             end
        end

        function [ as, mean_angle ]  = calc_angular_spreads(obj,ang,LOSang,wrap_angles)
            % CALC_ANGULAR_SPREADS Calculates the angular spread in degree
            valid_mask = ~isnan(ang);
            ang = ang(valid_mask);
            ang = ang(:).'*pi/180;
            pow = repmat((obj.Pn/obj.M).',1,obj.M);
            if ~isempty(obj.loss_blockage)
                pow = pow .* obj.loss_blockage;
            end
            pow = pow(valid_mask);
            pow = pow(:).';
            % Normalize powers
            pt = sum( pow,2 );
            pow = pow./pt;
            if obj.LOS && ~obj.O2I
                KF = 10^(obj.K/10);
                pow = pow*(1/(1+KF));
                ang(numel(pow)+1) = LOSang*pi/180;
                pow(numel(pow)+1) = KF/(1+KF);
            end
            if ~exist('wrap_angles','var')
                wrap_angles = true;
            end
            if wrap_angles
                mean_angle = angle( sum( pow.*exp( 1j*ang ) , 2 ) ); % [rad]
            else
                mean_angle = sum( pow.*ang,2 );
            end
            phi = ang - mean_angle;
            if wrap_angles
                phi = angle( exp( 1j*phi ) );
            end
            as = sqrt( sum(pow.*(phi.^2),2)); %  - sum( pow.*phi,2).^2
            mean_angle = mean_angle*180/pi;
            as = as*180/pi;
        end

        function [ ds, mean_delay ] = calc_delay_spread(obj)
            % CALC_DELAY_SPREAD Calculates the delay spread in [s]

            % Normalize powers
            pt = sum(obj.Pn);
            pow = obj.Pn./pt;
            taus = obj.tau_n;

            mean_delay = sum( pow.*taus );
            sort_idx = taus - mean_delay;
            ds = sqrt( sum(pow.*(sort_idx.^2))); %  - sum(pow.*sort_idx).^2
            % ds = sqrt( sum(pow.*((taus).^2)) - mean_delay.^2 ); %  - sum(pow.*sort_idx).^2
        end

    end

    methods(Access = private)
        function assertProcedureBSupported(obj)
            if ~ismember(obj.scenario.name, {'UrbanGrid', 'UMi', 'InH'})
                error('CommChannel:UnsupportedSpatialConsistencyProcedureB', ...
                    'Spatial consistency Procedure B is currently implemented only for UrbanGrid, UMi and InH scenarios.');
            end
        end

        function assertProcedureASupported(obj)
            if ~ismember(obj.scenario.name, {'UrbanGrid', 'UMi'})
                error('CommChannel:UnsupportedSpatialConsistencyProcedureA', ...
                    'Spatial consistency Procedure A is currently implemented only for UrbanGrid and UMi scenarios.');
            end
        end

        function limits = procedureBClusterLimits(obj)
            fc_GHz = obj.fc/1e9;
            if strcmp(obj.scenario.name, 'UrbanGrid')
                fc_eval = max(fc_GHz, 6);
                if obj.O2I
                    mu_lgDS = -6.62;
                    sigma_lgDS = 0.32;
                    mu_lgASD = 1.25;
                    sigma_lgASD = 0.42;
                    mu_lgASA = 1.76;
                    sigma_lgASA = 0.16;
                    mu_lgZSA = 1.01;
                    sigma_lgZSA = 0.43;
                elseif obj.LOS
                    mu_lgDS = -6.955 - 0.0963*log10(fc_eval);
                    sigma_lgDS = 0.66;
                    mu_lgASD = 1.06 + 0.1114*log10(fc_eval);
                    sigma_lgASD = 0.28;
                    mu_lgASA = 1.81;
                    sigma_lgASA = 0.20;
                    mu_lgZSA = 0.95;
                    sigma_lgZSA = 0.16;
                else
                    mu_lgDS = -6.28 - 0.204*log10(fc_eval);
                    sigma_lgDS = 0.39;
                    mu_lgASD = 1.5 - 0.1144*log10(fc_eval);
                    sigma_lgASD = 0.28;
                    mu_lgASA = 2.08 - 0.27*log10(fc_eval);
                    sigma_lgASA = 0.11;
                    mu_lgZSA = -0.3236*log10(fc_eval) + 1.512;
                    sigma_lgZSA = 0.16;
                end
            elseif strcmp(obj.scenario.name, 'UMi')
                fc_eval = max(fc_GHz, 2);
                if obj.O2I
                    mu_lgDS = -6.62;
                    sigma_lgDS = 0.32;
                    mu_lgASD = 1.25;
                    sigma_lgASD = 0.42;
                    mu_lgASA = 1.76;
                    sigma_lgASA = 0.16;
                    mu_lgZSA = 1.01;
                    sigma_lgZSA = 0.43;
                elseif obj.LOS
                    mu_lgDS = -0.24*log10(1 + fc_eval) - 7.14;
                    sigma_lgDS = 0.38;
                    mu_lgASD = -0.05*log10(1 + fc_eval) + 1.21;
                    sigma_lgASD = 0.41;
                    mu_lgASA = -0.08*log10(1 + fc_eval) + 1.73;
                    sigma_lgASA = 0.014*log10(1 + fc_eval) + 0.28;
                    mu_lgZSA = -0.1*log10(1 + fc_eval) + 0.73;
                    sigma_lgZSA = -0.04*log10(1 + fc_eval) + 0.34;
                else
                    mu_lgDS = -0.24*log10(1 + fc_eval) - 6.83;
                    sigma_lgDS = 0.16*log10(1 + fc_eval) + 0.28;
                    mu_lgASD = -0.23*log10(1 + fc_eval) + 1.53;
                    sigma_lgASD = 0.11*log10(1 + fc_eval) + 0.33;
                    mu_lgASA = -0.08*log10(1 + fc_eval) + 1.81;
                    sigma_lgASA = 0.05*log10(1 + fc_eval) + 0.3;
                    mu_lgZSA = -0.04*log10(1 + fc_eval) + 0.92;
                    sigma_lgZSA = -0.07*log10(1 + fc_eval) + 0.41;
                end
            elseif strcmp(obj.scenario.name, 'InH')
                fc_eval = max(fc_GHz, 6);
                if obj.LOS
                    mu_lgDS = -0.01*log10(1 + fc_eval) - 7.692;
                    sigma_lgDS = 0.18;
                    mu_lgASD = 1.60;
                    sigma_lgASD = 0.18;
                    mu_lgASA = -0.19*log10(1 + fc_eval) + 1.781;
                    sigma_lgASA = 0.12*log10(1 + fc_eval) + 0.119;
                    mu_lgZSA = -0.26*log10(1 + fc_eval) + 1.44;
                    sigma_lgZSA = -0.04*log10(1 + fc_eval) + 0.264;
                else
                    mu_lgDS = -0.28*log10(1 + fc_eval) - 7.173;
                    sigma_lgDS = 0.1*log10(1 + fc_eval) + 0.055;
                    mu_lgASD = 1.62;
                    sigma_lgASD = 0.25;
                    mu_lgASA = -0.11*log10(1 + fc_eval) + 1.863;
                    sigma_lgASA = 0.12*log10(1 + fc_eval) + 0.059;
                    mu_lgZSA = -0.15*log10(1 + fc_eval) + 1.387;
                    sigma_lgZSA = -0.09*log10(1 + fc_eval) + 0.746;
                end
            else
                obj.assertProcedureBSupported();
            end

            limits.delay = 2*10^(mu_lgDS + sigma_lgDS);
            limits.AOA = 2*10^(mu_lgASA + sigma_lgASA);
            limits.AOD = 2*10^(mu_lgASD + sigma_lgASD);
            limits.ZOA = 2*10^(mu_lgZSA + sigma_lgZSA);
            limits.ZOD = 2*10^(obj.mu_lgZSD + obj.sigma_lgZSD);
        end

        function uniform_values = drawProcedureBUniform(obj, field_name, count)
            uniform_values = rand(1, count);
            if ~obj.isSpatialConsistencyProcedureB()
                return;
            end

            [external_values, isAvailable] = obj.getExternalProcedureBRaw(field_name, count);
            if isAvailable
                uniform_values = external_values;
            else
                obj.warnSpatialConsistencyFallback();
            end
        end

        function gaussian_values = drawProcedureBGaussian(obj, field_name, count)
            gaussian_values = randn(1, count);
            if ~obj.isSpatialConsistencyProcedureB()
                return;
            end

            [external_values, isAvailable] = obj.getExternalProcedureBRaw(field_name, count);
            if isAvailable
                gaussian_values = external_values;
            else
                obj.warnSpatialConsistencyFallback();
            end
        end

        function [values, isAvailable] = getExternalProcedureBRaw(obj, field_name, count)
            values = [];
            isAvailable = false;
            if isempty(obj.TX) || isempty(obj.RX) || isempty(obj.RX.ID)
                return;
            end

            if ~isprop(obj.TX, 'SC_procB_raw') || isempty(obj.TX.SC_procB_raw)
                return;
            end

            rxId = obj.RX.ID;
            if rxId < 1 || numel(obj.TX.SC_procB_raw) < rxId
                return;
            end

            proc_raw = obj.TX.SC_procB_raw(rxId);
            if ~isfield(proc_raw, field_name)
                return;
            end

            raw_values = proc_raw.(field_name);
            if numel(raw_values) < count || any(~isfinite(raw_values(1:count)))
                return;
            end

            values = raw_values(1:count);
            isAvailable = true;
        end

        function raw_rand_vec = drawLspRawRandn(obj, vectorLength)
            raw_rand_vec = randn(vectorLength, 1);
            if ~obj.isSpatialConsistencyEnabled()
                return;
            end

            [external_lsp, isAvailable] = obj.getExternalLspRaw(vectorLength);
            if isAvailable
                raw_rand_vec = external_lsp;
            else
                obj.warnSpatialConsistencyFallback();
            end
        end

        function tf = isSpatialConsistencyEnabled(obj)
            tf = isprop(obj.scenario, 'spatial_consistency_enable') && ...
                obj.scenario.spatial_consistency_enable;
        end

        function tf = isSpatialConsistencyProcedureB(obj)
            tf = isprop(obj.scenario, 'spatial_consistency_enable') && ...
                obj.scenario.spatial_consistency_enable && ...
                isprop(obj.scenario, 'spatial_consistency_procedure') && ...
                strcmpi(obj.scenario.spatial_consistency_procedure, 'B');
        end

        function tf = isSpatialConsistencyProcedureA(obj)
            tf = isprop(obj.scenario, 'spatial_consistency_enable') && ...
                obj.scenario.spatial_consistency_enable && ...
                isprop(obj.scenario, 'spatial_consistency_procedure') && ...
                strcmpi(obj.scenario.spatial_consistency_procedure, 'A');
        end

        function tf = hasProcedureAState(obj)
            state = obj.getProcedureAState();
            tf = ~isempty(state) && isfield(state, 'time') && isfinite(state.time);
        end

        function state = getProcedureAState(obj)
            state = [];
            if isempty(obj.TX) || isempty(obj.RX) || isempty(obj.RX.ID) || ...
                    ~isprop(obj.TX, 'SC_procA_comm_state') || isempty(obj.TX.SC_procA_comm_state)
                return;
            end

            rxId = obj.RX.ID;
            if rxId < 1 || numel(obj.TX.SC_procA_comm_state) < rxId
                return;
            end
            candidate = obj.TX.SC_procA_comm_state(rxId);
            if isstruct(candidate) && isfield(candidate, 'initialized') && candidate.initialized
                state = candidate;
            end
        end

        function saveProcedureAState(obj)
            if isempty(obj.TX) || isempty(obj.RX) || isempty(obj.RX.ID) || ...
                    ~isprop(obj.TX, 'SC_procA_comm_state')
                return;
            end

            rxId = obj.RX.ID;
            state = obj.makeProcedureAState();
            if isempty(obj.TX.SC_procA_comm_state)
                state_table = repmat(obj.emptyProcedureAState(), 1, rxId);
            else
                state_table = obj.alignProcedureAStateTable(obj.TX.SC_procA_comm_state);
            end
            if numel(state_table) < rxId
                state_table(end+1:rxId) = repmat(obj.emptyProcedureAState(), 1, rxId - numel(state_table));
            end
            state_table(rxId) = state;
            obj.TX.SC_procA_comm_state = state_table;
            obj.TX.SC_time_nodes = unique([obj.TX.SC_time_nodes, obj.t]);
        end

        function state = makeProcedureAState(obj)
            state = obj.emptyProcedureAState();
            state.initialized = true;
            state.time = obj.t;
            state.tx_position = obj.TX.Position;
            state.rx_position = obj.RX.Position;
            state.tau_n = obj.tau_n;
            if isempty(obj.tau_absolute)
                state.tau_absolute = obj.d_3D/3e8 + obj.tau_n;
            else
                state.tau_absolute = obj.tau_absolute;
            end
            state.tau_prime = obj.tau_prime;
            state.tau_order = obj.tau_order;
            state.phi_prime_AOA = obj.procedureAPrimeFromRays(obj.phi_n_m_AOA, obj.phi_LOS_AOA, true);
            state.phi_prime_AOD = obj.procedureAPrimeFromRays(obj.phi_n_m_AOD, obj.phi_LOS_AOD, true);
            state.theta_prime_ZOA = obj.procedureAPrimeFromRays(obj.theta_n_m_ZOA, obj.theta_LOS_ZOA, false);
            state.theta_prime_ZOD = obj.procedureAPrimeFromRays(obj.theta_n_m_ZOD, obj.theta_LOS_ZOD, false);
            state.phi_LOS_AOA = obj.phi_LOS_AOA;
            state.phi_LOS_AOD = obj.phi_LOS_AOD;
            state.theta_LOS_ZOA = obj.theta_LOS_ZOA;
            state.theta_LOS_ZOD = obj.theta_LOS_ZOD;
            state.cluster_shadow_dB = obj.cluster_shadow_dB;
            previous_state = obj.getProcedureAState();
            state.procedureA_Xn = obj.procedureAXnFromState(previous_state, numel(state.tau_n));
            state.procedureA_Xn_correlation_distance = obj.procedureAXnCorrelationDistance();
        end

        function state = emptyProcedureAState(obj) %#ok<MANU>
            state = struct( ...
                'initialized', false, ...
                'time', nan, ...
                'tx_position', nan(1,3), ...
                'rx_position', nan(1,3), ...
                'tau_n', [], ...
                'tau_absolute', [], ...
                'tau_prime', [], ...
                'tau_order', [], ...
                'phi_prime_AOA', [], ...
                'phi_prime_AOD', [], ...
                'theta_prime_ZOA', [], ...
                'theta_prime_ZOD', [], ...
                'phi_LOS_AOA', nan, ...
                'phi_LOS_AOD', nan, ...
                'theta_LOS_ZOA', nan, ...
                'theta_LOS_ZOD', nan, ...
                'cluster_shadow_dB', [], ...
                'procedureA_Xn', [], ...
                'procedureA_Xn_correlation_distance', nan);
        end

        function state_table = alignProcedureAStateTable(obj, state_table_in)
            template = obj.emptyProcedureAState();
            if isempty(state_table_in)
                state_table = state_table_in;
                return;
            end

            state_table = repmat(template, 1, numel(state_table_in));
            field_names = fieldnames(template);
            for idx = 1:numel(state_table_in)
                for field_idx = 1:numel(field_names)
                    field_name = field_names{field_idx};
                    if isfield(state_table_in(idx), field_name)
                        state_table(idx).(field_name) = state_table_in(idx).(field_name);
                    end
                end
            end
        end

        function prime_angles = procedureAPrimeFromRays(obj, ray_angles, los_angle, is_azimuth)
            if isempty(ray_angles)
                prime_angles = [];
                return;
            end

            num_clusters = size(ray_angles, 1);
            prime_angles = nan(1, num_clusters);
            for cluster_idx = 1:num_clusters
                angles = ray_angles(cluster_idx, :);
                angles = angles(isfinite(angles));
                if isempty(angles)
                    continue;
                end
                if is_azimuth
                    cluster_angle = angle(sum(exp(1j * deg2rad(angles)))) * 180/pi;
                    prime_angles(cluster_idx) = obj.wrapAzimuthAngles(cluster_angle - los_angle);
                else
                    cluster_angle = mean(angles);
                    prime_angles(cluster_idx) = cluster_angle - los_angle;
                end
            end
        end

        function unit_vec = unitVectorFromAngles(obj, azimuth_deg, zenith_deg) %#ok<INUSL>
            unit_vec = [ ...
                sind(zenith_deg) * cosd(azimuth_deg), ...
                sind(zenith_deg) * sind(azimuth_deg), ...
                cosd(zenith_deg)];
        end

        function [azimuth_new, zenith_new] = updateProcedureAAnglePair(obj, azimuth_prev, zenith_prev, velocity_vec, delay_vec, delta_t)
            azimuth_new = azimuth_prev;
            zenith_new = zenith_prev;
            if delta_t <= 0 || all(velocity_vec(:) == 0)
                return;
            end

            for cluster_idx = 1:numel(azimuth_prev)
                delay_value = delay_vec(cluster_idx);
                if ~isfinite(delay_value) || delay_value <= 0
                    continue;
                end
                if size(velocity_vec, 1) == 1
                    cluster_velocity = velocity_vec;
                else
                    cluster_velocity = velocity_vec(cluster_idx, :);
                end
                if all(cluster_velocity == 0)
                    continue;
                end

                phi_rad = deg2rad(azimuth_prev(cluster_idx));
                theta_rad = deg2rad(zenith_prev(cluster_idx));
                sin_theta = max(abs(sin(theta_rad)), 1e-6);
                phi_hat = [-sin(phi_rad), cos(phi_rad), 0];
                theta_hat = [cos(theta_rad)*cos(phi_rad), cos(theta_rad)*sin(phi_rad), -sin(theta_rad)];

                d_phi = dot(cluster_velocity, phi_hat) * delta_t / (3e8 * delay_value * sin_theta);
                d_theta = dot(cluster_velocity, theta_hat) * delta_t / (3e8 * delay_value);
                azimuth_new(cluster_idx) = azimuth_prev(cluster_idx) + rad2deg(d_phi);
                zenith_new(cluster_idx) = zenith_prev(cluster_idx) + rad2deg(d_theta);
            end

            azimuth_new = obj.wrapAzimuthAngles(azimuth_new);
            zenith_new = obj.wrapZenithAngles(zenith_new);
        end

        function [rx_velocity_eff, tx_velocity_eff] = procedureAEffectiveVelocities(obj, phi_AOA, theta_ZOA, phi_AOD, theta_ZOD, tx_velocity, rx_velocity, x_n, is_los)
            num_clusters = numel(phi_AOA);
            rx_velocity_eff = zeros(num_clusters, 3);
            tx_velocity_eff = zeros(num_clusters, 3);
            for cluster_idx = 1:num_clusters
                if is_los
                    rx_velocity_eff(cluster_idx, :) = rx_velocity - tx_velocity;
                    tx_velocity_eff(cluster_idx, :) = tx_velocity - rx_velocity;
                    continue;
                end

                mirror_matrix = diag([1, x_n(cluster_idx), 1]);
                phi_aoa = deg2rad(phi_AOA(cluster_idx));
                phi_aod = deg2rad(phi_AOD(cluster_idx));
                theta_zoa = deg2rad(theta_ZOA(cluster_idx));
                theta_zod = deg2rad(theta_ZOD(cluster_idx));
                R_rx = obj.rotationZ(phi_aod + pi) * obj.rotationY(pi/2 - theta_zod) * ...
                    mirror_matrix * obj.rotationY(pi/2 - theta_zoa) * obj.rotationZ(-phi_aoa);
                R_tx = obj.rotationZ(-phi_aod) * obj.rotationY(pi/2 - theta_zod) * ...
                    mirror_matrix * obj.rotationY(pi/2 - theta_zoa) * obj.rotationZ(phi_aoa + pi);
                rx_velocity_eff(cluster_idx, :) = (R_rx * rx_velocity(:)).' - tx_velocity;
                tx_velocity_eff(cluster_idx, :) = (R_tx * tx_velocity(:)).' - rx_velocity;
            end
        end

        function x_n = procedureAXnFromState(obj, state, num_clusters)
            [external_xn, is_available] = obj.getExternalProcedureAXn(num_clusters);
            if is_available
                x_n = external_xn;
            elseif ~isempty(state) && isfield(state, 'procedureA_Xn') && numel(state.procedureA_Xn) >= num_clusters
                x_n = state.procedureA_Xn(1:num_clusters);
            else
                x_n = obj.drawProcedureAXn(num_clusters);
            end
        end

        function [x_n, is_available] = getExternalProcedureAXn(obj, num_clusters)
            x_n = [];
            is_available = false;
            if isempty(obj.TX) || isempty(obj.RX) || isempty(obj.RX.ID) || ...
                    ~isprop(obj.TX, 'SC_procA_comm_Xn') || isempty(obj.TX.SC_procA_comm_Xn)
                return;
            end

            rxId = obj.RX.ID;
            x_n_table = obj.TX.SC_procA_comm_Xn;
            if rxId < 1 || size(x_n_table, 1) < rxId || size(x_n_table, 2) < num_clusters
                return;
            end

            x_n = x_n_table(rxId, 1:num_clusters);
            is_available = isnumeric(x_n) && all(isfinite(x_n));
        end

        function x_n = drawProcedureAXn(obj, num_clusters) %#ok<INUSL>
            x_n = ones(1, num_clusters);
            x_n(rand(1, num_clusters) < 0.5) = -1;
        end

        function corr_dist = procedureAXnCorrelationDistance(obj)
            switch obj.scenario.name
                case 'RMa'
                    corr_dist = 60;
                case {'UMa','UrbanGrid'}
                    corr_dist = 50;
                case 'UMi'
                    corr_dist = 15;
                case {'InH','InF'}
                    corr_dist = 10;
                otherwise
                    corr_dist = nan;
            end
        end

        function R = rotationZ(obj, alpha) %#ok<INUSL>
            R = [cos(alpha), -sin(alpha), 0; ...
                sin(alpha), cos(alpha), 0; ...
                0, 0, 1];
        end

        function R = rotationY(obj, beta) %#ok<INUSL>
            R = [cos(beta), 0, sin(beta); ...
                0, 1, 0; ...
                -sin(beta), 0, cos(beta)];
        end

        function velocity_vec = velocityVector(obj, node) %#ok<INUSL>
            velocity_vec = [0, 0, 0];
            if isempty(node) || ~isprop(node, 'velocity') || isempty(node.velocity)
                return;
            end
            velocity_vec = node.velocity * [cosd(node.phi_v), sind(node.phi_v), cosd(node.theta_v)];
        end

        function [external_lsp, isAvailable] = getExternalLspRaw(obj, vectorLength)
            external_lsp = [];
            isAvailable = false;
            if isempty(obj.TX) || isempty(obj.RX) || isempty(obj.RX.ID)
                return;
            end

            field_name = obj.lspRawFieldName();
            if ~isprop(obj.TX, field_name) || isempty(obj.TX.(field_name))
                return;
            end

            rxId = obj.RX.ID;
            lsp_table = obj.TX.(field_name);
            if rxId < 1 || rxId > size(lsp_table, 2) || vectorLength > size(lsp_table, 1)
                return;
            end

            external_lsp = lsp_table(1:vectorLength, rxId);
            isAvailable = isnumeric(external_lsp) && all(isfinite(external_lsp));
        end

        function field_name = lspRawFieldName(obj)
            if obj.O2I
                field_name = 'LSP_raw_O2I';
            elseif obj.LOS
                field_name = 'LSP_raw_LOS';
            else
                field_name = 'LSP_raw_NLOS';
            end
        end

        function warnSpatialConsistencyFallback(~)
            persistent warned
            if isempty(warned)
                warning('CommChannel:SpatialConsistencyFallback', ...
                    'Spatial consistency LSPs were not applied successfully; falling back to random LSPs.');
                warned = true;
            end
        end

        function angles = wrapAzimuthAngles(~, angles)
            angles = mod(angles, 360);
            angles = angles - 360*floor(angles/180);
        end

        function angles = wrapZenithAngles(~, angles)
            wrapped = mod((angles + 360), 360);
            over_half_circle = ((wrapped <= 360) & (wrapped >= 180));
            angles = (~over_half_circle).*wrapped + over_half_circle.*(360 - wrapped);
        end

        function [rn1, rn2, rn3] = drawCouplingSeeds(obj)
            rn1 = rand(obj.N_new, obj.M);
            rn2 = rand(obj.N_new, obj.M);
            rn3 = rand(obj.N_new, obj.M);
        end

        function coupleCommunicationRays(obj, rn1, rn2, rn3, rayGroups)
            for clusterIdx = 1:obj.N_new
                if obj.hasSubClusters(clusterIdx)
                    for groupIdx = 1:numel(rayGroups)
                        obj.sortCoupledAngles(clusterIdx, rayGroups{groupIdx}, rn1, rn2, rn3);
                    end
                else
                    obj.sortCoupledAngles(clusterIdx, 1:obj.M, rn1, rn2, rn3);
                end
            end
        end

        function sortCoupledAngles(obj, clusterIdx, rayList, rn1, rn2, rn3)
            [~, aoaOrder] = sort(rn1(clusterIdx, rayList));
            obj.phi_n_m_AOA(clusterIdx, rayList) = obj.phi_n_m_AOA(clusterIdx, rayList(aoaOrder));

            [~, zoaOrder] = sort(rn2(clusterIdx, rayList));
            obj.theta_n_m_ZOA(clusterIdx, rayList) = obj.theta_n_m_ZOA(clusterIdx, rayList(zoaOrder));

            [~, aodOrder] = sort(rn3(clusterIdx, rayList));
            obj.phi_n_m_AOD(clusterIdx, rayList) = obj.phi_n_m_AOD(clusterIdx, rayList(aodOrder));
        end

        function coupleReflectionPointRays(obj, rayGroups)
            for clusterIdx = 1:obj.N_new
                if obj.hasSubClusters(clusterIdx)
                    for groupIdx = 1:numel(rayGroups)
                        obj.shuffleRpRayGroup(clusterIdx, rayGroups{groupIdx});
                    end
                else
                    obj.shuffleRpRayMask(clusterIdx, ~isnan(obj.phi_n_m_AOD(clusterIdx, :)));
                end
            end
        end

        function shuffleRpRayGroup(obj, clusterIdx, rayList)
            validRays = rayList(~isnan(obj.phi_n_m_AOD(clusterIdx, rayList)));
            if isempty(validRays)
                return;
            end
            obj.shuffleRpRayMask(clusterIdx, validRays);
        end

        function shuffleRpRayMask(obj, clusterIdx, raySelector)
            order = randperm(numel(raySelector));
            obj.phi_n_m_AOD(clusterIdx, raySelector) = obj.phi_n_m_AOD(clusterIdx, raySelector(order));
            obj.phi_n_m_AOA(clusterIdx, raySelector) = obj.phi_n_m_AOA(clusterIdx, raySelector(order));
        end

        function tf = hasSubClusters(obj, clusterIdx)
            tf = any(obj.strong_cluster_id == clusterIdx);
        end
    end

    methods(Static, Access = private)
        function rayGroups = strongClusterRayGroups()
            rayGroups = {
                [1, 2, 3, 4, 5, 6, 7, 8, 19, 20], ...
                [9, 10, 11, 12, 17, 18], ...
                [13, 14, 15, 16]
            };
        end
    end

end

% other function
function [phi, theta] = cart2sph_(x, y, z)
[phi,ele,~] = cart2sph(x,y,z);
phi         = phi*180/pi;          % [-180, 180]
theta       = 90 - ele*180/pi;     % [   0, 180]
end
