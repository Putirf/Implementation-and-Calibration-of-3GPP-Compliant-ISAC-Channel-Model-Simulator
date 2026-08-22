classdef Target_channel < handle
    % LINK
    properties
        fastfading_enable = false
        t
        static = 1
        TXRX
        sector
        BS_pos_wrap
        ST
        O2I = false
        scenario
        fc
        d_2D
        d_2D_in
        d_2D_out
        d_3D
        bistatic_angle
        option

        phi_LOS_AOD
        theta_LOS_ZOD
        phi_LOS_AOA
        theta_LOS_ZOA
        bLOS

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

        tau_order
        tau_absolute
        tau_n
        tau_n_LOS
        Pn
        Pn_LOS
        cluster_shadow_dB
        keep logical
        tau_order_keep
        tau_n_keep
        tau_n_LOS_keep
        Pn_keep
        Pn_LOS_keep

        N_new
        strong_cluster_id
        map_delay

        phi_n_m_AOA
        phi_n_m_AOD
        theta_n_m_ZOA
        theta_n_m_ZOD

        coupling_txrx_path_table
        path_num
        P_path
        keep_path
        P_path_keep
        coupling_txrx_path_table_keep
        path_num_keep
        RCS_sigma_S_n_m

        mu_lg_delta_tau
        sigma_lg_delta_tau
        path_delay

        XPR_txrx
        XPR_spst_path
        XPR_spst_n_m
        mu_XPR_sensing
        sigma_XPR_sensing
        PHI_n_m_txrx
        PHI_spst_path
        PHI_LOS
        tau_prime
        phi_prime_AOA
        phi_prime_AOD
        theta_prime_ZOA
        theta_prime_ZOD

        CPM_all
        phase_doppler_total

        phi_path_AOA
        phi_path_AOD
        theta_path_ZOA
        theta_path_ZOD

        r_tx_path
        r_rx_path
        r_tx_spst_path
        r_spst_rx_path

        Ftx
        Frx
        % phase_doppler



        loss_blockage

        sigma_RCS
        H_fastfading
        H_cali
        H_full
        RCS_sigma_M_dB

    end

    methods
        function obj = Target_channel(BSsector_TX,BSsector_RX, ST, scenario, fastfading_enable, t,varargin)
            obj.fastfading_enable = fastfading_enable;
            obj.t = t;
            if ~isempty(BSsector_RX)
                obj.TXRX = 2;
            else
                obj.TXRX = 1;
            end
            obj.sector = [BSsector_TX,BSsector_RX];
            if BSsector_TX.ID == BSsector_RX.ID
                obj.static = 1;
            else
                obj.static = 2;
            end
            obj.ST = ST;
            obj.scenario = scenario;
            obj.fc = scenario.frequency;
            obj.option  = scenario.option;
            % obj.fastfading_enable = fastfading_enable;
            % if strcmp(scenario.name,'RMa') ||strcmp(scenario.name,'3D-UMa')||strcmp(scenario.name,'3D-UMi')
            %     [~, wrapped_site_idx] = min(sum(((repmat(ST.pos, size(BSsector.Pos_wrapped,1), 1)-BSsector.Pos_wrapped).^2),2));
            %     obj.BS_pos_wrap = BSsector.Pos_wrapped(wrapped_site_idx,:);
            % elseif strcmp(scenario.name,'3D-InH') || strcmp(scenario.name,'InF')
            %     obj.BS_pos_wrap = BSsector.Position;
            % end
            obj.BS_pos_wrap = [BSsector_TX.Position;BSsector_RX.Position];


            %% step 1
            x = ST.Position(1)-obj.BS_pos_wrap(:,1);
            y = ST.Position(2)-obj.BS_pos_wrap(:,2);
            z = ST.height-obj.BS_pos_wrap(:,3);
            obj.d_2D = sqrt(x.^2+y.^2);
            obj.d_3D = sqrt(x.^2+y.^2+z.^2);
            if ismember(scenario.name,{'RMa','UMa','UMi','UrbanGrid'})
                obj.d_2D_in  = zeros(size(obj.d_2D));
                obj.d_2D_out = obj.d_2D - obj.d_2D_in;
            elseif ismember(scenario.name,{'InH','InF'})
                obj.d_2D_in  = obj.d_2D;
                obj.d_2D_out = zeros(size(obj.d_2D));
            end

            obj.bistatic_angle = acosd(max(-1, min(1, dot([x(1), y(1), z(1)] / norm([x(1), y(1), z(1)]), [x(2), y(2), z(2)] / norm([x(2), y(2), z(2)])))));
            % obj.sigma_RCS = obj.ST.type.RCS.sigma_M_dB -3*sind(obj.bistatic_angle/2)+ obj.ST.type.RCS.enable_small_scale*(obj.ST.type.RCS.sigma_D_dB+obj.ST.type.RCS_sigma_S);

            [obj.phi_LOS_AOD(1), obj.theta_LOS_ZOD(1)] = cart2sph_(x(1), y(1), z(1));
            [obj.phi_LOS_AOA(2), obj.theta_LOS_ZOA(2)] = cart2sph_(x(2), y(2), z(2));
            [obj.phi_LOS_AOA(1), obj.theta_LOS_ZOA(1)] = cart2sph_(-x(1), -y(1), -z(1));
            [obj.phi_LOS_AOD(2), obj.theta_LOS_ZOD(2)] = cart2sph_(-x(2), -y(2), -z(2));

            %% step 2
            obj.los_probability();

            %% step 3
            obj.pathloss();

            %% step 4
            obj.large_scale_para();
            % obj.ASA = min(obj.ASA, 104);
            % obj.ASD = min(obj.ASD, 104);
            % obj.ZSA = min(obj.ZSA, 52);
            % obj.ZSD = min(obj.ZSD, 52);

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
                    obj.cluster_power();
                    %% step 7: Generate arrival angles and departure angles for both azimuth and elevation
                    obj.AOA_calc();
                    obj.AOD_calc();
                    obj.ZOA_calc();
                    obj.ZOD_calc();
                end

                %% step 8: Coupling of rays within a cluster for both azimuth and elevation
                obj.RandomCouplingRays();

                %% step 9: Coupling of rays for a STX-SPST link and the corresponding SPST-SRX link of the same SPST.
                obj.Coupling_txrx();

                %% Step 10: Obtain the power for all generated paths
                obj.Path_power();

                %% Step 11: Obtain the absolute delay for each path in set R
                obj.Absolute_Delay()

                %% step 12: Generate XPRs
                obj.generate_XPRs();

                % The outcome of Steps 1-12 shall be identical for all the links from co-sited sectors to a STX/ST/SRX.
                %% step 13: Draw initial random phases
                obj.initial_random_phases();

                %% step 13.5: calulate part of step 14 parameter which is independent with antenna element
                obj.prepare_channel_parameter();                

                %% Step 14 generate channel
                obj.generate_channel()



            end
        end


        % caculator the LOS probability
        % The LOS probability is derived with assuming antenna heights of
        % 3m for indoor, 10m for UMi, and 25m for UMa
        function los_probability(obj)
            for txrx = 1:obj.TXRX
                d2D = obj.d_2D(txrx);
                d_2Din = obj.d_2D_in(txrx);
                height = obj.ST.height;
            switch obj.scenario.name
                case 'UMa'
                    if height>22.5
                        if strcmp(obj.scenario.subname,'AV')
                            if height<=100
                                d1 = max(460*log10(height)-700,18);
                                p1 = 4300*log10(height)-3800;
                                if d2D <= d1
                                    Pr_LOS = 1;
                                else
                                    Pr_LOS = d1/d2D + exp(-d2D/p1)*(1-d1/d2D);
                                end
                            elseif height<=300
                                Pr_LOS = 1;
                            else
                                error('The height is over limit of TR 36.777');
                            end
                        else
                            error('The height is not support in this scenario, change to case"AV"');
                        end
                    elseif d2D <= 18
                        Pr_LOS = 1;
                    else
                        C_tmp = (max((height-13)/10, 0)).^1.5;
                        Pr_LOS = (18./d2D + exp(-d2D/63).*(1-18./d2D)) .* (1+C_tmp*(5/4)*((d2D/100).^3).*exp(-d2D/150));
                    end
            case 'UMi'
                Pr_LOS = min(18./d2D,1).*(1-exp(-d2D./36))+exp(-d2D./36);
            case'RMa'
                Pr_LOS = min(exp(-(d2D-10)./1000),1);
            case'InH'
                    switch obj.scenario.subname
                        case 'A'
                            Pr_LOS = min(max(exp(-(d_2Din-18)./27),0.5),1);
                        case {'B','open_office'}
                            Pr_LOS = min(max(exp(-(d_2Din-5)./70.8),0.54*exp(-(d_2Din-49)./211.7)),1);
                        case 'mixed_office'
                            Pr_LOS = min(max(exp(-(d_2Din-1.2)/4.7),0.32*exp(-(d_2Din-6.5)/32.6)),1);
                        otherwise
                            error(['Case "' obj.scenario.subname '" for 3D-InH is not valid!']);
                    end
                case'InF'
                    switch obj.scenario.subname
                        case {'SL','DL'}
                            k_subsce = -obj.scenario.d_cluster/log(1-obj.scenario.r);
                        case {'SH','DH'}
                            k_subsce = -obj.scenario.d_cluster/log(1-obj.scenario.r)*((obj.sector(txrx).height-height)/(obj.scenario.hc - height));
                        case 'HH'
                            k_subsce = inf;
                        otherwise
                            error(['Case "' obj.scenario.subname '" for InF is not valid! ']);
                    end
                    if d_2Din < obj.scenario.d_subsce
                        Pr_LOS = obj.scenario.p_subsce;
                    else
                        Pr_LOS = obj.scenario.p_subsce*exp(-d_2Din/k_subsce);
                    end
                case 'UrbanGrid'
                    if obj.fc < 30e9
                        % TEMP UrbanGrid calibration workaround: FR1 follows TR 37.885 Urban V2X LOS probability.
                        % TR 37.885 itself does not define a frequency-dependent LOS split.
                        Pr_LOS = min(1, 1.05 * exp(-0.0114 * d2D));
                    else
                        % TEMP UrbanGrid calibration workaround: keep previous formula for FR2 because 30GHz benchmark matched.
                        % Original active formula before FR1 calibration:
                        if d2D <= 18
                            Pr_LOS = 1;
                        else
                            Pr_LOS = 18/d2D + exp(-d2D/63)*(1 - 18/d2D);
                        end
                    end
            end
            randlos = (txrx==1)* obj.ST.rand_LoS(obj.sector(1).ID) + (txrx==2)* obj.ST.rand_LoS(obj.sector(obj.TXRX).ID);
            obj.bLOS(txrx) = (randlos < Pr_LOS);
            if obj.static == 1
                obj.bLOS(2) = obj.bLOS(1);
                break;
            end
            end
        end

        % caculator the path loss distance depend
        function pathloss(obj)
            tx_height = obj.sector.height;
            rx_height = obj.ST.height;
            fc_GHz = obj.fc/1e9;    % to GHz
            shadow_sigma_dB = zeros(obj.TXRX,1);
            pathloss_dB = zeros(obj.TXRX,1);
            switch obj.scenario.name
                case 'UMa'
                for txrx=1:2
                    if rx_height <= 22.5
                        rx_breakpoint_height = rx_height;
                    else
                        rx_breakpoint_height = 1.5;
                    end

                    if rx_breakpoint_height <= 13
                        C = 0;
                    elseif (rx_breakpoint_height > 13) && (rx_breakpoint_height <= 23)
                        if obj.d_2D(txrx) <= 18
                            g = 0;
                        else
                            g = 5/4*((obj.d_2D(txrx)/100).^3)*exp(-obj.d_2D(txrx)/150);
                        end
                        C = (((rx_breakpoint_height-13)/10).^1.5)*g;
                    end
                    
                    if rand < 1/(1+C)
                        environment_height = 1;
                    else
                        environment_height = randsrc(1,1,[12,15,18,21]);
                    end
                    effective_tx_height = tx_height-environment_height;
                    effective_rx_height = rx_breakpoint_height-environment_height;
                    breakpoint_distance = 4*effective_tx_height*effective_rx_height*obj.fc/3e8;

                    if rx_height<=22.5 && rx_height>=1.5
                        if (obj.d_2D(txrx) >= 10) && (obj.d_2D(txrx) <= breakpoint_distance)
                            los_pathloss_dB = 22*log10(obj.d_3D(txrx))+28+20*log10(fc_GHz);
                        elseif (obj.d_2D(txrx) > breakpoint_distance) && (obj.d_2D(txrx) < 5000)
                            los_pathloss_dB = 40*log10(obj.d_3D(txrx))+28+20*log10(fc_GHz)-9*log10(breakpoint_distance.^2+(tx_height-rx_height).^2);
                        else
                            los_pathloss_dB = 22*log10(obj.d_3D(txrx))+28+20*log10(fc_GHz);
                        end
                        if obj.bLOS(txrx)
                            pathloss_dB(txrx) = los_pathloss_dB;
                            shadow_sigma_dB(txrx) = 4;
                        else
                            nlos_pathloss_dB = 13.54 + 39.08*log10(obj.d_3D(txrx)) + 20*log10(fc_GHz)-0.6*(rx_height-1.5);
                            pathloss_dB(txrx) = max(los_pathloss_dB, nlos_pathloss_dB);
                            shadow_sigma_dB(txrx) = 6;
                            % optional
                            % nlos_pathloss_dB = 32.4+20*log10(fc_GHz)+30*log10(obj.d_3D);
                            % shadow_sigma_dB = 7.8;
                        end
                    elseif strcmp(obj.scenario.subname,"AV")
                        if rx_height<=300 && rx_height>22.5
                            if obj.bLOS(txrx)
                                los_pathloss_dB = 22*log10(obj.d_3D(txrx))+28+20*log10(fc_GHz);
                                pathloss_dB(txrx) = los_pathloss_dB;
                                shadow_sigma_dB(txrx) = 4.64*exp(-0.0066*rx_height);
                            else
                                nlos_pathloss_dB = -17.5 + (46-7*log10(rx_height))*log10(obj.d_3D(txrx)) + 20*log10(40*pi*fc_GHz/3);
                                pathloss_dB(txrx) = nlos_pathloss_dB;
                                shadow_sigma_dB(txrx) = 6;
                            end
                        else
                            error("unsupport UT height in this scenario");
                        end
                    else
                        error("unsupport UT height in this scenario");
                    end

                    if obj.O2I
                        [penetration_loss_dB,penetration_sigma_dB] = obj.O2IpenetrationLoss(fc_GHz,obj.ST.O2IPL);
                        shadow_sigma_dB = 7;
                        pathloss_dB = pathloss_dB + penetration_loss_dB + 0.5*obj.d_2D_in + penetration_sigma_dB*obj.ST.O2Isigma;
                    end
                    if obj.static == 1
                        pathloss_dB(2)=pathloss_dB(1);
                        shadow_sigma_dB(2)=shadow_sigma_dB(1);
                        break;
                    end

                end
            case 'UMi'
                effective_tx_height = tx_height-1;
                effective_rx_height = rx_height-1;
                breakpoint_distance = 4*effective_tx_height*effective_rx_height*obj.fc/3e8;
                for txrx = 1:2
                    if (obj.d_2D(txrx) <= breakpoint_distance)% && (obj.d_2D >= 10)
                        los_pathloss_dB = 32.4 + 21*log10(obj.d_3D(txrx)) + 20*log10(fc_GHz); % TR38.901
                    elseif  (obj.d_2D(txrx) > breakpoint_distance) && (obj.d_2D(txrx) < 5000)
                        % los_pathloss_dB = 40*log10(obj.d_3D)+28+20*log10(fc_GHz)-9*log10(breakpoint_distance.^2+(tx_height-rx_height).^2); TR36.873
                        los_pathloss_dB = 32.4 + 40*log10(obj.d_3D(txrx)) + 20*log10(fc_GHz) - 9.5*log10(breakpoint_distance.^2+(tx_height-rx_height).^2); % TR38.901
                    end
                    if obj.bLOS(txrx)
                        pathloss_dB(txrx) = los_pathloss_dB;
                        shadow_sigma_dB(txrx) = 4; % TR38.901
                    else
                        nlos_pathloss_dB = 35.3*log10(obj.d_3D(txrx))+22.4+21.3*log10(fc_GHz)-0.3*(rx_height-1.5); % TR38.901
                        pathloss_dB(txrx) = max(los_pathloss_dB, nlos_pathloss_dB);
                        shadow_sigma_dB(txrx) = 7.82; % TR38.901
                        % optional
                        % nlos_pathloss_dB = 32.4+20*log10(fc_GHz)+31.9*log10(obj.d_3D);
                        % shadow_sigma_dB = 8.2;
                    end
                    if obj.O2I
                        [penetration_loss_dB,penetration_sigma_dB] = obj.O2IpenetrationLoss(fc_GHz,obj.ST.O2IPL);
                        shadow_sigma_dB = 7;
                        pathloss_dB = pathloss_dB + penetration_loss_dB + 0.5*obj.d_2D_in + penetration_sigma_dB*obj.ST.O2Isigma;
                    end
                    if obj.static == 1
                        pathloss_dB(2) = pathloss_dB(1);
                        shadow_sigma_dB(2) = shadow_sigma_dB(1);
                        break;
                    end
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
                if obj.bLOS
                    pathloss_dB = los_pathloss_dB;
                else
                    pathloss_dB = 161.04-7.1*log10(W)+7.5*log10(h)-(24.37-3.7*(h/tx_height)^2)*log10(tx_height)+ ...
                        (43.42-3.1*log10(tx_height))*(log10(obj.d_3D)-3)+20*log10(fc_GHz)-(3.2*(log10(11.75*rx_height))^2-4.97);
                    shadow_sigma_dB = 8;
                end
                if obj.O2I
                    [penetration_loss_dB,penetration_sigma_dB] = obj.O2IpenetrationLoss(fc_GHz,obj.ST.O2IPL);
                    shadow_sigma_dB = 8;
                    pathloss_dB = pathloss_dB + penetration_loss_dB + 0.5*obj.d_2D_in + penetration_sigma_dB*obj.ST.O2Isigma;
                end
            case 'InH'
               for txrx = 1:2  
                   pathloss_dB(txrx) = 32.4 + 17.3*log10(obj.d_3D(txrx)) + 20*log10(fc_GHz);
                   shadow_sigma_dB(txrx) = 3;
                   if ~obj.bLOS(txrx)
                       nlos_pathloss_dB = 17.3 + 38.3*log10(obj.d_3D(txrx)) + 24.9*log10(fc_GHz);
                       pathloss_dB(txrx) = max(pathloss_dB(txrx),nlos_pathloss_dB);
                       shadow_sigma_dB(txrx) = 8.03;
                       % optional
                       % nlos_pathloss_dB = 32.4+20*log10(fc_GHz)+31.9*log10(obj.d_3D);
                       % shadow_sigma_dB = 8.29;
                   end
               end
            case 'InF'
                for txrx = 1:2
                    pathloss_dB(txrx) = 31.84+21.5*log10(obj.d_3D(txrx))+19*log10(fc_GHz);
                    %                 shadow_sigma_dB = 4;
                    shadow_sigma_dB(txrx) = 4.3;
                    if ~obj.bLOS(txrx)
                        switch obj.scenario.subname
                            case 'SL'
                                pathloss_dB(txrx) = max(pathloss_dB(txrx),33+25.5*log10(obj.d_3D(txrx))+20*log10(fc_GHz));
                                shadow_sigma_dB(txrx) = 5.7;
                            case 'DL'
                                pathloss_dB(txrx) = max(pathloss_dB(txrx),18.6+35.7*log10(obj.d_3D(txrx))+20*log10(fc_GHz));
                                shadow_sigma_dB(txrx) = 7.2;
                            case 'SH'
                                pathloss_dB(txrx) = max(pathloss_dB(txrx),32.4+23*log10(obj.d_3D(txrx))+20*log10(fc_GHz));
                                shadow_sigma_dB(txrx) = 5.9;
                            case 'DH'
                                pathloss_dB(txrx) = max(pathloss_dB(txrx),33.63+21.9*log10(obj.d_3D(txrx))+20*log10(fc_GHz));
                                shadow_sigma_dB(txrx) = 4;
                        end
                    end
                    if obj.O2I
                        [penetration_loss_dB,penetration_sigma_dB] = obj.O2IpenetrationLoss(fc_GHz,obj.ST.O2IPL);
                        
                        switch obj.scenario.subname
                            case 'SL'
                                shadow_sigma_dB(txrx) = 5.7;
                            case 'DL'
                                shadow_sigma_dB(txrx) = 7.2;
                            case 'SH'
                                shadow_sigma_dB(txrx) = 5.9;
                            case 'DH'
                                shadow_sigma_dB(txrx) = 4;
                        end
                        pathloss_dB = pathloss_dB + penetration_loss_dB + 0.5*obj.d_2D_in + penetration_sigma_dB*obj.ST.O2Isigma;
                    end
                    if  obj.static == 1
                        pathloss_dB(2) = pathloss_dB(1);
                        shadow_sigma_dB(2) = shadow_sigma_dB(1);
                        break;
                    end
                end
            case 'UrbanGrid'
                for txrx=1:2
                    rx_breakpoint_height = rx_height;
                    environment_height = 0.25;
                    
                    effective_tx_height = tx_height-environment_height;
                    effective_rx_height = rx_breakpoint_height-environment_height;
                    breakpoint_distance = 4*effective_tx_height*effective_rx_height*obj.fc/3e8;

                    if (obj.d_2D(txrx) >= 10) && (obj.d_2D(txrx) <= breakpoint_distance)

                        los_pathloss_dB = 22*log10(obj.d_3D(txrx))+28+20*log10(fc_GHz);
                    elseif (obj.d_2D(txrx) > breakpoint_distance) && (obj.d_2D(txrx) < 5000)
                        los_pathloss_dB = 40*log10(obj.d_3D(txrx))+28+20*log10(fc_GHz)-9*log10(breakpoint_distance.^2+(tx_height-rx_height).^2);
                    else
                        los_pathloss_dB = 22*log10(obj.d_3D(txrx))+28+20*log10(fc_GHz);
                    end
                    if obj.bLOS(txrx)
                        pathloss_dB(txrx) = los_pathloss_dB;
                        shadow_sigma_dB(txrx) = 4;
                    else
                        nlos_pathloss_dB = 13.54 + 39.08*log10(obj.d_3D(txrx)) + 20*log10(fc_GHz)-0.6*(rx_height-1.5);
                        pathloss_dB(txrx) = max(los_pathloss_dB, nlos_pathloss_dB);
                        shadow_sigma_dB(txrx) = 6;
                        % optional
                        % nlos_pathloss_dB = 32.4+20*log10(fc_GHz)+30*log10(obj.d_3D);
                        % shadow_sigma_dB = 7.8;
                    end
                    if obj.static == 1
                        pathloss_dB(2)=pathloss_dB(1);
                        shadow_sigma_dB(2)=shadow_sigma_dB(1);
                        break;
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
                material_loss_dB = tools.get_MaterialLoss(fc);
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
            fc_GHz(fc_GHz<=6) = 6;
            for txrx = 1:2
                switch obj.scenario.name
                    case {'UMa','UrbanGrid'}
                        fc_GHz(fc_GHz<=6) = 6;
                        is_AV = strcmp(obj.scenario.subname,'AV');
                        if is_AV
                            av_alternative = obj.scenario.alternative;
                        else
                            av_alternative = [];
                        end
                        if ~is_AV || av_alternative == 3
                            if obj.O2I
                                if obj.bLOS
                                    obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D_out/1000)-0.01*(obj.ST.height-1.5)+0.75);
                                    obj.sigma_lgZSD = 0.4;
                                    obj.mu_offset_ZOD = 0;
                                else
                                    obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D_out/1000)-0.01*(obj.ST.height-1.5)+0.9);
                                    obj.sigma_lgZSD = 0.49;
                                    % obj.mu_offset_ZOD = -10^(-0.62*log10(max(10, obj.d_2D_out))+1.93-0.07*(obj.UE.height-1.5));
                                    obj.mu_offset_ZOD = (7.66*log10(fc_GHz)-5.96) -10^((0.208*log10(fc_GHz)-0.782)*log10(max((25), obj.d_2D_out))+(-0.13*log10(fc_GHz)+2.03)-0.07*(obj.ST.height-1.5));
                                end
                                raw_rand_vec = obj.drawLspRawRandn(6, txrx);
                                corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                                corr_rand_vec = obj.scenario.Cross_correlation_O2I*corr_rand_vec;
                                % SF/DS/ASD/ASA/ZSD/ZSA
                                lsp_std_vec = [obj.sigma_SF, 0.32, 0.7, 0.16, obj.sigma_lgZSD, 0.43]';
                                lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                                obj.SF = lsp_offset_vec(1);
                                obj.DS = 10^(lsp_offset_vec(2)-6.62);
                                obj.ASD = 10^(lsp_offset_vec(3)+0.58);
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
                                obj.c_ASD = 1.8;
                            else
                                if obj.bLOS(txrx) % 0 1
                                    obj.mu_lgZSD(txrx) = max(-0.5, -2.1*(obj.d_2D(txrx)/1000)-0.01*(obj.ST.height-1.5)+0.75);
                                    obj.sigma_lgZSD(txrx) = 0.4;
                                    obj.mu_offset_ZOD(txrx) = 0;
                                    raw_rand_vec = obj.drawLspRawRandn(7, txrx);
                                    corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                                    corr_rand_vec = obj.scenario.Cross_correlation_LOS*corr_rand_vec;
                                    % SF/K/DS/ASD/ASA/ZSD/ZSA
                                    % lsp_std_vec = [obj.sigma_SF, 5, 0.39, 0.43, 0.19, obj.sigma_lgZSD, 0.16]';
                                    lsp_std_vec = [obj.sigma_SF(txrx), 3.5, 0.57+0.026*log10(fc_GHz), 0.31, 0.19, obj.sigma_lgZSD(txrx), 0.15]';
                                    lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                                    obj.SF(txrx) = lsp_offset_vec(1);
                                    obj.K(txrx) = lsp_offset_vec(2)+9;
                                    if is_AV && av_alternative == 3
                                        obj.K(txrx) = 15;
                                    end
                                    obj.DS(txrx) = 10.^(lsp_offset_vec(3)+(-7.067-0.0794*log10(1+fc_GHz)));
                                    obj.ASD(txrx) = 10.^(lsp_offset_vec(4)+0.92);
                                    obj.ASA(txrx) = 10.^(lsp_offset_vec(5)+1.76);
                                    obj.ZSD(txrx) = 10.^(lsp_offset_vec(6)+obj.mu_lgZSD(txrx));
                                    obj.ZSA(txrx) = 10.^(lsp_offset_vec(7)+0.96);
                                    obj.r_tau(txrx) = 2.5;
                                    obj.mu_XPR(txrx) = 8;
                                    obj.sigma_XPR(txrx) = 4;
                                    obj.N(txrx) = 12;
                                    obj.M(txrx) = 20;
                                    obj.zeta(txrx) = 3;
                                    obj.c_DS(txrx)  = max(0.25, 6.5622 -3.4084*log10(fc_GHz));
                                    obj.c_ASA(txrx) = 11;
                                    obj.c_ZSA(txrx) = 7;
                                    obj.c_ASD(txrx) = 3.58;
                                else % NLOS
                                    obj.mu_lgZSD(txrx) = max(-0.5, -2.1*(obj.d_2D(txrx)/1000)-0.01*(obj.ST.height-1.5)+0.9);
                                    obj.sigma_lgZSD(txrx) = 0.49;
                                    obj.mu_offset_ZOD(txrx) = (7.66*log10(fc_GHz) - 5.96) ...
                                        - 10^((0.208*log10(fc_GHz) - 0.782) * log10(max(25, obj.d_2D(txrx)))) ...
                                        + (-0.13*log10(fc_GHz) + 2.03) ...
                                        - 0.07*(obj.ST.height - 1.5);
                                    raw_rand_vec = obj.drawLspRawRandn(6, txrx);
                                    corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                                    corr_rand_vec = obj.scenario.Cross_correlation_NLOS*corr_rand_vec;
                                    % SF/DS/ASD/ASA/ZSD/ZSA
                                    lsp_std_vec = [obj.sigma_SF(txrx), 0.39, 0.44, (0.17-0.03*log10(fc_GHz)), obj.sigma_lgZSD(txrx), 0.17]';
                                    lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                                    obj.SF(txrx) = lsp_offset_vec(1);
                                    obj.DS(txrx) = 10.^(lsp_offset_vec(2)+(-6.47-0.134*log10(fc_GHz)));
                                    obj.ASD(txrx) = 10.^(lsp_offset_vec(3)+1.09);
                                    obj.ASA(txrx) = 10.^(lsp_offset_vec(4)+(2.04-0.25*log10(fc_GHz)));
                                    obj.ZSD(txrx) = 10.^(lsp_offset_vec(5)+obj.mu_lgZSD(txrx));
                                    obj.ZSA(txrx) = 10.^(lsp_offset_vec(6)+(-0.2856*log10(fc_GHz)+1.445));
                                    obj.r_tau(txrx) = 2.3;
                                    obj.mu_XPR(txrx) = 7;
                                    obj.sigma_XPR(txrx) = 3;
                                    obj.N(txrx) = 20;
                                    obj.M(txrx) = 20;
                                    obj.zeta(txrx) = 3;
                                    obj.c_DS(txrx)  = max(0.25, 6.5622 -3.4084*log10(fc_GHz));
                                    obj.c_ASA(txrx) = 15;
                                    obj.c_ZSA(txrx) = 7;
                                    obj.c_ASD(txrx) = 1.8;
                                end
                            end
                        else
                            if av_alternative == 1
                                error('AV alternative 1 is not implemented yet.');
                            elseif av_alternative ~= 2
                                error('Unsupported AV alternative: %g.', av_alternative);
                            end
                            if obj.bLOS(txrx) % 0 1
                                obj.mu_lgZSD(txrx) = max(-0.5, -2.1*(obj.d_2D(txrx)/1000)-0.01*(obj.ST.height-1.5)+0.75);
                                obj.sigma_lgZSD(txrx) = 0.4;
                                obj.mu_offset_ZOD(txrx) = 0;
                                raw_rand_vec = obj.drawLspRawRandn(7, txrx);
                                corr_rand_vec = reshape(raw_rand_vec,[numel(raw_rand_vec), 1]);
                                corr_rand_vec = obj.scenario.Cross_correlation_LOS*corr_rand_vec;
                                % SF/K/DS/ASD/ASA/ZSD/ZSA
                                % lsp_std_vec = [obj.sigma_SF, 5, 0.39, 0.43, 0.19, obj.sigma_lgZSD, 0.16]';
                                K_std = 8.158 * exp(0.0046 *obj.ST.height );
                                DS_std = 0.7294 * exp(0.0014 *obj.ST.height );
                                ASD_std = 1.0188 * exp(-0.0001 *obj.ST.height );
                                ASA_std = 1.0389 * exp(0.0085 *obj.ST.height );                                
                                ZSD_std = 1.0757* exp(0.0059 *obj.ST.height );
                                ZSA_std = 0.9576 * exp(-0.0018 *obj.ST.height );
                                lsp_std_vec = [obj.sigma_SF(txrx), K_std, DS_std, ASD_std, ASA_std, ZSD_std, ZSA_std]';
                                lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                                obj.SF(txrx) = lsp_offset_vec(1);
                                K_mu = 4.217*log10(obj.ST.height)+ 5.787;
                                obj.K(txrx) = lsp_offset_vec(2)+K_mu;
                                DS_mu = -0.31*log10(obj.ST.height)+ -6.845;
                                obj.DS(txrx) = 10.^(lsp_offset_vec(3)+DS_mu);
                                ASD_mu = -0.0135*log10(obj.ST.height)+ 1.345;
                                obj.ASD(txrx) = 10.^(lsp_offset_vec(4)+ASD_mu);
                                ASA_mu = -2.4985*log10(obj.ST.height)-1.602;
                                obj.ASA(txrx) = 10.^(lsp_offset_vec(5)+ASA_mu);
                                ZSD_mu = -0.2975*log10(obj.ST.height)-0.5798;
                                obj.ZSD(txrx) = 10.^(lsp_offset_vec(6)+ZSD_mu);
                                ZSA_mu = -0.2895*log10(obj.ST.height)+ 0.225;
                                obj.ZSA(txrx) = 10.^(lsp_offset_vec(7)+ZSA_mu);

                                obj.r_tau(txrx) = 2.5;
                                obj.mu_XPR(txrx) = 8;
                                obj.sigma_XPR(txrx) = 4;
                                obj.N(txrx) = 12;
                                obj.M(txrx) = 20;
                                obj.zeta(txrx) = 3;
                                obj.c_DS(txrx)  = max(0.25, 6.5622 -3.4084*log10(fc_GHz));
                                obj.c_ASA(txrx) = 11;
                                obj.c_ZSA(txrx) = 7;
                                obj.c_ASD(txrx) = 3.58;
                            else % NLOS
                                obj.mu_lgZSD(txrx) = max(-0.5, -2.1*(obj.d_2D(txrx)/1000)-0.01*(obj.ST.height-1.5)+0.9);
                                obj.sigma_lgZSD(txrx) = 0.49;
                                obj.mu_offset_ZOD(txrx) = (7.66*log10(fc_GHz) - 5.96) ...
                                    - 10^((0.208*log10(fc_GHz) - 0.782) * log10(max(25, obj.d_2D(txrx)))) ...
                                    + (-0.13*log10(fc_GHz) + 2.03) ...
                                    - 0.07*(obj.ST.height - 1.5);
                                raw_rand_vec = obj.drawLspRawRandn(6, txrx);
                                corr_rand_vec = obj.scenario.Cross_correlation_NLOS*raw_rand_vec;
                                % SF/DS/ASD/ASA/ZSD/ZSA
                                DS_std = 0.9745 * exp(-0.0045 *obj.ST.height );
                                ASD_std = 1.2387 * exp(-0.0046 *obj.ST.height );
                                ASA_std = 1.022  * exp(0.009944 *obj.ST.height );                                
                                ZSD_std = 1.6421* exp(-0.0092 *obj.ST.height );
                                ZSA_std = 1.6237 * exp(-0.0076 *obj.ST.height );
                                lsp_std_vec = [obj.sigma_SF(txrx), DS_std, ASD_std, ASA_std, ZSD_std, ZSA_std]';

                                lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                                obj.SF(txrx) = lsp_offset_vec(1);
                                DS_mu = 0.0965*log10(obj.ST.height)+ -7.503;
                                obj.DS(txrx) = 10.^(lsp_offset_vec(2)+DS_mu);
                                ASD_mu = 1.17*log10(obj.ST.height)-0.665;
                                obj.ASD(txrx) = 10.^(lsp_offset_vec(3)+ASD_mu);
                                ASA_mu = -2.266*log10(obj.ST.height)- 2.666;
                                obj.ASA(txrx) = 10.^(lsp_offset_vec(4)+ASA_mu);
                                ZSD_mu = 0.925*log10(obj.ST.height)-2.725;
                                obj.ZSD(txrx) = 10.^(lsp_offset_vec(5)+ZSD_mu);
                                ZSA_mu = -0.0005*log10(obj.ST.height)- 0.4695;
                                obj.ZSA(txrx) = 10.^(lsp_offset_vec(6)+ZSA_mu);

                                obj.r_tau(txrx) = 2.3;
                                obj.mu_XPR(txrx) = 7;
                                obj.sigma_XPR(txrx) = 3;
                                obj.N(txrx) = 20;
                                obj.M(txrx) = 20;
                                obj.zeta(txrx) = 3;
                                obj.c_DS(txrx)  = max(0.25, 6.5622 -3.4084*log10(fc_GHz));
                                obj.c_ASA(txrx) = 15;
                                obj.c_ZSA(txrx) = 7;
                                obj.c_ASD(txrx) = 1.8;
                            end

                        end
                    case 'UMi'
                        fc_GHz(fc_GHz<=2) = 2;
                        if obj.O2I %
                            if obj.bLOS
                                % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D_out/1000)+0.01*abs(obj.UE.height-obj.sector.attached_BS.height)+0.75);
                                obj.mu_lgZSD = max(-0.21, -14.8*(obj.d_2D_out/1000)+0.01*abs(obj.ST.height-obj.sector.height)+0.83);
                                obj.sigma_lgZSD = 0.35; % 0.4
                                obj.mu_offset_ZOD = 0;
                            else
                                % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D_out/1000)+0.01*max(obj.UE.height-obj.sector.attached_BS.height, 0)+0.9);
                                obj.mu_lgZSD = max(-0.5, -3.1*(obj.d_2D_out/1000)+0.01*max(obj.ST.height-obj.sector.height, 0)+0.2);
                                obj.sigma_lgZSD = 0.35; %0.6
                                % obj.mu_offset_ZOD = -10^(-0.55*log10(max(10, obj.d_2D_out))+1.6);
                                obj.mu_offset_ZOD = -10^(-1.5*log10(max(10, obj.d_2D_out))+3.3);
                            end
                            raw_rand_vec = obj.drawLspRawRandn(6, txrx);
                            corr_rand_vec = obj.scenario.Cross_correlation_O2I*raw_rand_vec;
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
                        else
                            if obj.bLOS(txrx) % 0 1
                                % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)+0.01*abs(obj.UE.height-obj.sector.attached_BS.height)+0.75);
                                obj.mu_lgZSD(txrx) = max(-0.21, -14.8*(obj.d_2D(txrx)/1000)+0.01*abs(obj.ST.height-obj.sector(txrx).height)+0.83);
                                obj.sigma_lgZSD(txrx) = 0.35; % 0.4
                                obj.mu_offset_ZOD(txrx) = 0;
                                raw_rand_vec = obj.drawLspRawRandn(7, txrx);
                                corr_rand_vec = obj.scenario.Cross_correlation_LOS*raw_rand_vec;
                                % SF/K/DS/ASD/ASA/ZSD/ZSA
                                % lsp_std_vec = [obj.sigma_SF, 5, 0.4, 0.43, 0.19, obj.sigma_lgZSD, 0.16]';
                                lsp_std_vec = [obj.sigma_SF(txrx), 5, 0.38, 0.41, (0.014*log10(1+fc_GHz)+0.28), obj.sigma_lgZSD(txrx), (-0.04*log10(1+fc_GHz)+0.34)]';
                                lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                                obj.SF(txrx) = lsp_offset_vec(1);
                                obj.K(txrx) = lsp_offset_vec(2)+9;
                                obj.DS(txrx) = 10.^(lsp_offset_vec(3)+(-0.24*log10(1+fc_GHz)-7.14));
                                obj.ASD(txrx) = 10.^(lsp_offset_vec(4)+(-0.05*log10(1+fc_GHz)+1.21));
                                obj.ASA(txrx) = 10.^(lsp_offset_vec(5)+(-0.08*log10(1+fc_GHz)+1.73));
                                obj.ZSD(txrx) = 10.^(lsp_offset_vec(6)+obj.mu_lgZSD(txrx));
                                obj.ZSA(txrx) = 10.^(lsp_offset_vec(7)+(-0.1*log10(1+fc_GHz)+0.73));
                                obj.r_tau(txrx) = 3;
                                obj.mu_XPR(txrx) = 9;
                                obj.sigma_XPR(txrx) = 3;
                                obj.N(txrx) = 12;
                                obj.M(txrx) = 20;
                                obj.zeta(txrx) = 3;
                                obj.c_DS(txrx)  = 5;
                                obj.c_ASA(txrx) = 17;
                                obj.c_ZSA(txrx) = 7;
                                obj.c_ASD(txrx) = 3;
                            else % NLOS
                                % obj.mu_lgZSD = max(-0.5, -2.1*(obj.d_2D/1000)+0.01*max(obj.UE.height-obj.sector.attached_BS.height, 0)+0.9);
                                obj.mu_lgZSD(txrx) = max(-0.5, -3.1*(obj.d_2D(txrx)/1000)+0.01*max(obj.ST.height-obj.sector(txrx).height, 0)+0.2);
                                obj.sigma_lgZSD(txrx) = 0.35;
                                % obj.mu_offset_ZOD = -10^(-0.55*log10(max(10, obj.d_2D))+1.6);
                                obj.mu_offset_ZOD(txrx) = -10.^(-1.5*log10(max(10, obj.d_2D(txrx)))+3.3);
                                raw_rand_vec = obj.drawLspRawRandn(6, txrx);
                                corr_rand_vec = obj.scenario.Cross_correlation_NLOS*raw_rand_vec;
                                % SF/DS/ASD/ASA/ZSD/ZSA
                                lsp_std_vec = [obj.sigma_SF(txrx), (0.16*log10(1+fc_GHz)+0.28), (0.11*log10(1+fc_GHz)+0.33), (0.05*log10(1+fc_GHz)+0.3), obj.sigma_lgZSD(txrx), (-0.07*log10(1+fc_GHz)+0.41)]';
                                lsp_offset_vec = lsp_std_vec.*corr_rand_vec;
                                obj.SF(txrx) = lsp_offset_vec(1);
                                obj.DS(txrx) = 10.^(lsp_offset_vec(2)+(-0.24*log10(1+fc_GHz)-6.83));
                                obj.ASD(txrx) = 10.^(lsp_offset_vec(3)+(-0.23*log10(1+fc_GHz)+1.53));
                                obj.ASA(txrx) = 10.^(lsp_offset_vec(4)+(-0.08*log10(1+fc_GHz)+1.81));
                                obj.ZSD(txrx) = 10.^(lsp_offset_vec(5)+obj.mu_lgZSD(txrx));
                                obj.ZSA(txrx) = 10.^(lsp_offset_vec(6)+(-0.04*log10(1+fc_GHz)+0.92));
                                obj.r_tau(txrx) = 2.1;
                                obj.mu_XPR(txrx) = 8;
                                obj.sigma_XPR(txrx) = 3;
                                obj.N(txrx) = 19;
                                obj.M(txrx) = 20;
                                obj.zeta(txrx) = 3;
                                obj.c_DS(txrx)  = 11;
                                obj.c_ASA(txrx) = 22;
                                obj.c_ZSA(txrx) = 7;
                                obj.c_ASD(txrx) = 10;
                            end

                        end

                    case 'RMa'
                        if obj.O2I
                            obj.mu_lgZSD = max(-1, -0.19*(obj.d_2D_out/1000)-0.01*(obj.ST.height-1.5)+0.28);
                            obj.sigma_lgZSD = 0.30;
                            obj.mu_offset_ZOD = atand((35-3.5)/obj.d_2D_out)-atand((35-1.5)/obj.d_2D_out);
                            raw_rand_vec = obj.drawLspRawRandn(6, txrx);
                            corr_rand_vec = obj.scenario.Cross_correlation_O2I*raw_rand_vec;
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
                        elseif obj.bLOS
                            obj.mu_lgZSD = max(-1, -0.17*(obj.d_2D/1000)-0.01*(obj.ST.height-1.5)+0.22);
                            obj.sigma_lgZSD = 0.34;
                            obj.mu_offset_ZOD = 0;
                            raw_rand_vec = obj.drawLspRawRandn(7, txrx);
                            corr_rand_vec = obj.scenario.Cross_correlation_LOS*raw_rand_vec;
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
                            obj.mu_lgZSD = max(-1, -0.19*(obj.d_2D/1000)-0.01*(obj.ST.height-1.5)+0.28);
                            obj.sigma_lgZSD = 0.30;
                            obj.mu_offset_ZOD = atand((35-3.5)/obj.d_2D)-atand((35-1.5)/obj.d_2D);
                            raw_rand_vec = obj.drawLspRawRandn(6, txrx);
                            corr_rand_vec = obj.scenario.Cross_correlation_NLOS*raw_rand_vec;
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

                        if obj.bLOS(txrx)
                            % obj.mu_lgZSD = 1.02;
                            % obj.sigma_lgZSD = 0.41;
                            obj.mu_lgZSD(txrx) = -1.43*log10(1+fc_GHz)+2.228;
                            obj.sigma_lgZSD(txrx) = 0.13*log10(1+fc_GHz)+0.3;
                            obj.mu_offset_ZOD(txrx) = 0;
                                raw_rand_vec = obj.drawLspRawRandn(7, txrx);
                            corr_rand_vec = obj.scenario.Cross_correlation_LOS*raw_rand_vec;
                            obj.SF(txrx) = obj.sigma_SF(txrx)*corr_rand_vec(1);
                            obj.K(txrx) = 4*corr_rand_vec(2)+7;
                            obj.DS(txrx) = 10^(0.18*corr_rand_vec(3)+ (-0.01*log10(1+fc_GHz)-7.692));
                            obj.ASD(txrx) = 10^(0.18*corr_rand_vec(4)+1.60);
                            obj.ASA(txrx) = 10^((0.12*log10(1+fc_GHz)+0.119)*corr_rand_vec(5)+ (-0.19*log10(1+fc_GHz)+1.781));
                            obj.ZSD(txrx) = 10^(obj.sigma_lgZSD(txrx)*corr_rand_vec(6)+obj.mu_lgZSD(txrx));
                            obj.ZSA(txrx) = 10^((-0.04*log10(1+fc_GHz)+0.264)*corr_rand_vec(7)+ (-0.26*log10(1+fc_GHz)+1.44));
                            obj.r_tau(txrx) = 3.6;
                            obj.mu_XPR(txrx) = 11;
                            obj.sigma_XPR(txrx) = 4; % 3
                            obj.N(txrx) = 15;
                            obj.M(txrx) = 20;
                            obj.zeta(txrx) = 6;
                            obj.c_DS(txrx) = 0;
                            obj.c_ASA(txrx) = 8;
                            obj.c_ZSA(txrx) = 9;
                            obj.c_ASD(txrx) = 5;
                        else
                            obj.mu_lgZSD(txrx) = 1.08;
                            obj.sigma_lgZSD(txrx) = 0.36;
                            obj.mu_offset_ZOD(txrx) = 0;
                                raw_rand_vec = obj.drawLspRawRandn(6, txrx);
                            corr_rand_vec = obj.scenario.Cross_correlation_NLOS*raw_rand_vec;
                            obj.SF(txrx) = obj.sigma_SF(txrx)*corr_rand_vec(1);
                            obj.DS(txrx) = 10^((0.1*log10(1+fc_GHz)+0.055)*corr_rand_vec(2)+ (-0.28*log10(1+fc_GHz)-7.173));
                            obj.ASD(txrx) = 10^(0.25*corr_rand_vec(3)+1.62);
                            obj.ASA(txrx) = 10^((0.12*log10(1+fc_GHz)+0.059)*corr_rand_vec(4)+ (-0.11*log10(1+fc_GHz)+1.863));
                            obj.ZSD(txrx) = 10^(obj.sigma_lgZSD(txrx)*corr_rand_vec(5)+obj.mu_lgZSD(txrx));
                            obj.ZSA(txrx) = 10^((-0.09*log10(1+fc_GHz)+0.746)*corr_rand_vec(6)+ (-0.15*log10(1+fc_GHz)+1.387));
                            obj.r_tau(txrx) = 3;
                            obj.mu_XPR(txrx) = 10;
                            obj.sigma_XPR(txrx) = 4; % 3
                            obj.N(txrx) = 19;
                            obj.M(txrx) = 20;
                            obj.zeta(txrx) = 3;
                            obj.c_DS(txrx) = 0;
                            obj.c_ASA(txrx) = 11;
                            obj.c_ZSA(txrx) = 9;
                            obj.c_ASD(txrx) = 5;
                        end
                    case 'InF'
                        % V = hall volume in m^3, S = total surface area of hall in m^2 (walls+floor+ceiling)

                        if obj.bLOS(txrx)
                            obj.mu_lgZSD(txrx) = 1.35;
                            obj.sigma_lgZSD(txrx) = 0.35;
                            obj.mu_offset_ZOD(txrx) = 0;
                                raw_rand_vec = obj.drawLspRawRandn(7, txrx);
                            corr_rand_vec = obj.scenario.Cross_correlation_LOS*raw_rand_vec;
                            obj.SF(txrx) = obj.sigma_SF(txrx)*corr_rand_vec(1);
                            obj.K(txrx) = 8*corr_rand_vec(2)+7;
                            obj.DS(txrx) = 10^(0.15*corr_rand_vec(3)+ (log10(26*(obj.scenario.V/obj.scenario.S)+14)-9.35));
                            obj.ASD(txrx) = 10^(0.25*corr_rand_vec(4)+1.56);
                            obj.ASA(txrx) = 10^((0.12*log10(1+fc_GHz)+0.2)*corr_rand_vec(5)+ (-0.18*log10(1+fc_GHz)+1.78));
                            obj.ZSD(txrx) = 10^(obj.sigma_lgZSD*corr_rand_vec(6)+obj.mu_lgZSD);
                            obj.ZSA(txrx) = 10^(0.35*corr_rand_vec(7)+ (-0.2*log10(1+fc_GHz)+1.5));
                            obj.r_tau(txrx) = 2.7;
                            obj.mu_XPR(txrx) = 12;
                            obj.sigma_XPR(txrx) = 6;
                            obj.N(txrx) = 25;
                            obj.M(txrx) = 20;
                            obj.zeta(txrx) = 4;
                            obj.c_ASA(txrx) = 8;
                            obj.c_ZSA(txrx) = 9;
                            obj.c_ASD(txrx) = 5;
                        else
                            obj.mu_lgZSD(txrx) = 1.2;
                            obj.sigma_lgZSD(txrx) = 0.55;
                            obj.mu_offset_ZOD(txrx) = 0;
                                raw_rand_vec = obj.drawLspRawRandn(6, txrx);
                            corr_rand_vec = obj.scenario.Cross_correlation_NLOS*raw_rand_vec;
                            obj.SF(txrx) = obj.sigma_SF(txrx)*corr_rand_vec(1);
                            obj.DS(txrx) = 10^(0.19*corr_rand_vec(2)+ (log10(30*(obj.scenario.V/obj.scenario.S)+32)-9.44));
                            obj.ASD(txrx) = 10^(0.2*corr_rand_vec(3)+1.57);
                            obj.ASA(txrx) = 10^(0.3*corr_rand_vec(4)+1.72);
                            obj.ZSD(txrx) = 10^(obj.sigma_lgZSD(txrx)*corr_rand_vec(5)+obj.mu_lgZSD(txrx));
                            obj.ZSA(txrx) = 10^(0.45*corr_rand_vec(6)+ (-0.13*log10(1+fc_GHz)+1.45));
                            obj.r_tau(txrx) = 3;
                            obj.mu_XPR(txrx) = 11;
                            obj.sigma_XPR(txrx) = 6;
                            obj.N(txrx) = 25;
                            obj.M(txrx) = 20;
                            obj.zeta(txrx) = 3;
                            obj.c_ASA(txrx) = 8;
                            obj.c_ZSA(txrx) = 9;
                            obj.c_ASD(txrx) = 5;
                        end
                end
                obj.ASA(txrx) = min(obj.ASA(txrx), 104);
                obj.ASD(txrx) = min(obj.ASD(txrx), 104);
                obj.ZSA(txrx) = min(obj.ZSA(txrx), 52);
                obj.ZSD(txrx) = min(obj.ZSD(txrx), 52);
                if obj.static ==1
                    obj.mu_lgZSD(2) = obj.mu_lgZSD(1);
                    obj.sigma_lgZSD(2) = obj.sigma_lgZSD(1);
                    obj.mu_offset_ZOD(2) = obj.mu_offset_ZOD(1);
                    obj.SF(2) = obj.SF(1);
                    if obj.bLOS(1)
                        obj.K(2) = obj.K(1);
                    end
                    obj.DS(2) = obj.DS(1);
                    obj.ASD(2) = obj.ASD(1);
                    obj.ASA(2) = obj.ASA(1);
                    obj.ZSD(2) = obj.ZSD(1);
                    obj.ZSA(2) = obj.ZSA(1);
                    obj.r_tau(2) = obj.r_tau(1);
                    obj.mu_XPR(2) = obj.mu_XPR(1);
                    obj.sigma_XPR(2) = obj.sigma_XPR(1);
                    obj.N(2) = obj.N(1);
                    obj.M(2) = obj.M(1);
                    obj.zeta(2) = obj.zeta(1);
                    obj.c_ASA(2) = obj.c_ASA(1);
                    obj.c_ZSA(2) = obj.c_ZSA(1);
                    obj.c_ASD(2) = obj.c_ASD(1);                    
                    break;
                end
            end
        end

        % caculator the cluster delay parameters
        function cluster_delay(obj)
            max_N = max(obj.N);
            obj.tau_n = nan(2, max_N);
            obj.tau_order = nan(2, max_N);
            obj.tau_n_LOS = nan(2, max_N);
            for txrx = 1:2
                N_txrx = obj.N(txrx);
                    random_draw = rand(N_txrx,1);
                % end
                tau_n_tmp = -obj.r_tau(txrx)*obj.DS(txrx)*log(random_draw.');
                [tau_n_sorted, tau_order_sorted] = sort(tau_n_tmp - min(tau_n_tmp));
                obj.tau_n(txrx,1:N_txrx) = tau_n_sorted;
                obj.tau_order(txrx,1:N_txrx) = tau_order_sorted;
                obj.tau_n_LOS(txrx,1:N_txrx) = tau_n_sorted;
                if obj.bLOS(txrx) && ~obj.O2I
                    C_tau = 0.7705-0.0433*obj.K(txrx)+0.0002*obj.K(txrx)^2+0.000017*obj.K(txrx)^3;
                    obj.tau_n_LOS(txrx,1:N_txrx) = obj.tau_n(txrx,1:N_txrx)/C_tau;  % not to be used in cluster power generation
                end

                if obj.static == 1
                    obj.tau_n(2,:) = obj.tau_n(1,:);
                    obj.tau_order(2,:) = obj.tau_order(1,:);
                    obj.tau_n_LOS(2,:) = obj.tau_n_LOS(1,:);
                    break;
                end
            end
        end

        %% step 6: Generate cluster powers
        function cluster_power(obj)
            max_N = max(obj.N);
            max_M = max(obj.M);
            obj.Pn = nan(2, max_N);
            obj.Pn_LOS = nan(2, max_N);
            obj.cluster_shadow_dB = nan(2, max_N);
            obj.keep = false(2, max_N);
            obj.tau_n_keep = nan(2, max_N);
            obj.tau_n_LOS_keep = nan(2, max_N);
            obj.Pn_keep = nan(2, max_N);
            obj.Pn_LOS_keep = nan(2, max_N);
            obj.N_new = zeros(1, 2);
            obj.strong_cluster_id = nan(2, 2);
            obj.map_delay = zeros(2, max_M);
            for txrx = 1:2
                N_txrx = obj.N(txrx);
                    shadow_dB = randn(N_txrx,1)*obj.zeta(txrx);
                    obj.cluster_shadow_dB(txrx,1:N_txrx) = shadow_dB.';
                % end
                cluster_power_raw = exp(-obj.tau_n(txrx,1:N_txrx).*(obj.r_tau(txrx)-1)./obj.r_tau(txrx)./obj.DS(txrx)).*(10.^(-shadow_dB.'/10));
                obj.Pn(txrx,1:N_txrx) = cluster_power_raw/sum(cluster_power_raw);
                obj.Pn_LOS(txrx,1:N_txrx) = obj.Pn(txrx,1:N_txrx);
                if obj.bLOS(txrx) && ~obj.O2I
                    rice_linear = 10^(obj.K(txrx)/10);
                    P1_LOS = rice_linear/(rice_linear+1);
                    obj.Pn_LOS(txrx,1:N_txrx) = (1/(rice_linear+1)).*(cluster_power_raw./sum(cluster_power_raw));
                    obj.Pn_LOS(txrx,1) = obj.Pn_LOS(txrx,1) + P1_LOS;  % not be used in equation step 11(7.3-22)
                end
                    keep_tmp       = ((10*log10(obj.Pn(txrx,1:N_txrx)/max(obj.Pn(txrx,1:N_txrx)))) >= -25);
                % end

                obj.keep(txrx,1:N_txrx)  = keep_tmp;
                tau_absolute_valid = [];
                if ~isempty(obj.tau_absolute)
                    tau_absolute_valid = obj.tau_absolute(txrx,1:N_txrx);
                end
                tau_n_valid = obj.tau_n(txrx,1:N_txrx);
                tau_n_LOS_valid = obj.tau_n_LOS(txrx,1:N_txrx);
                Pn_valid = obj.Pn(txrx,1:N_txrx);
                Pn_LOS_valid = obj.Pn_LOS(txrx,1:N_txrx);
                shadow_valid = obj.cluster_shadow_dB(txrx,1:N_txrx);
                if ~isempty(tau_absolute_valid)
                    obj.tau_absolute(txrx,1:N_txrx) = [tau_absolute_valid(keep_tmp),nan(1,sum(~keep_tmp))];
                end
                obj.tau_n_keep(txrx,1:N_txrx)  = [tau_n_valid(keep_tmp),nan(1,sum(~keep_tmp))];
                obj.tau_n_LOS_keep(txrx,1:N_txrx)  = [tau_n_LOS_valid(keep_tmp),nan(1,sum(~keep_tmp))];
                obj.Pn_keep(txrx,1:N_txrx)     = [Pn_valid(keep_tmp),nan(1,sum(~keep_tmp))];
                obj.Pn_LOS_keep(txrx,1:N_txrx) = [Pn_LOS_valid(keep_tmp),nan(1,sum(~keep_tmp))];
                obj.cluster_shadow_dB(txrx,1:N_txrx) = [shadow_valid(keep_tmp),nan(1,sum(~keep_tmp))];
                obj.N_new(txrx)  = sum(obj.keep(txrx,1:N_txrx));
                [~,sort_idx]    = sort((1./obj.Pn(txrx,1:N_txrx)));
                obj.strong_cluster_id(txrx,:) = sort_idx(1:min(2,numel(sort_idx)));  % cluster number is 2
                if ~isempty(obj.c_DS)
                    obj.map_delay(txrx,[9,10,11,12,17,18]) = obj.c_DS(txrx)*1.28e-9;
                    obj.map_delay(txrx,13:16) = obj.c_DS(txrx)*2.56e-9;
                else
                    obj.map_delay(txrx,9:18) = 3.91e-9;
                end
                if obj.static == 1
                    obj.keep(2,:) = obj.keep(1,:);
                    if ~isempty(obj.tau_absolute)
                        obj.tau_absolute(2,:) = obj.tau_absolute(1,:);
                    end
                    obj.tau_n(2,:) = obj.tau_n(1,:);
                    obj.tau_n_LOS(2,:) = obj.tau_n_LOS(1,:);
                    obj.Pn(2,:) = obj.Pn(1,:);
                    obj.Pn_LOS(2,:) = obj.Pn_LOS(1,:);
                    obj.cluster_shadow_dB(2,:) = obj.cluster_shadow_dB(1,:);
                    obj.N_new(2) = obj.N_new(1);
                    obj.tau_n_keep(2,:)  = obj.tau_n_keep(1,:);
                    obj.tau_n_LOS_keep(2,:)  = obj.tau_n_LOS_keep(1,:);
                    obj.Pn_keep(2,:)     = obj.Pn_keep(1,:);
                    obj.Pn_LOS_keep(2,:) = obj.Pn_LOS_keep(1,:);
                    obj.strong_cluster_id(2,:) = obj.strong_cluster_id(1,:);
                    obj.map_delay(2,:) = obj.map_delay(1,:);
                    break;
                end
            end
        end

        function cluster_delay_procedureB(obj)
            max_N = max(obj.N);
            obj.tau_prime = nan(2, max_N);
            obj.tau_n = nan(2, max_N);
            obj.tau_order = nan(2, max_N);
            obj.tau_n_LOS = nan(2, max_N);
            for txrx = 1:2
                N_txrx = obj.N(txrx);
                limits = obj.procedureBClusterLimits(txrx);
                tau_prime_tmp = limits.delay*obj.drawProcedureBUniform(txrx, 'tau', N_txrx);
                [tau_sorted, tau_order_sorted] = sort(tau_prime_tmp - min(tau_prime_tmp));
                if obj.bLOS(txrx) && ~obj.O2I
                    tau_sorted(1) = 0;
                end
                obj.tau_prime(txrx,1:N_txrx) = tau_sorted;
                obj.tau_n(txrx,1:N_txrx) = tau_sorted;
                obj.tau_order(txrx,1:N_txrx) = tau_order_sorted;
                obj.tau_n_LOS(txrx,1:N_txrx) = tau_sorted;

                if obj.static == 1
                    obj.tau_prime(2,:) = obj.tau_prime(1,:);
                    obj.tau_n(2,:) = obj.tau_n(1,:);
                    obj.tau_order(2,:) = obj.tau_order(1,:);
                    obj.tau_n_LOS(2,:) = obj.tau_n_LOS(1,:);
                    break;
                end
            end
        end

        function cluster_angles_procedureB(obj)
            max_N = max(obj.N);
            obj.phi_prime_AOA = nan(2, max_N);
            obj.phi_prime_AOD = nan(2, max_N);
            obj.theta_prime_ZOA = nan(2, max_N);
            obj.theta_prime_ZOD = nan(2, max_N);
            for txrx = 1:2
                N_txrx = obj.N(txrx);
                order = obj.tau_order(txrx,1:N_txrx);
                limits = obj.procedureBClusterLimits(txrx);
                phi_prime_AOA_tmp = limits.AOA*(2*obj.drawProcedureBUniform(txrx, 'AOA', N_txrx)-1);
                phi_prime_AOD_tmp = limits.AOD*(2*obj.drawProcedureBUniform(txrx, 'AOD', N_txrx)-1);
                theta_prime_ZOA_tmp = limits.ZOA*(2*obj.drawProcedureBUniform(txrx, 'ZOA', N_txrx)-1);
                theta_prime_ZOD_tmp = limits.ZOD*(2*obj.drawProcedureBUniform(txrx, 'ZOD', N_txrx)-1);

                obj.phi_prime_AOA(txrx,1:N_txrx) = phi_prime_AOA_tmp(order);
                obj.phi_prime_AOD(txrx,1:N_txrx) = phi_prime_AOD_tmp(order);
                obj.theta_prime_ZOA(txrx,1:N_txrx) = theta_prime_ZOA_tmp(order);
                obj.theta_prime_ZOD(txrx,1:N_txrx) = theta_prime_ZOD_tmp(order);

                if obj.bLOS(txrx) && ~obj.O2I
                    obj.phi_prime_AOA(txrx,1) = 0;
                    obj.phi_prime_AOD(txrx,1) = 0;
                    obj.theta_prime_ZOA(txrx,1) = 0;
                    obj.theta_prime_ZOD(txrx,1) = 0;
                end

                if obj.static == 1
                    obj.phi_prime_AOA(2,:) = obj.phi_prime_AOD(1,:);
                    obj.phi_prime_AOD(2,:) = obj.phi_prime_AOA(1,:);
                    obj.theta_prime_ZOA(2,:) = obj.theta_prime_ZOD(1,:);
                    obj.theta_prime_ZOD(2,:) = obj.theta_prime_ZOA(1,:);
                    break;
                end
            end
        end

        function cluster_delay_procedureA(obj)
            max_N = max(obj.N);
            obj.tau_absolute = nan(2, max_N);
            obj.tau_prime = nan(2, max_N);
            obj.tau_n = nan(2, max_N);
            obj.tau_order = nan(2, max_N);
            obj.tau_n_LOS = nan(2, max_N);
            for txrx = 1:obj.TXRX
                state = obj.getProcedureAState(txrx);
                if isempty(state) || ~isfield(state, 'tau_absolute') || isempty(state.tau_absolute)
                    obj.cluster_delay_procedureA_initial();
                    return;
                end

                N_txrx = numel(state.tau_absolute);
                obj.N(txrx) = N_txrx;
                delta_t = max(obj.t - state.time, 0);
                [tx_velocity, rx_velocity] = obj.targetProcedureAVelocities(txrx);
                tx_delta = tx_velocity * delta_t;
                rx_delta = rx_velocity * delta_t;
                tau_update = zeros(1, N_txrx);
                for cluster_idx = 1:N_txrx
                    aoa_unit = obj.unitVectorFromAngles( ...
                        obj.phi_LOS_AOA(txrx) + state.phi_prime_AOA(cluster_idx), ...
                        obj.theta_LOS_ZOA(txrx) + state.theta_prime_ZOA(cluster_idx));
                    aod_unit = obj.unitVectorFromAngles( ...
                        obj.phi_LOS_AOD(txrx) + state.phi_prime_AOD(cluster_idx), ...
                        obj.theta_LOS_ZOD(txrx) + state.theta_prime_ZOD(cluster_idx));
                    tau_update(cluster_idx) = (-dot(rx_delta, aoa_unit) + dot(tx_delta, aod_unit)) / 3e8;
                end

                tau_abs = max(state.tau_absolute + tau_update, 0);
                tau_tmp = tau_abs - min(tau_abs);
                if obj.bLOS(txrx) && ~obj.O2I
                    tau_tmp(1) = 0;
                end
                obj.tau_absolute(txrx,1:N_txrx) = tau_abs;
                obj.tau_prime(txrx,1:N_txrx) = tau_tmp;
                obj.tau_n(txrx,1:N_txrx) = tau_tmp;
                obj.tau_order(txrx,1:N_txrx) = state.tau_order;
                obj.tau_n_LOS(txrx,1:N_txrx) = tau_tmp;

                if obj.static == 1
                    obj.tau_prime(2,:) = obj.tau_prime(1,:);
                    obj.tau_n(2,:) = obj.tau_n(1,:);
                    obj.tau_order(2,:) = obj.tau_order(1,:);
                    obj.tau_n_LOS(2,:) = obj.tau_n_LOS(1,:);
                    break;
                end
            end
        end

        function cluster_delay_procedureA_initial(obj)
            max_N = max(obj.N);
            obj.tau_absolute = nan(2, max_N);
            obj.tau_prime = nan(2, max_N);
            obj.tau_n = nan(2, max_N);
            obj.tau_order = nan(2, max_N);
            obj.tau_n_LOS = nan(2, max_N);
            for txrx = 1:obj.TXRX
                N_txrx = obj.N(txrx);
                tau_raw = -obj.r_tau(txrx)*obj.DS(txrx)*log(rand(1, N_txrx));
                tau_absolute_tmp = obj.d_3D(txrx)/3e8 + tau_raw;
                [tau_tmp, tau_order_tmp] = sort(tau_absolute_tmp - min(tau_absolute_tmp));
                tau_absolute_tmp = tau_absolute_tmp(tau_order_tmp);
                if obj.bLOS(txrx) && ~obj.O2I
                    tau_tmp(1) = 0;
                end
                obj.tau_absolute(txrx,1:N_txrx) = tau_absolute_tmp;
                obj.tau_prime(txrx,1:N_txrx) = tau_tmp;
                obj.tau_n(txrx,1:N_txrx) = tau_tmp;
                obj.tau_order(txrx,1:N_txrx) = tau_order_tmp;
                obj.tau_n_LOS(txrx,1:N_txrx) = tau_tmp;

                if obj.static == 1
                    obj.tau_absolute(2,:) = obj.tau_absolute(1,:);
                    obj.tau_prime(2,:) = obj.tau_prime(1,:);
                    obj.tau_n(2,:) = obj.tau_n(1,:);
                    obj.tau_order(2,:) = obj.tau_order(1,:);
                    obj.tau_n_LOS(2,:) = obj.tau_n_LOS(1,:);
                    break;
                end
            end
        end

        function cluster_angles_procedureA(obj)
            max_N = max(obj.N);
            obj.phi_prime_AOA = nan(2, max_N);
            obj.phi_prime_AOD = nan(2, max_N);
            obj.theta_prime_ZOA = nan(2, max_N);
            obj.theta_prime_ZOD = nan(2, max_N);
            for txrx = 1:obj.TXRX
                state = obj.getProcedureAState(txrx);
                N_txrx = obj.N(txrx);
                if isempty(state) || ~isfield(state, 'phi_prime_AOA') || numel(state.phi_prime_AOA) ~= N_txrx
                    obj.cluster_angles_procedureB();
                    return;
                end

                delta_t = max(obj.t - state.time, 0);
                [tx_velocity, rx_velocity] = obj.targetProcedureAVelocities(txrx);

                prev_phi_AOA = state.phi_LOS_AOA + state.phi_prime_AOA;
                prev_phi_AOD = state.phi_LOS_AOD + state.phi_prime_AOD;
                prev_theta_ZOA = state.theta_LOS_ZOA + state.theta_prime_ZOA;
                prev_theta_ZOD = state.theta_LOS_ZOD + state.theta_prime_ZOD;
                x_n = obj.procedureAXnFromState(state, N_txrx, txrx);
                [rx_velocity_eff, tx_velocity_eff] = obj.procedureAEffectiveVelocities( ...
                    prev_phi_AOA, prev_theta_ZOA, prev_phi_AOD, prev_theta_ZOD, ...
                    tx_velocity, rx_velocity, x_n, obj.bLOS(txrx) && ~obj.O2I);

                [phi_AOA, theta_ZOA] = obj.updateProcedureAAnglePair( ...
                    prev_phi_AOA, prev_theta_ZOA, rx_velocity_eff, obj.tau_n(txrx,1:N_txrx), delta_t);
                [phi_AOD, theta_ZOD] = obj.updateProcedureAAnglePair( ...
                    prev_phi_AOD, prev_theta_ZOD, tx_velocity_eff, obj.tau_n(txrx,1:N_txrx), delta_t);

                obj.phi_prime_AOA(txrx,1:N_txrx) = obj.wrapAzimuthAngles(phi_AOA - obj.phi_LOS_AOA(txrx));
                obj.phi_prime_AOD(txrx,1:N_txrx) = obj.wrapAzimuthAngles(phi_AOD - obj.phi_LOS_AOD(txrx));
                obj.theta_prime_ZOA(txrx,1:N_txrx) = theta_ZOA - obj.theta_LOS_ZOA(txrx);
                obj.theta_prime_ZOD(txrx,1:N_txrx) = theta_ZOD - obj.theta_LOS_ZOD(txrx);

                if obj.bLOS(txrx) && ~obj.O2I
                    obj.phi_prime_AOA(txrx,1) = 0;
                    obj.phi_prime_AOD(txrx,1) = 0;
                    obj.theta_prime_ZOA(txrx,1) = 0;
                    obj.theta_prime_ZOD(txrx,1) = 0;
                end

                if obj.static == 1
                    obj.phi_prime_AOA(2,:) = obj.phi_prime_AOD(1,:);
                    obj.phi_prime_AOD(2,:) = obj.phi_prime_AOA(1,:);
                    obj.theta_prime_ZOA(2,:) = obj.theta_prime_ZOD(1,:);
                    obj.theta_prime_ZOD(2,:) = obj.theta_prime_ZOA(1,:);
                    break;
                end
            end
        end

        function cluster_power_procedureB(obj)
            max_N = max(obj.N);
            max_M = max(obj.M);
            obj.Pn = nan(2, max_N);
            obj.Pn_LOS = nan(2, max_N);
            obj.cluster_shadow_dB = nan(2, max_N);
            obj.keep = false(2, max_N);
            obj.tau_n_keep = nan(2, max_N);
            obj.tau_n_LOS_keep = nan(2, max_N);
            obj.Pn_keep = nan(2, max_N);
            obj.Pn_LOS_keep = nan(2, max_N);
            obj.N_new = zeros(1, 2);
            obj.strong_cluster_id = nan(2, 2);
            obj.map_delay = zeros(2, max_M);
            for txrx = 1:2
                N_txrx = obj.N(txrx);
                shadow_dB = obj.drawProcedureBGaussian(txrx, 'shadow', N_txrx).' * obj.zeta(txrx);
                obj.cluster_shadow_dB(txrx,1:N_txrx) = shadow_dB.';
                ds_power = obj.DS(txrx);
                asa_power = obj.ASA(txrx);
                asd_power = obj.ASD(txrx);
                zsa_power = obj.ZSA(txrx);
                zsd_power = obj.ZSD(txrx);
                rice_linear = 0;

                if obj.bLOS(txrx) && ~obj.O2I
                    rice_linear = 10^(obj.K(txrx)/10);
                    ds_power = ds_power*(1 + rice_linear);
                    asa_power = asa_power*(1 + rice_linear);
                    asd_power = asd_power*(1 + rice_linear);
                    zsa_power = zsa_power*(1 + rice_linear);
                    zsd_power = zsd_power*(1 + rice_linear);
                end

                cluster_power_raw = exp(-obj.tau_prime(txrx,1:N_txrx)./ds_power) ...
                    .* exp(-sqrt(2)*abs(obj.phi_prime_AOA(txrx,1:N_txrx))./asa_power) ...
                    .* exp(-sqrt(2)*abs(obj.phi_prime_AOD(txrx,1:N_txrx))./asd_power) ...
                    .* exp(-sqrt(2)*abs(obj.theta_prime_ZOA(txrx,1:N_txrx))./zsa_power) ...
                    .* exp(-sqrt(2)*abs(obj.theta_prime_ZOD(txrx,1:N_txrx))./zsd_power) ...
                    .* 10.^(-shadow_dB.'/10);

                if obj.bLOS(txrx) && ~obj.O2I
                    obj.Pn(txrx,1:N_txrx) = (cluster_power_raw./sum(cluster_power_raw))/(1 + rice_linear);
                    obj.Pn(txrx,1) = obj.Pn(txrx,1) + rice_linear/(1 + rice_linear);
                else
                    obj.Pn(txrx,1:N_txrx) = cluster_power_raw./sum(cluster_power_raw);
                end
                obj.Pn_LOS(txrx,1:N_txrx) = obj.Pn(txrx,1:N_txrx);

                keep_tmp = true(1, N_txrx);
                obj.keep(txrx,1:N_txrx) = keep_tmp;
                tau_n_valid = obj.tau_n(txrx,1:N_txrx);
                tau_n_LOS_valid = obj.tau_n_LOS(txrx,1:N_txrx);
                Pn_valid = obj.Pn(txrx,1:N_txrx);
                Pn_LOS_valid = obj.Pn_LOS(txrx,1:N_txrx);
                phi_AOA_valid = obj.phi_prime_AOA(txrx,1:N_txrx);
                phi_AOD_valid = obj.phi_prime_AOD(txrx,1:N_txrx);
                theta_ZOA_valid = obj.theta_prime_ZOA(txrx,1:N_txrx);
                theta_ZOD_valid = obj.theta_prime_ZOD(txrx,1:N_txrx);
                shadow_valid = obj.cluster_shadow_dB(txrx,1:N_txrx);
                obj.tau_n_keep(txrx,1:N_txrx) = [tau_n_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.tau_n_LOS_keep(txrx,1:N_txrx) = [tau_n_LOS_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.Pn_keep(txrx,1:N_txrx) = [Pn_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.Pn_LOS_keep(txrx,1:N_txrx) = [Pn_LOS_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.phi_prime_AOA(txrx,1:N_txrx) = [phi_AOA_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.phi_prime_AOD(txrx,1:N_txrx) = [phi_AOD_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.theta_prime_ZOA(txrx,1:N_txrx) = [theta_ZOA_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.theta_prime_ZOD(txrx,1:N_txrx) = [theta_ZOD_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.cluster_shadow_dB(txrx,1:N_txrx) = [shadow_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.N_new(txrx) = sum(keep_tmp);
                [~, sort_idx] = sort((1./obj.Pn(txrx,1:N_txrx)));
                obj.strong_cluster_id(txrx,:) = sort_idx(1:min(2,numel(sort_idx)));
                if ~isempty(obj.c_DS)
                    obj.map_delay(txrx,[9,10,11,12,17,18]) = obj.c_DS(txrx)*1.28e-9;
                    obj.map_delay(txrx,13:16) = obj.c_DS(txrx)*2.56e-9;
                else
                    obj.map_delay(txrx,9:18) = 3.91e-9;
                end

                if obj.static == 1
                    obj.keep(2,:) = obj.keep(1,:);
                    obj.tau_n(2,:) = obj.tau_n(1,:);
                    obj.tau_n_LOS(2,:) = obj.tau_n_LOS(1,:);
                    obj.tau_prime(2,:) = obj.tau_prime(1,:);
                    obj.Pn(2,:) = obj.Pn(1,:);
                    obj.Pn_LOS(2,:) = obj.Pn_LOS(1,:);
                    obj.N_new(2) = obj.N_new(1);
                    obj.tau_n_keep(2,:) = obj.tau_n_keep(1,:);
                    obj.tau_n_LOS_keep(2,:) = obj.tau_n_LOS_keep(1,:);
                    obj.Pn_keep(2,:) = obj.Pn_keep(1,:);
                    obj.Pn_LOS_keep(2,:) = obj.Pn_LOS_keep(1,:);
                    obj.phi_prime_AOA(2,:) = obj.phi_prime_AOD(1,:);
                    obj.phi_prime_AOD(2,:) = obj.phi_prime_AOA(1,:);
                    obj.theta_prime_ZOA(2,:) = obj.theta_prime_ZOD(1,:);
                    obj.theta_prime_ZOD(2,:) = obj.theta_prime_ZOA(1,:);
                    obj.cluster_shadow_dB(2,:) = obj.cluster_shadow_dB(1,:);
                    obj.strong_cluster_id(2,:) = obj.strong_cluster_id(1,:);
                    obj.map_delay(2,:) = obj.map_delay(1,:);
                    break;
                end
            end
        end

        function cluster_power_procedureA(obj)
            max_N = max(obj.N);
            max_M = max(obj.M);
            obj.Pn = nan(2, max_N);
            obj.Pn_LOS = nan(2, max_N);
            obj.cluster_shadow_dB = nan(2, max_N);
            obj.keep = false(2, max_N);
            obj.tau_n_keep = nan(2, max_N);
            obj.tau_n_LOS_keep = nan(2, max_N);
            obj.Pn_keep = nan(2, max_N);
            obj.Pn_LOS_keep = nan(2, max_N);
            obj.N_new = zeros(1, 2);
            obj.strong_cluster_id = nan(2, 2);
            obj.map_delay = zeros(2, max_M);
            for txrx = 1:obj.TXRX
                N_txrx = obj.N(txrx);
                state = obj.getProcedureAState(txrx);
                if ~isempty(state) && isfield(state, 'cluster_shadow_dB') && ...
                        numel(state.cluster_shadow_dB) == N_txrx
                    obj.cluster_shadow_dB(txrx,1:N_txrx) = state.cluster_shadow_dB(:).';
                else
                    obj.cluster_shadow_dB(txrx,1:N_txrx) = (randn(N_txrx, 1) * obj.zeta(txrx)).';
                end

                rice_linear = 0;
                if obj.bLOS(txrx) && ~obj.O2I
                    rice_linear = 10^(obj.K(txrx)/10);
                end

                cluster_power_raw = exp(-obj.tau_n(txrx,1:N_txrx) ...
                    .*(obj.r_tau(txrx)-1)./obj.r_tau(txrx)./obj.DS(txrx)) ...
                    .* 10.^(-obj.cluster_shadow_dB(txrx,1:N_txrx)/10);

                if obj.bLOS(txrx) && ~obj.O2I
                    obj.Pn(txrx,1:N_txrx) = (cluster_power_raw./sum(cluster_power_raw))/(1 + rice_linear);
                    obj.Pn(txrx,1) = obj.Pn(txrx,1) + rice_linear/(1 + rice_linear);
                else
                    obj.Pn(txrx,1:N_txrx) = cluster_power_raw./sum(cluster_power_raw);
                end
                obj.Pn_LOS(txrx,1:N_txrx) = obj.Pn(txrx,1:N_txrx);

                keep_tmp = true(1, N_txrx);
                obj.keep(txrx,1:N_txrx) = keep_tmp;
                tau_n_valid = obj.tau_n(txrx,1:N_txrx);
                tau_n_LOS_valid = obj.tau_n_LOS(txrx,1:N_txrx);
                Pn_valid = obj.Pn(txrx,1:N_txrx);
                Pn_LOS_valid = obj.Pn_LOS(txrx,1:N_txrx);
                phi_AOA_valid = obj.phi_prime_AOA(txrx,1:N_txrx);
                phi_AOD_valid = obj.phi_prime_AOD(txrx,1:N_txrx);
                theta_ZOA_valid = obj.theta_prime_ZOA(txrx,1:N_txrx);
                theta_ZOD_valid = obj.theta_prime_ZOD(txrx,1:N_txrx);
                shadow_valid = obj.cluster_shadow_dB(txrx,1:N_txrx);
                obj.tau_n_keep(txrx,1:N_txrx) = [tau_n_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.tau_n_LOS_keep(txrx,1:N_txrx) = [tau_n_LOS_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.Pn_keep(txrx,1:N_txrx) = [Pn_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.Pn_LOS_keep(txrx,1:N_txrx) = [Pn_LOS_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.phi_prime_AOA(txrx,1:N_txrx) = [phi_AOA_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.phi_prime_AOD(txrx,1:N_txrx) = [phi_AOD_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.theta_prime_ZOA(txrx,1:N_txrx) = [theta_ZOA_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.theta_prime_ZOD(txrx,1:N_txrx) = [theta_ZOD_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.cluster_shadow_dB(txrx,1:N_txrx) = [shadow_valid(keep_tmp), nan(1,sum(~keep_tmp))];
                obj.N_new(txrx) = sum(keep_tmp);
                [~, sort_idx] = sort((1./obj.Pn(txrx,1:N_txrx)));
                obj.strong_cluster_id(txrx,:) = sort_idx(1:min(2,numel(sort_idx)));
                if ~isempty(obj.c_DS)
                    obj.map_delay(txrx,[9,10,11,12,17,18]) = obj.c_DS(txrx)*1.28e-9;
                    obj.map_delay(txrx,13:16) = obj.c_DS(txrx)*2.56e-9;
                else
                    obj.map_delay(txrx,9:18) = 3.91e-9;
                end

                if obj.static == 1
                    obj.keep(2,:) = obj.keep(1,:);
                    obj.tau_n(2,:) = obj.tau_n(1,:);
                    obj.tau_n_LOS(2,:) = obj.tau_n_LOS(1,:);
                    obj.tau_prime(2,:) = obj.tau_prime(1,:);
                    obj.Pn(2,:) = obj.Pn(1,:);
                    obj.Pn_LOS(2,:) = obj.Pn_LOS(1,:);
                    obj.cluster_shadow_dB(2,:) = obj.cluster_shadow_dB(1,:);
                    obj.N_new(2) = obj.N_new(1);
                    obj.tau_n_keep(2,:) = obj.tau_n_keep(1,:);
                    obj.tau_n_LOS_keep(2,:) = obj.tau_n_LOS_keep(1,:);
                    obj.Pn_keep(2,:) = obj.Pn_keep(1,:);
                    obj.Pn_LOS_keep(2,:) = obj.Pn_LOS_keep(1,:);
                    obj.phi_prime_AOA(2,:) = obj.phi_prime_AOD(1,:);
                    obj.phi_prime_AOD(2,:) = obj.phi_prime_AOA(1,:);
                    obj.theta_prime_ZOA(2,:) = obj.theta_prime_ZOD(1,:);
                    obj.theta_prime_ZOD(2,:) = obj.theta_prime_ZOA(1,:);
                    obj.strong_cluster_id(2,:) = obj.strong_cluster_id(1,:);
                    obj.map_delay(2,:) = obj.map_delay(1,:);
                    break;
                end
            end
        end

        function apply_angle_offsets_procedureB(obj)
            max_N_new = max(obj.N_new);
            max_M = max(obj.M);
            obj.phi_n_m_AOA = nan(2, max_N_new, max_M);
            obj.phi_n_m_AOD = nan(2, max_N_new, max_M);
            obj.theta_n_m_ZOA = nan(2, max_N_new, max_M);
            obj.theta_n_m_ZOD = nan(2, max_N_new, max_M);
            ray_angle_offset = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];

            for txrx = 1:2
                N_new_txrx = obj.N_new(txrx);
                phi_n_AOA = obj.phi_LOS_AOA(txrx) + obj.phi_prime_AOA(txrx,1:N_new_txrx);
                phi_n_AOD = obj.phi_LOS_AOD(txrx) + obj.phi_prime_AOD(txrx,1:N_new_txrx);
                ZOA = obj.theta_LOS_ZOA(txrx);
                ZOA(obj.O2I) = 90;
                theta_n_ZOA = ZOA + obj.theta_prime_ZOA(txrx,1:N_new_txrx);
                theta_n_ZOD = obj.theta_LOS_ZOD(txrx) + obj.mu_offset_ZOD(txrx) + obj.theta_prime_ZOD(txrx,1:N_new_txrx);

                obj.phi_n_m_AOA(txrx,1:N_new_txrx,1:obj.M(txrx)) = obj.wrapAzimuthAngles(repmat(phi_n_AOA.',1,obj.M(txrx)) + obj.c_ASA(txrx)*ones(N_new_txrx,1)*ray_angle_offset);
                obj.phi_n_m_AOD(txrx,1:N_new_txrx,1:obj.M(txrx)) = obj.wrapAzimuthAngles(repmat(phi_n_AOD.',1,obj.M(txrx)) + obj.c_ASD(txrx)*ones(N_new_txrx,1)*ray_angle_offset);
                obj.theta_n_m_ZOA(txrx,1:N_new_txrx,1:obj.M(txrx)) = obj.wrapZenithAngles(repmat(theta_n_ZOA.',1,obj.M(txrx)) + obj.c_ZSA(txrx)*ones(N_new_txrx,1)*ray_angle_offset);

                zsd_offset_scale = (3/8)*(10^obj.mu_lgZSD(txrx));
                obj.theta_n_m_ZOD(txrx,1:N_new_txrx,1:obj.M(txrx)) = obj.wrapZenithAngles(repmat(theta_n_ZOD.',1,obj.M(txrx)) + zsd_offset_scale*ones(N_new_txrx,1)*ray_angle_offset);

                if obj.static == 1
                    obj.phi_n_m_AOD(2,:,:) = obj.phi_n_m_AOA(1,:,:);
                    obj.phi_n_m_AOA(2,:,:) = obj.phi_n_m_AOD(1,:,:);
                    obj.theta_n_m_ZOD(2,:,:) = obj.theta_n_m_ZOA(1,:,:);
                    obj.theta_n_m_ZOA(2,:,:) = obj.theta_n_m_ZOD(1,:,:);
                    if isa(obj.ST,"elements.RP")
                        obj.phi_n_m_AOA = obj.phi_n_m_AOD;
                    end
                    break;
                end
            end
        end

        % caculator the AOA for every ray
        function AOA_calc(obj)
            for txrx = 1:2
                n_clusters     = [4 5 6 7 8 10 11 12 14 15 16 19 20 25];
                nlos_az_scale_table = [0.779, 0.860, 0.921, 0.973, 1.018, 1.090, 1.123, 1.146, 1.190, 1.211, 1.226, 1.273, 1.289, 1.358];% , 1.358
                nlos_az_scale     = nlos_az_scale_table(n_clusters == obj.N(txrx));
                if obj.bLOS(txrx) && ~obj.O2I
                    az_scale          = sum([nlos_az_scale ,nlos_az_scale*(0.1035-0.028*obj.K(txrx)-0.002*obj.K(txrx)^2+0.0001*obj.K(txrx)^3)]); % for NLOS and O2I, K = [].
                else
                    az_scale          = nlos_az_scale;
                end
                phi_n_AOA_tmp  = 2*(obj.ASA(txrx)/1.4)*sqrt(-log(obj.Pn_LOS_keep(txrx,1:obj.N_new(txrx))/max(obj.Pn_LOS_keep(txrx,1:obj.N_new(txrx)))))/az_scale;

                    random_draw = randsrc(1,obj.N_new(txrx),[-1, 1]);
                    angle_jitter = (obj.ASA(txrx)/7)*randn(1,obj.N_new(txrx));
                % end
                phi_n_AOA = random_draw.*phi_n_AOA_tmp+angle_jitter+obj.phi_LOS_AOA(txrx);
                if obj.bLOS(txrx) % && ~obj.O2I
                    phi_n_AOA = (random_draw.*phi_n_AOA_tmp+angle_jitter)-(random_draw(1)*phi_n_AOA_tmp(1)+angle_jitter(1)-obj.phi_LOS_AOA(txrx));
                end
                ray_angle_offset = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                    0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];
                phi_n_m_AOA_ = repmat(phi_n_AOA.',1,obj.M(txrx))+obj.c_ASA(txrx)*ones(obj.N_new(txrx),1)*ray_angle_offset;
                phi_n_m_AOA_ = mod(phi_n_m_AOA_,360);
                phi_n_m_AOA_ = phi_n_m_AOA_-360*floor(phi_n_m_AOA_/180);
                obj.phi_n_m_AOA(txrx,:,:) = phi_n_m_AOA_;
                if obj.static ==1
                    break;
                end
            end
        end

        % caculator the AOD for every ray
        function AOD_calc(obj)
            for txrx = 1:2
                n_clusters     = [4 5 6 7 8 10 11 12 14 15 16 19 20 25];
                nlos_az_scale_table = [0.779, 0.860, 0.921, 0.973, 1.018, 1.090, 1.123, 1.146, 1.190, 1.211, 1.226, 1.273, 1.289, 1.358];% , 1.358
                nlos_az_scale     = nlos_az_scale_table(n_clusters == obj.N(txrx));
                if obj.bLOS(txrx) && ~obj.O2I
                    az_scale          = sum([nlos_az_scale ,nlos_az_scale*(0.1035-0.028*obj.K(txrx)-0.002*obj.K(txrx)^2+0.0001*obj.K(txrx)^3)]); % for NLOS and O2I, K = [].
                else
                    az_scale          = nlos_az_scale;
                end
                phi_n_AOD_tmp  = 2*(obj.ASD(txrx)/1.4)*sqrt(-log(obj.Pn_LOS_keep(txrx,1:obj.N_new(txrx))/max(obj.Pn_LOS_keep(txrx,1:obj.N_new(txrx)))))/az_scale;

                    random_draw = randsrc(1,obj.N_new(txrx),[-1, 1]);
                    angle_jitter = (obj.ASD(txrx)/7)*randn(1,obj.N_new(txrx));
                % end
                phi_n_AOD = random_draw.*phi_n_AOD_tmp+angle_jitter+obj.phi_LOS_AOD(txrx);
                if obj.bLOS(txrx) % && ~obj.O2I
                    phi_n_AOD = (random_draw.*phi_n_AOD_tmp+angle_jitter)-(random_draw(1)*phi_n_AOD_tmp(1)+angle_jitter(1)-obj.phi_LOS_AOD(txrx));
                end
                ray_angle_offset = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                    0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];
                phi_n_m_AOD_ = repmat(phi_n_AOD.',1,obj.M(txrx))+obj.c_ASD(txrx)*ones(obj.N_new(txrx),1)*ray_angle_offset;
                phi_n_m_AOD_ = mod(phi_n_m_AOD_,360);
                phi_n_m_AOD_ = phi_n_m_AOD_-360*floor(phi_n_m_AOD_/180);
                obj.phi_n_m_AOD(txrx,:,:) = phi_n_m_AOD_;
                if obj.static ==1
                    obj.phi_n_m_AOD(2,:,:)= obj.phi_n_m_AOA(1,:,:);
                    obj.phi_n_m_AOA(2,:,:)= obj.phi_n_m_AOD(1,:,:);
                    if isa(obj.ST,"elements.RP")
                        obj.phi_n_m_AOA = obj.phi_n_m_AOD;
                    end
                    break;
                end
            end
        end

        % caculator the ZOA for every ray
        function ZOA_calc(obj)
            for txrx = 1:2
                n_clusters       = [6 7 8 10 11 12 14 15 16 19 20 25];
                nlos_zenith_scale_table = [0.788, 0.847, 0.889, 0.957, 1.031, 1.104, 1.1072, 1.1088, 1.1276, 1.184, 1.178, 1.282]; % 1.282
                nlos_zenith_scale     = nlos_zenith_scale_table(n_clusters == obj.N(txrx));
                if obj.bLOS(txrx) && ~obj.O2I
                    zenith_scale          = sum([nlos_zenith_scale ,nlos_zenith_scale*(0.3086+0.0339*obj.K(txrx)-0.0077*obj.K(txrx)^2+0.0002*obj.K(txrx)^3)]); % for NLOS and O2I, K = [].
                else
                    zenith_scale          = nlos_zenith_scale;
                end
                theta_n_ZOA_tmp  = -obj.ZSA(txrx)*log(obj.Pn_LOS_keep(txrx,1:obj.N_new(txrx))/max(obj.Pn_LOS_keep(txrx,1:obj.N_new(txrx))))/zenith_scale;

                    random_draw  = randsrc(1,obj.N_new(txrx),[-1, 1]);
                    angle_jitter  = (obj.ZSA(txrx)/7)*randn(1,obj.N_new(txrx));
                % end
                ZOA = obj.theta_LOS_ZOA(txrx);
                ZOA(obj.O2I) = 90;
                theta_n_ZOA = random_draw.*theta_n_ZOA_tmp+angle_jitter+ZOA;
                if obj.bLOS(txrx)  && ~obj.O2I
                    theta_n_ZOA = (random_draw.*theta_n_ZOA_tmp+angle_jitter)-(random_draw(1)*theta_n_ZOA_tmp(1)+angle_jitter(1)-obj.theta_LOS_ZOA(txrx));
                end
                ray_angle_offset = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                    0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];
                theta_n_m_ZOA_ = repmat(theta_n_ZOA.',1,obj.M(txrx)) +obj.c_ZSA(txrx)*ones(obj.N_new(txrx),1)*ray_angle_offset;
                wrapped_zoa_angle = mod((theta_n_m_ZOA_+360),360);  % 360
                over_half_circle = ((wrapped_zoa_angle <= 360)&(wrapped_zoa_angle >= 180));
                theta_n_m_ZOA_ = (~over_half_circle).*wrapped_zoa_angle + over_half_circle.*(360 - wrapped_zoa_angle);
                obj.theta_n_m_ZOA(txrx,:,:) = theta_n_m_ZOA_;
                if obj.static ==1
                    break;
                end
            end
        end

        % caculator the ZOD for every ray
        function ZOD_calc(obj)
            for txrx = 1:2
                n_clusters       = [6 7 8 10 11 12 14 15 16 19 20 25];
                nlos_zenith_scale_table = [0.788, 0.847, 0.889, 0.957, 1.031, 1.104, 1.1072, 1.1088, 1.1276, 1.184, 1.178, 1.282]; % 1.282
                nlos_zenith_scale     = nlos_zenith_scale_table(n_clusters == obj.N(txrx));
                if obj.bLOS(txrx) && ~obj.O2I
                    zenith_scale          = sum([nlos_zenith_scale ,nlos_zenith_scale*(0.3086+0.0339*obj.K(txrx)-0.0077*obj.K(txrx)^2+0.0002*obj.K(txrx)^3)]); % for NLOS and O2I, K = [].
                else
                    zenith_scale          = nlos_zenith_scale;
                end
                theta_n_ZOD_tmp  = -obj.ZSD(txrx)*log(obj.Pn_LOS_keep(txrx,1:obj.N_new(txrx))/max(obj.Pn_LOS_keep(txrx,1:obj.N_new(txrx))))/zenith_scale;

                    random_draw = randsrc(1,obj.N_new(txrx),[-1, 1]);
                    angle_jitter = (obj.ZSD(txrx)/7)*randn(1,obj.N_new(txrx));
                % end
                theta_n_ZOD = random_draw.*theta_n_ZOD_tmp+angle_jitter+obj.theta_LOS_ZOD(txrx)+obj.mu_offset_ZOD(txrx);
                if obj.bLOS(txrx) % && ~obj.O2I
                    theta_n_ZOD = (random_draw.*theta_n_ZOD_tmp+angle_jitter)-(random_draw(1)*theta_n_ZOD_tmp(1)+angle_jitter(1)-obj.theta_LOS_ZOD(txrx));
                end
                ray_angle_offset = [0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3715, 0.5129, -0.5129,...
                    0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551];
                zsd_offset_scale = (3/8)*(10^obj.mu_lgZSD(txrx));
                theta_n_m_ZOD_  = repmat(theta_n_ZOD.',1,obj.M(txrx)) +zsd_offset_scale*ones(obj.N_new(txrx),1)*ray_angle_offset;
                wrapped_zod_angle = mod((theta_n_m_ZOD_+360),360);   % 360
                over_half_circle = ((wrapped_zod_angle <= 360)&(wrapped_zod_angle >= 180));
                theta_n_m_ZOD_ = (~over_half_circle).*wrapped_zod_angle + over_half_circle.*(360 - wrapped_zod_angle);
                obj.theta_n_m_ZOD(txrx,:,:) = theta_n_m_ZOD_;
                if obj.static ==1
                    obj.theta_n_m_ZOD(2,:,:)= obj.theta_n_m_ZOA(1,:,:);
                    obj.theta_n_m_ZOA(2,:,:)= obj.theta_n_m_ZOD(1,:,:);
                    break;
                end
            end
        end

        function RandomCouplingRays(obj)
            rayGroups = obj.strongClusterRayGroups();
            for txrx = 1:2
                [rn1, rn2, rn3] = obj.drawTargetCouplingSeeds(txrx);
                obj.coupleTargetLinkRays(txrx, rn1, rn2, rn3, rayGroups);
                if obj.static ==1
                    obj.coupleTargetLinkRays(2, rn1, rn2, rn3, rayGroups);
                    break;
                end
            end
        end

        function Coupling_txrx(obj)
            coupling_txrx_n_m_table_ = [];
            if obj.bLOS(1) && obj.bLOS(2)
                coupling_txrx_n_m_table_= [coupling_txrx_n_m_table_;0,0,0,0];
            end
            if obj.bLOS(1)
                [coupling_los_nlos_rx_n, coupling_los_nlos_rx_m] = ndgrid(1:obj.N_new(2),1:obj.M(2));
                coupling_los_nlos = [zeros(numel(coupling_los_nlos_rx_n),2), coupling_los_nlos_rx_n(:), coupling_los_nlos_rx_m(:)];
                coupling_txrx_n_m_table_= [coupling_txrx_n_m_table_;coupling_los_nlos];
            end
            if obj.bLOS(2)
                [coupling_los_nlos_tx_n, coupling_los_nlos_tx_m] = ndgrid(1:obj.N_new(1),1:obj.M(1));
                coupling_los_nlos = [coupling_los_nlos_tx_n(:), coupling_los_nlos_tx_m(:),zeros(numel(coupling_los_nlos_tx_n),2)];
                coupling_txrx_n_m_table_= [coupling_txrx_n_m_table_;coupling_los_nlos];
            end
            if obj.option == 1
                [coupling_nlos_nlos_tx_n, coupling_nlos_nlos_tx_m] = ndgrid(1:obj.N_new(1),1:obj.M(1));
                [coupling_nlos_nlos_rx_n, coupling_nlos_nlos_rx_m] = ndgrid(1:obj.N_new(2),1:obj.M(2));
                coupling_nlos_nlos_tx = [coupling_nlos_nlos_tx_n(:), coupling_nlos_nlos_tx_m(:)];
                coupling_nlos_nlos_rx = [coupling_nlos_nlos_rx_n(:), coupling_nlos_nlos_rx_m(:)];
                coupling_nlos_nlos = [ repelem(coupling_nlos_nlos_tx, size(coupling_nlos_nlos_rx,1), 1),repmat(coupling_nlos_nlos_rx, size(coupling_nlos_nlos_rx,1), 1)];
                coupling_txrx_n_m_table_= [coupling_txrx_n_m_table_;coupling_nlos_nlos];
            else
                min_txrx_ray_num = min(obj.N_new(1)*obj.M(1), obj.N_new(2)*obj.M(2));
                randperm_tx = randperm(obj.N_new(1)*obj.M(1), min_txrx_ray_num);
                randperm_rx = randperm(obj.N_new(2)*obj.M(2), min_txrx_ray_num);
                [coupling_nlos_nlos_tx_n, coupling_nlos_nlos_tx_m] = ind2sub([obj.N_new(1),obj.M(1)], randperm_tx);
                [coupling_nlos_nlos_rx_n, coupling_nlos_nlos_rx_m]  = ind2sub([obj.N_new(2),obj.M(2)], randperm_rx);
                coupling_nlos_nlos = [coupling_nlos_nlos_tx_n(:), coupling_nlos_nlos_tx_m(:), coupling_nlos_nlos_rx_n(:), coupling_nlos_nlos_rx_m(:)];
                coupling_txrx_n_m_table_= [coupling_txrx_n_m_table_;coupling_nlos_nlos];
            end
            obj.coupling_txrx_path_table = coupling_txrx_n_m_table_;
            obj.path_num = size(obj.coupling_txrx_path_table,1);
        end

        
        function Path_power(obj)
            power_tx = zeros(obj.path_num,1);
            power_rx = zeros(obj.path_num,1);
            if obj.bLOS(1)
                KR_tx = 10^(obj.K(1)/10);

            else
                KR_tx = 0;
            end
            if obj.bLOS(2)
                KR_rx = 10^(obj.K(2)/10);
            else
                KR_rx = 0;
            end
            P_tx_LOS = KR_tx/(KR_tx+1);
            P_rx_LOS = KR_rx/(KR_rx+1);

            is_0_tx = (obj.coupling_txrx_path_table(:,1)==0 & obj.coupling_txrx_path_table(:,2)==0);
            power_tx(is_0_tx) = P_tx_LOS ;
            power_tx(~is_0_tx) = obj.Pn_keep(1,obj.coupling_txrx_path_table(~is_0_tx,1))/(KR_tx+1)/obj.M(1);

            is_0_rx = (obj.coupling_txrx_path_table(:,3)==0 & obj.coupling_txrx_path_table(:,4)==0);
            power_rx(is_0_rx) = P_rx_LOS;
            power_rx(~is_0_rx) = obj.Pn_keep(2,obj.coupling_txrx_path_table(~is_0_rx,3))/(KR_rx+1)/obj.M(2);
            if isa(obj.ST,"elements.Target")
                RCS_sigma_S_path = 10.^((randn(obj.path_num,1)*obj.ST.type.RCS.sigma_sigma_S_dB+obj.ST.type.RCS.mu_sigma_S_dB)/10);
                % if obj.ST.type.RCS_model == 2
                %     sigma_MD_dB = obj.RCS_MD_dB;
                %     sigma_D_dB = sigma_MD_dB - obj.ST.type.RCS.sigma_M_dB;
                %     sigma_D = db2pow(sigma_D_dB);
                % else 
                %     sigma_D = obj.ST.type.RCS.sigma_D;
                % end
                power_path = obj.ST.type.RCS.sigma_D * min(RCS_sigma_S_path,10^((obj.ST.type.RCS.mu_sigma_S_dB+3*obj.ST.type.RCS.sigma_sigma_S_dB)/10)).*power_rx.*power_tx;
                obj.RCS_sigma_S_n_m = RCS_sigma_S_path;
            end

            obj.P_path = power_path;
            keep_tmp = true(obj.path_num,1);
            if obj.option == 1
                keep_tmp(~(is_0_tx & is_0_rx)) = (10*log10(obj.P_path(~(is_0_tx & is_0_rx))/max(obj.P_path)) >= -25);
            else
                keep_tmp(~(is_0_tx & is_0_rx))  = (10*log10(obj.P_path(~(is_0_tx & is_0_rx))/max(obj.P_path)) >= -40);
            end
            keep_tmp(is_0_tx & is_0_rx) = true;

            obj.keep_path = keep_tmp;
            obj.P_path_keep  = obj.P_path(keep_tmp);
            obj.coupling_txrx_path_table_keep = obj.coupling_txrx_path_table(keep_tmp,:);
            obj.path_num_keep = size(obj.coupling_txrx_path_table_keep,1);

        end

        function Absolute_Delay(obj)
            tau_tx = zeros(obj.path_num_keep,1);
            tau_rx = zeros(obj.path_num_keep,1);
            % for sub-cluster 2
            ray_idlist_2 = [9, 10, 11, 12, 17, 18];
            % for sub-cluster 3
            ray_idlist_3 = [13, 14, 15, 16];

            isnot_0_tx = find((obj.coupling_txrx_path_table_keep(:,1) ~= 0)&(obj.coupling_txrx_path_table_keep(:,2) ~= 0));
            tau_tx(isnot_0_tx) = obj.tau_n_keep(1,obj.coupling_txrx_path_table_keep(isnot_0_tx, 1));
            idx_n = isnot_0_tx( ismember(obj.tau_n_keep(1, obj.coupling_txrx_path_table_keep(isnot_0_tx, 1)),obj.strong_cluster_id(1,:) ) );
            idx_m_2 = idx_n(ismember(obj.coupling_txrx_path_table_keep(idx_n,2) , ray_idlist_2));
            tau_tx(idx_m_2) = tau_tx(idx_m_2) + obj.map_delay(obj.coupling_txrx_path_table_keep(idx_m_2,2));
            idx_m_3 = idx_n(ismember(obj.coupling_txrx_path_table_keep(idx_n,2) , ray_idlist_3));
            tau_tx(idx_m_3) = tau_tx(idx_m_3) + obj.map_delay(obj.coupling_txrx_path_table_keep(idx_m_3,2));

            isnot_0_rx = find((obj.coupling_txrx_path_table_keep(:,3) ~= 0)&(obj.coupling_txrx_path_table_keep(:,4) ~= 0));
            tau_rx(isnot_0_rx) = obj.tau_n_keep(2,obj.coupling_txrx_path_table_keep(isnot_0_rx, 3));
            idx_n = isnot_0_rx( ismember(obj.tau_n_keep(2, obj.coupling_txrx_path_table_keep(isnot_0_rx, 3)),obj.strong_cluster_id(2,:) ) );
            idx_m_2 = idx_n(ismember(obj.coupling_txrx_path_table_keep(idx_n,4) , ray_idlist_2));
            tau_rx(idx_m_2) = tau_rx(idx_m_2) + obj.map_delay(obj.coupling_txrx_path_table_keep(idx_m_2,4));
            idx_m_3 = idx_n(ismember(obj.coupling_txrx_path_table_keep(idx_n,4) , ray_idlist_3));
            tau_rx(idx_m_3) = tau_rx(idx_m_3) + obj.map_delay(obj.coupling_txrx_path_table_keep(idx_m_3,4));
            
            switch obj.scenario.name
                case {'UMa','UrbanGrid'}
                obj.mu_lg_delta_tau = -7.4;
                obj.sigma_lg_delta_tau = 0.2;
                case 'UMi'
                obj.mu_lg_delta_tau = -7.5;
                obj.sigma_lg_delta_tau = 0.5;
                case 'InH'
                obj.mu_lg_delta_tau = -8.6;
                obj.sigma_lg_delta_tau = 0.1;
                case 'InF'
                obj.mu_lg_delta_tau = -7.5;
                obj.sigma_lg_delta_tau = 0.4;
            end

            delta_tx = zeros(obj.path_num_keep,1);
            delta_rx = zeros(obj.path_num_keep,1);
            % delta_tx(isnot_0_tx) = 10.^(randn(length(isnot_0_tx),1)*obj.sigma_lg_delta_tau + obj.mu_lg_delta_tau);
            delta_tx(isnot_0_tx) = 10.^(randn*obj.sigma_lg_delta_tau + obj.mu_lg_delta_tau);
            if obj.static == 1
                delta_rx = delta_tx;
            else
                % delta_rx(isnot_0_rx) = 10.^(randn(length(isnot_0_rx),1)*obj.sigma_lg_delta_tau + obj.mu_lg_delta_tau);
                delta_rx(isnot_0_rx) = 10.^(randn*obj.sigma_lg_delta_tau + obj.mu_lg_delta_tau);
            end

            tau_path_ = tau_rx + obj.d_3D(2)/3e8 + delta_rx + tau_tx + obj.d_3D(1)/3e8 + delta_tx;
            obj.path_delay = tau_path_;
        end
        % step 12: Generate XPRs
        function generate_XPRs(obj)
            % The two target-channel legs can have different numbers of
            % retained clusters. Allocate the full storage first so that
            % XPR generation does not depend on the first leg's size.
            max_N_new = max(obj.N_new);
            max_M = max(obj.M);
            obj.XPR_txrx = nan(2, max_N_new, max_M);

            for txrx = 1:2
                xpr_randn = randn(obj.N_new(txrx), obj.M(txrx));

                X_n_m = obj.sigma_XPR(txrx)*xpr_randn+obj.mu_XPR(txrx);
                obj.XPR_txrx(txrx, 1:obj.N_new(txrx), 1:obj.M(txrx)) = ...
                    10.^(X_n_m/10);

                if obj.static ==1
                    obj.XPR_txrx(2, 1:obj.N_new(2), 1:obj.M(2)) = ...
                        obj.XPR_txrx(1, 1:obj.N_new(1), 1:obj.M(1));
                    break;
                end
            end
            xpr_randn = randn(obj.path_num_keep,1);
            X_path_SPST = obj.ST.type.XPR.sigma_sensing *xpr_randn+obj.ST.type.XPR.mu_sensing;
            obj.XPR_spst_path = 10.^(X_path_SPST/10);           
        end

        % step 13: Draw initial random phases
         function initial_random_phases(obj)
             % Allocate phase storage for the largest target-channel leg.
             max_N_new = max(obj.N_new);
             max_M = max(obj.M);
             obj.PHI_n_m_txrx = nan(2, 2, 2, max_N_new, max_M);

             for txrx = 1:2
                 phase_tmp = (2*rand([2, 2, obj.N_new(txrx), obj.M(txrx)])-1)*pi;
                 obj.PHI_n_m_txrx(txrx, 1:2, 1:2, 1:obj.N_new(txrx), 1:obj.M(txrx)) = ...
                     reshape(phase_tmp, [1, 2, 2, obj.N_new(txrx), obj.M(txrx)]);
                 if obj.static == 1
                     obj.PHI_n_m_txrx(2, :, :, 1:obj.N_new(2), 1:obj.M(2)) = ...
                         permute(obj.PHI_n_m_txrx(1, :, :, 1:obj.N_new(1), 1:obj.M(1)), [1, 3, 2, 4, 5]);
                     break;
                 end
             end
             obj.PHI_spst_path(:,:,:) = (2*rand([2, 2, obj.path_num_keep])-1)*pi;
             return;

             for txrx = 1:2
                 obj.PHI_n_m_txrx(txrx,:,:,:,:) = (2*rand([2, 2, obj.N_new(txrx), obj.M(txrx)])-1)*pi;%（θθ, θϕ, ϕθ, ϕϕ）
                 % end
                 if obj.static ==1
                     obj.PHI_n_m_txrx(2,:,:,:,:) = permute(obj.PHI_n_m_txrx(1,:,:,:,:),[1,3,2,4,5]);
                     break;
                 end
             end
             obj.PHI_spst_path(:,:,:) = (2*rand([2, 2, obj.path_num_keep])-1)*pi;
         end

         % step 13.5: calulate part of step 14 parameter which is independent with antenna element
         function prepare_channel_parameter(obj)
             CPM_tx_nm = complex(zeros(2,2,obj.N_new(1),obj.M(1)));
             CPM_rx_nm = complex(zeros(2,2,obj.N_new(2),obj.M(2)));
             CPM_tx_nm(1,1,:,:) = reshape(exp(1j*obj.PHI_n_m_txrx(1,1,1,1:obj.N_new(1),1:obj.M(1))),[1,1,obj.N_new(1),obj.M(1)]);
            CPM_tx_nm(1,2,:,:) = reshape(sqrt(1./obj.XPR_txrx(1,1:obj.N_new(1),1:obj.M(1))),[1,1,obj.N_new(1),obj.M(1)]).*reshape(exp(1j*obj.PHI_n_m_txrx(1,1,2,1:obj.N_new(1),1:obj.M(1))),[1,1,obj.N_new(1),obj.M(1)]);
            CPM_tx_nm(2,1,:,:) = reshape(sqrt(1./obj.XPR_txrx(1,1:obj.N_new(1),1:obj.M(1))),[1,1,obj.N_new(1),obj.M(1)]).*reshape(exp(1j*obj.PHI_n_m_txrx(1,2,1,1:obj.N_new(1),1:obj.M(1))),[1,1,obj.N_new(1),obj.M(1)]);
             CPM_tx_nm(2,2,:,:) = reshape(exp(1j*obj.PHI_n_m_txrx(1,2,2,1:obj.N_new(1),1:obj.M(1))),[1,1,obj.N_new(1),obj.M(1)]);

             CPM_rx_nm(1,1,:,:) = reshape(exp(1j*obj.PHI_n_m_txrx(2,1,1,1:obj.N_new(2),1:obj.M(2))),[1,1,obj.N_new(2),obj.M(2)]);
            CPM_rx_nm(1,2,:,:) = reshape(sqrt(1./obj.XPR_txrx(2,1:obj.N_new(2),1:obj.M(2))),[1,1,obj.N_new(2),obj.M(2)]).*reshape(exp(1j*obj.PHI_n_m_txrx(2,1,2,1:obj.N_new(2),1:obj.M(2))),[1,1,obj.N_new(2),obj.M(2)]);
            CPM_rx_nm(2,1,:,:) = reshape(sqrt(1./obj.XPR_txrx(2,1:obj.N_new(2),1:obj.M(2))),[1,1,obj.N_new(2),obj.M(2)]).*reshape(exp(1j*obj.PHI_n_m_txrx(2,2,1,1:obj.N_new(2),1:obj.M(2))),[1,1,obj.N_new(2),obj.M(2)]);
             CPM_rx_nm(2,2,:,:) = reshape(exp(1j*obj.PHI_n_m_txrx(2,2,2,1:obj.N_new(2),1:obj.M(2))),[1,1,obj.N_new(2),obj.M(2)]);
             los_pol_coupling_mat = [1,0;0,-1];
             CPM_spst(1,1,:) = exp(1j*obj.PHI_spst_path(1,1,:));
             CPM_spst(1,2,:) = sqrt(1./(obj.XPR_spst_path)).*squeeze(exp(1j*obj.PHI_spst_path(1,2,:)));
             CPM_spst(2,1,:) = sqrt(1./(obj.XPR_spst_path)).*squeeze(exp(1j*obj.PHI_spst_path(2,1,:)));
             CPM_spst(2,2,:) = exp(1j*obj.PHI_spst_path(2,2,:));
             % ensure CPM_(n m n' m') = CPM_(n' m' n m)^T
             coupling_txrx_n_m_table_keep_swap = obj.coupling_txrx_path_table_keep(:, [3 4 1 2]);
             [tf_step_13_5, loc_step_13_5] = ismember(coupling_txrx_n_m_table_keep_swap, obj.coupling_txrx_path_table_keep, 'rows');
             idx_a_step_13_5 = find(tf_step_13_5);
             idx_b_step_13_5 = loc_step_13_5(idx_a_step_13_5);
             mask_step_13_5 = idx_a_step_13_5 < idx_b_step_13_5;
             idx_a_step_13_5 = idx_a_step_13_5(mask_step_13_5);
             idx_b_step_13_5 = idx_b_step_13_5(mask_step_13_5);
             CPM_spst(:,:,idx_a_step_13_5) = permute(CPM_spst(:,:,idx_b_step_13_5), [2 1 3]);

             % ========= 2) Build per-path CPM_tx_page and CPM_rx_page (each 2×2×A) =========
             % Flatten 2×2×N×M -> 2×2×(N*M) so we can over_half_circle with sub2ind
             N_tx_step_13_5 = size(CPM_tx_nm,3);
             M_tx_step_13_5 = size(CPM_tx_nm,4);
             N_rx_step_13_5 = size(CPM_rx_nm,3);
             M_rx_step_13_5 = size(CPM_rx_nm,4);
             CPM_tx_flat = reshape(CPM_tx_nm, 2, 2, []);
             CPM_rx_flat = reshape(CPM_rx_nm, 2, 2, []);

             % Indices for (n,m) and (n',m') in table
             nm_step_13_5  = obj.coupling_txrx_path_table_keep(:, 1:2);
             npm_step_13_5 = obj.coupling_txrx_path_table_keep(:, 3:4);

             is0_tx_step_13_5 = (nm_step_13_5(:,1)==0 & nm_step_13_5(:,2)==0);
             is0_rx_step_13_5 = (npm_step_13_5(:,1)==0 & npm_step_13_5(:,2)==0);

             % Initialize pages with LOS (so [0,0] cases are already handled)
             CPM_tx_page = repmat(los_pol_coupling_mat, 1, 1, obj.path_num_keep);   % 2×2×A
             CPM_rx_page = repmat(los_pol_coupling_mat, 1, 1, obj.path_num_keep);   % 2×2×A

             % Fill non-zero indices from CPM_tx / CPM_rx
             idx_tx_step_13_5 = sub2ind([N_tx_step_13_5 M_tx_step_13_5], nm_step_13_5(~is0_tx_step_13_5,1),  nm_step_13_5(~is0_tx_step_13_5,2));
             idx_rx_step_13_5 = sub2ind([N_rx_step_13_5 M_rx_step_13_5], npm_step_13_5(~is0_rx_step_13_5,1), npm_step_13_5(~is0_rx_step_13_5,2));

             CPM_tx_page(:,:,~is0_tx_step_13_5) = CPM_tx_flat(:,:,idx_tx_step_13_5);
             CPM_rx_page(:,:,~is0_rx_step_13_5) = CPM_rx_flat(:,:,idx_rx_step_13_5);

             % ========= 3) Compute CPM_all(:,:,a) = CPM_rx_page(:,:,a) * CPM_swap(:,:,a) * CPM_tx_page(:,:,a) =========
             % Uses page-wise matrix multiply (R2020b+)
             CPM_all_ = pagemtimes(CPM_rx_page, pagemtimes(CPM_spst, CPM_tx_page));  % 2×2×A
             d_AA_step_13_5 = squeeze(CPM_all_(1,1,:));   % A×1
             d_TT_step_13_5 = squeeze(CPM_all_(2,2,:));   % A×1
             normFactor_step_13_5 = sqrt( (abs(d_AA_step_13_5).^2 + abs(d_TT_step_13_5).^2) / 2 );   % A×1
             normFactor_step_13_5(normFactor_step_13_5 == 0) = 1;
             CPM_all_ = CPM_all_ ./ reshape(normFactor_step_13_5, 1, 1, []);

             phi_path_AOA_ = repmat(obj.phi_LOS_AOA, obj.path_num_keep, 1)';                 % default a for (0,0)
             phi_path_AOD_ = repmat(obj.phi_LOS_AOD, obj.path_num_keep, 1)';
             theta_path_ZOA_ = repmat(obj.theta_LOS_ZOA, obj.path_num_keep, 1)';
             theta_path_ZOD_ = repmat(obj.theta_LOS_ZOD, obj.path_num_keep, 1)';
             is0_tx_step_13_5 = (obj.coupling_txrx_path_table_keep(:,1)==0 & obj.coupling_txrx_path_table_keep(:,2)==0);      % rows where (n,m) == (0,0)
             is0_rx_step_13_5 = (obj.coupling_txrx_path_table_keep(:,3)==0 & obj.coupling_txrx_path_table_keep(:,4)==0);      % rows where (n,m) == (0,0)
             idx_tx_step_13_5  = sub2ind(size(obj.phi_n_m_AOA(1,:,:),2,3), obj.coupling_txrx_path_table_keep(~is0_tx_step_13_5,1),  obj.coupling_txrx_path_table_keep(~is0_tx_step_13_5,2));
             idx_rx_step_13_5  = sub2ind(size(obj.phi_n_m_AOA(2,:,:),2,3), obj.coupling_txrx_path_table_keep(~is0_rx_step_13_5,3),  obj.coupling_txrx_path_table_keep(~is0_rx_step_13_5,4));
             phi_path_AOA_(1,~is0_tx_step_13_5) = obj.phi_n_m_AOA(1,idx_tx_step_13_5);
             phi_path_AOD_(1,~is0_tx_step_13_5) = obj.phi_n_m_AOD(1,idx_tx_step_13_5);
             theta_path_ZOA_(1,~is0_tx_step_13_5) = obj.theta_n_m_ZOA(1,idx_tx_step_13_5);
             theta_path_ZOD_(1,~is0_tx_step_13_5) = obj.theta_n_m_ZOD(1,idx_tx_step_13_5);
             phi_path_AOA_(2,~is0_rx_step_13_5) = obj.phi_n_m_AOA(2,idx_rx_step_13_5);
             phi_path_AOD_(2,~is0_rx_step_13_5) = obj.phi_n_m_AOD(2,idx_rx_step_13_5);
             theta_path_ZOA_(2,~is0_rx_step_13_5) = obj.theta_n_m_ZOA(2,idx_rx_step_13_5);
             theta_path_ZOD_(2,~is0_rx_step_13_5) = obj.theta_n_m_ZOD(2,idx_rx_step_13_5);

             r_tx_path_(1,:,:) = sind(theta_path_ZOD_(1,:)).*cosd(phi_path_AOD_(1,:));
             r_tx_path_(2,:,:) = sind(theta_path_ZOD_(1,:)).*sind(phi_path_AOD_(1,:));
             r_tx_path_(3,:,:) = cosd(theta_path_ZOD_(1,:));

             r_rx_path_(1,:,:) = sind(theta_path_ZOA_(2,:)).*cosd(phi_path_AOA_(2,:)); % the spherical unit vector
             r_rx_path_(2,:,:) = sind(theta_path_ZOA_(2,:)).*sind(phi_path_AOA_(2,:));
             r_rx_path_(3,:,:) = cosd(theta_path_ZOA_(2,:));

             r_tx_spst_path_(1,:,:) = sind(theta_path_ZOA_(1,:)).*cosd(phi_path_AOA_(1,:)); % the spherical unit vector
             r_tx_spst_path_(2,:,:) = sind(theta_path_ZOA_(1,:)).*sind(phi_path_AOA_(1,:));
             r_tx_spst_path_(3,:,:) = cosd(theta_path_ZOA_(1,:));


             r_spst_rx_path_(1,:,:) = sind(theta_path_ZOD_(2,:)).*cosd(phi_path_AOD_(2,:));
             r_spst_rx_path_(2,:,:) = sind(theta_path_ZOD_(2,:)).*sind(phi_path_AOD_(2,:));
             r_spst_rx_path_(3,:,:) = cosd(theta_path_ZOD_(2,:));

             obj.CPM_all = CPM_all_;
             obj.r_tx_path = r_tx_path_;
             obj.r_rx_path = r_rx_path_;
             obj.r_tx_spst_path = r_tx_spst_path_;
             obj.r_spst_rx_path = r_spst_rx_path_;
             obj.phi_path_AOA = phi_path_AOA_;
             obj.phi_path_AOD = phi_path_AOD_;
             obj.theta_path_ZOA = theta_path_ZOA_;
             obj.theta_path_ZOD = theta_path_ZOD_;
         end

         % step 14
        
         function generate_channel(obj)
            lambda      = 3e8/obj.fc;
            if obj.fastfading_enable
                t0 = 0;
                if isempty(obj.sector(1).velocity)
                    v_tx = [0, 0, 0];
                else
                    v_tx = obj.sector(1).velocity * [cosd(obj.sector(1).phi_v),sind(obj.sector(1).phi_v),cosd(obj.sector(1).theta_v)];
                end
                if isempty(obj.sector(2).velocity)
                    v_rx = [0, 0, 0];
                else
                    v_rx = obj.sector(2).velocity * [cosd(obj.sector(2).phi_v),sind(obj.sector(2).phi_v),cosd(obj.sector(2).theta_v)];
                end
                v_spst = obj.ST.velocity * [cosd(obj.ST.phi_v),sind(obj.ST.phi_v),cosd(obj.ST.theta_v)]; % leak v_mi micro motion of each sp


                obj.sector(1).PHI_n_m= obj.PHI_n_m_txrx(1,:,:,:,:);
                obj.sector(2).PHI_n_m= obj.PHI_n_m_txrx(2,:,:,:,:);

                tx_panel_num = obj.sector(1).sector.antenna.num_panel;
                Ptx = obj.sector(1).antenna_params.panel.P;
                E_tx = obj.sector(1).sector.antenna.panel{1, 1}.num_element/Ptx;
                S = tx_panel_num*E_tx;

                rx_panel_num = obj.sector(2).sector.antenna.num_panel;
                Prx = obj.sector(2).antenna_params.panel.P;
                E_rx = obj.sector(2).sector.antenna.panel{1, 1}.num_element/Prx;
                U = rx_panel_num*E_rx;

                %% BS params
                pos_panelTX = reshape(permute(obj.sector(1).sector.antenna.pos_panel_LCS,[3,1,2]),3,[]);
                pos_panelRX = reshape(permute(obj.sector(2).sector.antenna.pos_panel_LCS,[3,1,2]),3,[]);

                dTX = zeros(3, S);
                for s = 1:S
                    panS = floor((s-1)/E_tx)+1;
                    idxS = mod((s-1),E_tx)+1;
                    tx_pos = reshape(permute(obj.sector(1).sector.antenna.panel{panS}.pos_element_LCS,[3,1,2]),3,[]);
                    dTX(:,s) = obj.sector(1).sector.antenna.R * (pos_panelTX(:,panS) + tx_pos(:,idxS)) * lambda;
                end


                dRX = zeros(3, U);
                for u = 1:U
                    panU = floor((u-1)/E_rx)+1;
                    idxU = mod((u-1),E_rx)+1;
                    rx_pos = reshape(permute(obj.sector(2).sector.antenna.panel{panU}.pos_element_LCS,[3,1,2]),3,[]);
                    dRX(:,u) = obj.sector(2).sector.antenna.R * (pos_panelRX(:,panU) + rx_pos(:,idxU)) * lambda;
                end

                % phase_tx(l,s) = exp(j*2*pi/lambda * r_tx(:,l)^T dTX(:,s))
                phase_tx = exp(1j*2*pi/lambda * squeeze(sum(dTX.*obj.r_tx_path,1)));   % (L x S)
                phase_rx = exp(1j*2*pi/lambda * squeeze(sum(dRX.*obj.r_rx_path,1)));    % (L x U)

                doppler_tx = (squeeze(sum(v_tx(:).*obj.r_tx_path,1))+squeeze(sum(v_spst(:).*obj.r_tx_spst_path,1)))/lambda;
                doppler_rx = (squeeze(sum(v_rx(:).*obj.r_rx_path,1))+squeeze(sum(v_spst(:).*obj.r_spst_rx_path,1)))/lambda;

                doppler_total_ = doppler_tx + doppler_rx;
                phase_doppler_ = exp(1j*2*pi*doppler_total_*(obj.t-t0)); % leak integrted

                tx_field_pattern = complex(zeros(2, S, obj.path_num_keep, Ptx));
                rx_field_pattern = complex(zeros(2, U, obj.path_num_keep, Prx));


                for ptx = 1:Ptx
                    % field_pattern should accept vector angles and return 2 x L (or 2 x 1 x L)
                    element_field_pattern = obj.sector(1).sector.antenna.panel{1}.element_list(ptx).field_pattern(obj.phi_path_AOD(1,:), obj.theta_path_ZOD(1,:)); % expect 2 x L
                    element_field_pattern = reshape(element_field_pattern, 2, 1, obj.path_num_keep);  % 2 x 1 x L
                    tx_field_pattern(:,:,:,ptx) = element_field_pattern;
                end

                for prx = 1:Prx
                    element_field_pattern = obj.sector(2).sector.antenna.panel{1}.element_list(prx).field_pattern(obj.phi_path_AOA(2,:), obj.theta_path_ZOA(2,:)); % 2 x L
                    element_field_pattern = reshape(element_field_pattern, 2, 1, obj.path_num_keep);
                    rx_field_pattern(:,:,:,prx) = element_field_pattern;
                end
                polarized_channel = complex(zeros(U, S, obj.path_num_keep, Prx, Ptx));

                for prx = 1:Prx
                    Frx_u2L = permute(rx_field_pattern(:,:,:,prx), [2 1 3]);   % 2 x U x L
                    for ptx = 1:Ptx
                        Ftx_2sL = tx_field_pattern(:,:,:,ptx);                 % 2 x S x L
                        polarized_channel(:,:,:,prx,ptx) =  pagemtimes(Frx_u2L, pagemtimes(obj.CPM_all, Ftx_2sL));  % U x S x L
                    end
                end
                %
                ray_power_sqrt = sqrt(obj.P_path_keep);
                ray_power_sqrt = reshape(ray_power_sqrt, 1,1,obj.path_num_keep);
                phase_doppler_ = reshape(phase_doppler_, 1,1,obj.path_num_keep);

                % ===== "P part same as original": noncoherent power sum over pols =====
                calibration_power = sum(abs(polarized_channel(:,:,:,1,1).* ray_power_sqrt).^2,"all");   % U x S
                obj.H_fastfading = sum(abs(polarized_channel.* ray_power_sqrt .*phase_doppler_.*reshape(phase_tx,1,1,[]).*reshape(phase_rx,1,1,[])).^2, [3 4 5]);   % U x S                
                obj.H_cali = pow2db(calibration_power);
                obj.RCS_sigma_M_dB = obj.ST.type.RCS.sigma_M_dB;
                obj.H_full = obj.H_fastfading * db2pow(-(sum(obj.PL) +  10*log10(lambda^2 /(4*pi)) -obj.RCS_MD_dB  + sum(obj.SF)));

                obj.Ftx = tx_field_pattern;
                obj.Frx = rx_field_pattern;
            else
                obj.H_full = db2pow(-(sum(obj.PL) +  10*log10(lambda^2 /(4*pi)) -obj.RCS_MD_dB  + sum(obj.SF)));
            end
         end

         function [couplingloss_ls,couplingloss_full, pathloss] = calc_coupling_loss(obj)
             lambda      = 3e8/obj.fc;
             couplingloss_ls = -(sum(obj.PL) +  10*log10(lambda^2 /(4*pi)) -obj.ST.type.RCS.sigma_M_dB + sum(obj.SF));
             pathloss = sum(obj.PL) + sum(obj.SF);
             if obj.fastfading_enable
                 couplingloss_full =-(sum(obj.PL) +  10*log10(lambda^2 /(4*pi)) -obj.ST.type.RCS.sigma_M_dB  + sum(obj.SF)) + obj.H_cali;
             else
                 couplingloss_full = nan;
             end
         end
       

        function [ as, mean_angle ]  = calc_angular_spreads(obj,ang,wrap_angles )
            % CALC_ANGULAR_SPREADS Calculates the angular spread in degree

            ang = ang(:)*pi/180;
            % Normalize powers
            pt = sum(obj.P_path_keep);
            pow = obj.P_path_keep./pt;
            if ~exist('wrap_angles','var')
                wrap_angles = true;
            end
            if wrap_angles
                mean_angle = angle( sum( pow.*exp( 1j*ang )) ); % [rad]
            else
                mean_angle = sum( pow.*ang);
            end
            phi = ang - mean_angle;
            if wrap_angles
                phi = angle( exp( 1j*phi ) );
            end
            %             as1 = sqrt(-2*log(abs(sum( pow.*exp( 1j*ang ),2))));
            as = sqrt( sum(pow.*(phi.^2))); %  - sum( pow.*phi,2).^2
            %             as  = min(as1,as2);
            mean_angle = mean_angle*180/pi;
            as = as*180/pi;
        end

        function [ ds, mean_delay ] = calc_delay_spread(obj)
            % CALC_DELAY_SPREAD Calculates the delay spread in [s]

            % Normalize powers
            pt = sum(obj.P_path_keep);
            pow = obj.P_path_keep./pt;
            taus = obj.path_delay;

            mean_delay = sum( pow.*taus );
            sort_idx = taus - mean_delay;
            ds = sqrt( sum(pow.*(sort_idx.^2))); %  - sum(pow.*sort_idx).^2
            % ds = sqrt( sum(pow.*((taus).^2)) - mean_delay.^2 ); %  - sum(pow.*sort_idx).^2
        end

        function sigma_MD_dB = RCS_MD_dB(obj)

            % ----- direction vectors -----
            vi = obj.ST.Position - squeeze(obj.BS_pos_wrap(1,:));     % incident
            vs = obj.ST.Position - squeeze(obj.BS_pos_wrap(2,:));     % scattering

            % ----- azimuth (degree) -----
            phi_i = atan2d(vi(2), vi(1));
            phi_s = atan2d(vs(2), vs(1));

            phi = abs(phi_i - phi_s);
            % phi = min(phi, 360 - phi);     % limit to [0,180]

            % ----- elevation (degree) -----
            theta_i = atan2d(vi(3), hypot(vi(1),vi(2)));
            theta_s = atan2d(vs(3), hypot(vs(1),vs(2)));

            theta = abs(theta_i - theta_s);
            theta = min(theta, 360 - theta);   % limit to [0,180]

            % ----- full 3D angle -----
            cosA = dot(vi,vs) / (norm(vi)*norm(vs));
            cosA = max(-1, min(1, cosA));   % clamp to [-1,1]
            beta = acosd(cosA);            % automatically [0,180]


            if obj.ST.type.RCS_model  == 1
                sigma_MD_dB = max(10*log10(obj.ST.type.RCS.sigma_M)-3*sind(beta/2),-Inf);
            elseif obj.ST.type.RCS_model  == 2
                if obj.ST.is_single_STSP
                    if (theta >= obj.ST.type.sigle_STSP.theta_range(obj.ST.SP,1)) && (theta <= obj.ST.type.sigle_STSP.theta_range(obj.ST.SP,2)) && ~isnan(obj.ST.type.sigle_STSP.theta_center(obj.ST.SP))
                        sigme_V_dB = -min(12*((theta-obj.ST.type.sigle_STSP.theta_center(obj.ST.SP))/obj.ST.type.sigle_STSP.theta_3dB(obj.ST.SP))^2,obj.ST.type.sigle_STSP.sigma_max(obj.ST.SP));
                    else
                        sigme_V_dB = 0;
                    end
                    if (phi >= obj.ST.type.sigle_STSP.phi_range(obj.ST.SP,1)) && (phi <= obj.ST.type.sigle_STSP.phi_range(obj.ST.SP,2)) && ~isnan(obj.ST.type.sigle_STSP.phi_center(obj.ST.SP))
                        sigme_H_dB = -min(12*((phi-obj.ST.type.sigle_STSP.phi_center(obj.ST.SP))/obj.ST.type.sigle_STSP.phi_3dB(obj.ST.SP))^2,obj.ST.type.sigle_STSP.sigma_max(obj.ST.SP));
                    else
                        sigme_H_dB = 0;
                    end
                    sigma_MD_dB = max([obj.ST.type.sigle_STSP.G_max(obj.ST.SP)-min(-(sigme_V_dB+sigme_H_dB),obj.ST.type.sigle_STSP.sigma_max(obj.ST.SP))-obj.ST.type.k1*sin(obj.ST.type.k2*beta/2)+5*log10(cos(beta/2)),...
                        obj.ST.type.sigle_STSP.G_max(obj.ST.SP)-obj.ST.type.sigle_STSP.sigma_max(obj.ST.SP),-Inf],[],"all");
                else
                    if all([theta >= obj.ST.type.multi_STSP.theta_range(obj.ST.SP,1),theta <= obj.ST.type.multi_STSP.theta_range(obj.ST.SP,2),~isnan(obj.ST.type.multi_STSP.theta_center(obj.ST.SP))])
                        sigme_V_dB = -min(12*((theta-obj.ST.type.multi_STSP.theta_center(obj.ST.SP))/obj.ST.type.multi_STSP.theta_3dB(obj.ST.SP))^2,obj.ST.type.multi_STSP.sigma_max(obj.ST.SP));
                    else
                        sigme_V_dB = 0;
                    end
                    if all([phi >= obj.ST.type.multi_STSP.phi_range(obj.ST.SP,1),phi <= obj.ST.type.multi_STSP.phi_range(obj.ST.SP,2) , ~isnan(obj.ST.type.multi_STSP.phi_center(obj.ST.SP))])
                        sigme_H_dB = -min(12*((phi-obj.ST.type.multi_STSP.phi_center(obj.ST.SP))/obj.ST.type.multi_STSP.phi_3dB(obj.ST.SP))^2,obj.ST.type.multi_STSP.sigma_max(obj.ST.SP));
                    else
                        sigme_H_dB = 0;
                    end
                    sigma_MD_dB = max([obj.ST.type.multi_STSP.G_max(obj.ST.SP)-min(-(sigme_V_dB+sigme_H_dB),obj.ST.type.multi_STSP.sigma_max(obj.ST.SP))-obj.ST.type.k1*sin(obj.ST.type.k2*beta/2)+5*log10(cos(beta/2)),...
                        obj.ST.type.multi_STSP.G_max(obj.ST.SP)-obj.ST.type.multi_STSP.sigma_max(obj.ST.SP),-Inf],[],"all");
                end
            end

        end

    end

    methods(Access = private)
        function assertProcedureBSupported(obj)
            if ~ismember(obj.scenario.name, {'UrbanGrid', 'UMi', 'InH'})
                error('TargetChannel:UnsupportedSpatialConsistencyProcedureB', ...
                    'Spatial consistency Procedure B is currently implemented only for UrbanGrid, UMi and InH scenarios.');
            end
        end

        function assertProcedureASupported(obj)
            if ~ismember(obj.scenario.name, {'UrbanGrid', 'UMi'})
                error('TargetChannel:UnsupportedSpatialConsistencyProcedureA', ...
                    'Spatial consistency Procedure A is currently implemented only for UrbanGrid and UMi scenarios.');
            end
        end

        function limits = procedureBClusterLimits(obj, txrx)
            fc_GHz = obj.fc/1e9;
            if strcmp(obj.scenario.name, 'UrbanGrid')
                fc_eval = max(fc_GHz, 6);
                if obj.O2I
                    mu_lgDS = -6.62;
                    sigma_lgDS = 0.32;
                    mu_lgASD = 0.58;
                    sigma_lgASD = 0.7;
                    mu_lgASA = 1.76;
                    sigma_lgASA = 0.16;
                    mu_lgZSA = 1.01;
                    sigma_lgZSA = 0.43;
                elseif obj.bLOS(txrx)
                    mu_lgDS = -7.067 - 0.0794*log10(1 + fc_eval);
                    sigma_lgDS = 0.57 + 0.026*log10(fc_eval);
                    mu_lgASD = 0.92;
                    sigma_lgASD = 0.31;
                    mu_lgASA = 1.76;
                    sigma_lgASA = 0.19;
                    mu_lgZSA = 0.96;
                    sigma_lgZSA = 0.15;
                else
                    mu_lgDS = -6.47 - 0.134*log10(fc_eval);
                    sigma_lgDS = 0.39;
                    mu_lgASD = 1.09;
                    sigma_lgASD = 0.44;
                    mu_lgASA = 2.04 - 0.25*log10(fc_eval);
                    sigma_lgASA = 0.17 - 0.03*log10(fc_eval);
                    mu_lgZSA = -0.2856*log10(fc_eval) + 1.445;
                    sigma_lgZSA = 0.17;
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
                elseif obj.bLOS(txrx)
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
                if obj.bLOS(txrx)
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
            limits.ZOD = 2*10^(obj.mu_lgZSD(txrx) + obj.sigma_lgZSD(txrx));
        end

        function uniform_values = drawProcedureBUniform(obj, txrx, field_name, count)
            uniform_values = rand(1, count);
            if ~obj.isSpatialConsistencyProcedureB()
                return;
            end
            if txrx == 1
                return;
            end

            [external_values, isAvailable] = obj.getExternalProcedureBRaw(txrx, field_name, count);
            if isAvailable
                uniform_values = external_values;
            else
                obj.warnSpatialConsistencyFallback();
            end
        end

        function gaussian_values = drawProcedureBGaussian(obj, txrx, field_name, count)
            gaussian_values = randn(1, count);
            if ~obj.isSpatialConsistencyProcedureB()
                return;
            end
            if txrx == 1
                return;
            end

            [external_values, isAvailable] = obj.getExternalProcedureBRaw(txrx, field_name, count);
            if isAvailable
                gaussian_values = external_values;
            else
                obj.warnSpatialConsistencyFallback();
            end
        end

        function [values, isAvailable] = getExternalProcedureBRaw(obj, txrx, field_name, count)
            values = [];
            isAvailable = false;

            if txrx > numel(obj.sector)
                return;
            end

            equipment = obj.sector(txrx);
            if ~isprop(equipment, 'SC_procB_raw') || isempty(equipment.SC_procB_raw)
                return;
            end

            st_id = obj.ST.ID;
            if st_id < 1 || numel(equipment.SC_procB_raw) < st_id
                return;
            end

            proc_raw = equipment.SC_procB_raw(st_id);
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

        function raw_rand_vec = drawLspRawRandn(obj, vectorLength, txrx)
            raw_rand_vec = randn(vectorLength, 1);
            if ~obj.isSpatialConsistencyEnabled()
                return;
            end

            [external_lsp, isAvailable] = obj.getExternalLspRaw(vectorLength, txrx);
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
            tf = true;
            max_txrx = obj.TXRX;
            if obj.static == 1
                max_txrx = 1;
            end
            for txrx = 1:max_txrx
                state = obj.getProcedureAState(txrx);
                tf = tf && ~isempty(state) && isfield(state, 'time') && isfinite(state.time);
            end
        end

        function state = getProcedureAState(obj, txrx)
            state = [];
            if isempty(obj.sector) || txrx < 1 || txrx > numel(obj.sector) || ...
                    isempty(obj.ST) || isempty(obj.ST.ID) || ...
                    ~isprop(obj.sector(txrx), 'SC_procA_target_state') || isempty(obj.sector(txrx).SC_procA_target_state)
                return;
            end

            stId = obj.ST.ID;
            if stId < 1 || numel(obj.sector(txrx).SC_procA_target_state) < stId
                return;
            end
            candidate = obj.sector(txrx).SC_procA_target_state(stId);
            if isstruct(candidate) && isfield(candidate, 'initialized') && candidate.initialized
                state = candidate;
            end
        end

        function saveProcedureAState(obj)
            max_txrx = obj.TXRX;
            if obj.static == 1
                max_txrx = 1;
            end
            for txrx = 1:max_txrx
                if isempty(obj.sector) || txrx > numel(obj.sector) || ...
                        isempty(obj.ST) || isempty(obj.ST.ID) || ...
                        ~isprop(obj.sector(txrx), 'SC_procA_target_state')
                    continue;
                end

                stId = obj.ST.ID;
                state = obj.makeProcedureAState(txrx);
                if isempty(obj.sector(txrx).SC_procA_target_state)
                    state_table = repmat(obj.emptyProcedureAState(), 1, stId);
                else
                    state_table = obj.alignProcedureAStateTable(obj.sector(txrx).SC_procA_target_state);
                end
                if numel(state_table) < stId
                    state_table(end+1:stId) = repmat(obj.emptyProcedureAState(), 1, stId - numel(state_table));
                end
                state_table(stId) = state;
                obj.sector(txrx).SC_procA_target_state = state_table;
                obj.sector(txrx).SC_time_nodes = unique([obj.sector(txrx).SC_time_nodes, obj.t]);
            end
        end

        function state = makeProcedureAState(obj, txrx)
            N_txrx = obj.N_new(txrx);
            state = obj.emptyProcedureAState();
            state.initialized = true;
            state.time = obj.t;
            state.sector_position = obj.sector(txrx).Position;
            state.st_position = obj.ST.Position;
            state.tau_n = obj.tau_n(txrx,1:N_txrx);
            if isempty(obj.tau_absolute)
                state.tau_absolute = obj.d_3D(txrx)/3e8 + state.tau_n;
            else
                state.tau_absolute = obj.tau_absolute(txrx,1:N_txrx);
            end
            state.tau_prime = obj.tau_prime(txrx,1:N_txrx);
            state.tau_order = obj.tau_order(txrx,1:N_txrx);
            state.phi_prime_AOA = obj.procedureAPrimeFromRays(obj.phi_n_m_AOA, txrx, obj.phi_LOS_AOA(txrx), true);
            state.phi_prime_AOD = obj.procedureAPrimeFromRays(obj.phi_n_m_AOD, txrx, obj.phi_LOS_AOD(txrx), true);
            state.theta_prime_ZOA = obj.procedureAPrimeFromRays(obj.theta_n_m_ZOA, txrx, obj.theta_LOS_ZOA(txrx), false);
            state.theta_prime_ZOD = obj.procedureAPrimeFromRays(obj.theta_n_m_ZOD, txrx, obj.theta_LOS_ZOD(txrx), false);
            state.phi_LOS_AOA = obj.phi_LOS_AOA(txrx);
            state.phi_LOS_AOD = obj.phi_LOS_AOD(txrx);
            state.theta_LOS_ZOA = obj.theta_LOS_ZOA(txrx);
            state.theta_LOS_ZOD = obj.theta_LOS_ZOD(txrx);
            state.cluster_shadow_dB = obj.cluster_shadow_dB(txrx,1:N_txrx);
            previous_state = obj.getProcedureAState(txrx);
            state.procedureA_Xn = obj.procedureAXnFromState(previous_state, N_txrx, txrx);
            state.procedureA_Xn_correlation_distance = obj.procedureAXnCorrelationDistance();
        end

        function state = emptyProcedureAState(obj) %#ok<MANU>
            state = struct( ...
                'initialized', false, ...
                'time', nan, ...
                'sector_position', nan(1,3), ...
                'st_position', nan(1,3), ...
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

        function prime_angles = procedureAPrimeFromRays(obj, ray_angles, txrx, los_angle, is_azimuth)
            if isempty(ray_angles) || size(ray_angles, 1) < txrx
                prime_angles = [];
                return;
            end

            num_clusters = size(ray_angles, 2);
            prime_angles = nan(1, num_clusters);
            for cluster_idx = 1:num_clusters
                angles = squeeze(ray_angles(txrx, cluster_idx, :));
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

        function x_n = procedureAXnFromState(obj, state, num_clusters, txrx)
            [external_xn, is_available] = obj.getExternalProcedureAXn(num_clusters, txrx);
            if is_available
                x_n = external_xn;
            elseif ~isempty(state) && isfield(state, 'procedureA_Xn') && numel(state.procedureA_Xn) >= num_clusters
                x_n = state.procedureA_Xn(1:num_clusters);
            else
                x_n = obj.drawProcedureAXn(num_clusters);
            end
        end

        function [x_n, is_available] = getExternalProcedureAXn(obj, num_clusters, txrx)
            x_n = [];
            is_available = false;
            if isempty(obj.sector) || txrx < 1 || txrx > numel(obj.sector) || ...
                    isempty(obj.ST) || isempty(obj.ST.ID)
                return;
            end

            if txrx == 1 && isprop(obj.sector(txrx), 'SC_procA_target_tx_Xn') && ...
                    ~isempty(obj.sector(txrx).SC_procA_target_tx_Xn)
                x_n_table = obj.sector(txrx).SC_procA_target_tx_Xn;
            elseif txrx == 2 && isprop(obj.sector(txrx), 'SC_procA_target_rx_Xn') && ...
                    ~isempty(obj.sector(txrx).SC_procA_target_rx_Xn)
                x_n_table = obj.sector(txrx).SC_procA_target_rx_Xn;
            elseif isprop(obj.sector(txrx), 'SC_procA_target_Xn') && ~isempty(obj.sector(txrx).SC_procA_target_Xn)
                x_n_table = obj.sector(txrx).SC_procA_target_Xn;
            else
                return;
            end

            stId = obj.ST.ID;
            if stId < 1 || size(x_n_table, 1) < stId || size(x_n_table, 2) < num_clusters
                return;
            end

            x_n = x_n_table(stId, 1:num_clusters);
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

        function [tx_velocity, rx_velocity] = targetProcedureAVelocities(obj, txrx)
            if txrx == 1
                tx_velocity = obj.velocityVector(obj.sector(txrx));
                rx_velocity = obj.velocityVector(obj.ST);
            else
                tx_velocity = obj.velocityVector(obj.ST);
                rx_velocity = obj.velocityVector(obj.sector(txrx));
            end
        end

        function velocity_vec = velocityVector(obj, node) %#ok<INUSL>
            velocity_vec = [0, 0, 0];
            if isempty(node) || ~isprop(node, 'velocity') || isempty(node.velocity)
                return;
            end
            velocity_vec = node.velocity * [cosd(node.phi_v), sind(node.phi_v), cosd(node.theta_v)];
        end

        function [external_lsp, isAvailable] = getExternalLspRaw(obj, vectorLength, txrx)
            external_lsp = [];
            isAvailable = false;
            if isempty(obj.sector) || txrx < 1 || txrx > numel(obj.sector) || ...
                    isempty(obj.ST) || isempty(obj.ST.ID)
                return;
            end

            field_name = obj.lspRawFieldName(txrx);
            if ~isprop(obj.sector(txrx), field_name) || isempty(obj.sector(txrx).(field_name))
                return;
            end

            stId = obj.ST.ID;
            lsp_table = obj.sector(txrx).(field_name);
            if stId < 1 || stId > size(lsp_table, 2) || vectorLength > size(lsp_table, 1)
                return;
            end

            external_lsp = lsp_table(1:vectorLength, stId);
            isAvailable = isnumeric(external_lsp) && all(isfinite(external_lsp));
        end

        function field_name = lspRawFieldName(obj, txrx)
            if obj.O2I
                field_name = 'LSP_raw_O2I';
            elseif obj.bLOS(txrx)
                field_name = 'LSP_raw_LOS';
            else
                field_name = 'LSP_raw_NLOS';
            end
        end

        function warnSpatialConsistencyFallback(~)
            persistent warned
            if isempty(warned)
                warning('TargetChannel:SpatialConsistencyFallback', ...
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


        function [rn1, rn2, rn3] = drawTargetCouplingSeeds(obj, txrx)
            rn1 = rand(obj.N_new(txrx), obj.M(txrx));
            rn2 = rand(obj.N_new(txrx), obj.M(txrx));
            rn3 = rand(obj.N_new(txrx), obj.M(txrx));
        end

        function coupleTargetLinkRays(obj, txrx, rn1, rn2, rn3, rayGroups)
            for clusterIdx = 1:obj.N_new(txrx)
                if obj.targetClusterHasSubpaths(txrx, clusterIdx)
                    for groupIdx = 1:numel(rayGroups)
                        obj.sortTargetCoupledAngles(txrx, clusterIdx, rayGroups{groupIdx}, rn1, rn2, rn3);
                    end
                else
                    obj.sortTargetCoupledAngles(txrx, clusterIdx, 1:obj.M(txrx), rn1, rn2, rn3);
                end
            end
        end

        function sortTargetCoupledAngles(obj, txrx, clusterIdx, rayList, rn1, rn2, rn3)
            [~, aoaOrder] = sort(rn1(clusterIdx, rayList));
            obj.phi_n_m_AOA(txrx, clusterIdx, rayList) = obj.phi_n_m_AOA(txrx, clusterIdx, rayList(aoaOrder));

            [~, zoaOrder] = sort(rn2(clusterIdx, rayList));
            obj.theta_n_m_ZOA(txrx, clusterIdx, rayList) = obj.theta_n_m_ZOA(txrx, clusterIdx, rayList(zoaOrder));

            [~, aodOrder] = sort(rn3(clusterIdx, rayList));
            obj.phi_n_m_AOD(txrx, clusterIdx, rayList) = obj.phi_n_m_AOD(txrx, clusterIdx, rayList(aodOrder));
        end

        function tf = targetClusterHasSubpaths(obj, txrx, clusterIdx)
            tf = any(obj.strong_cluster_id(txrx, :) == obj.tau_n_keep(txrx, clusterIdx));
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
% function [xyz] = sph2cart_(phi, theta, varargin)
% if isempty(varargin)
%     r = ones(size(phi));
% else
%     r = varargin{1};
% end
% x              = r.*sind(theta).*cosd(phi);
% y              = r.*sind(theta).*sind(phi);
% z              = r.*cosd(theta);
% [n, m]         = size(phi);
% xyz(1,1:n,1:m) = x;
% xyz(2,1:n,1:m) = y;
% xyz(3,1:n,1:m) = z;
% %     xyz = [x y z];
% end
