classdef UAV < sensing_types.Sensing_type
    % UAV sensing target (TR 38.901)
    methods
        function obj = UAV(scenario)
            obj.sensing_type  = 'UAV';
            obj.RCS_model  = 1;
            scenario.BS_ST_min_d = 0;
            scenario.UE_ST_min_d = 0;
            scenario.subname = 'AV';
            scenario.alternative = 3;
            obj.height = 200;
            
            
            % ----- RCS parameters (Table 7.9.6.1-2, outdoor) -----
            obj.RCS.sigma_M_dB        = -12.81;    % dBsm
            obj.RCS.sigma_M           = db2pow(obj.RCS.sigma_M_dB);    % dBsm
            obj.RCS.sigma_D_dB        = 0;        % human: deterministic
            obj.RCS.sigma_D           = db2pow(obj.RCS.sigma_D_dB);    % dBsm
            obj.RCS.sigma_sigma_S_dB  = 3.74;     % small-scale only
            obj.RCS.mu_sigma_S_dB     = -log(10) * obj.RCS.sigma_sigma_S_dB^2 / 20;

            % Large-scale calibration: DO NOT enable sigma_S
            obj.RCS.enable_small_scale = false;

            % ----- RCS parameters (Table 7.9.6.1-2, outdoor) model 2 -----
            % left back right front roof bottom 
            % sigle STSP
            obj.sigle_STSP.phi_center          = [90,180,270,0,nan,nan];
            obj.sigle_STSP.phi_3dB             = [7.13,10.09,7.13,14.19,nan,nan];
            obj.sigle_STSP.theta_center        = [90,90,90,90,180,0];
            obj.sigle_STSP.theta_3dB           = [8.68,11.43,8.68,16.53,4.93,4.93];
            obj.sigle_STSP.G_max               = [7.43,3.99,7.43,1.02,13.55,13.55];
            obj.sigle_STSP.sigma_max           = [14.30,10.86,14.30,7.89,20.42,20.42];
            obj.sigle_STSP.theta_range         = [45,135;45,135;45,135;45,135;135,180;0,45];
            obj.sigle_STSP.phi_range           = [45,135;135,225;225,315;-45,45;0,360;0,360];

            obj.k1                             = 6.05;
            obj.k2                             = 1.33;           

            % ----- Polarization -----
            obj.XPR.mu_sensing     = 13.75;          
            obj.XPR.sigma_sensing  = 7.07;

        end
        
        
    end

end
