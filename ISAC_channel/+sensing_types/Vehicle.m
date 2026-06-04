classdef Vehicle < sensing_types.Sensing_type
    % UAV sensing target (TR 38.901)
    
    methods
        function obj = Vehicle()
            obj.sensing_type  = 'Vehicle';
            obj.RCS_model  = 2;
            obj.height = 1.6;

            % ----- RCS parameters (Table 7.9.6.1-2, outdoor) model 2 -----
            % left back right front roof  
            % sigle STSP
            obj.RCS.sigma_M_dB        = 11.25;    % dBsm
            obj.RCS.sigma_M           = db2pow(obj.RCS.sigma_M_dB);    % dBsm
            obj.RCS.sigma_D_dB        = 0;        % human: deterministic
            obj.RCS.sigma_D           = db2pow(obj.RCS.sigma_D_dB);    % dBsm
            obj.RCS.sigma_sigma_S_dB  = 3.41;     % small-scale only
            obj.RCS.mu_sigma_S_dB     = -log(10) * obj.RCS.sigma_sigma_S_dB^2 / 20;

            obj.sigle_STSP.phi_center          = [90,180,270,0,nan];
            obj.sigle_STSP.phi_3dB             = [26.90,36.32,26.90,40.54,nan];
            obj.sigle_STSP.theta_center        = [79.70,79.65,79.70,71.75,0];
            obj.sigle_STSP.theta_3dB           = [44.42,36.73,44.42,29.13,18.13];
            obj.sigle_STSP.G_max               = [20.75,14.56,20.75,15.52,21.26];
            obj.sigle_STSP.sigma_max           = [13.68,7.5,13.68,8.45,14.19];
            obj.sigle_STSP.theta_range         = [30,180;30,180;30,180;30,180;0,30];
            obj.sigle_STSP.phi_range           = [45,135;135,225;225,315;-45,45;0,360];

            % multi STSP
            obj.multi_STSP.phi_center          = [90,180,270,0,nan];
            obj.multi_STSP.phi_3dB             = [26.90,36.32,26.90,40.54,nan];
            obj.multi_STSP.theta_center        = [79.70,79.65,79.70,71.75,0];
            obj.multi_STSP.theta_3dB           = [44.42,36.73,44.42,29.13,18.13];
            obj.multi_STSP.G_max               = [20.6,13.9,20.6,14.99,21.12];
            obj.multi_STSP.sigma_max           = [20.52,13.82,20.52,14.91,21.05];
            obj.multi_STSP.theta_range         = [0,180;0,180;0,180;0,180;0,180;];
            obj.multi_STSP.phi_range           = [0,360;0,360;0,360;0,360;0,360;];

            obj.k1                             = 6;
            obj.k2                             = 1.65;           

            % ----- Polarization -----
            obj.XPR.mu_sensing     = 21.12;          
            obj.XPR.sigma_sensing  = 6.88;
        end
        
        
    end

end
