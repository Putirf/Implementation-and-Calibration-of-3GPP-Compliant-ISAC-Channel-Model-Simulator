classdef AGV < sensing_types.Sensing_type
    % UAV sensing target (TR 38.901)
       
    methods
        function obj = AGV()
            obj.sensing_type  = 'AGV';
            obj.RCS_model  = 2;

            % ----- RCS parameters (Table 7.9.6.1-2, outdoor) model 2 -----
            % left back right front roof  
            % sigle STSP
            obj.RCS.sigma_M_dB        = -4.25;    % dBsm
            obj.RCS.sigma_M           = db2pow(obj.RCS.sigma_M_dB);    % dBsm
            obj.RCS.sigma_D_dB        = 0;        % human: deterministic
            obj.RCS.sigma_D           = db2pow(obj.RCS.sigma_D_dB);    % dBsm
            obj.RCS.sigma_sigma_S_dB  = 2.51;     % small-scale only
            obj.RCS.mu_sigma_S_dB     = -log(10) * obj.RCS.sigma_sigma_S_dB^2 / 20;

            obj.sigle_STSP.phi_center          = [90,180,270,0,nan];
            obj.sigle_STSP.phi_3dB             = [15.53,12.49,15.53,13.68,nan];
            obj.sigle_STSP.theta_center        = [75,90,75,90,0];
            obj.sigle_STSP.theta_3dB           = [20.03,11.89,20.03,13.68,11.44];
            obj.sigle_STSP.G_max               = [7.33,11.01,7.33,13.02,11.79];
            obj.sigle_STSP.sigma_max           = [17.6,21.28,17.6,23.29,22.06];
            obj.sigle_STSP.theta_range         = [30,180;30,180;30,180;30,180;0,30];
            obj.sigle_STSP.phi_range           = [-45,45;45,135;135,225;225,315;0,360];

            % multi STSP
            obj.multi_STSP.phi_center          = [90,180,270,0,nan];
            obj.multi_STSP.phi_3dB             = [15.53,12.49,15.53,13.68,nan];
            obj.multi_STSP.theta_center        = [75,90,75,90,0];
            obj.multi_STSP.theta_3dB           = [20.03,11.89,20.03,13.68,11.44];
            obj.multi_STSP.G_max               = [7.27,10.98,7.27,13,11.77];
            obj.multi_STSP.sigma_max           = [24.53,28.24,24.53,30.26,29.03];
            obj.multi_STSP.theta_range         = [0,180;0,180;0,180;0,180;0,180;];
            obj.multi_STSP.phi_range           = [0,360;0,360;0,360;0,360;0,360;];

            obj.k1                             = 12;
            obj.k2                             = 1.45;           

            % ----- Polarization -----
            obj.XPR.mu_sensing     = 9.6;          
            obj.XPR.sigma_sensing  = 6.85;
        end

        
        
    end

end
