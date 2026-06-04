classdef Human < sensing_types.Sensing_type
    % Human sensing target (TR 38.901)


    methods
        function obj = Human()
            obj.sensing_type  = 'Human';
            obj.RCS_model  = 1;

            % ----- RCS parameters (Table 7.9.6.1-2, outdoor) model 1 -----
            obj.RCS.sigma_M_dB        = -1.37;    % dBsm
            obj.RCS.sigma_M           = db2pow(obj.RCS.sigma_M_dB);    % dBsm
            obj.RCS.sigma_D_dB        = 0;        % human: deterministic
            obj.RCS.sigma_D           = db2pow(obj.RCS.sigma_D_dB);    % dBsm
            obj.RCS.sigma_sigma_S_dB  = 3.94;     % small-scale only
            obj.RCS.mu_sigma_S_dB     = -log(10) * obj.RCS.sigma_sigma_S_dB^2 / 20;

            % ----- RCS parameters (Table 7.9.6.1-2, outdoor) model 2 -----
            % front back
            % sigle STSP
            obj.sigle_STSP.phi_center          = [0,180];
            obj.sigle_STSP.phi_3dB             = [216.65,216.65];
            obj.sigle_STSP.theta_center        = [90,90];
            obj.sigle_STSP.theta_3dB           = [55.7,55.7];
            obj.sigle_STSP.G_max               = [2.14,2.14];
            obj.sigle_STSP.sigma_max           = [7.7,7.7];
            obj.sigle_STSP.theta_range         = [0,180;0,180];
            obj.sigle_STSP.phi_range           = [-90,90;90,270];

            obj.k1                             = 0.5714;
            obj.k2                             = 0.1;


            % Large-scale calibration: DO NOT enable sigma_S
            obj.RCS.enable_small_scale = false;

            % ----- Polarization -----
            obj.XPR.mu_sensing     = 19.81;           % typical human
            obj.XPR.sigma_sensing  = 4.25;
        end
        
    end

end
