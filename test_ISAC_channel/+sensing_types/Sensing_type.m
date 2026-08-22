classdef Sensing_type < handle
    % Base class for ISAC sensing target types (TR 38.901 aligned)

    properties
        sensing_type                % string
        sensing_sub_type
        RCS_model                   % integer
        RCS
        XPR
        sigle_STSP
        multi_STSP
        height = 1.5

        k1                          % parameter of multi_STSP RCS
        k2                          % parameter of multi_STSP RCS
    end

    methods
        function obj = Sensing_type()
            % empty base constructor
        end

        % function sigma_MD_dB = RCS_MD_dB(obj,TX_pos,RX_pos,ST_pos,varargin)
        % 
        %     if ~isempty(varargin) 
        %         SP = varargin{1};
        %         if length(varargin)>1
        %             is_single_STSP = varargin{2};
        %         else
        %             is_single_STSP = true;
        %         end
        %     end
        % 
        %     % ----- direction vectors -----
        %     vi = ST_pos - TX_pos;     % incident
        %     vs = RX_pos - ST_pos;     % scattering
        % 
        %     % ----- azimuth (degree) -----
        %     phi_i = atan2d(vi(2), vi(1));
        %     phi_s = atan2d(vs(2), vs(1));
        % 
        %     phi = abs(phi_i - phi_s);
        %     % phi = min(phi, 360 - phi);     % limit to [0,180]
        % 
        %     % ----- elevation (degree) -----
        %     theta_i = atan2d(vi(3), hypot(vi(1),vi(2)));
        %     theta_s = atan2d(vs(3), hypot(vs(1),vs(2)));
        % 
        %     theta = abs(theta_i - theta_s);
        %     theta = min(theta, 360 - theta);   % limit to [0,180]
        % 
        %     % ----- full 3D angle -----
        %     cosA = dot(vi,vs) / (norm(vi)*norm(vs));
        %     cosA = max(-1, min(1, cosA));   % clamp to [-1,1]
        %     beta = acosd(cosA);            % automatically [0,180]
        % 
        % 
        %     if obj.RCS_model  == 1
        %         sigma_MD_dB = max(10*log10(obj.RCS.sigma_M_dB)-3*sin(beta/2),-Inf);
        %     elseif obj.RCS_model  == 2
        %         if is_single_STSP
        %             if theta > obj.sigle_STSP.theta_range(SP,1) && theta < obj.sigle_STSP.theta_range(SP,2) && ~isnan(obj.sigle_STSP.theta_center(SP))
        %                 sigme_V_dB = -min(12*((theta-obj.sigle_STSP.theta_center(SP))/obj.sigle_STSP.theta_3dB(SP))^2,obj.sigle_STSP.sigma_max(SP));
        %             else
        %                 sigme_V_dB = 0;
        %             end
        %             if phi > obj.sigle_STSP.phi_range(SP,1) && phi < obj.sigle_STSP.phi_range(SP,2) && ~isnan(obj.sigle_STSP.phi_center(SP))
        %                 sigme_H_dB = -min(12*((phi-obj.sigle_STSP.phi_center(SP))/obj.sigle_STSP.phi_3dB(SP))^2,obj.sigle_STSP.sigma_max(SP));
        %             else
        %                 sigme_H_dB = 0;
        %             end                
        %             sigma_MD_dB = max(obj.sigle_STSP.G_max(SP)-min(-(sigme_V_dB+sigme_H_dB),obj.sigle_STSP.sigma_max(SP))-obj.k1*sin(obj.k2*beta/2)+5*log10(cos(beta/2)),...
        %                 obj.sigle_STSP.G_max(SP)-obj.sigle_STSP.sigma_max(SP),-Inf);
        %         else
        %             if theta > obj.multi_STSP.theta_range(SP,1) && theta < obj.multi_STSP.theta_range(SP,2) && ~isnan(obj.multi_STSP.theta_center(SP))
        %                 sigme_V_dB = -min(12*((theta-obj.multi_STSP.theta_center(SP))/obj.multi_STSP.theta_3dB(SP))^2,obj.multi_STSP.sigma_max(SP));
        %             else
        %                 sigme_V_dB = 0;
        %             end
        %             if phi > obj.multi_STSP.phi_range(SP,1) && phi < obj.multi_STSP.phi_range(SP,2) && ~isnan(obj.multi_STSP.phi_center(SP))
        %                 sigme_H_dB = -min(12*((phi-obj.multi_STSP.phi_center(SP))/obj.multi_STSP.phi_3dB(SP))^2,obj.multi_STSP.sigma_max(SP));
        %             else
        %                 sigme_H_dB = 0;
        %             end                
        %             sigma_MD_dB = max(obj.multi_STSP.G_max(SP)-min(-(sigme_V_dB+sigme_H_dB),obj.multi_STSP.sigma_max(SP))-obj.k1*sin(obj.k2*beta/2)+5*log10(cos(beta/2)),...
        %                 obj.multi_STSP.G_max(SP)-obj.multi_STSP.sigma_max(SP),-Inf);
        %         end
        %     end
        % 
        % end
        % 

    end
end
