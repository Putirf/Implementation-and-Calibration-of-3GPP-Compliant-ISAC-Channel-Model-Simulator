classdef TRP_TRP_background_channel < channel.Comm_channel
    %TRP_TRP_BACKGROUND_CHANNEL Direct bistatic background channel.
    %   Implements the TRP-to-TRP Case-1 path used by TR 38.901
    %   Clause 7.9.4.2, with the gNB-gNB updates from TR 38.858 Annex A.3.

    methods
        function obj = TRP_TRP_background_channel(STX, SRX, scenario, fastfading_enable, t)
            if STX.ID == SRX.ID
                error('TrpTrpBackground:MonostaticPair', ...
                    'TRP-TRP bistatic background requires two different TRPs.');
            end

            obj@channel.Comm_channel( ...
                STX, SRX, scenario, fastfading_enable, t, 'TRP-TRP');

            % Clause 7.9.4.2 applies the absolute time of arrival from
            % Clause 7.6.9 to the direct STX-SRX background channel.
            if ~isempty(obj.tau_n_LOS)
                obj.tau_absolute = obj.d_3D/3e8 + obj.delta_tau + obj.tau_n_LOS;
            end
            if fastfading_enable && ~isempty(obj.tau_channel)
                obj.tau_channel = obj.tau_channel + obj.d_3D/3e8 + obj.delta_tau;
            end
        end
    end
end
