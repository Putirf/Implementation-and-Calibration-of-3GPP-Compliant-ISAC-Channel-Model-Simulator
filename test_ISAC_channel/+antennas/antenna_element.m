classdef antenna_element < handle
    properties
        attachedPanel
        slant_angle
        pol_model = 'model-2'
        antenna_model = 'isotropic'
    end

    properties(Dependent)
        max_gain
        attachedType
        alpha
        beta
        gamma
    end

    methods
        function obj = antenna_element(varargin)
            if nargin == 0
                obj.slant_angle = 0;
            else
                obj.slant_angle = varargin{1};
            end
        end

        function lcsgain = lcs_gain(obj, phi, theta)
            type = obj(1).attachedType;
            if isempty(type) || strcmp(type, 'BS')
                horizontal = -min(12 * (phi / 65).^2, 30);
                vertical = -min(12 * ((theta - 90) / 65).^2, 30);
                attenuation = -min(-(horizontal + vertical), 30);
                lcsgain = 10.^((attenuation + obj.max_gain) / 10);
                return;
            end

            if any(strcmp(type, {'UE', 'BSInF', 'BS3D-InH'}))
                lcsgain = ones(size(phi));
                return;
            end
        end

        function gcsgain = gcs_gain(obj, phi, theta)
            [localPhi, localTheta] = obj.toLocalAngles(phi, theta);
            gcsgain = obj.lcs_gain(localPhi, localTheta);
        end

        function F_theta_phi = field_pattern(obj, phi, theta)
            localField = obj.fieldInLocalBasis(phi, theta);
            rot = obj.gcsPolarizationRotation(phi, theta);
            F_theta_phi = multiprod(rot, localField, [1 2], [1 2]);
        end

        function plt(obj)
            figure;
            subplot(221);
            zenithGain = 10 * log10(obj.lcs_gain(0, 0:1:180));
            plot(0:1:180, zenithGain, 'linewIDth', 1.5);
            grid on; xlim([0 180]);
            xlabel('Zenith angle \theta (^o)'); ylabel('Gain (dB)');

            subplot(222);
            azGain = 10 * log10(obj.lcs_gain(-180:1:180, 90));
            plot(-180:1:180, azGain, 'linewIDth', 1.5);
            grid on; xlim([-180 180]);
            xlabel('Azimuth angle \phi (^o)'); ylabel('Gain (dB)');

            theta = linspace(0, pi, 200);
            phi = linspace(0, 2*pi - 0.0001, 200);
            [phiGrid, thetaGrid] = meshgrid(phi, theta);
            gain = obj.lcs_gain(phiGrid*180/pi - 360*floor(phiGrid/pi), thetaGrid*180/pi);
            [x, y, z] = sph2cart((phiGrid*180/pi - 360*floor((phiGrid - eps)/pi))*pi/180, pi/2 - thetaGrid, gain);
            subplot(2, 2, [3 4]);
            mesh(x, y, z);
            title('Pattern of TR36.873 3D Antenna element');
        end

        function out = get.max_gain(obj)
            type = obj.attachedType;
            if isempty(type) || strcmp(type(1:min(2, end)), 'BS')
                out = 8;
            elseif strcmp(type, 'UE')
                out = 1;
            else
                error('Unknown attached type of device.');
            end
        end

        function out = get.attachedType(obj)
            out = obj.attachedPanel.attachedType;
        end

        function out = get.alpha(obj)
            out = obj.attachedPanel.alpha;
        end

        function out = get.beta(obj)
            out = obj.attachedPanel.beta;
        end

        function out = get.gamma(obj)
            out = obj.attachedPanel.gamma;
        end
    end

    methods(Access = private)
        function field = fieldInLocalBasis(obj, phi, theta)
            switch obj.pol_model
                case 'model-2'
                    field = obj.model2Field(phi, theta);
                case 'model-1'
                    field = obj.model1Field(phi, theta);
                otherwise
                    error('only support model-1 and model-2');
            end
        end

        function field = model2Field(obj, phi, theta)
            if strcmp(obj.antenna_model, 'isotropic')
                gain = 1;
            else
                gain = obj.gcs_gain(phi, theta);
            end
            field = obj.makeFieldVector(sqrt(gain) * cosd(obj.slant_angle), ...
                                        sqrt(gain) * sind(obj.slant_angle));
        end

        function field = model1Field(obj, phi, theta)
            [localPhi, localTheta] = obj.toLocalAngles(phi, theta);
            basisField = obj.makeFieldVector(sqrt(obj.gcs_gain(phi, theta)), zeros(size(phi)));
            slantRot = obj.slantRotation(localPhi, localTheta, obj.slant_angle);
            field = multiprod(slantRot, basisField, [1 2], [1 2]);
        end

        function [phiLocal, thetaLocal] = toLocalAngles(obj, phi, theta)
            xyz = obj.sphericalToCartesian(phi, theta);
            local = multiprod(obj.rotationMatrix(), xyz, [1 2], 1);
            [phiLocal, thetaLocal] = obj.cartesianToSpherical(local(1, :, :), local(2, :, :), local(3, :, :));
        end

        function rot = rotationMatrix(obj)
            rot = obj.eulerRotation(obj.alpha, obj.beta, obj.gamma);
        end

        function rot = gcsPolarizationRotation(obj, phi, theta)
            a = obj.alpha;
            b = obj.beta;
            c = obj.gamma;
            denom = sqrt(1 - (cosd(b)*cosd(c).*cosd(theta) + ...
                (sind(b)*cosd(c).*cosd(phi-a) - sind(c).*sind(phi-a)).*sind(theta)).^2);
            cosPsi = (cosd(b)*cosd(c).*sind(theta) - ...
                (sind(b)*cosd(c).*cosd(phi-a) - sind(c).*sind(phi-a)).*cosd(theta)) ./ denom;
            sinPsi = (sind(b)*cosd(c).*sind(phi-a) + sind(c).*cosd(phi-a)) ./ denom;
            rot = obj.rotationTensor(cosPsi, sinPsi);
        end
    end

    methods(Static, Access = private)
        function field = makeFieldVector(thetaComponent, phiComponent)
            field(1, :, :, :) = thetaComponent;
            field(2, :, :, :) = phiComponent;
        end

        function rot = slantRotation(phi, theta, slant)
            denom = sqrt(1 - (cosd(slant).*cosd(theta) - sind(slant).*sind(phi).*sind(theta)).^2);
            cosPsi = (cosd(slant).*sind(theta) + sind(slant).*sind(phi).*cosd(theta)) ./ denom;
            sinPsi = sind(slant).*cosd(phi) ./ denom;
            rot = antennas.antenna_element.rotationTensor(cosPsi, sinPsi);
        end

        function rot = rotationTensor(cosPsi, sinPsi)
            rot(1, 1, :, :) = cosPsi;
            rot(1, 2, :, :) = -sinPsi;
            rot(2, 1, :, :) = sinPsi;
            rot(2, 2, :, :) = cosPsi;
        end

        function xyz = sphericalToCartesian(phi, theta)
            xyz(1, 1:size(phi, 1), 1:size(phi, 2)) = sind(theta).*cosd(phi);
            xyz(2, 1:size(phi, 1), 1:size(phi, 2)) = sind(theta).*sind(phi);
            xyz(3, 1:size(phi, 1), 1:size(phi, 2)) = cosd(theta);
        end

        function [phi, theta] = cartesianToSpherical(x, y, z)
            [phi, ele] = cart2sph(x, y, z);
            phi = phi * 180 / pi;
            theta = 90 - ele * 180 / pi;
        end

        function rot = eulerRotation(alpha, beta, gamma)
            ca = cosd(alpha); sa = sind(alpha);
            cb = cosd(beta);  sb = sind(beta);
            cg = cosd(gamma); sg = sind(gamma);

            rot = [ca*cb,                  sa*cb,                  -sb;
                   ca*sb*sg - sa*cg,       sa*sb*sg + ca*cg,       cb*sg;
                   ca*sb*cg + sa*sg,       sa*sb*cg - ca*sg,       cb*cg];
        end
    end
end
