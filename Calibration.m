clc;close all;
% load("Calibration_parameter.mat");
% load("Calibration_parameter_6GHZ.mat");
% load("Calibration_parameter_30GHZ.mat");

load("Calibration_parameter_6GHz_Background.mat");
load("Calibration_parameter_6GHz_Target.mat");
load("Calibration_parameter_30GHZ_Target.mat");
load("Calibration_parameter_30GHZ_Background.mat");

%% 6GHZ
lambda = 3e8/6e9;
couplingloss_target_ls_6GHZ = -(cellfun(@sum, Cali_GHZ6_Target.PL) +  10*log10(lambda^2 /(4*pi)) -cell2mat(Cali_GHZ6_Target.RCS) + cellfun(@sum, Cali_GHZ6_Target.SF));
couplingloss_target_full_6GHZ = couplingloss_target_ls_6GHZ + cell2mat(Cali_GHZ6_Target.H);
DS_target_6GHZ  = calc_delay_spread(Cali_GHZ6_Target.p_path,Cali_GHZ6_Target.tau);
ASD_target_6GHZ = calc_angular_spreads(Cali_GHZ6_Target.p_path,cellfun(@(x) x(1,:), Cali_GHZ6_Target.AOD, 'UniformOutput', false));
ASA_target_6GHZ = calc_angular_spreads(Cali_GHZ6_Target.p_path,cellfun(@(x) x(2,:), Cali_GHZ6_Target.AOA, 'UniformOutput', false));
ZSD_target_6GHZ = calc_angular_spreads(Cali_GHZ6_Target.p_path,cellfun(@(x) x(1,:), Cali_GHZ6_Target.ZOD, 'UniformOutput', false));
ZSA_target_6GHZ = calc_angular_spreads(Cali_GHZ6_Target.p_path,cellfun(@(x) x(2,:), Cali_GHZ6_Target.ZOA, 'UniformOutput', false));

couplingloss_background_ls_6GHZ = -cell2mat(Cali_GHZ6_Background.PL) + cell2mat(Cali_GHZ6_Background.SF);
couplingloss_background_full_6GHZ_tmp = couplingloss_background_ls_6GHZ + cell2mat(Cali_GHZ6_Background.H);
couplingloss_background_full_6GHZ  = pow2db(sum(reshape(db2pow(couplingloss_background_full_6GHZ_tmp), 3, []), 1));
DS_background_6GHZ = calc_delay_spread_bg(Cali_GHZ6_Background.tau,Cali_GHZ6_Background.p_path,Cali_GHZ6_Background.map_delay,Cali_GHZ6_Background.strong_cluster_id );
ASD_background_6GHZ = calc_angular_spreads_bg(Cali_GHZ6_Background.AOD,Cali_GHZ6_Background.p_path);
ASA_background_6GHZ = calc_angular_spreads_bg(Cali_GHZ6_Background.AOA,Cali_GHZ6_Background.p_path);
ZSD_background_6GHZ = calc_angular_spreads_bg(Cali_GHZ6_Background.ZOD,Cali_GHZ6_Background.p_path);
ZSA_background_6GHZ = calc_angular_spreads_bg(Cali_GHZ6_Background.ZOA,Cali_GHZ6_Background.p_path);


%% 30GHZ
lambda = 3e8/30e9;
couplingloss_target_ls_30GHZ = -(cellfun(@sum, Cali_GHZ30_Target.PL) +  10*log10(lambda^2 /(4*pi)) -cell2mat(Cali_GHZ30_Target.RCS) + cellfun(@sum, Cali_GHZ30_Target.SF));
couplingloss_target_full_30GHZ = couplingloss_target_ls_30GHZ + cell2mat(Cali_GHZ30_Target.H);
DS_target_30GHZ  = calc_delay_spread(Cali_GHZ30_Target.p_path,Cali_GHZ30_Target.tau);
ASD_target_30GHZ = calc_angular_spreads(Cali_GHZ30_Target.p_path,cellfun(@(x) x(1,:), Cali_GHZ30_Target.AOD, 'UniformOutput', false));
ASA_target_30GHZ = calc_angular_spreads(Cali_GHZ30_Target.p_path,cellfun(@(x) x(2,:), Cali_GHZ30_Target.AOA, 'UniformOutput', false));
ZSD_target_30GHZ = calc_angular_spreads(Cali_GHZ30_Target.p_path,cellfun(@(x) x(1,:), Cali_GHZ30_Target.ZOD, 'UniformOutput', false));
ZSA_target_30GHZ = calc_angular_spreads(Cali_GHZ30_Target.p_path,cellfun(@(x) x(2,:), Cali_GHZ30_Target.ZOA, 'UniformOutput', false));

couplingloss_background_ls_30GHZ = -cell2mat(Cali_GHZ30_Background.PL) + cell2mat(Cali_GHZ30_Background.SF);
couplingloss_background_full_30GHZ_tmp = couplingloss_background_ls_30GHZ + cell2mat(Cali_GHZ30_Background.H);
couplingloss_background_full_30GHZ  = pow2db(sum(reshape(db2pow(couplingloss_background_full_30GHZ_tmp), 3, []), 1));
DS_background_30GHZ = calc_delay_spread_bg(Cali_GHZ30_Background.tau,Cali_GHZ30_Background.p_path,Cali_GHZ30_Background.map_delay,Cali_GHZ30_Background.strong_cluster_id );
ASD_background_30GHZ = calc_angular_spreads_bg(Cali_GHZ30_Background.AOD,Cali_GHZ30_Background.p_path);
ASA_background_30GHZ = calc_angular_spreads_bg(Cali_GHZ30_Background.AOA,Cali_GHZ30_Background.p_path);
ZSD_background_30GHZ = calc_angular_spreads_bg(Cali_GHZ30_Background.ZOD,Cali_GHZ30_Background.p_path);
ZSA_background_30GHZ = calc_angular_spreads_bg(Cali_GHZ30_Background.ZOA,Cali_GHZ30_Background.p_path);
%% plot
plot_figure(couplingloss_target_ls_6GHZ,couplingloss_target_ls_30GHZ,"couplingloss target LS (dB)");
plot_figure(couplingloss_target_full_6GHZ,couplingloss_target_full_30GHZ,"couplingloss target Full (dB)");
plot_figure(DS_target_6GHZ,DS_target_30GHZ,"DS target (ns)");
plot_figure(ASD_target_6GHZ,ASD_target_30GHZ,"ASD target (degree)");
plot_figure(ASA_target_6GHZ,ASA_target_30GHZ,"ASA target (degree)");
plot_figure(ZSD_target_6GHZ,ZSD_target_30GHZ,"ZSD target (degree)");
plot_figure(ZSA_target_6GHZ,ZSA_target_30GHZ,"ZSA target (degree)");

plot_figure(couplingloss_background_ls_6GHZ,couplingloss_background_ls_30GHZ,"couplingloss background LS (dB)");
plot_figure(couplingloss_background_full_6GHZ,couplingloss_background_full_30GHZ,"couplingloss background Full (dB)");
plot_figure(DS_background_6GHZ,DS_background_30GHZ,"DS background");
plot_figure(ASD_background_6GHZ,ASD_background_30GHZ,"ASD background (degree)");
plot_figure(ASA_background_6GHZ,ASA_background_30GHZ,"ASA background (degree)");
plot_figure(ZSD_background_6GHZ,ZSD_background_30GHZ,"ZSD background (degree)");
plot_figure(ZSA_background_6GHZ,ZSA_background_30GHZ,"ZSA background (degree)");

%%
function [ds] = calc_delay_spread(pow,tau)
pt = cellfun(@sum, pow);
pow = cellfun(@(x,p) x./p, pow, num2cell(pt), 'UniformOutput', false);

mean_delay = cellfun(@(p,t) sum(p.*t), pow, tau);

tmp = cellfun(@(t,m) t - m, tau, num2cell(mean_delay), ...
    'UniformOutput', false);

ds = sqrt( cellfun(@(p,x) sum(p.*(x.^2)), pow, tmp) )*1e9;
end

% function [as]  = calc_angular_spreads(pow,ang)
% ang = ang(:)*pi/180;
% pt = sum(pow);
% pow = pow./pt;
% mean_angle = angle( sum( pow.*exp( 1j*ang )) ); % [rad]
% phi = ang - mean_angle;
% phi = angle( exp( 1j*phi) );
% as = sqrt( sum(pow.*(phi.^2)));
% as = as*180/pi;
% end

function as = calc_angular_spreads(pow, ang)

as = cellfun(@calc_one_as, pow, ang);

end

% ========================
function as = calc_one_as(pow, ang)

% --- Step 1: deg → rad ---
ang = ang(:) * pi/180;

% --- Step 2: normalize power ---
pt = sum(pow);
pow = pow(:) ./ pt;

% --- Step 3: circular mean ---
mean_angle = angle( sum( pow .* exp(1j*ang) ) );

% --- Step 4: wrap angle ---
phi = ang - mean_angle;
phi = angle( exp(1j*phi) );

% --- Step 5: angular spread ---
as = sqrt( sum( pow .* (phi.^2) ) );

% --- Step 6: rad → deg ---
as = as * 180/pi;

end

% function [ds] = calc_delay_spread_bg(tau, pow,map_delay,strong_cluster_id )
% taus = repmat( tau',1,20 ); 
% pt = sum( pow,2 );
% pow = pow./pt;
% pow = repmat( (pow/20).',1,20 );
% taus(strong_cluster_id,:) = taus(strong_cluster_id,:)+repmat(map_delay,numel(strong_cluster_id),1);
% pow = pow(:);
% taus = taus(:);
% pow = pow/sum(pow);
% mean_delay = sum( pow.*taus );
% ds = sqrt( sum(pow.*((taus).^2)) - mean_delay.^2 );
% end

function ds = calc_delay_spread_bg(tau, pow, map_delay, strong_cluster_id)
ds = cellfun(@calc_one_ds, tau, pow, map_delay, strong_cluster_id)*1e9;
end
% ========================
function ds = calc_one_ds(tau, pow, map_delay, strong_cluster_id)
% --- Step 1: tau expand ---
taus = repmat(tau(:), 1, 20);
% --- Step 2: power normalize ---
pt = sum(pow);
pow = pow / pt;

% --- Step 3: power expand ---
pow = repmat((pow(:)/20), 1, 20);

% --- Step 4: strong cluster delay ---
if ~isempty(strong_cluster_id)
    taus(strong_cluster_id,:) = ...
        taus(strong_cluster_id,:) + ...
        repmat(map_delay(:).', numel(strong_cluster_id),1);
end

% --- Step 5: flatten ---
pow = pow(:);
taus = taus(:);

% --- Step 6: normalize again ---
pow = pow / sum(pow);

% --- Step 7: mean delay ---
mean_delay = sum(pow .* taus);

% --- Step 8: DS ---
ds = sqrt( sum(pow .* (taus.^2)) - mean_delay.^2 );

end


% function [as]  = calc_angular_spreads_bg(ang,pow)
% ang = ang(:).'*pi/180;
% pow = repmat((pow/20).',1,20).*ones(size(ang));
% pow = pow(:).';
% pt = sum(pow,2);
% pow = pow./pt;
% mean_angle = angle( sum( pow.*exp( 1j*ang ) , 2 ) ); % [rad]
% phi = ang - mean_angle;
% phi = angle( exp( 1j*phi ) );
% as = sqrt( sum(pow.*(phi.^2),2));
% as = as*180/pi;
% end


function as = calc_angular_spreads_bg(ang, pow)

as = cellfun(@calc_one_as_bg, ang, pow);

end

function as = calc_one_as_bg(ang, pow)

ang = ang(:).';              % row vector
ang = ang * pi/180;          % deg -> rad

pow = repmat((pow(:)/20).', 1, 20) .* ones(size(ang));
pow = pow(:).';
pow = pow / sum(pow);

mean_angle = angle(sum(pow .* exp(1j*ang), 2));   % [rad]

phi = ang - mean_angle;
phi = angle(exp(1j*phi));

as = sqrt(sum(pow .* (phi.^2), 2));
as = as * 180/pi;

end


function plot_figure(oursData_6GHz,oursData_30GHz,metricName)
figure();
clf; hold on; grid on;
set(gca,'FontSize',12);

styleMap.GHZ6_simColor   = [0.10 0.45 0.85];
styleMap.GHZ30_simColor   = [0.90 0.45 0.10];
oursData_6GHz = oursData_6GHz(~isnan(oursData_6GHz));
x_ours_6GHz = sort(oursData_6GHz);
N = numel(x_ours_6GHz);
f_ours_6GHz = (1:N).' / N*100;

plot(x_ours_6GHz, f_ours_6GHz, '-','Color', styleMap.GHZ6_simColor,'LineWidth', 1.8);

oursData_30GHz = oursData_30GHz(~isnan(oursData_30GHz));
x_ours_30GHz = sort(oursData_30GHz);
N = numel(x_ours_30GHz);
f_ours_30GHz = (1:N).' / N*100;

plot(x_ours_30GHz, f_ours_30GHz, '-','Color', styleMap.GHZ30_simColor,'LineWidth', 1.8);

xlabel(metricName);
ylabel('CDF (%)');
ylim([0 100]);
yticks(0:10:100);
legend("6GHz","30GHz", 'Location', 'best');
end