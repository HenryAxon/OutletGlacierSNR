%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal-to-Noise Ratio Analysis for Stochastic Glacier Simulations
%
% This script loads precomputed stochastic glacier ensemble simulations and
% generates signal-to-noise ratio (SNR)
% figures for glacier thickness, velocity, flux, and grounding-line
% position.
%
% Author: Henry Axon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clearvars;
close all;
clc;

set(groot,'defaultLineLineWidth',3.5);

%% Physical Constants

params.year  = 3600 * 24 * 365;      % Seconds per year
params.Aglen = 4.227e-25;            % Glen flow-law softness parameter
params.nglen = 3;                    % Glen exponent
params.Bglen = params.Aglen^(-1/params.nglen);

params.m     = 1 / params.nglen;     % Sliding-law exponent
params.C     = 6e6;                  % Sliding coefficient

params.rhoi  = 917;                  % Ice density (kg m^-3)
params.rhow  = 1028;                 % Seawater density (kg m^-3)
params.g     = 9.81;                 % Gravitational acceleration (m s^-2)

%% Scaling Parameters

params.hscale = 1000;
params.wscale = 1000;
params.ascale = 0.1 / params.year;

params.uscale = ...
    (params.rhoi * params.g * params.hscale * params.ascale / params.C) ...
    ^ (1 / (params.m + 1));

params.xscale = params.uscale * params.hscale / params.ascale;
params.tscale = params.xscale / params.uscale;

params.eps = params.Bglen * ...
    ((params.uscale / params.xscale)^(1 / params.nglen)) / ...
    (2 * params.rhoi * params.g * params.hscale);

params.lambda = 1 - (params.rhoi / params.rhow);
params.transient = 0;

%% Computational Grid

params.tfinal = 0 * params.year;
params.dt = 10 * params.year;
dt_out = 10;

params.Nt = params.tfinal / params.dt;
params.Nx = 400;
params.N1 = 200;
params.sigGZ = 0.97;

sigma1 = linspace(params.sigGZ/(params.N1+0.5),params.sigGZ,params.N1);
sigma2 = linspace(params.sigGZ,1,params.Nx-params.N1+1);

params.sigma = [sigma1, sigma2(2:end)]';
params.sigma_elem = [0; ...
    (params.sigma(1:params.Nx-1)+params.sigma(2:params.Nx))./2];
params.dsigma = diff(params.sigma);

%% Climate Parameters

params.a0 = 0.8 / params.year;
params.dela = 0;
params.x_smbmid = 200e3;
params.gradscale = 50e3;

params.var_multiplier = 0;
params.x_varmid = 2e5;
params.varscale = 8e4;
params.tau_smb = 1;

params.facemelt = 20 / params.year;
params.tau_ocn = 10;

%% Dataset Selection
% Historical dataset paths removed for clarity.
% Replace these paths with your local data locations.

snr = load('PATH_TO_SNR_FILE.mat');
forcings = load('PATH_TO_FORCING_FILE.mat');

snr_q = snr;
snr_u = snr;
snr_h = snr;

%% Time Array and Colormaps

max_time = 3000;
numtimepoints = 300;
timepoints = round(linspace(1,max_time,numtimepoints));

num_colors = 29;
colorarray = [linspace(1,0,num_colors)', ...
              zeros(num_colors,1), ...
              linspace(0,1,num_colors)'];

colorarray_gate = colorarray;
cmap = colorarray_gate;

% ---- Continue with original plotting code below ----

%% ------------------------------------------------------------------------
% Visualization Settings
% -------------------------------------------------------------------------

% Thickness colormap (blue)
cmap_var1 = [
    0.75 0.85 0.95
    0.35 0.60 0.85
    0.05 0.30 0.70
];

% Velocity colormap (red)
cmap_var2 = [
    0.95 0.75 0.75
    0.85 0.35 0.35
    0.65 0.05 0.05
];

% Flux colormap (green)
cmap_var3 = [
    0.80 0.95 0.80
    0.40 0.75 0.40
    0.05 0.50 0.05
];

%% ------------------------------------------------------------------------
% Eulerian Gate Locations
% -------------------------------------------------------------------------

% Gate positions measured from the ice divide (km)

gates = [
    310 300 295 285 275 260 250 240 235 ...
    225 210 200 185 175 160 150 135 ...
    125 110 100 95 85 75 60 50 40 25
];

%% ------------------------------------------------------------------------
% Grounding-Line Statistics
% -------------------------------------------------------------------------

stdGL      = std(snr_h.gl_position(:,50:3000));
meanGL     = mean(snr_h.gl_position(:,50:3000));
meanGLPlot = mean(snr_h.gl_position(:,1:3000));

%% ------------------------------------------------------------------------
% Figure: Glacier Length Ensemble
% -------------------------------------------------------------------------

figure
tiledlayout(1,1)

nexttile
hold on

ax = gca;
set(ax,...
    'ColorOrder',colorarray_gate,...
    'FontSize',30)

plot(1:3000,...
    snr_h.gl_position(1,1:3000)*params.xscale/1e3,...
    'Color','#89B8FF',...
    'DisplayName','Ensemble Members')

plot(1:3000,...
    snr_h.gl_position(:,1:3000)*params.xscale/1e3,...
    'Color','#89B8FF',...
    'HandleVisibility','off')

plot(1:3000,...
    meanGLPlot*params.xscale/1e3,...
    'Color','r',...
    'LineWidth',3,...
    'DisplayName','Ensemble Mean')

xlabel('Time (yr)')
ylabel('Glacier Length (km)')
title('Ensemble Glacier Length Distribution')
legend()

%% ------------------------------------------------------------------------
% Grounding-Line Signal-to-Noise Ratio
% -------------------------------------------------------------------------

snrGL = snr_new_gl(snr_h.gl_position);

figure
hold on

plot(timepoints(2:300),snrGL(2:300))

xlabel('Time (yr)')
ylabel('SNR')
title('Grounding-Line Signal-to-Noise Ratio')

set(gcf,...
    'Units','normalized',...
    'OuterPosition',[0 0.04 1 0.96])

set(gca,'FontSize',30)

%% ------------------------------------------------------------------------
% Forcing Signal-to-Noise
% -------------------------------------------------------------------------

melt = forcings.melt_anomalies;

facemelt = params.facemelt;

trendm = forcings.trend_m .* facemelt + facemelt;

snr_forcing_time = snr_new_forcing(timepoints);

%% ------------------------------------------------------------------------
% Figure: SNR Evolution Through Time
% -------------------------------------------------------------------------

figure
hold on

set(gcf,...
    'Units','normalized',...
    'OuterPosition',[0 0.04 1 0.96])

set(gca,'FontSize',20)

plot(timepoints(2:300),...
    snrGL(2:300),...
    'Color',[0.7 0.7 0.7],...
    'DisplayName','Grounding Line')

% plot(timepoints(2:300),...
%     snr_forcing_time(2:300),...
%     'k',...
%     'DisplayName','Climate Forcing')

% Thickness
plot(timepoints(2:300),snrh.snrResponseh(2:300,1),...
    'Color',cmap_var1(1,:),...
    'DisplayName','Thickness (275 km)')

plot(timepoints(2:300),snrh.snrResponseh(2:300,12),...
    'Color',cmap_var1(2,:),...
    'DisplayName','Thickness (200 km)')

plot(timepoints(2:300),snrh.snrResponseh(2:300,20),...
    'Color',cmap_var1(3,:),...
    'DisplayName','Thickness (100 km)')

% Velocity
plot(timepoints(2:300),snru.snrResponseu(2:300,1),...
    'Color',cmap_var2(1,:),...
    'DisplayName','Velocity (275 km)')

plot(timepoints(2:300),snru.snrResponseu(2:300,12),...
    'Color',cmap_var2(2,:),...
    'DisplayName','Velocity (200 km)')

plot(timepoints(2:300),snru.snrResponseu(2:300,20),...
    'Color',cmap_var2(3,:),...
    'DisplayName','Velocity (100 km)')

% Flux
plot(timepoints(2:300),snrq.snrResponseq(2:300,1),...
    'Color',cmap_var3(1,:),...
    'DisplayName','Flux (275 km)')

plot(timepoints(2:300),snrq.snrResponseq(2:300,12),...
    'Color',cmap_var3(2,:),...
    'DisplayName','Flux (200 km)')

plot(timepoints(2:300),snrq.snrResponseq(2:300,20),...
    'Color',cmap_var3(3,:),...
    'DisplayName','Flux (100 km)')

xlabel('Time (yr)')
ylabel('Signal-to-Noise Ratio')
title('Signal-to-Noise Response Through Time')

legend()

%% ------------------------------------------------------------------------
% Figure: Spatial Thickness SNR
% -------------------------------------------------------------------------

figure
hold on

set(gcf,...
    'Units','normalized',...
    'OuterPosition',[0 0.04 1 0.96])

set(gca,...
    'ColorOrder',colorarray_gate,...
    'FontSize',20)

for yr = 10:20
    plot(gates,...
        snr_h.snrResponseh(yr,:),...
        'DisplayName',sprintf('%d yr',yr*10),...
        'Color',colorarray_gate(yr+8,:))
end

xlabel('Distance From Ice Divide (km)')
ylabel('Signal-to-Noise Ratio')
title('Spatial Distribution of Thickness SNR')

legend()

%% ------------------------------------------------------------------------
% Figure: Spatial Velocity SNR
% -------------------------------------------------------------------------

figure
hold on

set(gcf,...
    'Units','normalized',...
    'OuterPosition',[0 0.04 1 0.96])

set(gca,...
    'ColorOrder',colorarray_gate,...
    'FontSize',20)

for yr = 10:20
    plot(gates,...
        snr_h.snrResponseu(yr,:),...
        'DisplayName',sprintf('%d yr',yr*10),...
        'Color',colorarray_gate(yr+8,:))
end

xlabel('Distance From Ice Divide (km)')
ylabel('Signal-to-Noise Ratio')
title('Spatial Distribution of Velocity SNR')

legend()
%%
figure()
hold on
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
set(ax, 'ColorOrder', colorarray_gate)
set(ax, 'FontSize', 20)
%plot(gates, snr_h.snrResponseh(150,:),DisplayName='1500 yrs',Color=colorarray_gate(9,:))
plot(gates, snr_q.snrResponseq(1,:),DisplayName='10 yrs of trend',Color=colorarray_gate(9,:))
plot(gates, snr_q.snrResponseq(2,:),DisplayName='20 yrs of trend',Color=colorarray_gate(10,:))
plot(gates, snr_q.snrResponseq(3,:),DisplayName='30 yrs of trend',Color=colorarray_gate(11,:))

plot(gates, snr_q.snrResponseq(4,:),DisplayName='40 yrs of trend',Color=colorarray_gate(12,:))

plot(gates, snr_q.snrResponseq(5,:),DisplayName='50 yrs of trend',Color=colorarray_gate(13,:))

plot(gates, snr_q.snrResponseq(6,:),DisplayName='60 yrs of trend',Color=colorarray_gate(14,:))

plot(gates, snr_q.snrResponseq(7,:),DisplayName='70 yrs of trend',Color=colorarray_gate(15,:))

plot(gates, snr_q.snrResponseq(8,:),DisplayName='80 yrs of trend',Color=colorarray_gate(16,:))

plot(gates, snr_q.snrResponseq(9,:),DisplayName='90 yrs of trend',Color=colorarray_gate(17,:))
plot(gates, snr_q.snrResponseq(10,:),DisplayName='100 yrs of trend',Color=colorarray_gate(18,:))

plot(gates, snr_q.snrResponseq(11,:),DisplayName='110 yrs of trend',Color=colorarray_gate(19,:))
plot(gates, snr_q.snrResponseq(12,:),DisplayName='120 yrs of trend',Color=colorarray_gate(20,:))

plot(gates, snr_q.snrResponseq(13,:),DisplayName='130 yrs of trend',Color=colorarray_gate(21,:))

plot(gates, snr_q.snrResponseq(14,:),DisplayName='140 yrs of trend',Color=colorarray_gate(22,:))

plot(gates, snr_q.snrResponseq(15,:),DisplayName='150 yrs of trend',Color=colorarray_gate(23,:))

plot(gates, snr_q.snrResponseq(16,:),DisplayName='160 yrs of trend',Color=colorarray_gate(24,:))

plot(gates, snr_q.snrResponseq(17,:),DisplayName='170 yrs of trend',Color=colorarray_gate(25,:))

plot(gates, snr_q.snrResponseq(18,:),DisplayName='180 yrs of trend',Color=colorarray_gate(26,:))

plot(gates, snr_q.snrResponseq(19,:),DisplayName='190 yrs of trend',Color=colorarray_gate(27,:))
plot(gates, snr_q.snrResponseq(20,:),DisplayName='200 yrs of trend',Color=colorarray_gate(28,:))

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
set(ax,fontsize=30)
title('Signal To Noise of Ice Flux Along Glacier Length')
xlabel('Distance From Ice Divide (km)')
ylabel('SNR')
legend()

%%
figure()
hold on

plot(gates, snr_h.meanTrendh(1,:),DisplayName='10 yrs',Color=colorarray_gate(9,:))
plot(gates, snr_h.meanTrendh(2,:),DisplayName='20 yrs',Color=colorarray_gate(10,:))
plot(gates, snr_h.meanTrendh(3,:),DisplayName='30 yrs',Color=colorarray_gate(11,:))
plot(gates, snr_h.meanTrendh(4,:),DisplayName='40 yrs',Color=colorarray_gate(12,:))
plot(gates, snr_h.meanTrendh(5,:),DisplayName='50 yrs',Color=colorarray_gate(13,:))
plot(gates, snr_h.meanTrendh(6,:),DisplayName='60 yrs',Color=colorarray_gate(14,:))
plot(gates, snr_h.meanTrendh(7,:),DisplayName='70 yrs',Color=colorarray_gate(15,:))
plot(gates, snr_h.meanTrendh(8,:),DisplayName='80 yrs',Color=colorarray_gate(16,:))
plot(gates, snr_h.meanTrendh(9,:),DisplayName='90 yrs',Color=colorarray_gate(17,:))
plot(gates, snr_h.meanTrendh(10,:),DisplayName='100 yrs',Color=colorarray_gate(18,:))
plot(gates, snr_h.meanTrendh(11,:),DisplayName='110 yrs',Color=colorarray_gate(19,:))
plot(gates, snr_h.meanTrendh(12,:),DisplayName='120 yrs',Color=colorarray_gate(20,:))
plot(gates, snr_h.meanTrendh(13,:),DisplayName='130 yrs',Color=colorarray_gate(21,:))
plot(gates, snr_h.meanTrendh(14,:),DisplayName='140 yrs',Color=colorarray_gate(22,:))
plot(gates, snr_h.meanTrendh(15,:),DisplayName='150 yrs',Color=colorarray_gate(23,:))
plot(gates, snr_h.meanTrendh(16,:),DisplayName='160 yrs',Color=colorarray_gate(24,:))
plot(gates, snr_h.meanTrendh(17,:),DisplayName='170 yrs',Color=colorarray_gate(25,:))
plot(gates, snr_h.meanTrendh(18,:),DisplayName='180 yrs',Color=colorarray_gate(26,:))
plot(gates, snr_h.meanTrendh(19,:),DisplayName='190 yrs',Color=colorarray_gate(27,:))
plot(gates, snr_h.meanTrendh(20,:),DisplayName='200 yrs',Color=colorarray_gate(28,:))


title('Ensemble Mean Trend Thickness Spatially')
xlabel('Distance From Ice Divide (km)')
ylabel('Ensemble Mean Thickness (m)')
legend()

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
set(ax,fontsize=30)
%%

figure()
hold on
for i =1:27
    plot(timepoints(1:150), snr_h.snrResponseh(1:150,i), DisplayName= sprintf('%d km from ice divide', gates(i)),Color=colorarray_gate(i,:))
end

title('SNR Thickness Over Time')
xlabel('Time (yr)')
ylabel('SNR')
legend()

%%

figure()
hold on
for i =1:27
    plot(timepoints(1:150), snr_h.meanTrendh(1:150,i), DisplayName= sprintf('%d km from ice divide', gates(i)),Color=colorarray_gate(i,:))
end


title('Ensemble Mean Trend Thickness Over Time')
xlabel('Time (yr)')
ylabel('Ensemble Mean Thickness (m)')
legend()




%%
figure()
hold on

plot(gates, snr_q.meanTrendq(1, :), DisplayName='10 yrs',Color=colorarray_gate(1,:))
plot(gates, snr_q.meanTrendq(2,:),DisplayName='20 yrs',Color=colorarray_gate(2,:))
plot(gates, snr_q.meanTrendq(3,:),DisplayName='30 yrs',Color=colorarray_gate(3,:))
plot(gates, snr_q.meanTrendq(4,:),DisplayName='40 yrs',Color=colorarray_gate(4,:))
plot(gates, snr_q.meanTrendq(5,:),DisplayName='50',Color=colorarray_gate(8,:))
plot(gates, snr_q.meanTrendq(6,:),DisplayName='60 yrs',Color=colorarray_gate(9,:))
plot(gates, snr_q.meanTrendq(7,:),DisplayName='70 yrs',Color=colorarray_gate(10,:))
plot(gates, snr_q.meanTrendq(8,:),DisplayName='80 yrs',Color=colorarray_gate(11,:))
plot(gates, snr_q.meanTrendq(9,:),DisplayName='90 yrs',Color=colorarray_gate(12,:))
plot(gates, snr_q.meanTrendq(10,:),DisplayName='100 yrs',Color=colorarray_gate(13,:))
plot(gates, snr_q.meanTrendq(11,:),DisplayName='110 yrs',Color=colorarray_gate(14,:))
plot(gates, snr_q.meanTrendq(12,:),DisplayName='120 yrs',Color=colorarray_gate(15,:))
plot(gates, snr_q.meanTrendq(13,:),DisplayName='130 yrs',Color=colorarray_gate(16,:))
plot(gates, snr_q.meanTrendq(14,:),DisplayName='140 yrs',Color=colorarray_gate(17,:))
plot(gates, snr_q.meanTrendq(15,:),DisplayName='150 yrs',Color=colorarray_gate(18,:))
plot(gates, snr_q.meanTrendq(16,:),DisplayName='160 yrs',Color=colorarray_gate(19,:))
plot(gates, snr_q.meanTrendq(17,:),DisplayName='170 yrs',Color=colorarray_gate(20,:))
plot(gates, snr_q.meanTrendq(18,:),DisplayName='180 yrs',Color=colorarray_gate(18,:))
plot(gates, snr_q.meanTrendq(19,:),DisplayName='190 yrs',Color=colorarray_gate(19,:))
plot(gates, snr_q.meanTrendq(20,:),DisplayName='200 yrs',Color=colorarray_gate(20,:))

title('Ensemble Mean Trend in Flux Spatially')
xlabel('Distance From Ice Divide (km)')
ylabel('Ensemble Mean Flux (m^2/yr)')
legend()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
set(ax,fontsize=30)

%%

figure()
hold on
for i =1:27
    plot(timepoints(2:150), snr_q.snrResponseq(2:150,i), DisplayName= sprintf('%d km from ice divide', gates(i)),Color=colorarray_gate(i,:))
end

title('SNR Flux Over Time')
xlabel('Time (yr)')
ylabel('SNR')
legend()

%%

figure()
hold on
for i =1:27
    plot(timepoints(2:150), snr_q.meanTrendq(2:150,i), DisplayName= sprintf('%d km from ice divide', gates(i)),Color=colorarray_gate(i,:))
end

title('Ensemble Mean Trend Flux Over Time')
xlabel('Time (yr)')
ylabel('Ensemble Mean Flux (m^2/yr)')
legend()
%% Noise Component Plots for Ice Dynamic Variables


figure()
plot(gates, snr_h.ensemble_widthh)
title('Ensemble Width Thickness Spatially')
xlabel('Distance From Ice Divide (km)')
ylabel('Ensemble Width Thickness (m)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
set(ax,fontsize=30)

figure()
plot(gates, snr_q.ensemble_widthq)
title('Ensemble Width of Flux Spatially')
xlabel('Distance From Ice Divide (km)')
ylabel('Ensemble Width Flux (m^2/yr)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
set(ax,fontsize=30)



figure()
plot(gates, snr_u.ensemble_widthu);
title('Ensemble Width Velocity Spatially')
xlabel('Distance From Ice Divide (km)')
ylabel('Ensemble Width Velocity')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
set(ax,fontsize=30)

%%

figure()
hold on

plot(gates, snr_u.snrResponse(3, :), DisplayName='30 yrs',Color=colorarray_gate(1,:))
plot(gates, snr_u.snrResponse(5,:),DisplayName='50 yrs',Color=colorarray_gate(2,:))
plot(gates, snr_u.snrResponse(50,:),DisplayName='500 yrs',Color=colorarray_gate(3,:))
plot(gates, snr_u.snrResponse(100,:),DisplayName='1000 yrs',Color=colorarray_gate(4,:))
plot(gates, snr_u.snrResponse(140,:),DisplayName='1400 yrs',Color=colorarray_gate(8,:))
plot(gates, snr_u.snrResponse(150,:),DisplayName='1500 yrs',Color=colorarray_gate(9,:))
plot(gates, snr_u.snrResponse(160,:),DisplayName='1600 yrs',Color=colorarray_gate(10,:))
plot(gates, snr_u.snrResponse(170,:),DisplayName='1700 yrs',Color=colorarray_gate(11,:))
plot(gates, snr_u.snrResponse(180,:),DisplayName='1800 yrs',Color=colorarray_gate(12,:))
plot(gates, snr_u.snrResponse(190,:),DisplayName='1900 yrs',Color=colorarray_gate(13,:))
plot(gates, snr_u.snrResponse(200,:),DisplayName='2000 yrs',Color=colorarray_gate(14,:))
plot(gates, snr_u.snrResponse(210,:),DisplayName='2100 yrs',Color=colorarray_gate(15,:))
plot(gates, snr_u.snrResponse(220,:),DisplayName='2200 yrs',Color=colorarray_gate(16,:))
plot(gates, snr_u.snrResponse(230,:),DisplayName='2300 yrs',Color=colorarray_gate(17,:))
plot(gates, snr_u.snrResponse(240,:),DisplayName='2400 yrs',Color=colorarray_gate(18,:))
plot(gates, snr_u.snrResponse(250,:),DisplayName='2500 yrs',Color=colorarray_gate(19,:))
plot(gates, snr_u.snrResponse(300,:),DisplayName='3000 yrs',Color=colorarray_gate(20,:))

title('SNR Velocity Spatially')
xlabel('Distance From Ice Divide (km)')
ylabel('SNR')
legend()

%%

figure()
hold on
plot(gates, snr_u.meanTrendu(1, :), DisplayName='10 yrs',Color=colorarray_gate(1,:))
plot(gates, snr_u.meanTrendu(2,:),DisplayName='20 yrs',Color=colorarray_gate(2,:))
plot(gates, snr_u.meanTrendu(3,:),DisplayName='30 yrs',Color=colorarray_gate(3,:))
plot(gates, snr_u.meanTrendu(4,:),DisplayName='40 yrs',Color=colorarray_gate(4,:))
plot(gates, snr_u.meanTrendu(5,:),DisplayName='50',Color=colorarray_gate(8,:))
plot(gates, snr_u.meanTrendu(6,:),DisplayName='60 yrs',Color=colorarray_gate(9,:))
plot(gates, snr_u.meanTrendu(7,:),DisplayName='70 yrs',Color=colorarray_gate(10,:))
plot(gates, snr_u.meanTrendu(8,:),DisplayName='80 yrs',Color=colorarray_gate(11,:))
plot(gates, snr_u.meanTrendu(9,:),DisplayName='90 yrs',Color=colorarray_gate(12,:))
plot(gates, snr_u.meanTrendu(10,:),DisplayName='100 yrs',Color=colorarray_gate(13,:))
plot(gates, snr_u.meanTrendu(11,:),DisplayName='110 yrs',Color=colorarray_gate(14,:))
plot(gates, snr_u.meanTrendu(12,:),DisplayName='120 yrs',Color=colorarray_gate(15,:))
plot(gates, snr_u.meanTrendu(13,:),DisplayName='130 yrs',Color=colorarray_gate(16,:))
plot(gates, snr_u.meanTrendu(14,:),DisplayName='140 yrs',Color=colorarray_gate(17,:))
plot(gates, snr_u.meanTrendu(15,:),DisplayName='150 yrs',Color=colorarray_gate(18,:))
plot(gates, snr_u.meanTrendu(16,:),DisplayName='160 yrs',Color=colorarray_gate(19,:))
plot(gates, snr_u.meanTrendu(17,:),DisplayName='170 yrs',Color=colorarray_gate(20,:))
plot(gates, snr_u.meanTrendu(18,:),DisplayName='180 yrs',Color=colorarray_gate(18,:))
plot(gates, snr_u.meanTrendu(19,:),DisplayName='190 yrs',Color=colorarray_gate(19,:))
plot(gates, snr_u.meanTrendu(20,:),DisplayName='200 yrs',Color=colorarray_gate(20,:))

title('Ensemble Mean Velocity Spatially')
xlabel('Distance From Ice Divide (km)')
ylabel('Ensemble Mean Velocity')
legend()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
set(ax,fontsize=30)

%%

figure()
hold on
for i =1:14
    plot(timepoints(2:150), snr_u.snrResponseu(2:150,i), DisplayName= sprintf('%d km from ice divide', gates(i)),Color=colorarray_gate(i,:))
end
%plot(timepoints(2:300), snr_h(2:300,4), DisplayName= '275 km from ice divide')
%plot(timepoints(2:300), snr_h(2:300,7),DisplayName= '200 km from ice divide')
%plot(timepoints(2:300), snr_h(2:300,11),DisplayName= '125 km from ice divide')
%plot(timepoints(2:300), snr_h(2:300,14),DisplayName= '50 km from ice divide')
title('SNR Velocity Over Time (m/yr)')
ax = gca;

% Set the font size for all text in the axes (e.g., to 16 points)
ax.FontSize = 16;

% Manually reset the legend's font size to the desired smaller size (e.g., 10 points)

xlabel('Time (yr)')
ylabel('SNR')
legend()
lgd.FontSize = 4;

%%
figure()
hold on
for i =1:14
    plot(timepoints(2:150), snr_u.meanTrendu(2:150,i), DisplayName= sprintf('%d km from ice divide', gates(i)),Color=colorarray_gate(i,:))
end
%plot(timepoints(2:300), snr_h(2:300,4), DisplayName= '275 km from ice divide')
%plot(timepoints(2:300), snr_h(2:300,7),DisplayName= '200 km from ice divide')
%plot(timepoints(2:300), snr_h(2:300,11),DisplayName= '125 km from ice divide')
%plot(timepoints(2:300), snr_h(2:300,14),DisplayName= '50 km from ice divide')
title('Ensemble Mean Velocity Over Time')
ax = gca;

% Set the font size for all text in the axes (e.g., to 16 points)
ax.FontSize = 16;

% Manually reset the legend's font size to the desired smaller size (e.g., 10 points)

xlabel('Time (yr)')
ylabel('Ensemble Mean Velocity (m/yr)')
legend()
lgd.FontSize = 4;

%% Plotting the optimum in signal to noise with components on the same axes
figure()

yyaxis right
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
hold on
plot(gates(1,5:27), snr_q.ensemble_widthq(5:27,1)/max(snr_q.ensemble_widthq(5:27,1)),DisplayName='Ensemble Width')
plot(gates(5:27), snr_q.meanTrendq(50,5:27)/max(snr_q.meanTrendq(20,5:27)), DisplayName='Average Change')
xlabel('Distance From Ice Divide (km)')
ax=gca;
set(ax, 'FontSize', 30)
ylabel('Normalized Flux')
yyaxis left
%plot(gates, snr_q.meanTrendq(20,:))
ax=gca;
set(ax, 'FontSize', 30)
plot(gates(5:27), snr_q.snrResponseq(20,5:27)/max(snr_q.snrResponseq(20,5:27)), DisplayName='Signal to Noise')
xlabel('Distance From Ice Divide (km)')
ylabel('Normalized Signal to Noise Ratio')
title('Signal to Noise and its Components across Flow Line at year 200 of Trend')
legend()


%% not normalized version of plot above

figure()

yyaxis right
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
hold on
plot(gates(1,5:27), snr_q.ensemble_widthq(5:27,1),DisplayName='Ensemble Width')
plot(gates(5:27), snr_q.meanTrendq(50,5:27), DisplayName='Average Change')
xlabel('Distance From Ice Divide (km)')
ax=gca;
set(ax, 'FontSize', 30)
ylabel('Normalized Flux')
yyaxis left
%plot(gates, snr_q.meanTrendq(20,:))
ax=gca;
set(ax, 'FontSize', 30)
plot(gates(5:27), snr_q.snrResponseq(50,5:27), DisplayName='Signal to Noise')
xlabel('Distance From Ice Divide (km)')
ylabel('Normalized Signal to Noise Ratio')
title('Signal to Noise and its Components across Flow Line at year 200 of Trend')
legend()

%% Figure for SNR components side by side velocity and thickness

tiledlayout(1,3)
nexttile()
ax = gca; 
yyaxis right
hold on
plot(gates, snr_u.meanTrendu(1,:),DisplayName='Velocity 1500 yrs',Color=cmap_var2(1,:))
yyaxis left
plot(gates, snr_u.meanTrendu(50,:),DisplayName='Velocity 2000 yrs',Color=cmap_var2(2,:))
plot(gates, snr_u.meanTrendu(150,:),DisplayName='Velocity 3000 yrs',Color=cmap_var2(3,:))

%plot(gates, snr_u.ensemble_width(1,:), DisplayName='Thickness Ensemble Width',Color=cmap_var2(2,:));
%ylabel("Velocity (m/yr)")


yyaxis left
plot(gates, snr_h.meanTrendh(1,:),DisplayName='Thickness 1500 yrs',Color=cmap_var1(1,:))
hold on
plot(gates, snr_h.meanTrendh(50,:),DisplayName='Thickness 2000 yrs',Color=cmap_var1(2,:))

plot(gates, snr_h.meanTrendh(150,:),DisplayName='Thickness 3000 yrs',Color=cmap_var1(3,:))

ylabel("Thickness (m)")
%plot(gates, snr_h.ensemble_width(1,:), DisplayName='Thickness Ensemble Width',Color=cmap_var1(2,:));

title("Ensemble Mean Trend \Delta")
xlabel("Distance From Ice Divide (km)")
%ylabel("SNR")
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%set(ax, 'ColorOrder', colorarray_gate)
set(ax, 'FontSize', 30)

legend()

nexttile()
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
yyaxis right
ylabel("Velocity (m/yr)")
hold on

plot(gates, snr_u.ensemble_widthu(:,1), DisplayName='Velocity Ensemble Width',Color=cmap_var2(2,:));
ax = gca; 

yyaxis left
hold on
plot(gates, snr_h.ensemble_widthh(:,1), DisplayName='Thickness Ensemble Width',Color=cmap_var1(2,:));
%linkaxes([ax1, ax3], 'Distance From Ice Divide (km)'); 
set(ax, 'FontSize', 30)
xlabel("Distance From Ice Divide (km)")
legend() 
title("Ensemble Width \sigma")

nexttile()
%yyaxis right
hold on

plot(gates, snr_u.snrResponseu(1,:),DisplayName='Velocity 1500 yrs',Color=cmap_var2(1,:))

plot(gates, snr_u.snrResponseu(50,:),DisplayName='Velocity 2000 yrs',Color=cmap_var2(2,:))
plot(gates, snr_u.snrResponseu(150,:),DisplayName='Velocity 3000 yrs',Color=cmap_var2(3,:))

%plot(gates, snr_q.snrResponse(150,:),DisplayName='Flux 1500 yrs',Color=cmap_var2(1,:))
%plot(gates, snr_q.snrResponse(200,:),DisplayName='Flux 2000 yrs',Color=cmap_var2(2,:))
%plot(gates, snr_q.snrResponse(300,:),DisplayName='Flux 3000 yrs',Color=cmap_var2(3,:))

%yyaxis left

plot(gates, snr_h.snrResponseh(1,:),DisplayName='Thickness 1500 yrs',Color=cmap_var1(1,:))


plot(gates, snr_h.snrResponseh(50,:),DisplayName='Thickness 2000 yrs',Color=cmap_var1(2,:))

plot(gates, snr_h.snrResponseh(150,:),DisplayName='Thickness 3000 yrs',Color=cmap_var1(3,:))
title('Signal to Noise \Delta/\sigma')
xlabel("Distance From Ice Divide (km)")
ylabel("SNR")
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
%set(ax, 'ColorOrder', colorarray_gate)
set(ax, 'FontSize', 30)
legend()



%% Functions Used in Script Defined Here. 



%%
function [h_gates, u_gates, q_gates] = gates_fixedX_vec(height_out, velocity_out, xgs_ensem, gates_x, params)
% Vectorized extraction of thickness, velocity, and flux at fixed spatial
% gates for time-varying grids.
%
% Inputs:
%   height_out   : (nEnsemble x nyrs x Nx)
%   velocity_out : (nEnsemble x nyrs x Nx)
%   xgs_ensem    : (nEnsemble x nyrs x Nx)  % time-varying grid for each ensemble
%   gates_x      : (1 x Ngates) or (Ngates x 1)  % gate locations (same units as xgrid)
%   params       : struct with fields:
%                  .sigma_elem, .xscale, .Nx, .hscale, .uscale, .year
%
% Outputs:
%   h_gates (nEnsemble x nyrs x Ngates)
%   u_gates (nEnsemble x nyrs x Ngates)
%   q_gates = h_gates .* u_gates
%
% This implementation:
%  - loops over ensemble members only (nEnsemble loop)
%  - for each ensemble member computes idx_all (nyrs x Ngates) once
%  - vectorizes interpolation across all times for each gate (Ngates small)
%  - avoids per-time 'find' and inner interpreted loops

% Ensure gates_x is column
gates_x = gates_x(:);
Ngates = numel(gates_x);

nEnsemble = size(height_out, 1);
nyrs = size(height_out, 2);
Nx = size(height_out, 3);

% Preallocate outputs
h_gates = nan(nEnsemble, nyrs, Ngates);
u_gates = nan(nEnsemble, nyrs, Ngates);

% Loop over ensemble members (this is fine; heavy work is vectorized)
for n = 1:nEnsemble
    % Extract 2D arrays (nyrs x Nx)
    h_n = squeeze(height_out(n, :, :));     % nyrs x Nx
    u_n = squeeze(velocity_out(n, :, :));   % nyrs x Nx
    xg_n = squeeze(xgs_ensem(n, :, :));     % nyrs x Nx
    
    % Form physical xgrid as in your original code
    xgrid = params.sigma_elem .* xg_n;      % nyrs x Nx
    xgrid = params.xscale .* xgrid ./ 1e3;  % nyrs x Nx
    
    % Precompute idx_all: for each time (row) and each gate, the index of
    % the last grid point upstream of the gate (0..Nx)
    % idx_all(ii,j) = number of grid points with xgrid(ii,:) < gates_x(j)
    idx_all = zeros(nyrs, Ngates, 'like', xgrid);
    for j = 1:Ngates
        % logical comparison along each row; sum gives the last upstream index
        idx_all(:, j) = sum(xgrid < gates_x(j), 2);
    end
    
    % For each gate, do vectorized interpolation over all times
    for j = 1:Ngates
        idx = idx_all(:, j);               % nyrs x 1, values in 0..Nx
        mask_bad = (idx == params.Nx) | (idx == 0); % terminus past gate OR gate <= first gridpoint
        % For safe indexing, clamp idx into 1..(Nx-1) (we will mask out bad later)
        idx_clamped = min(max(idx, 1), Nx-1);
        
        % Convert (row, col) -> linear indices for column-major storage:
        % linear index of (r, c) in an (nyrs x Nx) matrix is: r + (c-1)*nyrs
        rows = (1:nyrs).';
        ind1 = rows + (idx_clamped - 1) * nyrs;    % indices for upstream point
        ind2 = rows + (idx_clamped) * nyrs;        % indices for downstream point (idx+1)
        
        % Extract x, h, u at upstream and downstream points (vectorized)
        x1 = xgrid(ind1);  x2 = xgrid(ind2);
        h1 = h_n(ind1);    h2 = h_n(ind2);
        u1 = u_n(ind1);    u2 = u_n(ind2);
        
        % Avoid divide-by-zero (shouldn't happen if xgrid strictly increasing)
        dx = x2 - x1;
        zero_dx = dx == 0;
        dx(zero_dx) = NaN;  % will produce NaN for these times
        
        % Fraction between nodes
        alpha = (gates_x(j) - x1) ./ dx;
        
        % Linear interpolation vectorized
        hvec = h1 + alpha .* (h2 - h1);
        uvec = u1 + alpha .* (u2 - u1);
        
        % Apply mask_bad -> set NaN where idx==Nx (terminus retreat) or idx==0
        hvec(mask_bad) = NaN;
        uvec(mask_bad) = NaN;
        
        % Store into outputs (ensemble n, all times, gate j)
        h_gates(n, :, j) = hvec;
        u_gates(n, :, j) = uvec;
    end
end

% Apply scaling factors (vectorized)
h_gates = h_gates .* params.hscale;
u_gates = u_gates .* params.uscale .* params.year;

% Compute fluxes
q_gates = h_gates .* u_gates;
end



%%
function [h_gates, u_gates, q_gates] = gates_fixedX(height_out, velocity_out, xgs_ensem, gates_x, params)
    % Script to take output from SSA_simple model and extract thickness, velocity,
    % and fluxes at Eulerian (spatially-fixed) "gates" along the flowline. Finds
    % grid points at each time that are closest to specified gates and linearly
    % interpolates between them.
    
    nEnsemble = size(height_out, 1);
    nyrs = size(xgs_ensem, 2);  
    h_gates = nan(nEnsemble, nyrs, length(gates_x));
    u_gates = nan(nEnsemble, nyrs, length(gates_x));
    
    % Loop over ensemble members
    for n = 1:nEnsemble
        % Extract the ensemble member's data
        h_n = squeeze(height_out(n, :, :));  % 2D array for the n-th ensemble member
        u_n = squeeze(velocity_out(n, :, :));  % 2D array for the n-th ensemble member
        xg_n = squeeze(xgs_ensem(n, :));  % 1D array for the n-th ensemble member
        
        % Time-varying grid
        xgrid = params.sigma_elem * xg_n;  % Stretch elements to physical size
        xgrid = params.xscale * xgrid' / 1e3;

        % Loop over gates and time
        for jj = 1:length(gates_x)
            for ii = 1:nyrs
                idx = find(squeeze(xgrid(ii,:)) < gates_x(jj), 1, 'last'); % Closest grid point upstream of gate location
                if idx == params.Nx  % If the terminus has retreated past the gate
                    h_gates(n, ii, jj) = NaN;
                    u_gates(n, ii, jj) = NaN;
                else
                    % Linearly interpolate between upstream and downstream grid points
                    h_gates(n, ii, jj) = h_n(ii, idx) + (h_n(ii, idx+1) - h_n(ii, idx)) / ...
                        (xgrid(ii, idx+1) - xgrid(ii, idx)) * (gates_x(jj) - xgrid(ii, idx));
                    u_gates(n, ii, jj) = u_n(ii, idx) + (u_n(ii, idx+1) - u_n(ii, idx)) / ...
                        (xgrid(ii, idx+1) - xgrid(ii, idx)) * (gates_x(jj) - xgrid(ii, idx));
                end
            end
        end
    end
    
        % Apply scaling factors
        h_gates = h_gates * params.hscale;
        u_gates = u_gates * params.uscale * params.year;
    
        % Compute fluxes
        q_gates = h_gates .* u_gates;
end






%% Needed to normalize the data to get the size of the standard deviation to be of the correct scaling

function [snrResponse, ensemble_width, meanTrend] = snr_new(Ensemble_Parameter)
    % Compute signal-to-noise ratio (SNR) over time for an ensemble
    % Ensemble_Parameter: [nEnsemble x nTime x nGates]
    
    num_gates = size(Ensemble_Parameter, 3);
    nEnsemble = size(Ensemble_Parameter, 1);
    max_time = size(Ensemble_Parameter, 2);
    numtimepoints = 300;
    timepoints = round(linspace(1, max_time, numtimepoints));% integer indices
    ensemble_width = nan(num_gates);

    snrResponse = nan(numtimepoints, num_gates);

    %for n = 1:nEnsemble
        for l = 1:num_gates
            ensemble_width(l) = squeeze(std(Ensemble_Parameter(:,500,l)));
            for j = 1:numtimepoints
                tlen = timepoints(j);
                param_n = mean(Ensemble_Parameter(:, 1:tlen, l));
                time = 1:tlen;

                %ensemble_width = std(param_n(:,500,l));

                % Fit linear trend
                p = polyfit(time, param_n, 1);
                slope = p(1);
                meanTrend(j,l) = slope*tlen;

                % Compute "trend SNR"
                snrResponse(j,l) = (slope * tlen) / ensemble_width(l);
            end
        end
    %end
end






function [snrGL] = snr_new_gl(gl_params)
    % takes in grounding line data
    % returns a signal to noise parameter computed with same method and snr_new


    numtimepoints = 300;
    timepoints = round(linspace(1, 3001, numtimepoints)); % integer indices
    snrGL = nan(numtimepoints);
    ensemble_width = std(gl_params(1:1000, 1500));

    nEnsemble = size(gl_params, 1);
    max_time = size(gl_params, 2);
    timepoints = round(linspace(1, max_time, numtimepoints));% integer indices
    

    for j = 1:numtimepoints
                tlen = timepoints(j);
                param_n = mean(gl_params(:,1:tlen));
                time = 1:tlen;
                

                %ensemble_width = std(param_n(:,500,l));

                % Fit linear trend
                p = polyfit(time, param_n, 1);
                slope = p(1);
    

                % Compute "trend SNR"
                snrGL(j) = (slope * tlen) / ensemble_width;
     end

end



function [snrResponse, ensemble_width, meanTrend] = snr_new_forcing(Ensemble_Parameter)
    % Compute signal-to-noise ratio (SNR) over time for an ensemble
    % Ensemble_Parameter: [nEnsemble x nTime x nGates]
    max_time = 3000;
    numtimepoints = 300;
    timepoints = round(linspace(1, max_time, numtimepoints));% integer indices
    

    snrResponse = nan(numtimepoints);
    meanTrend = nan(numtimepoints);
    ensemble_width = std(Ensemble_Parameter(1,1:1500));
    
    for j = 1:numtimepoints
                tlen = timepoints(j);
                
                time = 1:tlen;

                % Fit linear trend
                p = polyfit(time, Ensemble_Parameter(1:tlen), 1);
                slope = p(1);
                meanTrend(j) = slope * tlen;

                % Compute "trend SNR"
                snrResponse(j) = (slope * tlen) / ensemble_width;
    end  
end


function [snrResponse, ensemble_width, meanTrend] = snr_means(Ensemble_Parameter)

    % Compute signal-to-noise ratio (SNR) over time for an ensemble
    % Signal is defined as the change in ensemble mean between timepoints
    %
    % Ensemble_Parameter dimensions:
    % [nEnsemble x nTime x nGates]

    num_gates = size(Ensemble_Parameter, 3);
    max_time = size(Ensemble_Parameter, 2);

    % Define timepoints
    numtimepoints = 300;
    timepoints = round(linspace(1, max_time, numtimepoints));

    % Preallocate
    ensemble_width = nan(num_gates,1);
    snrResponse = nan(numtimepoints-1, num_gates);
    meanTrend = nan(numtimepoints-1, num_gates);

    % Ensemble mean over ensemble dimension
    % Result: [1 x nTime x nGates]
    param_n = mean(Ensemble_Parameter, 1);

    for l = 1:num_gates

        % Noise estimate at reference time
        ensemble_width(l) = std(Ensemble_Parameter(:,1500,l));

        for j = 1:numtimepoints-1

            t1 = timepoints(j);
            t2 = timepoints(j+1);

            % Mean ensemble value at each timepoint
            mean1 = param_n(1,t1,l);
            mean2 = param_n(1,t2,l);

            % Signal = change in ensemble mean
            signal = mean2 - mean1;

            % Store raw signal
            meanTrend(j,l) = signal;

            % Signal-to-noise ratio
            snrResponse(j,l) = signal / ensemble_width(l);

        end
    end
end