%% ========================================================================
%  Amplification-Factor Analysis
% =========================================================================
%
%  DESCRIPTION
%    Analyzes the raw forcing output from a single ensemble run of a
%    flowline glacier model. The script:
%      1. Loads a forcing time series (ocean melt anomaly + linear trend)
%         and the corresponding "gate" response data (flux q, thickness h,
%         velocity u) at a set of distances from the ice divide.
%      2. Computes the forcing signal-to-noise ratio (SNR) through time.
%      3. Computes the "amplification factor" at each gate/time point,
%         defined as (response SNR) / (forcing SNR).

%    Inputs:
%    Forcings file (struct with at least):
%      melt_anomalies   - ocean melt anomaly time series
%      trend_m          - imposed linear trend component (melt)
%      accum_anomalies  - accumulation anomaly time series (optional,
%                          only needed for the accumulation plots)
%      trend_a          - imposed linear trend component (accumulation,
%                          optional)
%
%    SNR/gate file (struct with at least one of):
%      snrResponseq, snrResponseh, snrResponseu
%      each sized [nGates x nTimepoints], matching `gates` and
%      `timepoints` defined below.
%
% =========================================================================

%% CONFIGURATION -- edit these paths for your local environment
% Directory containing the model output .mat files used in this run.
dataDir = '<path-to-local-data-directory>';   % EDIT ME

% Active input files for this analysis (ocean forcing, shallow-bump
% glacier geometry, far-side / no-ablation-zone run).
snrFile      = fullfile(dataDir, 'GlacierGeometries', 'BumpShallow', ...
    '3kyr-1k_OM1stoch_A01Stoch_500yrOM2.0_Noablation_farside_shallowhump_SNR.mat');
forcingsFile = fullfile(dataDir, 'GlacierGeometries', 'BumpShallow', ...
    'OMFOrcing1.0stoch2.0Trend.mat');

% ---------------------------------------------------------------------
% Other input files previously used with this script, retained here for
% reference/reproducibility. Swap in as needed by editing snrFile /
% forcingsFile above.
% ---------------------------------------------------------------------
%   Trends_Workspace.mat
%   AccumForcingOutput.mat
%   3kyr-1k_OM2stoch_A01Stoch_500yrOM2.0_GreenlandGlacierNoAblationFarSide_SNR_{u,q,h}gates.mat
%   ForcingDataA01.0TrendHillGL.mat
%   3kyr-1k_OM2stoch_A01Stoch_500yrA01.0_GreenlandGlacierNoAblationFarSide_SNR_{h,u,q}gates.mat
%   Forcing_OM20Stoch2.0_NoHillNoAblation.mat
%   3kyr-1k_OM2stoch_A01Stoch_500yrOM2.0_GreenlandGlacierNoAblation_nohump_SNR.mat
%   ForcingOM2.0.mat / MeltTrendOm2.0Trend1.0OMstoch.mat

%% PLOT COLORMAP SETUP
% Two custom colormaps used to color the many gate/time-series lines
% below. Both currently interpolate from green->blue and are identical;
% kept as separate variables (colorarray, colorarray_gate) for
% compatibility with how they are referenced elsewhere in this script.
num_colors = 20;
green_rgb = [1, 0, 0];
blue_rgb  = [0, 0, 1];

green_channel = linspace(green_rgb(1), blue_rgb(1), num_colors)';
red_channel   = zeros(num_colors, 1);          % red channel stays at 0
blue_channel  = linspace(green_rgb(3), blue_rgb(3), num_colors)';
colorarray = [green_channel, red_channel, blue_channel];

num_colors_gate = 20;
green_channel_gate = linspace(green_rgb(1), blue_rgb(1), num_colors_gate)';
red_channel_gate   = zeros(num_colors_gate, 1);
blue_channel_gate  = linspace(green_rgb(3), blue_rgb(3), num_colors_gate)';
colorarray_gate = [green_channel_gate, red_channel_gate, blue_channel_gate];
cmap = colorarray_gate;

set(groot, 'defaultLineLineWidth', 3.5)

%% LOAD DATA
% snr/forcings loaded from the paths configured above. The same `snr`
% struct is reused for q, h, and u responses (all three fields are
% expected to live in the one SNR/gate file for this run).
snr      = load(snrFile);
forcings = load(forcingsFile);

snr_q = snr;
snr_u = snr;
snr_h = snr;

% If your SNR file instead stores the three response types as separate
% top-level arrays rather than as nested fields, use this form instead:
% snr_u = snr.snrResponseu;
% snr_h = snr.snrResponseh;
% snr_q = snr.snrResponseq;

%% GATE DISTANCES (km from ice divide)
% Distance from ice divide, in km, for each response "gate" in the
% SNR/gate data (index-matched to columns of snrResponseq/h/u).
gates = [300,295,285,275,260,250,240,235,225,210,200,185,175,160,150,135, ...
         125,110,100,95,85,75,60,50,40,25,20];

% Alternative gate spacing used in some earlier runs (kept for reference):
% gates = [275,270,265,260,255,250,245,240,235,225,210,200,185,175,160, ...
%          150,135,125,110,100,95,85,75,60,50,40,25];

%% BUILD FORCING TIME SERIES AND FORCING SNR
melt = forcings.melt_anomalies;
% accum = forcings.accum_anomalies;   % NOTE: needed for accumulation plots below

trendm = forcings.trend_m;

snr_forcing = trendm / std(melt);

melt = trendm + melt;                 % melt anomaly + imposed trend
snr_forcing_time = snr_forcing(timepoints); %#ok<NODEF> -- see next cell for `timepoints`

% NOTE: `total_trend = truncation(facemelt, melt);` and related fitting
% calls (fitting_functions on melt/accum pre/during/post windows) are
% intentionally left commented out from the original workflow. Uncomment
% and supply `facemelt`/`accum` if those outputs are needed.

%% SIGNAL-TO-NOISE RATIO OF THE FORCING
[snrResponse, ensemble_width, meanTrend] = snr_new_forcing(melt);

max_time = 3000;
numtimepoints = 300;
timepoints = round(linspace(1, max_time, numtimepoints));  % integer indices

%% AMPLIFICATION FACTORS (response SNR / forcing SNR) FOR q, h, u
amplification_factorq = amp_factor(snr_q.snrResponseq, snr_forcing_time);
amplification_factorh = amp_factor(snr_h.snrResponseh, snrResponse);
amplification_factoru = amp_factor(snr_u.snrResponseu, snrResponse);

%% PLOT: Forcing SNR vs. time
time = linspace(1, 3000, 3000);
figure()
plot(timepoints, snrResponse)
title('Forcing Signal to Noise Ratio')
xlabel('Years')
ylabel('SNR Forcing')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
set(ax, 'ColorOrder', colorarray_gate)
set(ax, 'FontSize', 20)

%% PLOT: Flux amplification factor vs. time, one line per gate
figure()
hold on
plot(timepoints, amplification_factorq(1,:),  DisplayName='310 km from ice divide')
plot(timepoints, amplification_factorq(2,:),  DisplayName='300 km from ice divide')
plot(timepoints, amplification_factorq(3,:),  DisplayName='295 km from ice divide')
plot(timepoints, amplification_factorq(5,:),  DisplayName='275 km from ice divide')
plot(timepoints, amplification_factorq(7,:),  DisplayName='250 km from ice divide')
plot(timepoints, amplification_factorq(9,:),  DisplayName='225 km from ice divide')
plot(timepoints, amplification_factorq(11,:), DisplayName='200 km from ice divide')
plot(timepoints, amplification_factorq(13,:), DisplayName='175 km from ice divide')
plot(timepoints, amplification_factorq(15,:), DisplayName='150 km from ice divide')
plot(timepoints, amplification_factorq(17,:), DisplayName='125 km from ice divide')
plot(timepoints, amplification_factorq(19,:), DisplayName='100 km from ice divide')
plot(timepoints, amplification_factorq(21,:), DisplayName='75 km from ice divide')
plot(timepoints, amplification_factorq(23,:), DisplayName='50 km from ice divide')
plot(timepoints, amplification_factorq(25,:), DisplayName='25 km from ice divide')

legend()
title('Amplification Factor For Flux in Time')
xlabel('Years')
ylabel('SNR Forcing')

%% PLOT: Flux amplification factor vs. distance from ice divide, per trend year
figure()
hold on
plot(gates, abs(amplification_factorq(1:27,160)), DisplayName='Year 100 of Trend')
plot(gates, abs(amplification_factorq(1:27,161)), DisplayName='Year 110 of Trend')
plot(gates, abs(amplification_factorq(1:27,162)), DisplayName='Year 120 of Trend')
plot(gates, abs(amplification_factorq(1:27,163)), DisplayName='Year 130 of Trend')
plot(gates, abs(amplification_factorq(1:27,164)), DisplayName='Year 140 of Trend')
plot(gates, abs(amplification_factorq(1:27,165)), DisplayName='Year 150 of Trend')
plot(gates, abs(amplification_factorq(1:27,166)), DisplayName='Year 160 of Trend')
plot(gates, abs(amplification_factorq(1:27,167)), DisplayName='Year 170 of Trend')
plot(gates, abs(amplification_factorq(1:27,168)), DisplayName='Year 180 of Trend')
plot(gates, abs(amplification_factorq(1:27,169)), DisplayName='Year 190 of Trend')
plot(gates, abs(amplification_factorq(1:27,170)), DisplayName='Year 200 of Trend')
title('Amplification of Forcing in Ice Flux')
xlabel('Distance From Ice Divide (km)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
set(ax, 'ColorOrder', colorarray_gate)
set(ax, 'FontSize', 30)
legend()
ylabel('Response SNR / Forcing SNR')

%% PLOT: Velocity amplification factor vs. distance from ice divide, per trend year
figure()
hold on
plot(gates, abs(amplification_factoru(1:27,160)), DisplayName='Year 100 of Trend')
plot(gates, abs(amplification_factoru(1:27,161)), DisplayName='Year 110 of Trend')
plot(gates, abs(amplification_factoru(1:27,162)), DisplayName='Year 120 of Trend')
plot(gates, abs(amplification_factoru(1:27,163)), DisplayName='Year 130 of Trend')
plot(gates, abs(amplification_factoru(1:27,164)), DisplayName='Year 140 of Trend')
plot(gates, abs(amplification_factoru(1:27,165)), DisplayName='Year 150 of Trend')
plot(gates, abs(amplification_factoru(1:27,166)), DisplayName='Year 160 of Trend')
plot(gates, abs(amplification_factoru(1:27,167)), DisplayName='Year 170 of Trend')
plot(gates, abs(amplification_factoru(1:27,168)), DisplayName='Year 180 of Trend')
plot(gates, abs(amplification_factoru(1:27,169)), DisplayName='Year 190 of Trend')
plot(gates, abs(amplification_factoru(1:27,170)), DisplayName='Year 200 of Trend')
title('Amplification of Forcing in Velocity')
xlabel('Distance From Ice Divide (km)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
set(ax, 'ColorOrder', colorarray_gate)
set(ax, 'FontSize', 20)
legend()
ylabel('Response SNR / Forcing SNR')

%% PLOT: Thickness amplification factor vs. distance from ice divide, per trend year
figure()
hold on
plot(gates, abs(amplification_factorh(1:27,160)), DisplayName='Year 100 of Trend')
plot(gates, abs(amplification_factorh(1:27,161)), DisplayName='Year 110 of Trend')
plot(gates, abs(amplification_factorh(1:27,162)), DisplayName='Year 120 of Trend')
plot(gates, abs(amplification_factorh(1:27,163)), DisplayName='Year 130 of Trend')
plot(gates, abs(amplification_factorh(1:27,164)), DisplayName='Year 140 of Trend')
plot(gates, abs(amplification_factorh(1:27,165)), DisplayName='Year 150 of Trend')
plot(gates, abs(amplification_factorh(1:27,166)), DisplayName='Year 160 of Trend')
plot(gates, abs(amplification_factorh(1:27,167)), DisplayName='Year 170 of Trend')
plot(gates, abs(amplification_factorh(1:27,168)), DisplayName='Year 180 of Trend')
plot(gates, abs(amplification_factorh(1:27,169)), DisplayName='Year 190 of Trend')
plot(gates, abs(amplification_factorh(1:27,170)), DisplayName='Year 200 of Trend')
title('Amplification of Forcing in Thickness')
xlabel('Distance From Ice Divide (km)')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
ax = gca;
set(ax, 'ColorOrder', colorarray_gate)
set(ax, 'FontSize', 30)
legend()
ylabel('Response SNR / Forcing SNR')

%% PLOT: Thickness amplification factor vs. time, one line per gate
figure()
hold on
plot(timepoints, amplification_factorh(1,:),  DisplayName='310 km from ice divide')
plot(timepoints, amplification_factorh(2,:),  DisplayName='300 km from ice divide')
plot(timepoints, amplification_factorh(3,:),  DisplayName='295 km from ice divide')
plot(timepoints, amplification_factorh(5,:),  DisplayName='275 km from ice divide')
plot(timepoints, amplification_factorh(7,:),  DisplayName='250 km from ice divide')
plot(timepoints, amplification_factorh(9,:),  DisplayName='225 km from ice divide')
plot(timepoints, amplification_factorh(11,:), DisplayName='200 km from ice divide')
plot(timepoints, amplification_factorh(13,:), DisplayName='175 km from ice divide')
plot(timepoints, amplification_factorh(15,:), DisplayName='150 km from ice divide')
plot(timepoints, amplification_factorh(17,:), DisplayName='125 km from ice divide')
plot(timepoints, amplification_factorh(19,:), DisplayName='100 km from ice divide')
plot(timepoints, amplification_factorh(21,:), DisplayName='75 km from ice divide')
plot(timepoints, amplification_factorh(23,:), DisplayName='50 km from ice divide')
plot(timepoints, amplification_factorh(25,:), DisplayName='25 km from ice divide')

legend()
title('Amplification Factor For Thickness in Time')
xlabel('Years')
ylabel('SNR Response / SNR Forcing')

%% PLOT: Velocity amplification factor vs. time, one line per gate
figure()
hold on
plot(timepoints, amplification_factoru(1,:),  DisplayName='310 km from ice divide')
plot(timepoints, amplification_factoru(2,:),  DisplayName='300 km from ice divide')
plot(timepoints, amplification_factoru(3,:),  DisplayName='295 km from ice divide')
plot(timepoints, amplification_factoru(5,:),  DisplayName='275 km from ice divide')
plot(timepoints, amplification_factoru(7,:),  DisplayName='250 km from ice divide')
plot(timepoints, amplification_factoru(9,:),  DisplayName='225 km from ice divide')
plot(timepoints, amplification_factoru(11,:), DisplayName='200 km from ice divide')
plot(timepoints, amplification_factoru(13,:), DisplayName='175 km from ice divide')
plot(timepoints, amplification_factoru(15,:), DisplayName='150 km from ice divide')
plot(timepoints, amplification_factoru(17,:), DisplayName='125 km from ice divide')
plot(timepoints, amplification_factoru(19,:), DisplayName='100 km from ice divide')
plot(timepoints, amplification_factoru(21,:), DisplayName='75 km from ice divide')
plot(timepoints, amplification_factoru(23,:), DisplayName='50 km from ice divide')
plot(timepoints, amplification_factoru(25,:), DisplayName='25 km from ice divide')

legend()
title('Amplification Factor For Velocity in Time')
xlabel('Years')
ylabel('SNR Response / SNR Forcing')

%% PLOT: Raw melt forcing time series
time = linspace(1, length(melt), length(melt));

figure()
plot(time, melt, DisplayName='Melt Anomaly')
hold on
ax = gca();
set(ax, 'FontSize', 30)
xlabel('Time (yrs)')
ylabel('Ocean Melt (km^3 / yr)')
title('External Climate Trend Timeseries')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%% PLOT: Raw accumulation forcing time series
% NOTE: requires `accum = forcings.accum_anomalies;` to be uncommented
% above.
figure()
plot(time, accum, DisplayName='Accumulation Anomaly')
hold on
xlabel('Time (yrs)')
ylabel('Forcing in Accumulation')
title('Model Forcing Output of Single Ensemble Member')

%% PLOT: Melt trend ramp only
figure()
plot(time, trendm, DisplayName='Melt Anomaly')
hold on
xlabel('Time (yrs)')
ylabel('Forcing in Ocean Melt Ramp')
title('Model Forcing Ramp Output of Single Ensemble Member')

%% PLOT: Accumulation trend ramp only
% NOTE: requires `trenda = forcings.trend_a;` to be uncommented above.
figure()
plot(time, trenda, DisplayName='Accumulation Anomaly Ramp')
hold on
xlabel('Time (yrs)')
ylabel('Forcing in Accumulation Ramp')
title('Model Forcing Ramp Output of Single Ensemble Member')

%% PLOT: Truncated total trend (melt + facemelt, floored at zero)
% NOTE: requires `total_trend = truncation(facemelt, melt);` to be
% uncommented above (and `facemelt` to be defined) before this cell will
% run without error.
figure()
plot(time, total_trend, DisplayName='Melt Anomaly')
hold on
xlabel('Time (yrs)')
ylabel('Forcing in Ocean Melt Ramp')
title('Model Forcing Output of Single Ensemble Member Truncated at Zero')

%% PLOT: Combined tiled summary of all forcing time series
% NOTE: this cell requires accum, trenda, and total_trend (see NOTEs
% above) to all be defined.
figure()
tiledlayout(3,2)

ax1 = nexttile;
plot(time, trendm, DisplayName='Melt Anomaly')
hold on
xlabel('Time (yrs)')
ylabel('Forcing in Ocean Melt Ramp')
title('Model Forcing Ramp Output of Single Ensemble Member')

ax2 = nexttile;
plot(time, trenda, DisplayName='Accumulation Anomaly Ramp')
hold on
xlabel('Time (yrs)')
ylabel('Forcing in Accumulation Ramp')
title('Model Forcing Ramp Output of Single Ensemble Member')

ax3 = nexttile;
plot(time, melt, DisplayName='Melt Anomaly')
hold on
xlabel('Time (yrs)')
ylabel('Forcing in Ocean Melt')
title('Model Forcing Output of Single Ensemble Member')

ax4 = nexttile;
plot(time, accum, DisplayName='Accumulation Anomaly')
hold on
xlabel('Time (yrs)')
ylabel('Forcing in Accumulation')
title('Model Forcing Output of Single Ensemble Member')

ax5 = nexttile;
plot(time, total_trend, DisplayName='Melt Anomaly')
hold on
xlabel('Time (yrs)')
ylabel('Forcing in Ocean Melt')
title('Model Forcing Output of Single Ensemble Member Truncated at Zero')


%% ========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function [fit1, fit2, fit3] = fitting_functions(forcing_timeseries, trend_start, trend_end, simlength)
    % FITTING_FUNCTIONS  Fit linear trends to a forcing time series over
    % three windows: pre-trend, during-trend, and post-trend.
    %
    %   [fit1, fit2, fit3] = fitting_functions(forcing_timeseries, ...
    %       trend_start, trend_end, simlength)
    %
    %   Inputs:
    %     forcing_timeseries - vector of forcing values over time
    %     trend_start        - time index marking the start of the trend
    %     trend_end          - time index marking the end of the trend
    %     simlength          - total simulation length (time index)
    %
    %   Outputs:
    %     fit1 - polyfit coefficients [slope, intercept] for the
    %            pre-trend window (1:trend_start)
    %     fit2 - polyfit coefficients for the during-trend window
    %            (trend_start:trend_end)
    %     fit3 - polyfit coefficients for the post-trend window
    %            (trend_end:simlength)
    time = linspace(1, length(forcing_timeseries), length(forcing_timeseries));
    fit1 = polyfit(time(1:trend_start), forcing_timeseries(1:trend_start), 1);
    fit2 = polyfit(time(trend_start:trend_end), forcing_timeseries(trend_start:trend_end), 1);
    fit3 = polyfit(time(trend_end:simlength), forcing_timeseries(trend_end:simlength), 1);
end


function a = smb(x, ap, params)
    % SMB  Surface mass balance as a function of along-flow position x.
    %
    %   a = smb(x, ap, params)
    %
    %   Inputs:
    %     x      - along-flow position(s)
    %     ap     - accumulation-pattern amplitude/parameter
    %     params - struct with fields: var_multiplier, x_varmid, varscale,
    %              a0, dela, x_smbmid, gradscale
    %
    %   Output:
    %     a - surface mass balance value(s) at x
    var_pattern = params.var_multiplier + (1 - params.var_multiplier) * 0.5 * ...
        (1 + erf((x - params.x_varmid) / params.varscale));

    a = ap * var_pattern + params.a0 + 0.5 * params.dela * ...
        (1 + erf((x - params.x_smbmid) / params.gradscale));
end


function facemelt_trun = truncation(inputs, trend)
    % TRUNCATION  Add a trend to a baseline input series and floor the
    % result at zero (melt/ablation cannot be negative).
    %
    %   facemelt_trun = truncation(inputs, trend)
    %
    %   Inputs:
    %     inputs - baseline (e.g. facemelt) time series
    %     trend  - trend/forcing time series to add
    %
    %   Output:
    %     facemelt_trun - inputs + trend, truncated so no value is < 0
    trend_applied = inputs + trend;
    facemelt_trun = max(trend_applied, 0);
end


function [snrResponse, ensemble_width, meanTrend] = snr_new_forcing(Ensemble_Parameter)
    % SNR_NEW_FORCING  Compute the signal-to-noise ratio (SNR) of a linear
    % trend fit to a forcing/response time series, evaluated at a fixed
    % set of sample points ("timepoints") between year 1 and year
    % max_time.
    %
    %   [snrResponse, ensemble_width, meanTrend] = snr_new_forcing(Ensemble_Parameter)
    %
    %   Input:
    %     Ensemble_Parameter - forcing or response time series (single
    %                          ensemble member)
    %
    %   Outputs:
    %     snrResponse    - SNR of the fitted trend at each sample time
    %                      point: (slope * elapsed_time) / ensemble_width
    %     ensemble_width - noise estimate, taken as the std. dev. of the
    %                      series over its first 1500 time steps
    %     meanTrend      - fitted trend magnitude (slope * elapsed_time)
    %                      at each sample time point
    %
    %   NOTE: `snrResponse` and `meanTrend` are preallocated with
    %   nan(numtimepoints), which (per standard MATLAB semantics) creates
    %   a numtimepoints x numtimepoints matrix rather than a 1 x
    %   numtimepoints vector. Only the first row is ever written to /
    %   used downstream. Left unchanged from the original implementation.
    max_time = 3000;
    numtimepoints = 300;
    timepoints = round(linspace(1, max_time, numtimepoints));  % integer indices

    snrResponse = nan(numtimepoints);
    meanTrend = nan(numtimepoints);
    ensemble_width = std(Ensemble_Parameter(1, 1:1500));

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


function [amplification_factor] = amp_factor(SNR_Data, SNR_Forcing)
    % AMP_FACTOR  Compute the amplification factor (response SNR /
    % forcing SNR) at every gate and time point.
    %
    %   amplification_factor = amp_factor(SNR_Data, SNR_Forcing)
    %
    %   Inputs:
    %     SNR_Data    - response SNR, indexed [time x gates] (rows =
    %                   time points, columns = gates)
    %     SNR_Forcing - forcing SNR time series, one value per time point
    %
    %   Output:
    %     amplification_factor - [gates x time] matrix where
    %       amplification_factor(:, i) = SNR_Data(i, :) ./ SNR_Forcing(i)
    time_points = length(SNR_Forcing);
    gates = length(SNR_Data(1, :));
    amplification_factor = nan(gates, time_points);

    for i = 1:length(SNR_Forcing)
        amplification_factor(:, i) = SNR_Data(i, :) ./ SNR_Forcing(i);
    end
end
