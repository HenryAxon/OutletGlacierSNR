
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