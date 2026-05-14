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
            ensemble_width(l) = squeeze(std(Ensemble_Parameter(:,1500,l)));
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