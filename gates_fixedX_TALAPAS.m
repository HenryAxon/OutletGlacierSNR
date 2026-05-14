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