function SSA_SNR_Control

%%
clearvars

target_ensemble = load('') %local talapas filepath

directory_name = strcat('./gate_output_', date, '_', num2str(params.ttrendm), 'yr');
mkdir(directory_name)
save(strcat(directory_name, '/', 'control_params.mat'), 'params')

% Call the gates function with params
[h_gates, u_gates, q_gates] = gates_fixedX(height_out, velocity_out, xgs_ensem, [340, 330,325,310,300,285,275,260, 250,235, 225,210,200,185, 175, 160,150,135 ,125,110,100,85,75,60,50,40,25], params);

filename = strcat('3kyr-1k_OM2stoch_A01Stoch_500yrOM.5_GreenlandGlacierwithAblation_OnOverdeepening_gates.mat'); %desired Save name
save(strcat(directory_name, '/', filename),'-v7.3')
disp('done')