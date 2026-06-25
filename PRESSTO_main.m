function PRESSTO_main(experiment_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRESSTO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREdictive Switching for STOchastic Control of Gene Regulatory Networks %
%
% This function runs the Predictive-Switching Control (PSC) algorithm for 
% a specified gene regulatory network configuration.
%
% Author: Christian Fernández Pérez
% Date:    2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ========================= LOAD EXPERIMENT =============================
% Add experiments folder to path
experiments_path = fullfile(pwd, 'EXPERIMENTS');
addpath(experiments_path); 

% Check if the requested experiment file exists
if nargin < 1 || isempty(experiment_name)
    error('PRESSTO:NoInput', 'You must provide an experiment name.');
end

if ~exist(experiment_name, 'file')
    error('PRESSTO:ExperimentNotFound', 'The experiment file "%s" was not found in the EXPERIMENTS folder.', experiment_name);
end

% Call the experiment function dynamically using str2func
experiment_fcn = str2func(experiment_name);
config = experiment_fcn();

grn    = config.grn;
meshxt = config.meshxt;
ic     = config.ic;
psc    = config.psc;

%% =============================== PSC ====================================
% Run PSC
fprintf('Running PSC for experiment: %s...\n', config.name);
Results = PSC(grn, meshxt, ic, psc);

%% ============================= SAVE RESULTS =============================
% Create output folder if it doesn't exist
output_folder = fullfile(pwd, 'SIMULATIONS', config.name);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Save results and configuration
fprintf('Saving results to %s...\n', output_folder);
save(fullfile(output_folder,'Results.mat'),'Results','-v7.3');
save(fullfile(output_folder,'config.mat'),'config','-v7.3');

%% ============================= PLOT RESULTS =============================
% Plot results
fprintf('Generating plots...\n');
PLOT(output_folder);

end