function config = Example_2D_1I()
%% ========================= OUTPUT FOLDER ================================
config.name = mfilename;

%% ===================== GENE REGULATORY NETWORK ==========================
% Define a network of interacting genes with transcription, translation,
% degradation, and regulatory interactions.

% The dynamics are modeled with:
%   - Transcription rates (k_m)
%   - Translation rates (k_x)
%   - mRNA  degradation (g_m)
%   - Protein degradation (g_x)
%   - Leakage rates (eps)
%   - Inducer functions F
%   - Regulatory functions rho
%   - Input functions c (net effect of regulation + leakage)

% Number of genes
grn.n_gene = 2;

% --- Rates ---
grn.k_m = [10 10];         % transcription rates
grn.k_x = [100 100];       % translation rates
grn.g_m = [10 10];         % mRNA degradation
grn.g_x = [1 1];           % protein degradation
grn.eps = [0.1 0.1];       % leakage rates

% --- Inducer functions (F) ---
% F(I) represents the effect of the inducer on transcription
grn.F = {};
grn.theta = 0.1;
grn.mu    = 2;
grn.F{1} = @(I) 1 ./ (1 + (I ./ grn.theta(1)).^grn.mu(1));

% --- Regulatory functions (rho) ---
% Gene inhibition functions (nonlinear response, e.g., Hill-type)
grn.rho = cell(1, grn.n_gene);
grn.K = [40 40];
grn.H    = [4 4];
grn.rho{1} = @(x,u) grn.K(1)^grn.H(1) ./ ...
    (grn.K(1)^grn.H(1) + (x{2}.^grn.H(1)));
grn.rho{2} = @(x,u) grn.K(2)^grn.H(2) ./ ...
    (grn.K(2)^grn.H(2) + (x{1}.^grn.H(2)) .* grn.F{1}(u(1)));

% --- Net input functions (c) ---
% Combines regulation and leakage
grn.c = cell(1, grn.n_gene);
grn.c{1} = @(x,u) grn.rho{1}(x,u) + grn.eps(1)*(1-grn.rho{1}(x,u));
grn.c{2} = @(x,u) grn.rho{2}(x,u) + grn.eps(2)*(1-grn.rho{2}(x,u));

config.grn = grn;

%% ============================= MESH =====================================
% Spatial discretization (protein copy numbers) and temporal mesh
meshxt.Prot_mesh = [0 300 300; 0 300 300]; % [min, max, number of spatial intervals] per gene
meshxt.Time_mesh = [0 15 3000];            % [t0, t_end, number of time intervals]
meshxt.w = 1;                              % control window size

config.meshxt = meshxt;

%% ======================== INITIAL CONDITION =============================
ic.mode = 1;  % 0 = Gaussian, 1 = load external distribution

switch ic.mode
    case 1
        % Load external distribution
        PX0_folder = fullfile(pwd, 'SIMULATIONS', config.name);
        ic.PX0 = load(fullfile(PX0_folder, 'PX0.mat'), 'PX').PX;
end

config.ic = ic;

%% ======================== PSC PARAMETERS ================================
% Activation level = 1 - alpha_i (one per inducer)
psc.alpha = 0.01;

% Cost function 
target_point = [6 96];
psc.cost_f = @(PX) PX(target_point(1),target_point(2))-PX(target_point(2),target_point(1));

% Optimization mode: 0 = minimize, 1 = maximize 
psc.opt = 1;

config.psc = psc;

end
