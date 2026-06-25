function config = Example_1D()
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
grn.n_gene = 1;

% --- Rates ---
grn.k_m = 20;         % transcription rates
grn.k_x = 40;         % translation rates
grn.g_m = 15;         % mRNA degradation
grn.g_x = 1;          % protein degradation
grn.eps = 0.15;       % leakage rates

% --- Inducer functions (F) ---
% F(I) represents the effect of the inducer on transcription
grn.F = cell(1, grn.n_gene);
grn.theta = 0.1;
grn.mu   = 4;
grn.F{1} = @(I) 1 ./ (1 + (I ./ grn.theta(1)).^grn.mu(1));

% --- Regulatory functions (rho) ---
% Gene inhibition functions (nonlinear response, e.g., Hill-type)
grn.rho = cell(1, grn.n_gene);
grn.K = 45;
grn.H = 4;
grn.rho{1} = @(x,u) grn.K(1)^grn.H(1) ./ ...
    (grn.K(1)^grn.H(1) + (x{1}.^grn.H(1)).*(grn.F{1}(u(1))));

% --- Net input functions (c) ---
% Combines regulation and leakage
grn.c = cell(1, grn.n_gene);
grn.c{1} = @(x,u) grn.rho{1}(x,u) + grn.eps(1)*(1-grn.rho{1}(x,u));

config.grn = grn;

%% ============================= MESH =====================================
% Spatial discretization (protein copy numbers) and temporal mesh
meshxt.Prot_mesh = [0 150 600]; % [min, max, number of spatial intervals] per gene
meshxt.Time_mesh = [0 2 400];   % [t0, t_end, number of time intervals]
meshxt.w = 10;                  % Control horizon window size

config.meshxt = meshxt;

%% ======================= INITIAL CONDITION ==============================
ic.mode = 0;  % 0 = Gaussian, 1 = load external distribution

switch ic.mode
    case 0
        % Gaussian IC
        ic.mean = 100;
        ic.var  = 5;
end

config.ic = ic;

%% ======================== PSC PARAMETERS ================================
% Activation level = 1 - alpha_i (one per inducer)
psc.alpha = 0.05;
 
% Cost function
target_mean = 54;
x = linspace(meshxt.Prot_mesh(1), meshxt.Prot_mesh(2), meshxt.Prot_mesh(3) + 1);

psc.cost_f = @(PX) (sum(x(:) .* PX(:)) / sum(PX(:)) - target_mean).^2;

% Optimization mode: 0 = minimize, 1 = maximize
psc.opt = 0;

config.psc = psc;

end