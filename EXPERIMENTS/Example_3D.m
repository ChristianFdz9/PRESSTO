function config = Example_3D()
%% ========================= OUTPUT FOLDER ================================
config.name = mfilename;

%% ===================== GENE REGULATORY NETWORK ==========================
% Define a network of interacting genes with transcription, translation,
% degradation, and regulatory interactions.
%
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
grn.n_gene = 3;

% --- Rates ---
grn.k_m = [121 96 108];                 % transcription rates
grn.k_x = [75 60 110];                  % translation rates
grn.g_m = [15 10 20];                   % mRNA degradation
grn.g_x = [1 1.2 1.15];                 % protein degradation
grn.eps = [0.12 0.15 0.13];             % leakage rates

% --- Inducer functions (F) ---
% F(I) represents the effect of the inducer on transcription
grn.F = cell(1, grn.n_gene);
grn.theta = [1.4 3.2 0.89];
grn.mu    = [3 4 2];
grn.F{1} = @(I) 1 ./ (1 + (I ./ grn.theta(1)).^grn.mu(1));
grn.F{2} = @(I) 1 ./ (1 + (I ./ grn.theta(2)).^grn.mu(2));
grn.F{3} = @(I) 1 ./ (1 + (I ./ grn.theta(3)).^grn.mu(3));

% --- Regulatory functions (rho) ---
% Gene inhibition functions (nonlinear response, e.g., Hill-type)
grn.rho = cell(1, grn.n_gene);
grn.K = [210 202 198];
grn.H    = [3 4 2];
grn.rho{1} = @(x,I) grn.K(1)^grn.H(1) ./ ...
    (grn.K(1)^grn.H(1) + x{2}.^grn.H(1) .* grn.F{1}(I(1)));
grn.rho{2} = @(x,I) grn.K(2)^grn.H(2) ./ ...
    (grn.K(2)^grn.H(2) + x{3}.^grn.H(2) .* grn.F{2}(I(2)));
grn.rho{3} = @(x,I) grn.K(3)^grn.H(3) ./ ...
    (grn.K(3)^grn.H(3) + x{1}.^grn.H(3) .* grn.F{3}(I(3)));

% --- Net input functions (c) ---
% Combines regulation and leakage
grn.c = cell(1, grn.n_gene);
grn.c{1} = @(x,I) grn.rho{1}(x,I) + grn.eps(1).*(1-grn.rho{1}(x,I));
grn.c{2} = @(x,I) grn.rho{2}(x,I) + grn.eps(2).*(1-grn.rho{2}(x,I));
grn.c{3} = @(x,I) grn.rho{3}(x,I) + grn.eps(3).*(1-grn.rho{3}(x,I));

config.grn = grn;

%% ============================= MESH =====================================
% Spatial discretization (protein copy numbers) and temporal mesh
meshxt.Prot_mesh = [0 700 175; 0 700 175; 0 700 175]; % [min, max, number of spatial intervals] per gene
meshxt.Time_mesh = [0 0.5 100];                       % [t0, t_end, number of time intervals]
meshxt.w = 1;                                         % control window size

config.meshxt = meshxt;

%% ======================= INITIAL CONDITION ==============================
ic.mode = 0; % 0 = Gaussian, 1 = load external distribution

switch ic.mode
    case 0
        % Gaussian IC
        ic.mean = [260 260 260];
        ic.var  = [100 100 100];
end

config.ic = ic;

%% ======================= PSC PARAMETERS =================================
% Activation level = 1 - alpha_i (one per inducer)
psc.alpha = [0.005 0.01 0.015];

% Cost function
target_point = [260, 260, 260];
dx = (meshxt.Prot_mesh(:,2)' - meshxt.Prot_mesh(:,1)') ./ meshxt.Prot_mesh(:,3)'; 
target_index = round((target_point - meshxt.Prot_mesh(:,1)') ./ dx) + 1;

psc.cost_f = @(PX) PX(target_index(1), target_index(2), target_index(3)) / max(PX(:));

% Optimization mode: 0 = minimize, 1 = maximize 
psc.opt = 1;

config.psc = psc;

end